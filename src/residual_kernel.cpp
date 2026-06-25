// dMod shared residual kernel implementation.
// Math reference: residual_kernel.h header.

#include "residual_kernel.h"

#include <Rcpp.h>     // brings in Rmath via R::dnorm / R::pnorm

#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

extern "C" {
  void dgemv_(const char* trans,
              const int* m, const int* n, const double* alpha,
              const double* a, const int* lda,
              const double* x, const int* incx,
              const double* beta, double* y, const int* incy);
}

namespace dmod {

namespace {

// Convenience wrappers around R math: identical numerics to R's stats::*.
inline double phi_log(double x) { return R::dnorm(x, 0.0, 1.0, /*log=*/1); }
inline double Phi(double x)     { return R::pnorm(x, 0.0, 1.0, /*lower=*/1, /*log=*/0); }
inline double Phi_log(double x) { return R::pnorm(x, 0.0, 1.0, /*lower=*/1, /*log=*/1); }

// G(w1, w2) = phi(w1) / Phi(w2). Computed in log space then exp'd for
// numerical stability when Phi(w2) underflows.
inline double G_by_Phi(double w1, double w2) {
  return std::exp(phi_log(w1) - Phi_log(w2));
}

// Symmetric add y += A^T x where A is a per-row [n_par, n_par] block stored
// row-major (i.e. col-major [n_par*n_par, n_obs] when viewed as a flat array).
// `W` is the per-row weight; the contracted result is added into col-major
// hess_acc[k1 + k2*n_par]. Since the per-row block is symmetric in (k1, k2),
// the row-major (k1, k2) index aligns with col-major hess by symmetry.
//
// Uses dgemv when n_obs * n_par * n_par is large enough to pay the BLAS
// dispatch cost (~32 doubles); otherwise a plain loop is faster and avoids
// the workspace allocation. Reused for both d2pred and d2sigma contractions.
void contract_d2_block(int n_obs, int n_par,
                       const double* d2_block, const double* W,
                       double* hess_acc, std::vector<double>& scratch) {
  const int N2 = n_par * n_par;
  if (N2 == 0 || n_obs == 0) return;
  scratch.assign(N2, 0.0);

  if (n_obs * N2 >= 64) {
    // dgemv: y = alpha A x + beta y, A col-major [N2, n_obs], lda = N2.
    const int incx = 1;
    const double alpha = 1.0, beta = 0.0;
    dgemv_("N", &N2, &n_obs, &alpha,
           d2_block, &N2, W, &incx, &beta, scratch.data(), &incx);
  } else {
    for (int i = 0; i < n_obs; ++i) {
      const double w = W[i];
      if (w == 0.0) continue;
      const double* row = d2_block + i * N2;
      for (int j = 0; j < N2; ++j) scratch[j] += w * row[j];
    }
  }

  for (int k1 = 0; k1 < n_par; ++k1) {
    for (int k2 = 0; k2 < n_par; ++k2) {
      hess_acc[k1 + k2 * n_par] += scratch[k1 * n_par + k2];
    }
  }
}

}  // namespace


void accumulate_aloq_residual(
    int n_obs,
    int n_par,
    const double* pred,
    const double* dpred,
    const double* d2pred,
    const double* y_data,
    const double* sigma,
    const double* dsigma,
    const double* d2sigma,
    const double* lloq,
    const AccumOpts& opts,
    double& value_acc,
    double* grad_acc,
    double* hess_acc) {

  if (n_obs == 0) return;
  if (n_par < 0) throw std::runtime_error("accumulate_aloq: n_par must be >= 0");

  const bool   has_dsig    = (dsigma != nullptr) && opts.sigma_depends_on_par;
  const bool   has_d2sig   = opts.use_deriv2_exact && has_dsig
                             && (d2sigma != nullptr);
  const bool   m4beal_aloq = (opts.bloq_mode == BloqMode::M4BEAL);
  const double bessel      = opts.bessel;
  const double LOG_2PI     = std::log(2.0 * M_PI);

  // Per-row scratch (small, local).
  std::vector<double> dwr(n_par), dw0(n_par);
  std::vector<double> dlogs(n_par);

  // Per-row deriv2-exact weights (accumulated across rows).
  std::vector<double> W_d2(opts.use_deriv2_exact ? n_obs : 0, 0.0);
  std::vector<double> W_d2sig(has_d2sig ? n_obs : 0, 0.0);
  std::vector<double> scratch_d2;

  for (int i = 0; i < n_obs; ++i) {
    const double sig    = sigma[i];
    const double inv_s  = 1.0 / sig;
    const double inv_s2 = inv_s * inv_s;
    const double pred_i = pred[i];
    const double y_i    = y_data[i];

    // Unscaled wr/w0; bessel applied to both for downstream math.
    const double wr0_unscaled = (pred_i - y_i) * inv_s;
    const double wr  = bessel * wr0_unscaled;

    // Value: ALOQ likelihood contribution.
    value_acc += wr * wr + LOG_2PI + 2.0 * std::log(sig);

    // M4BEAL ALOQ-side correction term: +2 log Phi(w0).
    double w0 = 0.0, G_w0 = 0.0;
    if (m4beal_aloq) {
      const double pred_val_for_w0 = pred_i;  // w0 uses the prediction, not pred - lloq
      w0 = bessel * pred_val_for_w0 * inv_s;
      value_acc += 2.0 * Phi_log(w0);
      G_w0 = G_by_Phi(w0, w0);
    }

    // dpred/dpar and dsigma/dpar rows.
    const double* dpred_i  = dpred  + (std::size_t) i * n_par;
    const double* dsigma_i = has_dsig ? (dsigma + (std::size_t) i * n_par) : nullptr;

    // Build dwr[k] = bessel * inv_s * dpred[k] - wr * inv_s * dsigma[k];
    //       dw0[k] = bessel * inv_s * dpred[k] - w0 * inv_s * dsigma[k] (M4BEAL or deriv2-exact)
    //       dlogs[k] = inv_s * dsigma[k].
    // The bessel correction scales the dpred coefficient only; the wr factor
    // multiplying dsigma is already bessel-scaled.
    if (has_dsig) {
      for (int k = 0; k < n_par; ++k) {
        const double dx = dpred_i[k];
        const double ds = dsigma_i[k];
        dwr[k]   = bessel * inv_s * dx - wr * inv_s * ds;
        dlogs[k] = inv_s * ds;
        if (m4beal_aloq) dw0[k] = bessel * inv_s * dx - w0 * inv_s * ds;
      }
    } else {
      for (int k = 0; k < n_par; ++k) {
        dwr[k]   = bessel * inv_s * dpred_i[k];
        dlogs[k] = 0.0;
        if (m4beal_aloq) dw0[k] = dwr[k];
      }
    }

    // Gradient: 2 * (wr * dwr + dlogs); M4BEAL adds 2 * G(w0) * dw0.
    for (int k = 0; k < n_par; ++k) {
      grad_acc[k] += 2.0 * (wr * dwr[k] + dlogs[k]);
    }
    if (m4beal_aloq) {
      for (int k = 0; k < n_par; ++k) grad_acc[k] += 2.0 * G_w0 * dw0[k];
    }

    // Combine all per-row Hessian additions into a single (k1, k2) loop with
    // pre-merged coefficients. Each term contributes a fixed bilinear form
    // (dwr⊗dwr, dpred⊗dsigma symmetric, dsigma⊗dsigma, or dw0⊗dw0); the
    // per-form coefficients sum across the Hessian Parts and the optional
    // M4BEAL corrections.
    //
    //   coef_dwr2:   dwr ⊗ dwr      (always; ALOQ Part0)
    //   coef_cross:  dpred ⊗ dsigma + dsigma ⊗ dpred
    //                (ALOQ Part1 if has_dsig, plus M4BEAL cross if has_dsig)
    //   coef_dsig2:  dsigma ⊗ dsigma
    //                (ALOQ Part2 + Part3, plus M4BEAL Part1-gated tail)
    //   coef_dw02:   dw0 ⊗ dw0      (M4BEAL only)
    const double coef_dwr2 = 2.0;
    double coef_cross = 0.0;
    double coef_dsig2 = 0.0;
    double coef_dw02  = 0.0;
    if (has_dsig) {
      if (opts.aloq_part1) coef_cross += -2.0 * wr * inv_s2;
      if (opts.aloq_part2) coef_dsig2 += 4.0 * wr * wr * inv_s2;
      if (opts.aloq_part3) coef_dsig2 += -2.0 * inv_s2;
    }
    if (m4beal_aloq) {
      coef_dw02 = 2.0 * (-w0 * G_w0 - G_w0 * G_w0);
      if (has_dsig) {
        coef_cross += -2.0 * G_w0 * inv_s2;
        coef_dsig2 += 4.0 * G_w0 * w0 * inv_s2;
      }
    }

    // Single per-row Hessian double-loop. Branchless on the inner k1 loop,
    // since each coefficient defaults to 0 when its contribution is absent.
    for (int k2 = 0; k2 < n_par; ++k2) {
      const double dwr_k2 = dwr[k2];
      const double dx_k2  = has_dsig ? dpred_i[k2]  : 0.0;
      const double ds_k2  = has_dsig ? dsigma_i[k2] : 0.0;
      const double dw0_k2 = m4beal_aloq ? dw0[k2] : 0.0;
      double* h_col = hess_acc + (std::size_t) k2 * n_par;
      for (int k1 = 0; k1 < n_par; ++k1) {
        const double dwr_k1 = dwr[k1];
        const double dx_k1  = has_dsig ? dpred_i[k1]  : 0.0;
        const double ds_k1  = has_dsig ? dsigma_i[k1] : 0.0;
        const double dw0_k1 = m4beal_aloq ? dw0[k1] : 0.0;
        h_col[k1] +=
            coef_dwr2 * dwr_k1 * dwr_k2
          + coef_cross * (dx_k1 * ds_k2 + ds_k1 * dx_k2)
          + coef_dsig2 * ds_k1 * ds_k2
          + coef_dw02  * dw0_k1 * dw0_k2;
      }
    }

    // d2pred exact weight: 2 * wr/sigma from the base term; M4BEAL adds the
    // exact second-order part of its +2 log Phi(w0) correction, 2 * G(w0)/sigma.
    if (opts.use_deriv2_exact) {
      W_d2[i] = 2.0 * wr * inv_s;
      if (m4beal_aloq) W_d2[i] += 2.0 * G_w0 * inv_s;
    }
    // d2sigma exact weight: (2/sigma)(1 - wr^2) from the base term; M4BEAL adds
    // -2 * w0 * G(w0)/sigma.
    if (has_d2sig) {
      W_d2sig[i] = 2.0 * inv_s * (1.0 - wr * wr);
      if (m4beal_aloq) W_d2sig[i] += -2.0 * w0 * G_w0 * inv_s;
    }
  }  // for each row

  if (opts.use_deriv2_exact && d2pred != nullptr) {
    contract_d2_block(n_obs, n_par, d2pred, W_d2.data(), hess_acc, scratch_d2);
  }
  if (has_d2sig) {
    contract_d2_block(n_obs, n_par, d2sigma, W_d2sig.data(), hess_acc,
                      scratch_d2);
  }
}


void accumulate_bloq_residual(
    int n_obs,
    int n_par,
    const double* pred,
    const double* dpred,
    const double* d2pred,
    const double* y_data,
    const double* sigma,
    const double* dsigma,
    const double* d2sigma,
    const double* lloq,
    const AccumOpts& opts,
    double& value_acc,
    double* grad_acc,
    double* hess_acc) {

  if (n_obs == 0) return;
  if (opts.bloq_mode == BloqMode::NONE || opts.bloq_mode == BloqMode::M1) {
    return;
  }
  if (n_par < 0) throw std::runtime_error("accumulate_bloq: n_par must be >= 0");
  if (lloq == nullptr) {
    throw std::runtime_error(
        "accumulate_bloq: lloq must be non-null for BLOQ rows");
  }

  // M4 modes (M4NM/M4BEAL) require non-negative lloq.
  if (opts.bloq_mode == BloqMode::M4NM || opts.bloq_mode == BloqMode::M4BEAL) {
    for (int i = 0; i < n_obs; ++i) {
      if (y_data[i] < 0.0) {
        throw std::runtime_error(
            "accumulate_bloq: M4 method cannot handle LLOQ < 0; "
            "use M3 or exponentiate log-transformed DV.");
      }
    }
  }

  const bool   has_dsig  = (dsigma != nullptr) && opts.sigma_depends_on_par;
  const bool   has_d2sig = opts.use_deriv2_exact && has_dsig
                           && (d2sigma != nullptr);
  const bool   is_m3     = (opts.bloq_mode == BloqMode::M3);
  const double bessel    = opts.bessel;

  std::vector<double> dwr(n_par), dw0(n_par);

  // Per-row deriv2-exact weights, computed alongside the per-row math.
  std::vector<double> W_d2(opts.use_deriv2_exact ? n_obs : 0, 0.0);
  std::vector<double> W_d2sig(has_d2sig ? n_obs : 0, 0.0);
  std::vector<double> scratch_d2;

  for (int i = 0; i < n_obs; ++i) {
    const double sig    = sigma[i];
    const double inv_s  = 1.0 / sig;
    const double inv_s2 = inv_s * inv_s;
    const double pred_i = pred[i];
    const double val_i  = y_data[i];   // LOQ-substituted value; equals lloq on BLOQ rows

    // wr and w0 (bessel-scaled), matching res() output for BLOQ rows
    //   wr = (pred - val) / sigma, w0 = pred / sigma; both * bessel.
    const double wr = bessel * (pred_i - val_i) * inv_s;
    const double w0 = bessel * pred_i           * inv_s;

    // ---- Value contribution ----
    if (is_m3) {
      value_acc += -2.0 * Phi_log(-wr);
    } else {
      // M4NM/M4BEAL: -2 * log(1 - Phi(wr)/Phi(w0)) with stability fallback.
      const double obj_raw = -2.0 * std::log(1.0 - Phi(wr) / Phi(w0));
      double obj_i = obj_raw;
      if (!std::isfinite(obj_i)) {
        // Stability fallback expansion when the closed form overflows.
        const double diff_w = w0 - wr;
        const double logd   = std::log(diff_w);
        const double intercept = (logd > 0.0) ? 1.8 : (-1.9 * logd + 0.9);
        const double lin       = (logd > 0.0) ? 0.9 : 0.5;
        obj_i = intercept + lin * w0 + 0.95 * w0 * w0;
      }
      value_acc += obj_i;
    }

    // ---- dwr, dw0 ----
    const double* dpred_i  = dpred  + (std::size_t) i * n_par;
    const double* dsigma_i = has_dsig ? (dsigma + (std::size_t) i * n_par) : nullptr;

    if (has_dsig) {
      for (int k = 0; k < n_par; ++k) {
        const double dx = dpred_i[k];
        const double ds = dsigma_i[k];
        dwr[k] = bessel * inv_s * dx - wr * inv_s * ds;
        dw0[k] = bessel * inv_s * dx - w0 * inv_s * ds;
      }
    } else {
      for (int k = 0; k < n_par; ++k) {
        const double dx = bessel * inv_s * dpred_i[k];
        dwr[k] = dx;
        dw0[k] = dx;
      }
    }

    // ---- Gradient ----
    // M3:  grad += 2 G(-wr) dwr
    // M4*: grad += 2 (c1 dwr - c2 dw0 + c3 dw0)
    //   with c1 = phi(wr) / (Phi(w0) - Phi(wr))   ("stable" form for numerical safety)
    //        c2 = phi(w0) / (Phi(w0) - Phi(wr))
    //        c3 = G(w0, w0)
    double w_deriv2 = 0.0;  // per-row weight for the d2pred exact contribution
    double w_deriv2_sig = 0.0;  // per-row weight for the d2sigma exact contribution
    if (is_m3) {
      const double G_neg_wr = G_by_Phi(-wr, -wr);
      for (int k = 0; k < n_par; ++k) {
        grad_acc[k] += 2.0 * G_neg_wr * dwr[k];
      }
      if (opts.use_deriv2_exact) {
        w_deriv2 = 2.0 * G_neg_wr * inv_s;
        // d2sigma weight: 2*G(-wr) * d2wr propagates -wr/sigma onto d2sigma.
        if (has_d2sig) w_deriv2_sig = -2.0 * wr * G_neg_wr * inv_s;
      }
    } else {
      // c1 = 1 / (1/G(wr,w0) - 1/G(wr,wr))
      // c2 = 1 / (1/G(w0,w0) - 1/G(w0,wr))
      // c3 = G(w0,w0)
      // Closed form: c1 = phi(wr) / (Phi(w0) - Phi(wr)); c2 = phi(w0) / same.
      // Stable closed form, avoiding the indirect "1/(1/G - 1/G)" expression.
      const double dP   = Phi(w0) - Phi(wr);
      const double phi_wr = std::exp(phi_log(wr));
      const double phi_w0 = std::exp(phi_log(w0));
      const double c1 = phi_wr / dP;
      const double c2 = phi_w0 / dP;
      const double c3 = G_by_Phi(w0, w0);
      for (int k = 0; k < n_par; ++k) {
        grad_acc[k] += 2.0 * (c1 * dwr[k] + (c3 - c2) * dw0[k]);
      }
      if (opts.use_deriv2_exact) {
        // Both dwr and dw0 contribute inv_s * dpred to their d2pred-via-d2wr term;
        // sum the coefficients accordingly: c1 + (c3 - c2).
        w_deriv2 = 2.0 * inv_s * (c1 + c3 - c2);
        // d2sigma weight: 2*(c1*d2wr + (c3-c2)*d2w0) propagates -wr/sigma and
        // -w0/sigma respectively onto d2sigma. Net coefficient on d2sigma:
        // -2/sigma * (c1*wr + (c3 - c2)*w0).
        if (has_d2sig) {
          w_deriv2_sig = -2.0 * inv_s * (c1 * wr + (c3 - c2) * w0);
        }
      }
    }

    // ---- Hessian ----
    if (is_m3) {
      // M3 Hessian: combined single (k1,k2) loop over Parts 1/2/3.
      //   coef_dwr2:  2 * (-wr*G + G^2) * dwr ⊗ dwr        (BLOQ_part1)
      //   coef_cross: -2 * G_neg_wr * inv_s^2 * symm-cross (BLOQ_part2)
      //   coef_dsig2: +4 * G_neg_wr * wr * inv_s^2 * dsigma⊗dsigma (BLOQ_part3)
      const double G_neg_wr = G_by_Phi(-wr, -wr);
      const double coef_dwr2_m3  = opts.bloq_part1
          ? 2.0 * (-wr * G_neg_wr + G_neg_wr * G_neg_wr)
          : 0.0;
      double coef_cross_m3 = 0.0;
      double coef_dsig2_m3 = 0.0;
      if (has_dsig) {
        if (opts.bloq_part2) coef_cross_m3 = -2.0 * G_neg_wr * inv_s2;
        if (opts.bloq_part3) coef_dsig2_m3 =  4.0 * G_neg_wr * wr * inv_s2;
      }
      for (int k2 = 0; k2 < n_par; ++k2) {
        const double dwr_k2 = dwr[k2];
        const double dx_k2  = has_dsig ? dpred_i[k2]  : 0.0;
        const double ds_k2  = has_dsig ? dsigma_i[k2] : 0.0;
        double* h_col = hess_acc + (std::size_t) k2 * n_par;
        for (int k1 = 0; k1 < n_par; ++k1) {
          const double dwr_k1 = dwr[k1];
          const double dx_k1  = has_dsig ? dpred_i[k1]  : 0.0;
          const double ds_k1  = has_dsig ? dsigma_i[k1] : 0.0;
          h_col[k1] +=
              coef_dwr2_m3  * dwr_k1 * dwr_k2
            + coef_cross_m3 * (dx_k1 * ds_k2 + ds_k1 * dx_k2)
            + coef_dsig2_m3 * ds_k1 * ds_k2;
        }
      }
    } else {
      // M4* branch: exact Hessian of f = -2 log(Phi(w0) - Phi(wr)) + 2 log Phi(w0).
      // First derivatives: df = a dwr + b dw0 with a = 2 phi(wr)/D and
      // b = 2 phi(w0) (1/Phi(w0) - 1/D), D = Phi(w0) - Phi(wr). The Hessian is
      // the exact d(a dwr + b dw0): the dwr/dw0 quadratic forms plus the
      // a*d2wr + b*d2w0 chain-rule terms (whose d2pred / d2sigma parts are the
      // w_deriv2 / w_deriv2_sig weights computed above).
      const double dP     = Phi(w0) - Phi(wr);
      const double phi_wr = std::exp(phi_log(wr));
      const double phi_w0 = std::exp(phi_log(w0));
      double swr = phi_wr / dP;   // phi(wr) / D
      double sw0 = phi_w0 / dP;   // phi(w0) / D
      if (!std::isfinite(swr)) swr = 1.0 / (w0 - wr) + wr;
      if (!std::isfinite(sw0)) sw0 = 0.0;
      const double G_w0 = G_by_Phi(w0, w0);   // phi(w0) / Phi(w0)

      // Half the first-derivative coefficients: a = 2 half_a, b = 2 half_b.
      const double half_a = swr;
      const double half_b = G_w0 - sw0;

      // Exact second-derivative coefficients of the dwr/dw0 quadratic forms.
      const double C_wrwr  = 2.0 * (-wr * swr + swr * swr);
      const double C_w0w0  = 2.0 * (-w0 * G_w0 - G_w0 * G_w0 + w0 * sw0 + sw0 * sw0);
      const double C_cross = -2.0 * swr * sw0;

      if (opts.bloq_part1) {
        for (int k2 = 0; k2 < n_par; ++k2) {
          const double dwr_k2 = dwr[k2];
          const double dw0_k2 = dw0[k2];
          for (int k1 = 0; k1 < n_par; ++k1) {
            hess_acc[k1 + k2 * n_par] +=
                C_wrwr  * dwr[k1] * dwr_k2
              + C_w0w0  * dw0[k1] * dw0_k2
              + C_cross * (dwr[k1] * dw0_k2 + dw0[k1] * dwr_k2);
          }
        }
      }

      // dpred/dsigma cross term from the a*d2wr + b*d2w0 chain rule.
      if (opts.bloq_part2 && has_dsig) {
        const double c = -2.0 * (half_a + half_b) * inv_s2;   // -(a + b) inv_s^2
        for (int k2 = 0; k2 < n_par; ++k2) {
          const double dxk2 = dpred_i[k2];
          const double dsk2 = dsigma_i[k2];
          for (int k1 = 0; k1 < n_par; ++k1) {
            hess_acc[k1 + k2 * n_par] +=
                c * (dpred_i[k1] * dsk2 + dsigma_i[k1] * dxk2);
          }
        }
      }

      // dsigma^2 term from the a*d2wr + b*d2w0 chain rule.
      if (opts.bloq_part3 && has_dsig) {
        const double c = 4.0 * inv_s2 * (half_a * wr + half_b * w0);  // 2(a wr + b w0) inv_s^2
        for (int k2 = 0; k2 < n_par; ++k2) {
          const double dsk2 = dsigma_i[k2];
          for (int k1 = 0; k1 < n_par; ++k1) {
            hess_acc[k1 + k2 * n_par] += c * dsigma_i[k1] * dsk2;
          }
        }
      }
    }

    if (opts.use_deriv2_exact) W_d2[i] = w_deriv2;
    if (has_d2sig) W_d2sig[i] = w_deriv2_sig;
  }  // for each BLOQ row

  if (opts.use_deriv2_exact && d2pred != nullptr) {
    contract_d2_block(n_obs, n_par, d2pred, W_d2.data(), hess_acc, scratch_d2);
  }
  if (has_d2sig) {
    contract_d2_block(n_obs, n_par, d2sigma, W_d2sig.data(), hess_acc,
                      scratch_d2);
  }
}

}  // namespace dmod
