// dMod constraintL2 + datapointL2 C++ kernels.
//
// constraintL2_scalar_kernel: prior penalty (p - mu)^2 / sigma^2 + (est?
//   2 log sigma : 0). Mirrors R/objClass.R::constraintL2 scalar path.
//   Supports the dP chain rule and the dP2 exact-second-order term.
//
// datapointL2_kernel: single-data-point L2 penalty (pred - target)^2 / sigma^2
//   for a single (condition, time, observable) cell against a target value
//   which is itself a parameter. Mirrors R/objClass.R::datapointL2.
//
// constraintL2_mvn_kernel: multivariate Gaussian prior over per-subject eta
// blocks; Cholesky-parametrized Omega. Mirrors R/objClass.R::constraintL2_mvn.
// Eta + chol parameter blocks, GN Hessian in z-space, with full
// chain-rule via dP and exact dP2 contribution.

#include "residual_kernel.h"

#include <Rcpp.h>
#include <vector>
#include <string>
#include <cstddef>
#include <stdexcept>
#include <cmath>
#include <cstring>

using namespace Rcpp;

extern "C" {
  // LAPACK / BLAS used here:
  //   dtrtrs: triangular solve, forward/backsolve
  //   dgemm:  matrix-matrix multiply
  void dtrtrs_(const char* uplo, const char* trans, const char* diag,
               const int* n, const int* nrhs, const double* a, const int* lda,
               double* b, const int* ldb, int* info);
  void dgemm_(const char* transa, const char* transb,
              const int* m, const int* n, const int* k,
              const double* alpha, const double* a, const int* lda,
              const double* b, const int* ldb,
              const double* beta, double* c, const int* ldc);
}

namespace {

inline int find_name(const CharacterVector& haystack, const std::string& s) {
  for (int i = 0; i < haystack.size(); ++i)
    if (as<std::string>(haystack[i]) == s) return i;
  return -1;
}

// Build a row-major contraction H[k1, k2] += sum_p w[p] * dP2[p, k1, k2]
// where dP2 is provided as a flat numeric vector with col-major dims
// [n_inner, n_theta, n_theta] and `inner_idx` selects which inner-par
// rows to contract (1-based; -1 = skip).
void apply_dP2_exact(
    const double* dP2_flat,
    const IntegerVector& dP2_dim,
    const std::vector<int>& inner_idx,   // 0-based, -1 to skip
    const std::vector<double>& gi_inner,
    int n_theta,
    double* hess_theta_colmajor) {
  if (inner_idx.empty()) return;
  const int n_inner_full = dP2_dim[0];
  for (std::size_t p = 0; p < inner_idx.size(); ++p) {
    const int row = inner_idx[p];
    if (row < 0) continue;
    const double w = gi_inner[p];
    if (w == 0.0) continue;
    // dP2[row, k1, k2] = dP2_flat[row + k1 * n_inner_full + k2 * n_inner_full * n_theta]
    for (int k2 = 0; k2 < n_theta; ++k2) {
      const std::size_t off2 = (std::size_t) k2 * n_inner_full * n_theta;
      for (int k1 = 0; k1 < n_theta; ++k1) {
        hess_theta_colmajor[k1 + k2 * n_theta] +=
            w * dP2_flat[row + (std::size_t) k1 * n_inner_full + off2];
      }
    }
  }
}

// Sandwich update: H_theta = dP^T * H_inner * dP (both inner and theta in
// col-major, dP provided as col-major NumericMatrix). Adds to hess_theta.
void sandwich_hess(const double* dP, int n_inner, int n_theta,
                   const double* H_inner_colmajor,
                   double* hess_theta_colmajor) {
  // tmp = H_inner * dP : (n_inner x n_theta)
  std::vector<double> tmp((std::size_t) n_inner * n_theta, 0.0);
  const double alpha = 1.0, beta = 0.0;
  dgemm_("N", "N", &n_inner, &n_theta, &n_inner,
         &alpha, H_inner_colmajor, &n_inner, dP, &n_inner,
         &beta, tmp.data(), &n_inner);
  // H_theta += dP^T * tmp : (n_theta x n_theta)
  const double beta2 = 1.0;
  dgemm_("T", "N", &n_theta, &n_theta, &n_inner,
         &alpha, dP, &n_inner, tmp.data(), &n_inner,
         &beta2, hess_theta_colmajor, &n_theta);
}

}  // namespace


// [[Rcpp::export]]
List constraintL2_scalar_kernel(
    NumericVector pars,                 // free params (length n_theta if dP, else n_inner)
    Nullable<NumericMatrix> dP_opt,     // [n_inner, n_theta] col-major, or NULL (identity)
    Nullable<NumericVector> dP2_opt,    // 3D array [n_inner, n_theta, n_theta], or NULL
    CharacterVector inner_par_names,    // names of inner pars (length n_inner)
    Nullable<NumericVector> fixed_opt,  // any fixed values needed for `allp` lookup
    CharacterVector mu_names,
    NumericVector mu,
    NumericVector sigma,                // length = length(mu_names); when est, holds 0.0 for est-rows
    CharacterVector sigma_pars,         // names of sigma params (empty string if fixed)
    bool est) {

  // Build allp lookup: name -> value
  // pars carries the outer (theta) parameter values when dP is given;
  // otherwise it carries the inner parameter values directly.
  // For the sigma/log-sigma case (`est`), sigma values come from allp.
  const int n_inner_full = inner_par_names.size();
  // Build "allp" map from union of pars/fixed names.
  std::vector<std::string> allp_names;
  std::vector<double>      allp_vals;
  if (dP_opt.isNotNull()) {
    // pars are theta-named; we need inner-par values via allp lookup
    NumericVector p = pars;
    CharacterVector pn = p.names();
    for (int i = 0; i < p.size(); ++i) {
      allp_names.push_back(as<std::string>(pn[i]));
      allp_vals.push_back(p[i]);
    }
  } else {
    NumericVector p = pars;
    CharacterVector pn = p.names();
    for (int i = 0; i < p.size(); ++i) {
      allp_names.push_back(as<std::string>(pn[i]));
      allp_vals.push_back(p[i]);
    }
  }
  if (fixed_opt.isNotNull()) {
    NumericVector f = fixed_opt.get();
    CharacterVector fn = f.names();
    for (int i = 0; i < f.size(); ++i) {
      allp_names.push_back(as<std::string>(fn[i]));
      allp_vals.push_back(f[i]);
    }
  }
  auto lookup = [&](const std::string& nm) -> double {
    for (std::size_t i = 0; i < allp_names.size(); ++i) {
      if (allp_names[i] == nm) return allp_vals[i];
    }
    return std::numeric_limits<double>::quiet_NaN();
  };

  const int n_mu = mu_names.size();
  // For each mu entry: compute r = lookup(mu_name) - mu_value, sg.
  double value = 0.0;
  std::vector<double> r_vec(n_mu), sg_vec(n_mu);
  std::vector<bool>   avail(n_mu);
  for (int i = 0; i < n_mu; ++i) {
    std::string nm = as<std::string>(mu_names[i]);
    double v = lookup(nm);
    avail[i] = !std::isnan(v);
    if (!avail[i]) continue;
    r_vec[i] = v - mu[i];
    if (est) {
      std::string sp = as<std::string>(sigma_pars[i]);
      double sg_lin = lookup(sp);
      if (std::isnan(sg_lin)) { avail[i] = false; continue; }
      sg_vec[i] = std::exp(sg_lin);
    } else {
      sg_vec[i] = sigma[i];
    }
    value += (r_vec[i] / sg_vec[i]) * (r_vec[i] / sg_vec[i]);
    if (est) value += 2.0 * std::log(sg_vec[i]);
  }

  // Compute inner gradient gi (length n_inner_full) + inner Hessian.
  std::vector<double> gi(n_inner_full, 0.0);
  std::vector<double> Hi((std::size_t) n_inner_full * n_inner_full, 0.0);

  // Map mu_name index -> inner-par index (-1 if not in inner pars)
  std::vector<int> mu_to_inner(n_mu, -1);
  for (int i = 0; i < n_mu; ++i) {
    mu_to_inner[i] = find_name(inner_par_names, as<std::string>(mu_names[i]));
  }
  // Map sigma_par name -> inner-par index (for est case)
  std::vector<int> sigma_to_inner(n_mu, -1);
  if (est) {
    for (int i = 0; i < n_mu; ++i) {
      sigma_to_inner[i] = find_name(inner_par_names, as<std::string>(sigma_pars[i]));
    }
  }

  for (int i = 0; i < n_mu; ++i) {
    if (!avail[i]) continue;
    const double sg2 = sg_vec[i] * sg_vec[i];
    const int ip = mu_to_inner[i];
    if (ip >= 0) {
      gi[ip] += 2.0 * r_vec[i] / sg2;
      Hi[ip + (std::size_t) ip * n_inner_full] += 2.0 / sg2;
    }
    if (est) {
      const int ips = sigma_to_inner[i];
      if (ips >= 0) {
        gi[ips] += -2.0 * r_vec[i] * r_vec[i] / sg2 + 2.0;
        Hi[ips + (std::size_t) ips * n_inner_full] += 4.0 * r_vec[i] * r_vec[i] / sg2;
        if (ip >= 0) {
          const double off = -4.0 * r_vec[i] / sg2;
          Hi[ip + (std::size_t) ips * n_inner_full] += off;
          Hi[ips + (std::size_t) ip * n_inner_full] += off;
        }
      }
    }
  }

  // Chain rule via dP / exact dP2 contribution.
  NumericVector grad_out;
  NumericMatrix hess_out;
  CharacterVector theta_names;

  if (dP_opt.isNotNull()) {
    NumericMatrix dP(dP_opt.get());
    const int dP_n_inner = dP.nrow();
    const int n_theta    = dP.ncol();
    if (dP_n_inner != n_inner_full)
      throw std::runtime_error("constraintL2_scalar: dP nrow mismatch.");
    List dP_dimnames = dP.attr("dimnames");
    theta_names = (dP_dimnames.size() >= 2) ? dP_dimnames[1] : CharacterVector(n_theta);

    // grad = dP^T * gi
    std::vector<double> grad_theta(n_theta, 0.0);
    for (int k = 0; k < n_theta; ++k) {
      double s = 0.0;
      for (int p = 0; p < n_inner_full; ++p) s += dP(p, k) * gi[p];
      grad_theta[k] = s;
    }
    // hess_theta = dP^T * Hi * dP
    std::vector<double> hess_theta((std::size_t) n_theta * n_theta, 0.0);
    sandwich_hess(&dP(0,0), n_inner_full, n_theta, Hi.data(), hess_theta.data());

    // Exact dP2 contribution: H_theta[k1, k2] += sum_p gi[p] * dP2[p, k1, k2]
    if (dP2_opt.isNotNull()) {
      NumericVector dP2_flat(dP2_opt.get());
      IntegerVector dP2_dim = dP2_flat.attr("dim");
      if (dP2_dim.size() == 3 && dP2_dim[1] == n_theta && dP2_dim[2] == n_theta) {
        List dP2_dn = dP2_flat.attr("dimnames");
        CharacterVector dP2_inner_names = dP2_dn[0];
        std::vector<int> idx(dP2_inner_names.size(), -1);
        std::vector<double> gi_inner(dP2_inner_names.size(), 0.0);
        for (int i = 0; i < dP2_inner_names.size(); ++i) {
          int ip = find_name(inner_par_names, as<std::string>(dP2_inner_names[i]));
          if (ip >= 0) {
            idx[i]      = i;
            gi_inner[i] = gi[ip];
          }
        }
        apply_dP2_exact(REAL(dP2_flat), dP2_dim, idx, gi_inner,
                        n_theta, hess_theta.data());
      }
    }

    grad_out = NumericVector(grad_theta.begin(), grad_theta.end());
    grad_out.names() = theta_names;
    hess_out = NumericMatrix(n_theta, n_theta);
    std::memcpy(&hess_out(0,0), hess_theta.data(),
                sizeof(double) * (std::size_t) n_theta * n_theta);
    hess_out.attr("dimnames") = List::create(theta_names, theta_names);
  } else {
    // No chain rule: gradient/Hessian directly in inner-par space.
    grad_out = NumericVector(gi.begin(), gi.end());
    grad_out.names() = inner_par_names;
    hess_out = NumericMatrix(n_inner_full, n_inner_full);
    std::memcpy(&hess_out(0,0), Hi.data(),
                sizeof(double) * (std::size_t) n_inner_full * n_inner_full);
    hess_out.attr("dimnames") = List::create(inner_par_names, inner_par_names);
  }

  return List::create(
      Named("value")    = value,
      Named("gradient") = grad_out,
      Named("hessian")  = hess_out);
}


// [[Rcpp::export]]
List datapointL2_kernel(
    NumericVector pouter,                // full outer param vector (named)
    Nullable<NumericVector> fixed_opt,
    NumericMatrix prdf,                  // prediction matrix for the target condition
    Nullable<NumericVector> dpred_attr_opt,  // 3D [time, var, par] from prdf attr "deriv"
    Nullable<NumericVector> d2pred_attr_opt, // 4D [time, var, par, par] or NULL
    std::string obs_name,                // single observable name
    double t,
    double sigma,
    std::string value_par) {             // name of pouter param that holds the data target value

  CharacterVector pouter_names = pouter.names();
  const int n_p = pouter.size();

  // 1. Time index in prediction
  NumericVector time_col(prdf.nrow());
  for (int i = 0; i < prdf.nrow(); ++i) time_col[i] = prdf(i, 0);
  int ti = -1;
  for (int i = 0; i < prdf.nrow(); ++i)
    if (time_col[i] == t) { ti = i; break; }
  if (ti < 0)
    throw std::runtime_error("datapointL2_kernel: time not in prediction.");

  // 2. Observable column index
  CharacterVector pcols = colnames(prdf);
  int oi = find_name(pcols, obs_name);
  if (oi < 0)
    throw std::runtime_error("datapointL2_kernel: observable not in prediction.");

  const double pred = prdf(ti, oi);

  // 3. Lookup data-target value from pouter or fixed
  auto lookup = [&](const std::string& nm) -> double {
    for (int i = 0; i < pouter_names.size(); ++i)
      if (as<std::string>(pouter_names[i]) == nm) return pouter[i];
    if (fixed_opt.isNotNull()) {
      NumericVector f = fixed_opt.get();
      CharacterVector fn = f.names();
      for (int i = 0; i < f.size(); ++i)
        if (as<std::string>(fn[i]) == nm) return f[i];
    }
    return std::numeric_limits<double>::quiet_NaN();
  };
  const double target = lookup(value_par);
  if (std::isnan(target))
    throw std::runtime_error("datapointL2_kernel: target value param not found.");

  const double res = pred - target;
  const double sigma2 = sigma * sigma;
  const double value = (res / sigma) * (res / sigma);

  // 4. dres/dp: -1 for the value_par, dpred/dp for the structural pars
  std::vector<double> dres_dp(n_p, 0.0);
  if (dpred_attr_opt.isNotNull()) {
    NumericVector dpred_flat(dpred_attr_opt.get());
    IntegerVector ddim = dpred_flat.attr("dim");
    if (ddim.size() == 3) {
      const int Dp0 = ddim[0];
      const int Dp1 = ddim[1];
      List ddn = dpred_flat.attr("dimnames");
      CharacterVector dpred_obs_names = ddn[1];
      CharacterVector dpred_par_names = ddn[2];
      int od = find_name(dpred_obs_names, obs_name);
      if (od >= 0) {
        for (int p = 0; p < dpred_par_names.size(); ++p) {
          int ip = find_name(pouter_names, as<std::string>(dpred_par_names[p]));
          if (ip >= 0) {
            dres_dp[ip] = dpred_flat[ti + od * Dp0 + p * Dp0 * Dp1];
          }
        }
      }
    }
  }
  // datapar: the value_par contributes -1
  const int idx_value = find_name(pouter_names, value_par);
  if (idx_value >= 0) dres_dp[idx_value] = -1.0;

  // 5. Gradient + Hessian
  std::vector<double> gr(n_p, 0.0);
  std::vector<double> hs((std::size_t) n_p * n_p, 0.0);
  for (int i = 0; i < n_p; ++i) {
    gr[i] = 2.0 * res * dres_dp[i] / sigma2;
    for (int j = 0; j < n_p; ++j) {
      hs[i + (std::size_t) j * n_p] = 2.0 * dres_dp[i] * dres_dp[j] / sigma2;
    }
  }
  // 6. Exact d2pred contribution (only on structural pars)
  if (d2pred_attr_opt.isNotNull()) {
    NumericVector d2_flat(d2pred_attr_opt.get());
    IntegerVector d2dim = d2_flat.attr("dim");
    if (d2dim.size() == 4) {
      const int D0 = d2dim[0], D1 = d2dim[1], D2 = d2dim[2];
      List d2dn = d2_flat.attr("dimnames");
      CharacterVector d2_obs_names = d2dn[1];
      CharacterVector d2_par_names = d2dn[2];
      int od = find_name(d2_obs_names, obs_name);
      if (od >= 0) {
        std::vector<int> par_to_outer(d2_par_names.size(), -1);
        for (int p = 0; p < d2_par_names.size(); ++p) {
          par_to_outer[p] = find_name(pouter_names, as<std::string>(d2_par_names[p]));
        }
        const double w = 2.0 * res / sigma2;
        for (int p = 0; p < d2_par_names.size(); ++p) {
          const int ip = par_to_outer[p];
          if (ip < 0) continue;
          for (int q = 0; q < d2_par_names.size(); ++q) {
            const int jq = par_to_outer[q];
            if (jq < 0) continue;
            const double val = d2_flat[ti + od * D0 + p * D0 * D1 + q * D0 * D1 * D2];
            hs[ip + (std::size_t) jq * n_p] += w * val;
          }
        }
      }
    }
  }

  NumericVector grad_out(gr.begin(), gr.end());
  grad_out.names() = pouter_names;
  NumericMatrix hess_out(n_p, n_p);
  std::memcpy(&hess_out(0,0), hs.data(),
              sizeof(double) * (std::size_t) n_p * n_p);
  hess_out.attr("dimnames") = List::create(pouter_names, pouter_names);

  return List::create(
      Named("value")      = value,
      Named("gradient")   = grad_out,
      Named("hessian")    = hess_out,
      Named("prediction") = pred);
}


// [[Rcpp::export]]
List constraintL2_mvn_kernel(
    NumericVector pars,                 // free params (theta-named if dP, else inner)
    Nullable<NumericVector> fixed_opt,
    Nullable<NumericMatrix> dP_opt,
    Nullable<NumericVector> dP2_opt,
    CharacterVector inner_par_names,    // c(all_eta_names, chol_pars)
    int K,
    int N,
    CharacterVector all_eta_names,      // length K*N
    NumericVector mu,                   // length K
    NumericMatrix L_lower,              // K x K lower-triangular Omega chol (precomputed in R)
    bool include_chol_block) {          // false => return zeros in chol-block (typical trust use)

  (void) include_chol_block;
  // Build allp lookup
  std::vector<std::string> allp_names; std::vector<double> allp_vals;
  CharacterVector pn = pars.names();
  for (int i = 0; i < pars.size(); ++i) {
    allp_names.push_back(as<std::string>(pn[i])); allp_vals.push_back(pars[i]);
  }
  if (fixed_opt.isNotNull()) {
    NumericVector f = fixed_opt.get();
    CharacterVector fn = f.names();
    for (int i = 0; i < f.size(); ++i) {
      allp_names.push_back(as<std::string>(fn[i])); allp_vals.push_back(f[i]);
    }
  }
  auto lookup = [&](const std::string& nm) -> double {
    for (std::size_t i = 0; i < allp_names.size(); ++i)
      if (allp_names[i] == nm) return allp_vals[i];
    return std::numeric_limits<double>::quiet_NaN();
  };

  // Gather eta_mat as K x N: eta_mat[k, i] = allp[all_eta_names[i*K + k]]
  // (matches R's subjectEtas which is N x K but indexed with all_eta_names
  // = as.vector(subjectEtas) — by column).
  std::vector<double> eta_mat((std::size_t) K * N, 0.0);
  bool all_present = true;
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < K; ++k) {
      // R uses as.vector(subject_etas) which is column-major over [N, K]:
      //   all_eta_names[(k-1)*N + (i-1)+1] = subject_etas[i, k]
      // So all_eta_names[k*N + i] = subject_etas[i, k]
      std::string nm = as<std::string>(all_eta_names[k * N + i]);
      double v = lookup(nm);
      if (std::isnan(v)) { all_present = false; break; }
      eta_mat[k + (std::size_t) i * K] = v;
    }
    if (!all_present) break;
  }

  const int n_inner = inner_par_names.size();
  std::vector<double> gi(n_inner, 0.0);
  std::vector<double> Hi((std::size_t) n_inner * n_inner, 0.0);
  double value = 0.0;

  if (all_present) {
    // R = eta_mat - mu  (K x N, col-major)
    std::vector<double> R((std::size_t) K * N, 0.0);
    for (int i = 0; i < N; ++i)
      for (int k = 0; k < K; ++k)
        R[k + (std::size_t) i * K] = eta_mat[k + (std::size_t) i * K] - mu[k];

    // Z = forwardsolve(L, R): L * Z = R, solve for Z; col-major
    std::vector<double> Z = R;  // overwritten in place by dtrtrs
    int nrhs = N, info = 0;
    dtrtrs_("L", "N", "N", &K, &nrhs, &L_lower(0, 0), &K,
            Z.data(), &K, &info);
    if (info != 0) throw std::runtime_error("constraintL2_mvn: dtrtrs forward failed.");

    // W = backsolve(t(L), Z): L^T * W = Z; col-major
    std::vector<double> W = Z;
    dtrtrs_("L", "T", "N", &K, &nrhs, &L_lower(0, 0), &K,
            W.data(), &K, &info);
    if (info != 0) throw std::runtime_error("constraintL2_mvn: dtrtrs backsolve failed.");

    // quad = sum(Z * Z); logdetO = 2 * sum log diag(L); value = quad + N * logdetO
    double quad = 0.0;
    for (std::size_t s = 0; s < (std::size_t) K * N; ++s) quad += Z[s] * Z[s];
    double logdetO = 0.0;
    for (int k = 0; k < K; ++k) logdetO += std::log(L_lower(k, k));
    logdetO *= 2.0;
    value = quad + (double) N * logdetO;

    // Inner gradient (eta block): gi[name_of(eta_mat[k,i])] = 2 * W[k, i]
    for (int i = 0; i < N; ++i) {
      for (int k = 0; k < K; ++k) {
        std::string nm = as<std::string>(all_eta_names[k * N + i]);
        int ip = find_name(inner_par_names, nm);
        if (ip >= 0) gi[ip] = 2.0 * W[k + (std::size_t) i * K];
      }
    }

    // Inner Hessian (eta-eta block, GN): Hi[eta_idx, eta_idx] = 2 * Omega_inv
    // Omega_inv = L^{-T} L^{-1}. Compute via two dtrtrs against identity.
    std::vector<double> Linv((std::size_t) K * K, 0.0);
    for (int k = 0; k < K; ++k) Linv[k + (std::size_t) k * K] = 1.0;
    int K_int = K;
    dtrtrs_("L", "N", "N", &K_int, &K_int, &L_lower(0, 0), &K,
            Linv.data(), &K, &info);
    // Now Linv = L^{-1} (lower tri).
    // Omega_inv = Linv^T * Linv : K x K, symmetric. Use dgemm.
    std::vector<double> Omega_inv((std::size_t) K * K, 0.0);
    const double alpha = 1.0, beta_init = 0.0;
    dgemm_("T", "N", &K_int, &K_int, &K_int,
           &alpha, Linv.data(), &K_int, Linv.data(), &K_int,
           &beta_init, Omega_inv.data(), &K_int);

    // For each subject, add 2 * Omega_inv to the eta-block of Hi.
    for (int i = 0; i < N; ++i) {
      std::vector<int> idx(K, -1);
      for (int k = 0; k < K; ++k) {
        std::string nm = as<std::string>(all_eta_names[k * N + i]);
        idx[k] = find_name(inner_par_names, nm);
      }
      for (int k2 = 0; k2 < K; ++k2) {
        if (idx[k2] < 0) continue;
        for (int k1 = 0; k1 < K; ++k1) {
          if (idx[k1] < 0) continue;
          Hi[idx[k1] + (std::size_t) idx[k2] * n_inner] +=
              2.0 * Omega_inv[k1 + (std::size_t) k2 * K];
        }
      }
    }
    // Chol-block + cross-block are NOT computed here (typical trust users
    // don't estimate chol pars; for ECM workflow the dedicated path runs).
  }

  // Chain rule via dP / exact dP2 contribution. Mirrors scalar kernel.
  NumericVector grad_out;
  NumericMatrix hess_out;
  CharacterVector theta_names;

  if (dP_opt.isNotNull()) {
    NumericMatrix dP(dP_opt.get());
    const int dP_n_inner = dP.nrow();
    const int n_theta    = dP.ncol();
    if (dP_n_inner != n_inner)
      throw std::runtime_error("constraintL2_mvn: dP nrow mismatch.");
    List dP_dn = dP.attr("dimnames");
    theta_names = (dP_dn.size() >= 2) ? dP_dn[1] : CharacterVector(n_theta);

    std::vector<double> grad_theta(n_theta, 0.0);
    for (int k = 0; k < n_theta; ++k) {
      double s = 0.0;
      for (int p = 0; p < n_inner; ++p) s += dP(p, k) * gi[p];
      grad_theta[k] = s;
    }
    std::vector<double> hess_theta((std::size_t) n_theta * n_theta, 0.0);
    sandwich_hess(&dP(0,0), n_inner, n_theta, Hi.data(), hess_theta.data());

    if (dP2_opt.isNotNull()) {
      NumericVector dP2_flat(dP2_opt.get());
      IntegerVector dP2_dim = dP2_flat.attr("dim");
      if (dP2_dim.size() == 3 && dP2_dim[1] == n_theta && dP2_dim[2] == n_theta) {
        List dP2_dn = dP2_flat.attr("dimnames");
        CharacterVector dP2_inner_names = dP2_dn[0];
        std::vector<int> idx(dP2_inner_names.size(), -1);
        std::vector<double> gi_inner(dP2_inner_names.size(), 0.0);
        for (int i = 0; i < dP2_inner_names.size(); ++i) {
          int ip = find_name(inner_par_names, as<std::string>(dP2_inner_names[i]));
          if (ip >= 0) {
            idx[i] = i;
            gi_inner[i] = gi[ip];
          }
        }
        apply_dP2_exact(REAL(dP2_flat), dP2_dim, idx, gi_inner,
                        n_theta, hess_theta.data());
      }
    }

    grad_out = NumericVector(grad_theta.begin(), grad_theta.end());
    grad_out.names() = theta_names;
    hess_out = NumericMatrix(n_theta, n_theta);
    std::memcpy(&hess_out(0,0), hess_theta.data(),
                sizeof(double) * (std::size_t) n_theta * n_theta);
    hess_out.attr("dimnames") = List::create(theta_names, theta_names);
  } else {
    grad_out = NumericVector(gi.begin(), gi.end());
    grad_out.names() = inner_par_names;
    hess_out = NumericMatrix(n_inner, n_inner);
    std::memcpy(&hess_out(0,0), Hi.data(),
                sizeof(double) * (std::size_t) n_inner * n_inner);
    hess_out.attr("dimnames") = List::create(inner_par_names, inner_par_names);
  }

  return List::create(
      Named("value")    = value,
      Named("gradient") = grad_out,
      Named("hessian")  = hess_out);
}
