// Analytical kernel for the LKJ + half-Normal / half-Cauchy prior on a
// Cholesky-parametrised Omega (dMod's omegaSpec convention: L_kk =
// exp(omega_kk), L_kl = omega_kl for k > l). Returns -2 log p and its
// gradient + Hessian over the Cholesky parameters.
//
// The log-prior factorises by row of L, so the Hessian over the Cholesky
// parameters is block-diagonal by row. Per row k of size n_k we compute:
//   value contribution:
//     (eta_k + 1) * omega_kk - 0.5 * eta_k * log(s_k) + scaleTerm(s_k)
//   gradient (per row entry m at column l_m):
//     diag (l_m == k):
//       (eta_k + 1) - eta_k * L_corr_kk^2 + dScale/d(omega_kk)
//     off-diag (l_m != k):
//       - eta_k * L_kl / s_k + dScale/d(omega_kl)
//   Hessian block per row: see formulas in priorOmega.R header.
//
// eta_k := K - k + 2 * lkjEta - 2 (standard LKJ exponent with the K-k
// term arising from the L_corr Jacobian; see Stan's
// lkj_corr_cholesky_lpdf).
//
// `kindFlag` selects the marginal-scale density:
//   0 = LKJHalfNormal:  scaleTerm(s) = log(2) - 0.5 log(2 pi tau2)
//                                       - s / (2 tau2)
//   1 = LKJHalfCauchy:  scaleTerm(s) = log(2) - log(pi) - log(tau)
//                                       - log(1 + s/tau2)

#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;


// Compute -2 log p, gradient (length P) and Hessian (P x P) over the
// Cholesky parameters. `cholLoc` is the (P x 2) 1-based (row, col) matrix
// from omegaSpec$cholLoc. `omegaVec` is the named numeric vector indexed
// in the same order. `isDiag` is a logical P-vector with TRUE for the
// log-diagonal entries.
//
// [[Rcpp::export(name = "priorOmegaKernel")]]
List prior_omega_kernel(NumericVector omegaVec,
                        IntegerMatrix cholLoc,
                        LogicalVector isDiag,
                        int K,
                        double lkjEta,
                        double scaleSD,
                        int kindFlag) {

  const int P = omegaVec.size();
  if (cholLoc.nrow() != P || isDiag.size() != P)
    stop("priorOmegaKernel: cholLoc and isDiag must have length P.");

  const double tau2 = scaleSD * scaleSD;

  // --- Build L (lower-triangular K x K) and the row index sets ---------
  std::vector<double> L(K * K, 0.0);
  std::vector<std::vector<int>> rowIdx(K);   // 0-based param indices per row
  for (int m = 0; m < P; ++m) {
    const int k = cholLoc(m, 0) - 1;          // 0-based
    const int l = cholLoc(m, 1) - 1;
    const double v = omegaVec[m];
    L[k + K * l] = isDiag[m] ? std::exp(v) : v;
    rowIdx[k].push_back(m);
  }

  // sigma_k^2 = sum_l L_kl^2 (l = 0..k)
  std::vector<double> s_vec(K, 0.0);
  for (int k = 0; k < K; ++k) {
    double s = 0.0;
    for (int l = 0; l <= k; ++l) {
      const double Lk = L[k + K * l];
      s += Lk * Lk;
    }
    s_vec[k] = s;
  }

  // L_corr_kl = L_kl / sigma_k
  std::vector<double> L_corr(K * K, 0.0);
  for (int k = 0; k < K; ++k) {
    const double sigma_k = std::sqrt(std::max(s_vec[k], std::numeric_limits<double>::min()));
    for (int l = 0; l <= k; ++l)
      L_corr[k + K * l] = L[k + K * l] / sigma_k;
  }

  // --- Value -----------------------------------------------------------
  double log_p = 0.0;
  for (int k = 0; k < K; ++k) {
    const double eta_k = (double)(K - 1 - k) + 2.0 * lkjEta - 2.0;
    const double Lkk   = L[k + K * k];
    const double sk    = s_vec[k];

    const double log_lkj_jac = (eta_k + 1.0) * std::log(Lkk)
                              - 0.5 * eta_k * std::log(sk);

    double log_scale;
    if (kindFlag == 0) { // half-Normal
      log_scale = std::log(2.0) - 0.5 * std::log(2.0 * M_PI * tau2)
                  - 0.5 * sk / tau2;
    } else { // half-Cauchy
      log_scale = std::log(2.0) - std::log(M_PI) - std::log(scaleSD)
                  - std::log1p(sk / tau2);
    }
    log_p += log_lkj_jac + log_scale;
  }
  const double value = -2.0 * log_p;

  // --- Gradient + Hessian (block-diagonal by row of L) ----------------
  NumericVector gradient(P);
  NumericMatrix H(P, P);

  for (int k = 0; k < K; ++k) {
    const std::vector<int>& idx_k = rowIdx[k];
    const int n_k = (int)idx_k.size();
    if (n_k == 0) continue;

    const double eta_k    = (double)(K - 1 - k) + 2.0 * lkjEta - 2.0;
    const double sk       = s_vec[k];
    const double Lkk      = L[k + K * k];
    const double Lcorr_kk = L_corr[k + K * k];
    const double sigma_k  = std::sqrt(std::max(sk, std::numeric_limits<double>::min()));
    const double q        = tau2 + sk;

    // Per-row entry data
    std::vector<int>    col_a(n_k);
    std::vector<bool>   is_diag_a(n_k);
    for (int a = 0; a < n_k; ++a) {
      const int m = idx_k[a];
      col_a[a]    = cholLoc(m, 1) - 1;
      is_diag_a[a] = (col_a[a] == k);
    }

    // --- gradient ----------------------------------------------------
    for (int a = 0; a < n_k; ++a) {
      const int m = idx_k[a];
      const int l_a = col_a[a];
      double g_lkj, g_scl;

      if (is_diag_a[a]) {
        g_lkj = (eta_k + 1.0) - eta_k * Lcorr_kk * Lcorr_kk;
        if (kindFlag == 0) g_scl = -Lkk * Lkk / tau2;
        else               g_scl = -2.0 * Lkk * Lkk / q;
      } else {
        const double Lkla = L[k + K * l_a];
        g_lkj = -eta_k * Lkla / sk;
        if (kindFlag == 0) g_scl = -Lkla / tau2;
        else               g_scl = -2.0 * Lkla / q;
      }
      // -2 log p convention
      gradient[m] = -2.0 * (g_lkj + g_scl);
    }

    // --- Hessian block (symmetric) -----------------------------------
    for (int a = 0; a < n_k; ++a) {
      for (int b = 0; b <= a; ++b) {
        const int m_a = idx_k[a]; const int m_b = idx_k[b];
        const int l_a = col_a[a]; const int l_b = col_a[b];

        // LKJ Hessian
        double h_lkj;
        if (is_diag_a[a] && is_diag_a[b]) {
          h_lkj = -2.0 * eta_k * Lcorr_kk * Lcorr_kk
                       * (1.0 - Lcorr_kk * Lcorr_kk);
        } else if (is_diag_a[a] && !is_diag_a[b]) {
          h_lkj = 2.0 * eta_k * Lcorr_kk * Lcorr_kk
                      * L_corr[k + K * l_b] / sigma_k;
        } else if (!is_diag_a[a] && is_diag_a[b]) {
          h_lkj = 2.0 * eta_k * Lcorr_kk * Lcorr_kk
                      * L_corr[k + K * l_a] / sigma_k;
        } else if (l_a == l_b) {
          const double Lcorr_kla = L_corr[k + K * l_a];
          h_lkj = -eta_k * (1.0 - 2.0 * Lcorr_kla * Lcorr_kla) / sk;
        } else {
          h_lkj = 2.0 * eta_k * L_corr[k + K * l_a] * L_corr[k + K * l_b] / sk;
        }

        // Scale Hessian
        double h_scl;
        if (kindFlag == 0) {
          if (is_diag_a[a] && is_diag_a[b])         h_scl = -2.0 * Lkk * Lkk / tau2;
          else if (a == b)                          h_scl = -1.0 / tau2;
          else                                      h_scl = 0.0;
        } else {
          const double ds_a = is_diag_a[a] ? (2.0 * Lkk * Lkk)
                                            : (2.0 * L[k + K * l_a]);
          const double ds_b = is_diag_a[b] ? (2.0 * Lkk * Lkk)
                                            : (2.0 * L[k + K * l_b]);
          double d2s = 0.0;
          if (is_diag_a[a] && is_diag_a[b])      d2s = 4.0 * Lkk * Lkk;
          else if (a == b)                       d2s = 2.0;
          h_scl = ds_a * ds_b / (q * q) - d2s / q;
        }

        const double val = -2.0 * (h_lkj + h_scl);
        H(m_a, m_b) = val;
        if (m_a != m_b) H(m_b, m_a) = val;
      }
    }
  }

  return List::create(_["value"]    = value,
                      _["gradient"] = gradient,
                      _["hessian"]  = H);
}
