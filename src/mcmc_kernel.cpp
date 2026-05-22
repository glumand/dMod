// MCMC inner-loop primitives. See mcmc_kernel.h for sign / metric / threading.

#include "mcmc_kernel.h"

#include <Rmath.h>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <limits>

extern "C" {
  double norm_rand(void);
  // LAPACK
  void dpotrf_(const char* uplo, const int* n, double* a, const int* lda,
               int* info);
  void dpotrs_(const char* uplo, const int* n, const int* nrhs,
               const double* a, const int* lda,
               double* b, const int* ldb, int* info);
  void dpotri_(const char* uplo, const int* n, double* a, const int* lda,
               int* info);
  // BLAS
  void dtrsv_(const char* uplo, const char* trans, const char* diag,
              const int* n, const double* a, const int* lda,
              double* x, const int* incx);
  void dtrmv_(const char* uplo, const char* trans, const char* diag,
              const int* n, const double* a, const int* lda,
              double* x, const int* incx);
  void dsymv_(const char* uplo, const int* n,
              const double* alpha, const double* a, const int* lda,
              const double* x, const int* incx,
              const double* beta, double* y, const int* incy);
  void dgemm_(const char* transa, const char* transb,
              const int* m, const int* n, const int* k,
              const double* alpha, const double* a, const int* lda,
              const double* b, const int* ldb,
              const double* beta, double* c, const int* ldc);
}

namespace dmod {


void da_init(DAState& da, double eps0, double target,
             double gamma, double t0, double kappa) {
  da.eps     = eps0;
  da.eps_bar = 1.0;
  da.Hbar    = 0.0;
  da.mu      = std::log(10.0 * eps0);
  da.gamma   = gamma;
  da.t0      = t0;
  da.kappa   = kappa;
  da.target  = target;
}


void da_update(DAState& da, double alpha, int t) {
  if (!std::isfinite(alpha)) alpha = 0.0;
  const double tt = static_cast<double>(t);
  const double w  = 1.0 / (tt + da.t0);
  da.Hbar = (1.0 - w) * da.Hbar + w * (da.target - alpha);
  const double log_eps = da.mu - std::sqrt(tt) / da.gamma * da.Hbar;
  const double eta = std::pow(tt, -da.kappa);
  const double log_eps_bar =
      eta * log_eps + (1.0 - eta) * std::log(da.eps_bar);
  da.eps     = std::exp(log_eps);
  da.eps_bar = std::exp(log_eps_bar);
}


bool chol_with_ridge(const double* G, int K, double ridge0, int max_iter,
                     double* L_out, double* ridge_out) {
  double ridge = ridge0;
  const int n = K;
  for (int it = 0; it < max_iter; ++it) {
    // Copy G + ridge*I into L_out (column-major).
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        L_out[i + j * n] = G[i + j * n];
      }
      L_out[j + j * n] += ridge;
    }
    int info = 0;
    const char uplo = 'U';
    dpotrf_(&uplo, &n, L_out, &n, &info);
    if (info == 0) {
      // Zero the strictly lower triangle so downstream code can safely
      // touch the full matrix.
      for (int j = 0; j < n; ++j)
        for (int i = j + 1; i < n; ++i)
          L_out[i + j * n] = 0.0;
      if (ridge_out) *ridge_out = ridge;
      return true;
    }
    ridge = std::max(2.0 * ridge,
                     std::sqrt(std::numeric_limits<double>::epsilon()));
  }
  return false;
}


// Solves G x = b via dpotrs given upper Cholesky factor L. b is overwritten
// with x. K = n.
static void chol_solve(const double* L_upper, int K, double* b) {
  const int n = K;
  const int nrhs = 1;
  int info = 0;
  const char uplo = 'U';
  // dpotrs declares its matrix arg as const-correct in some headers and not
  // in others; pass a writable copy.
  std::vector<double> Lcopy(static_cast<size_t>(n) * static_cast<size_t>(n));
  std::memcpy(Lcopy.data(), L_upper,
              sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(n));
  dpotrs_(&uplo, &n, &nrhs, Lcopy.data(), &n, b, &n, &info);
}


// In-place upper-triangular solve U x = b.
static void utri_N_solve(const double* U_upper, int K, double* b) {
  const int n = K;
  const char uplo = 'U';
  const char trans = 'N';
  const char diag = 'N';
  const int incx = 1;
  dtrsv_(&uplo, &trans, &diag, &n, U_upper, &n, b, &incx);
}


void mala_drift(const double* grad, const double* L_upper, int K,
                const double* dG_row_major, double* drift_out) {
  // Stage A/B: drift = G^{-1} grad.
  std::vector<double> b(grad, grad + K);
  chol_solve(L_upper, K, b.data());
  for (int i = 0; i < K; ++i) drift_out[i] = b[i];

  if (dG_row_major == nullptr) return;

  // Stage C (Xifara 2014 corrected drift): drift += -2 * termB + termA.
  // Build G^{-1} as a dense K x K matrix once via dpotri on a copy.
  std::vector<double> Ginv(static_cast<size_t>(K) * static_cast<size_t>(K));
  std::memcpy(Ginv.data(), L_upper, sizeof(double) * Ginv.size());
  int info = 0;
  const char uplo = 'U';
  const int n = K;
  dpotri_(&uplo, &n, Ginv.data(), &n, &info);
  // Symmetrise: dpotri writes only the upper triangle.
  for (int j = 0; j < K; ++j)
    for (int i = j + 1; i < K; ++i)
      Ginv[i + j * K] = Ginv[j + i * K];

  // v_k = tr(G^{-1} dG_k); dG is row-major [K, K, K], indexed (a, b, k) ->
  // [a*K*K + b*K + k]. Reference R uses dG[, , k] which means slice by k
  // last dimension.
  std::vector<double> v(K, 0.0);
  for (int k = 0; k < K; ++k) {
    double tr = 0.0;
    for (int a = 0; a < K; ++a) {
      for (int j = 0; j < K; ++j) {
        const double dG_aj_k = dG_row_major[a * K * K + j * K + k];
        tr += Ginv[a + j * K] * dG_aj_k;
      }
    }
    v[k] = tr;
  }

  // termA = Ginv %*% v
  std::vector<double> termA(K, 0.0);
  for (int i = 0; i < K; ++i) {
    double s = 0.0;
    for (int j = 0; j < K; ++j) s += Ginv[i + j * K] * v[j];
    termA[i] = s;
  }

  // termB_k = sum_j (Ginv * dG_k * Ginv)[j, k]
  std::vector<double> M(static_cast<size_t>(K) * static_cast<size_t>(K));
  std::vector<double> tmp(static_cast<size_t>(K) * static_cast<size_t>(K));
  std::vector<double> dGk(static_cast<size_t>(K) * static_cast<size_t>(K));
  std::vector<double> termB(K, 0.0);

  const double one_ = 1.0, zero_ = 0.0;
  const char transN_ = 'N';
  const int nK = K;

  for (int k = 0; k < K; ++k) {
    for (int a = 0; a < K; ++a)
      for (int j = 0; j < K; ++j)
        dGk[a + j * K] = dG_row_major[a * K * K + j * K + k];

    // tmp = Ginv * dG_k
    dgemm_(&transN_, &transN_, &nK, &nK, &nK,
           &one_, Ginv.data(), &nK,
           dGk.data(), &nK,
           &zero_, tmp.data(), &nK);
    // M = tmp * Ginv
    dgemm_(&transN_, &transN_, &nK, &nK, &nK,
           &one_, tmp.data(), &nK,
           Ginv.data(), &nK,
           &zero_, M.data(), &nK);

    double s = 0.0;
    for (int j = 0; j < K; ++j) s += M[j + k * K];
    termB[k] = s;
  }

  for (int i = 0; i < K; ++i)
    drift_out[i] += -2.0 * termB[i] + termA[i];
}


void mala_propose(const double* theta_old, const double* drift_old,
                  const double* L_upper, int K, double eps,
                  double* z_buf, double* theta_new) {
  for (int i = 0; i < K; ++i) z_buf[i] = ::norm_rand();
  // noise = L^{-1} z   (Cov(noise) = G^{-1} when G = L^T L)
  utri_N_solve(L_upper, K, z_buf);
  const double sq_eps = std::sqrt(eps);
  for (int i = 0; i < K; ++i)
    theta_new[i] = theta_old[i] + 0.5 * eps * drift_old[i] + sq_eps * z_buf[i];
}


double mala_log_q(const double* theta_new, const double* theta_old,
                  const double* drift_old, const double* G_old,
                  int K, double eps, double logdet_old) {
  std::vector<double> d(K);
  for (int i = 0; i < K; ++i)
    d[i] = theta_new[i] - theta_old[i] - 0.5 * eps * drift_old[i];

  std::vector<double> Gd(K, 0.0);
  const double alpha_ = 1.0, beta_ = 0.0;
  const char uplo_ = 'U';
  const int incx_ = 1, incy_ = 1;
  const int n = K;
  dsymv_(&uplo_, &n, &alpha_, G_old, &n,
         d.data(), &incx_, &beta_, Gd.data(), &incy_);

  double q = 0.0;
  for (int i = 0; i < K; ++i) q += d[i] * Gd[i];

  return -0.5 / eps * q + 0.5 * logdet_old;
}


void mh_propose(const double* theta_old, const double* L_upper, int K,
                double eps, double* z_buf, double* theta_new) {
  for (int i = 0; i < K; ++i) z_buf[i] = ::norm_rand();
  utri_N_solve(L_upper, K, z_buf);
  const double sq_eps = std::sqrt(eps);
  for (int i = 0; i < K; ++i)
    theta_new[i] = theta_old[i] + sq_eps * z_buf[i];
}


double kinetic_energy(const double* p, const double* Minv_chol_upper, int K) {
  // 0.5 * p^T Minv p; given Minv = U^T U, p^T Minv p = || U p ||^2.
  std::vector<double> Up(p, p + K);
  const int n = K, incx = 1;
  const char uplo = 'U', trans = 'N', diag = 'N';
  dtrmv_(&uplo, &trans, &diag, &n, Minv_chol_upper, &n, Up.data(), &incx);
  double s = 0.0;
  for (int i = 0; i < K; ++i) s += Up[i] * Up[i];
  return 0.5 * s;
}


void sample_momentum(const double* LM_upper, int K, double* z_buf,
                     double* p_out) {
  for (int i = 0; i < K; ++i) z_buf[i] = ::norm_rand();
  // p = LM^T z  where LM is the upper Cholesky of M (M = LM^T LM).
  for (int i = 0; i < K; ++i) p_out[i] = z_buf[i];
  const int n = K, incx = 1;
  const char uplo = 'U', trans = 'T', diag = 'N';
  dtrmv_(&uplo, &trans, &diag, &n, LM_upper, &n, p_out, &incx);
}


bool leapfrog(double* theta, double* p, double* grad_lp,
              const double* Minv_chol_upper,
              int K, double eps, int n_steps,
              std::function<bool(const double*, double&, double*)> grad_cb) {
  // Half-step in p
  for (int i = 0; i < K; ++i) p[i] += 0.5 * eps * grad_lp[i];

  for (int s = 0; s < n_steps; ++s) {
    // Full step in theta: theta += eps * Minv * p
    // Minv * p = U^T (U * p) when Minv = U^T U
    std::vector<double> Mp(p, p + K);
    const int n = K, incx = 1;
    const char uplo = 'U', diag = 'N';
    const char transN_ = 'N', transT_ = 'T';
    dtrmv_(&uplo, &transN_, &diag, &n, Minv_chol_upper, &n, Mp.data(), &incx);
    dtrmv_(&uplo, &transT_, &diag, &n, Minv_chol_upper, &n, Mp.data(), &incx);
    for (int i = 0; i < K; ++i) theta[i] += eps * Mp[i];

    // Refresh gradient at the new theta
    double logp_dummy = 0.0;
    if (!grad_cb(theta, logp_dummy, grad_lp)) return false;

    if (s < n_steps - 1) {
      for (int i = 0; i < K; ++i) p[i] += eps * grad_lp[i];
    } else {
      // Final half-step
      for (int i = 0; i < K; ++i) p[i] += 0.5 * eps * grad_lp[i];
    }
  }
  return true;
}


}  // namespace dmod
