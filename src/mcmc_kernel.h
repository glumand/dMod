// MCMC inner-loop primitives.
//
// Sign convention. objfn returns -2 log L; the kernel maps this to
// log-posterior space:  logp = -value / 2,  grad_lp = -gradient / 2,
// G = hessian / 2.
//
// Metric convention. All Cholesky-using primitives below take the upper
// factor L with G = L^T L (matches R's chol() and LAPACK dpotrf(uplo='U')).
// Drawing u = L^{-1} z with z ~ N(0, I) gives Cov(u) = G^{-1}; solving
// G x = b means L^T L x = b, done by dpotrs.
//
// Threading. Every function here is reentrant: callers manage their own
// buffers. Rcpp::Function callbacks (the only R touchpoints) must happen
// on the orchestrator thread; OpenMP regions must not re-enter R.

#ifndef DMOD_MCMC_KERNEL_H
#define DMOD_MCMC_KERNEL_H

#include <Rcpp.h>
#include <functional>
#include <vector>

namespace dmod {

// Metric / preconditioner selection inside evaluate_target.
enum class PreCond {
  LOCAL    = 0,    // G_local = hessian / 2 from objfn
  FIXED    = 1,    // G = GFixed / 2 (GFixed is in -2 log L units)
  IDENTITY = 2     // G = I
};

// What state the kernel returns from one objfn call. Sized to the parameter
// dimension K; grad / G / L all carry par names as Rcpp dimnames so the
// caller can re-index by name.
struct LogPState {
  double logp;                  // log target density at theta
  Rcpp::NumericVector grad;     // d logp / d theta, length K
  Rcpp::NumericMatrix G;        // metric tensor G(theta), K x K
  Rcpp::NumericMatrix L;        // upper-triangular U with G = U^T U
  double logdetG;               // 2 * sum log diag(L)
  double ridge_used;            // ridge actually added to G before chol
  bool ok;                      // false if objfn was non-finite or chol failed
};

// Dual-averaging state (Hoffman-Gelman 2014 Algorithm 5).
struct DAState {
  double eps;        // current step size
  double eps_bar;    // averaged step size, returned at end of warmup
  double Hbar;       // running statistic
  double mu;         // shrinkage target = log(10 * eps0)
  double gamma;      // adaptation rate (default 0.05)
  double t0;         // adaptation delay (default 10)
  double kappa;      // averaging power (default 0.75)
  double target;     // target acceptance probability
};

void da_init(DAState& da, double eps0, double target,
             double gamma = 0.05, double t0 = 10.0,
             double kappa = 0.75);

void da_update(DAState& da, double alpha, int t);

// Cholesky with adaptive ridge: tries G + ridge*I, doubling ridge up to
// max_iter times if dpotrf fails. Writes the upper factor into L_out
// (column-major K*K), and the actual ridge used into ridge_out.
// Returns true on success.
bool chol_with_ridge(const double* G, int K, double ridge0, int max_iter,
                     double* L_out, double* ridge_out);

// Natural-gradient drift d = G^{-1} grad. With non-null dG (K*K*K
// row-major), adds the Xifara-corrected geodesic terms:
//   d = G^{-1} grad - 2 * termB + termA,
//   termA_j = (G^{-1} v)_j,  v_k = tr(G^{-1} dG_k),
//   termB_k = sum_j (G^{-1} dG_k G^{-1})_{j,k}.
// The eps/2 multiplier is applied by the caller.
void mala_drift(const double* grad, const double* L_upper, int K,
                const double* dG_row_major,        // K*K*K or nullptr
                double* drift_out);

// MALA proposal: theta_new = theta_old + 0.5 * eps * drift_old +
// sqrt(eps) * L_upper^{-1} z, with z drawn into z_buf.
void mala_propose(const double* theta_old, const double* drift_old,
                  const double* L_upper, int K, double eps,
                  double* z_buf, double* theta_new);

// Log proposal density q(theta_new | theta_old) up to additive constants
// that cancel in the MH ratio:
//   -0.5 / eps * (theta_new - mu)^T G_old (theta_new - mu) + 0.5 logdet_old
double mala_log_q(const double* theta_new, const double* theta_old,
                  const double* drift_old, const double* G_old,
                  int K, double eps, double logdet_old);

// Random-walk MH proposal: theta_new = theta_old + sqrt(eps) * L^{-1} z.
// Same draw as MALA but without the gradient drift; symmetric proposal so
// log q drops out of the MH ratio.
void mh_propose(const double* theta_old, const double* L_upper, int K,
                double eps, double* z_buf, double* theta_new);

// Leapfrog integrator for HMC. Mass matrix is given through its inverse
// Cholesky factor (M = LM^T LM, Minv = LM^{-1} LM^{-T}). The kernel
// invokes grad_cb once per half-step (twice per full step except at the
// trajectory endpoints, where the half-steps share). p and theta are
// updated in place. grad_lp is filled with the gradient at the final
// theta. Returns true on success.
//
// grad_cb(theta) must return a LogPState with at least logp + grad
// populated (G and L need not be touched).
bool leapfrog(double* theta, double* p, double* grad_lp,
              const double* Minv_chol_upper,   // K x K, upper U with Minv = U^T U
              int K, double eps, int n_steps,
              std::function<bool(const double*, double&, double*)> grad_cb);

// Inverse mass matrix log-density evaluation: 0.5 * p^T Minv p. Used in
// the HMC / NUTS Hamiltonian.
double kinetic_energy(const double* p, const double* Minv_chol_upper, int K);

// Sample momentum p ~ N(0, M) into p_out. Receives the *upper* Cholesky
// of the mass matrix M = LM^T LM (i.e. p = LM^T z for z ~ N(0, I)).
void sample_momentum(const double* LM_upper, int K, double* z_buf,
                     double* p_out);

}  // namespace dmod

#endif  // DMOD_MCMC_KERNEL_H
