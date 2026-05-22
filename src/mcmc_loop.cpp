// mcmcChainRun: single-chain mh / langevin / hmc / nuts driver.
// mcmcSmcReweight: per-level SMC reweight + log-evidence increment.

#include "mcmc_kernel.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>

extern "C" {
  void dpotri_(const char* uplo, const int* n, double* a, const int* lda,
               int* info);
  void dtrmv_(const char* uplo, const char* trans, const char* diag,
              const int* n, const double* a, const int* lda,
              double* x, const int* incx);
  double unif_rand(void);
  double norm_rand(void);
}

using namespace Rcpp;


// Decode an objlist into log-posterior space, reordering grad and hessian to
// the canonical par_names order. dMod's value/gradient/hessian are in
// -2 log L units; we divide by 2 here.
static bool decode_objout(const List& out,
                          const CharacterVector& par_names,
                          double& logp_out,
                          std::vector<double>& grad_out,
                          std::vector<double>& G_out) {
  const int K = par_names.size();
  if (!out.containsElementNamed("value") || Rf_isNull(out["value"]))
    return false;
  const double v = as<double>(out["value"]);
  if (!std::isfinite(v)) {
    logp_out = -std::numeric_limits<double>::infinity();
    return false;
  }
  logp_out = -0.5 * v;

  std::fill(grad_out.begin(), grad_out.end(), 0.0);
  if (out.containsElementNamed("gradient") && !Rf_isNull(out["gradient"])) {
    NumericVector g = as<NumericVector>(out["gradient"]);
    CharacterVector gn;
    if (g.hasAttribute("names")) gn = g.names();
    if (gn.size() == g.size()) {
      for (int j = 0; j < K; ++j) {
        const std::string nm = as<std::string>(par_names[j]);
        for (int gi = 0; gi < g.size(); ++gi) {
          if (as<std::string>(gn[gi]) == nm) {
            grad_out[j] = -0.5 * g[gi];
            break;
          }
        }
      }
    } else if (g.size() == K) {
      for (int j = 0; j < K; ++j) grad_out[j] = -0.5 * g[j];
    }
  }

  std::fill(G_out.begin(), G_out.end(), 0.0);
  if (out.containsElementNamed("hessian") && !Rf_isNull(out["hessian"])) {
    NumericMatrix H = as<NumericMatrix>(out["hessian"]);
    if (H.nrow() == K && H.ncol() == K) {
      List dn = H.hasAttribute("dimnames")
                  ? as<List>(H.attr("dimnames")) : List();
      if (dn.size() == 2 && !Rf_isNull(dn[0]) && !Rf_isNull(dn[1])) {
        CharacterVector rn = as<CharacterVector>(dn[0]);
        CharacterVector cn = as<CharacterVector>(dn[1]);
        std::vector<int> r_perm(K, -1), c_perm(K, -1);
        for (int j = 0; j < K; ++j) {
          const std::string nm = as<std::string>(par_names[j]);
          for (int i = 0; i < rn.size(); ++i)
            if (as<std::string>(rn[i]) == nm) { r_perm[j] = i; break; }
          for (int i = 0; i < cn.size(); ++i)
            if (as<std::string>(cn[i]) == nm) { c_perm[j] = i; break; }
        }
        for (int j = 0; j < K; ++j)
          for (int i = 0; i < K; ++i)
            if (r_perm[i] >= 0 && c_perm[j] >= 0)
              G_out[i + j * K] = 0.5 * H(r_perm[i], c_perm[j]);
      } else {
        for (int j = 0; j < K; ++j)
          for (int i = 0; i < K; ++i)
            G_out[i + j * K] = 0.5 * H(i, j);
      }
    }
  }
  return true;
}


// Invoke objfun(pars = theta, deriv = TRUE, deriv2 = FALSE). Any user
// "dots" must be baked into the objfun closure at the R level before it
// reaches this driver (cleaner than reconstructing variadic R calls in
// C++). Returns false on R-side error or non-finite value.
static bool eval_objfun(Function& objfun,
                        const NumericVector& theta,
                        const CharacterVector& par_names,
                        double& logp_out,
                        std::vector<double>& grad_out,
                        std::vector<double>& G_out) {
  try {
    SEXP res = objfun(_["pars"] = theta, _["deriv"] = true,
                      _["deriv2"] = false);
    List out = as<List>(res);
    return decode_objout(out, par_names, logp_out, grad_out, G_out);
  } catch (...) {
    logp_out = -std::numeric_limits<double>::infinity();
    std::fill(grad_out.begin(), grad_out.end(), 0.0);
    std::fill(G_out.begin(),    G_out.end(),    0.0);
    return false;
  }
}


// Build the metric G and its upper Cholesky L given the preconditioner
// code: 0 = LOCAL (use objfun hessian, already in G_buf), 1 = FIXED
// (overwrite G with GFixed/2), 2 = IDENTITY. Returns the success flag of
// the Cholesky decomposition and writes ridge_used.
static bool build_metric(int precond_code, int K,
                         std::vector<double>& G_buf,
                         const double* GFixed_half,
                         double ridge0,
                         std::vector<double>& L_buf,
                         double& ridge_used,
                         double& logdetG_out) {
  if (precond_code == 1) {
    std::memcpy(G_buf.data(), GFixed_half,
                sizeof(double) * static_cast<size_t>(K) * K);
  } else if (precond_code == 2) {
    std::fill(G_buf.begin(), G_buf.end(), 0.0);
    for (int i = 0; i < K; ++i) G_buf[i + i * K] = 1.0;
  }
  if (!dmod::chol_with_ridge(G_buf.data(), K, ridge0, 6,
                             L_buf.data(), &ridge_used))
    return false;
  double s = 0.0;
  for (int j = 0; j < K; ++j) s += std::log(L_buf[j + j * K]);
  logdetG_out = 2.0 * s;
  return true;
}


// Decode dG callback (returns a K x K x K array, column-major).
// Returns true on success.
static bool eval_dG_cb(Function& dG_cb, const NumericVector& theta,
                       int K, std::vector<double>& dG_buf) {
  try {
    SEXP res = dG_cb(theta);
    NumericVector arr = as<NumericVector>(res);
    if (arr.size() != K * K * K) return false;
    // R column-major arr[i + j*K + k*K*K] = dG[i, j, k] (where the
    // last dimension is the differentiation index). mala_drift expects
    // row-major dG_buf[a*K*K + b*K + c] = dG[a, b, c].
    for (int k = 0; k < K; ++k)
      for (int j = 0; j < K; ++j)
        for (int i = 0; i < K; ++i)
          dG_buf[i * K * K + j * K + k] = arr[i + j * K + k * K * K];
    return true;
  } catch (...) {
    return false;
  }
}


// Bound check.
static bool inside_bounds(const double* theta, const double* lower,
                          const double* upper, int K) {
  for (int i = 0; i < K; ++i) {
    if (!std::isfinite(theta[i])) return false;
    if (theta[i] > upper[i] || theta[i] < lower[i]) return false;
  }
  return true;
}


// Heuristic initial step size for preconditioned MALA. Since the proposal
// covariance is eps * G^{-1}, the "natural" scale is already in the
// preconditioning; Roberts-Rosenthal (1998) suggests eps ~ K^{-1/3} for
// the optimal acceptance regime in this preconditioned coordinate. Using
// a trace(G^{-1})/K factor (the previous heuristic) makes eps too small
// when G has large eigenvalues, causing dual averaging to overshoot
// during warmup. The L argument is kept for the signature but no longer
// consulted.
static double heuristic_eps(int K, const std::vector<double>& /*L*/) {
  double eps = 0.5 * std::pow(static_cast<double>(K), -1.0 / 3.0);
  if (!std::isfinite(eps) || eps <= 0) eps = 0.1;
  return std::max(eps, 1e-6);
}


// Naive ESS (first-zero autocovariance truncation).
static NumericVector naive_ess(const NumericMatrix& samples) {
  const int N = samples.nrow();
  const int K = samples.ncol();
  NumericVector ess(K);
  for (int j = 0; j < K; ++j) {
    double m = 0.0;
    for (int i = 0; i < N; ++i) m += samples(i, j);
    m /= N;
    double v = 0.0;
    std::vector<double> xc(N);
    for (int i = 0; i < N; ++i) { xc[i] = samples(i, j) - m; v += xc[i] * xc[i]; }
    v /= N;
    if (v <= 0) { ess[j] = N; continue; }
    const int max_lag = std::max(1, N / 3);
    double tau = 0.0;
    for (int lag = 1; lag <= max_lag; ++lag) {
      double r = 0.0;
      for (int i = 0; i < N - lag; ++i) r += xc[i] * xc[i + lag];
      r /= (v * N);
      const bool below = r < 0.05;
      if (r > 0) tau += r;
      if (below) break;
    }
    ess[j] = static_cast<double>(N) / (1.0 + 2.0 * tau);
  }
  CharacterVector cn = colnames(samples);
  if (cn.size() == K) ess.attr("names") = cn;
  return ess;
}


struct NutsState {
  std::vector<double> theta;
  std::vector<double> p;
  std::vector<double> grad_lp;
  double logp;
};

struct BuildTreeResult {
  NutsState minus, plus;
  std::vector<double> theta_prime;
  double logp_prime;
  double n_prime;
  bool   s_prime;
  double alpha;
  int    n_alpha;
};


// One leapfrog step; returns false on objfn failure.
static bool leapfrog_step(Function& objfun, NutsState& s,
                          const std::vector<double>& Minv_chol_upper,
                          int K, double eps,
                          const CharacterVector& par_names,
                          std::vector<double>& dummy_G) {
  for (int i = 0; i < K; ++i) s.p[i] += 0.5 * eps * s.grad_lp[i];

  std::vector<double> Mp(s.p);
  const int n = K, incx = 1;
  const char uplo = 'U', diag = 'N';
  const char transN = 'N', transT = 'T';
  dtrmv_(&uplo, &transN, &diag, &n, Minv_chol_upper.data(), &n,
         Mp.data(), &incx);
  dtrmv_(&uplo, &transT, &diag, &n, Minv_chol_upper.data(), &n,
         Mp.data(), &incx);
  for (int i = 0; i < K; ++i) s.theta[i] += eps * Mp[i];

  NumericVector th_nv(s.theta.begin(), s.theta.end());
  th_nv.attr("names") = par_names;
  std::vector<double> grad_buf(K);
  if (!eval_objfun(objfun, th_nv, par_names, s.logp, grad_buf, dummy_G))
    return false;
  s.grad_lp = grad_buf;

  for (int i = 0; i < K; ++i) s.p[i] += 0.5 * eps * s.grad_lp[i];
  return true;
}


static BuildTreeResult build_tree(Function& objfun, NutsState s,
                                   double u, int v, int j,
                                   double eps, double H0,
                                   const std::vector<double>& Minv_chol_upper,
                                   int K, double delta_max,
                                   const CharacterVector& par_names,
                                   std::vector<double>& dummy_G) {
  BuildTreeResult res;
  if (j == 0) {
    const bool ok = leapfrog_step(objfun, s, Minv_chol_upper, K,
                                   static_cast<double>(v) * eps,
                                   par_names, dummy_G);
    res.minus = res.plus = s;
    res.theta_prime.assign(s.theta.begin(), s.theta.end());
    if (!ok) {
      res.logp_prime = -std::numeric_limits<double>::infinity();
      res.n_prime = 0.0;
      res.s_prime = false;
      res.alpha   = 0.0;
      res.n_alpha = 1;
      return res;
    }
    const double H = -s.logp + dmod::kinetic_energy(s.p.data(),
                                                     Minv_chol_upper.data(), K);
    res.logp_prime = s.logp;
    res.n_prime = std::exp(-H);
    res.s_prime = (res.n_prime > u) && std::isfinite(H) &&
                  (H - H0 < delta_max);
    res.alpha   = std::min(1.0, std::exp(H0 - H));
    if (!std::isfinite(res.alpha)) res.alpha = 0.0;
    res.n_alpha = 1;
    return res;
  }

  BuildTreeResult L = build_tree(objfun, s, u, v, j - 1, eps, H0,
                                  Minv_chol_upper, K, delta_max,
                                  par_names, dummy_G);
  if (!L.s_prime) return L;
  BuildTreeResult R;
  if (v == -1)
    R = build_tree(objfun, L.minus, u, v, j - 1, eps, H0,
                   Minv_chol_upper, K, delta_max, par_names, dummy_G);
  else
    R = build_tree(objfun, L.plus,  u, v, j - 1, eps, H0,
                   Minv_chol_upper, K, delta_max, par_names, dummy_G);

  const double total = L.n_prime + R.n_prime;
  BuildTreeResult out;
  if (total > 0.0 && ::unif_rand() < R.n_prime / total) {
    out.theta_prime = R.theta_prime;
    out.logp_prime  = R.logp_prime;
  } else {
    out.theta_prime = L.theta_prime;
    out.logp_prime  = L.logp_prime;
  }

  NutsState minus_out = (v == -1) ? R.minus : L.minus;
  NutsState plus_out  = (v == -1) ? L.plus  : R.plus;
  std::vector<double> dtheta(K);
  for (int i = 0; i < K; ++i)
    dtheta[i] = plus_out.theta[i] - minus_out.theta[i];
  double dot_minus = 0.0, dot_plus = 0.0;
  for (int i = 0; i < K; ++i) {
    dot_minus += dtheta[i] * minus_out.p[i];
    dot_plus  += dtheta[i] * plus_out.p[i];
  }

  out.minus = minus_out;
  out.plus  = plus_out;
  out.n_prime = total;
  out.s_prime = L.s_prime && R.s_prime &&
                (dot_minus >= 0.0) && (dot_plus >= 0.0);
  out.alpha   = L.alpha + R.alpha;
  out.n_alpha = L.n_alpha + R.n_alpha;
  return out;
}


// [[Rcpp::export(name = "mcmcChainRun")]]
Rcpp::List mcmc_chain_run(Rcpp::Function objfun,
                          Rcpp::NumericVector parinit,
                          int n, int warmup,
                          int moveType,
                          Rcpp::List control,
                          Rcpp::List bounds,
                          Rcpp::NumericVector parscale_,
                          Rcpp::Nullable<Rcpp::Function> dG_cb_opt) {
  (void)parscale_;
  const int K = parinit.size();
  Rcpp::CharacterVector par_names = parinit.names();
  if (par_names.size() != K)
    Rcpp::stop("mcmcChainRun: parinit must have names.");

  Rcpp::NumericVector upper = bounds["upper"];
  Rcpp::NumericVector lower = bounds["lower"];

  const double ridge0       = Rcpp::as<double>(control["ridge"]);
  const double daGamma      = Rcpp::as<double>(control["daGamma"]);
  const double daT0         = Rcpp::as<double>(control["daT0"]);
  const double daKappa      = Rcpp::as<double>(control["daKappa"]);
  const double acceptTarget = Rcpp::as<double>(control["acceptTarget"]);
  const bool has_stepsize = !Rf_isNull(control["stepsize"]);
  const double init_stepsize = has_stepsize
      ? Rcpp::as<double>(control["stepsize"]) : -1.0;

  const int precond_code = Rcpp::as<int>(control["preconditioner"]);
  const bool correction_full = (moveType == 1) &&
      Rcpp::as<bool>(control["correctionFull"]) && !dG_cb_opt.isNull();

  std::vector<double> GFixed_half;
  if (precond_code == 1) {
    Rcpp::NumericMatrix Gf = Rcpp::as<Rcpp::NumericMatrix>(control["GFixed"]);
    if (Gf.nrow() != K || Gf.ncol() != K)
      Rcpp::stop("mcmcChainRun: GFixed must be K x K matching parinit.");
    // GFixed is in -2 log L units (matches trust()$hessian); the metric G
    // is the Fisher info = hessian / 2.
    GFixed_half.assign(static_cast<size_t>(K) * K, 0.0);
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i)
        GFixed_half[i + j * K] = 0.5 * Gf(i, j);
  }

  const int leapfrogSteps = (moveType == 2)
                              ? Rcpp::as<int>(control["leapfrogSteps"]) : 0;
  const int maxTreeDepth  = (moveType == 3)
                              ? Rcpp::as<int>(control["maxTreeDepth"]) : 10;
  const double deltaMax   = (moveType == 3)
                              ? Rcpp::as<double>(control["deltaMax"]) : 1000.0;

  std::vector<double> M_chol_upper(static_cast<size_t>(K) * K, 0.0);
  std::vector<double> Minv_chol_upper(static_cast<size_t>(K) * K, 0.0);
  if (moveType == 2 || moveType == 3) {
    Rcpp::NumericMatrix Mc  = Rcpp::as<Rcpp::NumericMatrix>(control["MCholUpper"]);
    Rcpp::NumericMatrix MiC = Rcpp::as<Rcpp::NumericMatrix>(control["MinvCholUpper"]);
    if (Mc.nrow() != K || MiC.nrow() != K)
      Rcpp::stop("mcmcChainRun: M[Inv]CholUpper must be K x K.");
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i) {
        M_chol_upper[i + j * K]    = Mc(i, j);
        Minv_chol_upper[i + j * K] = MiC(i, j);
      }
  }

  Rcpp::RNGScope rng_scope;

  std::vector<double> theta_cur(parinit.begin(), parinit.end());
  std::vector<double> theta_new(K, 0.0);
  std::vector<double> grad_cur(K, 0.0);
  std::vector<double> grad_new(K, 0.0);
  std::vector<double> G_cur(static_cast<size_t>(K) * K, 0.0);
  std::vector<double> G_new(static_cast<size_t>(K) * K, 0.0);
  std::vector<double> L_cur(static_cast<size_t>(K) * K, 0.0);
  std::vector<double> L_new(static_cast<size_t>(K) * K, 0.0);
  std::vector<double> drift_cur(K, 0.0);
  std::vector<double> drift_new(K, 0.0);
  std::vector<double> z_buf(K, 0.0);
  std::vector<double> dG_cur(correction_full ?
                              (size_t)K * K * K : 0, 0.0);
  std::vector<double> dG_new(correction_full ?
                              (size_t)K * K * K : 0, 0.0);
  std::vector<double> dummy_G_unused(static_cast<size_t>(K) * K, 0.0);

  double logp_cur = -std::numeric_limits<double>::infinity();
  {
    Rcpp::NumericVector th_init(theta_cur.begin(), theta_cur.end());
    th_init.attr("names") = par_names;
    if (!eval_objfun(objfun, th_init, par_names, logp_cur, grad_cur, G_cur)
        || !std::isfinite(logp_cur))
      Rcpp::stop("Initial log-posterior is non-finite at parinit.");
  }

  double ridge_used_cur = 0.0, logdetG_cur = 0.0;

  Rcpp::Function dG_cb = correction_full
      ? Rcpp::as<Rcpp::Function>(dG_cb_opt)
      : Rcpp::Function("identity");

  if (moveType == 0 || moveType == 1) {
    if (!build_metric(precond_code, K, G_cur, GFixed_half.data(), ridge0,
                       L_cur, ridge_used_cur, logdetG_cur))
      Rcpp::stop("Initial metric Cholesky failed.");
    if (moveType == 1) {
      if (correction_full) {
        Rcpp::NumericVector th_nv(theta_cur.begin(), theta_cur.end());
        th_nv.attr("names") = par_names;
        if (!eval_dG_cb(dG_cb, th_nv, K, dG_cur))
          Rcpp::stop("Initial dG callback failed for correction='full'.");
      }
      dmod::mala_drift(grad_cur.data(), L_cur.data(), K,
                       correction_full ? dG_cur.data() : nullptr,
                       drift_cur.data());
    }
  }

  double eps = init_stepsize > 0.0
      ? init_stepsize
      : ((moveType == 0 || moveType == 1)
           ? heuristic_eps(K, L_cur) : 0.1);

  dmod::DAState da;
  dmod::da_init(da, eps, acceptTarget, daGamma, daT0, daKappa);

  const int total = warmup + n;
  Rcpp::NumericMatrix samples_full(total, K);
  Rcpp::colnames(samples_full) = par_names;
  Rcpp::NumericVector logp_full(total);
  Rcpp::LogicalVector accept_full(total);
  Rcpp::NumericVector eps_full(total);
  Rcpp::IntegerVector treedepth_full(total);

  for (int it = 0; it < total; ++it) {
    const bool is_warmup = (it < warmup);
    bool accepted = false;
    double alpha = 0.0;

    if (moveType == 0) {
      // ---- MH ----
      dmod::mh_propose(theta_cur.data(), L_cur.data(), K, eps,
                       z_buf.data(), theta_new.data());
      if (!inside_bounds(theta_new.data(), lower.begin(), upper.begin(), K)) {
        alpha = 0.0;
      } else {
        Rcpp::NumericVector th_nv(theta_new.begin(), theta_new.end());
        th_nv.attr("names") = par_names;
        double logp_prop = -std::numeric_limits<double>::infinity();
        const bool ok = eval_objfun(objfun, th_nv, par_names,
                                    logp_prop, grad_new, G_new);
        if (!ok || !std::isfinite(logp_prop)) {
          alpha = 0.0;
        } else {
          double log_alpha = logp_prop - logp_cur;
          alpha = std::min(1.0, std::exp(log_alpha));
          if (!std::isfinite(alpha)) alpha = 0.0;
          if (::unif_rand() < alpha) {
            theta_cur = theta_new;
            logp_cur  = logp_prop;
            grad_cur  = grad_new;
            G_cur     = G_new;
            accepted  = true;
          }
        }
      }
    } else if (moveType == 1) {
      // ---- Langevin (MALA) ----
      dmod::mala_propose(theta_cur.data(), drift_cur.data(),
                         L_cur.data(), K, eps, z_buf.data(),
                         theta_new.data());
      double logp_prop = -std::numeric_limits<double>::infinity();
      if (!inside_bounds(theta_new.data(), lower.begin(), upper.begin(), K)) {
        alpha = 0.0;
      } else {
        Rcpp::NumericVector th_nv(theta_new.begin(), theta_new.end());
        th_nv.attr("names") = par_names;
        const bool ok = eval_objfun(objfun, th_nv, par_names,
                                    logp_prop, grad_new, G_new);
        if (!ok || !std::isfinite(logp_prop)) {
          alpha = 0.0;
        } else {
          double ridge_used_new = 0.0, logdetG_new = 0.0;
          bool metric_ok = build_metric(precond_code, K, G_new,
                                         GFixed_half.data(), ridge0, L_new,
                                         ridge_used_new, logdetG_new);
          if (!metric_ok) {
            alpha = 0.0;
          } else {
            const double* dG_new_ptr = nullptr;
            if (correction_full) {
              if (!eval_dG_cb(dG_cb, th_nv, K, dG_new)) {
                alpha = 0.0;
                metric_ok = false;
              } else {
                dG_new_ptr = dG_new.data();
              }
            }
            if (metric_ok) {
              dmod::mala_drift(grad_new.data(), L_new.data(), K, dG_new_ptr,
                               drift_new.data());
              const double logq_fwd = dmod::mala_log_q(
                  theta_new.data(), theta_cur.data(), drift_cur.data(),
                  G_cur.data(), K, eps, logdetG_cur);
              const double logq_rev = dmod::mala_log_q(
                  theta_cur.data(), theta_new.data(), drift_new.data(),
                  G_new.data(), K, eps, logdetG_new);
              const double log_alpha = logp_prop - logp_cur +
                                       logq_rev - logq_fwd;
              alpha = std::min(1.0, std::exp(log_alpha));
              if (!std::isfinite(alpha)) alpha = 0.0;
              if (::unif_rand() < alpha) {
                theta_cur = theta_new;
                logp_cur  = logp_prop;
                grad_cur  = grad_new;
                G_cur     = G_new;
                L_cur     = L_new;
                drift_cur = drift_new;
                if (correction_full) dG_cur = dG_new;
                logdetG_cur = logdetG_new;
                accepted = true;
              }
            }
          }
        }
      }
    } else if (moveType == 2) {
      // ---- HMC ----
      std::vector<double> p_cur(K, 0.0);
      dmod::sample_momentum(M_chol_upper.data(), K, z_buf.data(),
                            p_cur.data());
      const double H0 = -logp_cur + dmod::kinetic_energy(
          p_cur.data(), Minv_chol_upper.data(), K);

      NutsState s;
      s.theta = theta_cur; s.p = p_cur; s.grad_lp = grad_cur; s.logp = logp_cur;
      bool diverged = false;
      for (int lf = 0; lf < leapfrogSteps; ++lf) {
        if (!leapfrog_step(objfun, s, Minv_chol_upper, K, eps,
                            par_names, dummy_G_unused)) {
          diverged = true; break;
        }
        if (!std::isfinite(s.logp)) { diverged = true; break; }
      }
      if (diverged ||
          !inside_bounds(s.theta.data(), lower.begin(), upper.begin(), K)) {
        alpha = 0.0;
      } else {
        const double H1 = -s.logp + dmod::kinetic_energy(
            s.p.data(), Minv_chol_upper.data(), K);
        const double log_alpha = H0 - H1;
        alpha = std::min(1.0, std::exp(log_alpha));
        if (!std::isfinite(alpha)) alpha = 0.0;
        if (::unif_rand() < alpha) {
          theta_cur = s.theta;
          logp_cur  = s.logp;
          grad_cur  = s.grad_lp;
          accepted  = true;
        }
      }
    } else if (moveType == 3) {
      // ---- NUTS (efficient algorithm 6) ----
      std::vector<double> p_cur(K, 0.0);
      dmod::sample_momentum(M_chol_upper.data(), K, z_buf.data(),
                            p_cur.data());
      const double H0 = -logp_cur + dmod::kinetic_energy(
          p_cur.data(), Minv_chol_upper.data(), K);
      const double u = ::unif_rand() * std::exp(-H0);

      NutsState init_state;
      init_state.theta   = theta_cur;
      init_state.p       = p_cur;
      init_state.grad_lp = grad_cur;
      init_state.logp    = logp_cur;
      NutsState minus = init_state, plus = init_state;

      std::vector<double> chosen_theta = theta_cur;
      double chosen_logp = logp_cur;
      double n_total = 1.0;
      bool   s_total = true;
      int    j = 0;
      double alpha_sum = 0.0;
      int    n_alpha = 0;

      while (s_total && j < maxTreeDepth) {
        const int v = (::unif_rand() < 0.5) ? -1 : 1;
        NutsState base = (v == -1) ? minus : plus;
        BuildTreeResult br = build_tree(objfun, base, u, v, j, eps, H0,
                                         Minv_chol_upper, K, deltaMax,
                                         par_names, dummy_G_unused);
        if (v == -1) minus = br.minus; else plus = br.plus;
        alpha_sum += br.alpha;
        n_alpha   += br.n_alpha;
        if (br.s_prime) {
          const double accept_prob = (n_total > 0.0)
              ? std::min(1.0, br.n_prime / n_total) : 1.0;
          if (::unif_rand() < accept_prob) {
            chosen_theta = br.theta_prime;
            chosen_logp  = br.logp_prime;
          }
        }
        n_total += br.n_prime;
        std::vector<double> dtheta(K);
        for (int i = 0; i < K; ++i)
          dtheta[i] = plus.theta[i] - minus.theta[i];
        double dot_m = 0.0, dot_p = 0.0;
        for (int i = 0; i < K; ++i) {
          dot_m += dtheta[i] * minus.p[i];
          dot_p += dtheta[i] * plus.p[i];
        }
        s_total = br.s_prime && (dot_m >= 0.0) && (dot_p >= 0.0);
        ++j;
      }

      if (inside_bounds(chosen_theta.data(), lower.begin(), upper.begin(), K)
          && chosen_logp > -std::numeric_limits<double>::infinity()) {
        Rcpp::NumericVector th_nv(chosen_theta.begin(), chosen_theta.end());
        th_nv.attr("names") = par_names;
        std::vector<double> grad_chosen(K, 0.0);
        double logp_chosen = chosen_logp;
        if (eval_objfun(objfun, th_nv, par_names,
                        logp_chosen, grad_chosen, dummy_G_unused)) {
          theta_cur = chosen_theta;
          logp_cur  = logp_chosen;
          grad_cur  = grad_chosen;
          accepted  = true;
        }
      }
      alpha = (n_alpha > 0) ? (alpha_sum / n_alpha) : 0.0;
      treedepth_full[it] = j;
    }

    for (int i = 0; i < K; ++i) samples_full(it, i) = theta_cur[i];
    logp_full[it]   = logp_cur;
    accept_full[it] = accepted;
    eps_full[it]    = eps;

    if (is_warmup) {
      dmod::da_update(da, alpha, it + 1);
      eps = da.eps;
    } else if (it == warmup && warmup > 0) {
      eps = da.eps_bar;
    }

    if ((it & 1023) == 0) Rcpp::checkUserInterrupt();
  }

  Rcpp::NumericMatrix samples_post(n, K);
  Rcpp::colnames(samples_post) = par_names;
  Rcpp::NumericVector logp_post(n);
  Rcpp::LogicalVector accept_post(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < K; ++j) samples_post(i, j) = samples_full(warmup + i, j);
    logp_post[i]   = logp_full[warmup + i];
    accept_post[i] = accept_full[warmup + i];
  }

  Rcpp::NumericVector ess = (n > 1) ? naive_ess(samples_post)
                                    : Rcpp::NumericVector(K, NA_REAL);

  double acc_rate = NA_REAL;
  if (n > 0) {
    int cnt = 0;
    for (int i = 0; i < n; ++i) if (accept_post[i]) ++cnt;
    acc_rate = static_cast<double>(cnt) / n;
  }

  return Rcpp::List::create(
      Rcpp::_["samples"]    = samples_post,
      Rcpp::_["logp"]       = logp_post,
      Rcpp::_["accept"]     = accept_full,
      Rcpp::_["acceptRate"] = acc_rate,
      Rcpp::_["stepsize"]   = eps_full,
      Rcpp::_["finalStep"]  = eps,
      Rcpp::_["ess"]        = ess,
      Rcpp::_["treedepth"]  = treedepth_full);
}


// [[Rcpp::export(name = "mcmcSmcReweight")]]
Rcpp::List mcmc_smc_reweight(const Rcpp::NumericVector& logL,
                              const Rcpp::NumericVector& logwPrev,
                              double betaOld, double betaNew) {
  const int N = logL.size();
  Rcpp::NumericVector logw_new(N);
  std::vector<double> logwUn(N);
  double m_un = -std::numeric_limits<double>::infinity();
  double m_prev = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < N; ++i) {
    logwUn[i] = logwPrev[i] + (betaNew - betaOld) * logL[i];
    if (logwUn[i]   > m_un)   m_un   = logwUn[i];
    if (logwPrev[i] > m_prev) m_prev = logwPrev[i];
  }
  double s_un = 0.0, s_prev = 0.0;
  for (int i = 0; i < N; ++i) {
    s_un   += std::exp(logwUn[i] - m_un);
    s_prev += std::exp(logwPrev[i] - m_prev);
  }
  const double lse_un   = m_un   + std::log(s_un);
  const double lse_prev = m_prev + std::log(s_prev);
  const double logZinc  = lse_un - lse_prev;
  for (int i = 0; i < N; ++i) logw_new[i] = logwUn[i] - lse_un;

  return Rcpp::List::create(Rcpp::_["logw"]    = logw_new,
                            Rcpp::_["logZinc"] = logZinc);
}
