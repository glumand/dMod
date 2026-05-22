// SMC sampler helper kernels: stable log-sum-exp, ESS, systematic
// resampling, and adaptive beta bisection. Each function is a thin Rcpp
// wrapper around a numerically robust C++ implementation; the SMC R
// orchestrator (R/sampleSMC.R) calls these to keep per-iteration overhead
// out of the R interpreter even for nParticles in the thousands.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Rcpp;


// Stable log-sum-exp over a numeric vector.
//
// [[Rcpp::export(name = "smcLogSumExp")]]
double smc_log_sum_exp(const NumericVector& x) {
  const int n = x.size();
  if (n == 0L) return -std::numeric_limits<double>::infinity();
  double m = x[0];
  for (int i = 1; i < n; ++i) if (x[i] > m) m = x[i];
  if (!std::isfinite(m)) return m;
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += std::exp(x[i] - m);
  return m + std::log(s);
}


// Effective sample size from log-weights.
// ESS = (sum w)^2 / sum w^2 = exp(2*LSE(logw) - LSE(2*logw))
//
// [[Rcpp::export(name = "smcESS")]]
double smc_ess(const NumericVector& logw) {
  const int n = logw.size();
  if (n == 0L) return 0.0;
  const double lse1 = smc_log_sum_exp(logw);
  NumericVector two = 2.0 * logw;
  const double lse2 = smc_log_sum_exp(two);
  if (!std::isfinite(lse1) || !std::isfinite(lse2)) return 0.0;
  return std::exp(2.0 * lse1 - lse2);
}


// O(N) systematic resampling. `weights` must be normalised. `u` is a single
// U[0,1) variate; offsets u/N, (u+1)/N, ... pick indices via the inverse
// CDF (cumulative sum). Returns 1-based indices for direct R use.
//
// [[Rcpp::export(name = "smcSystematicResample")]]
IntegerVector smc_systematic_resample(const NumericVector& weights, double u) {
  const int N = weights.size();
  IntegerVector idx(N);
  if (N == 0L) return idx;
  double c = weights[0];
  int j = 0;
  for (int i = 0; i < N; ++i) {
    const double up = (i + u) / static_cast<double>(N);
    while (c < up && j < N - 1) {
      ++j;
      c += weights[j];
    }
    idx[i] = j + 1; // 1-based for R
  }
  return idx;
}


// Adaptive beta step (Jasra et al. 2011): bisection on
//   ESS(logw_prev + (beta_new - beta_old) * (-value_lik / 2)) = targetESS
// where `logL = -value_lik / 2`. Returns delta = beta_new - beta_old,
// clamped to [0, 1 - beta_old]. If even delta = 1 - beta_old keeps ESS
// above the target, returns 1 - beta_old.
//
// [[Rcpp::export(name = "smcBetaBisect")]]
double smc_beta_bisect(const NumericVector& logL,
                       const NumericVector& logwPrev,
                       double betaOld,
                       double targetESS,
                       double tol = 1e-6,
                       int maxIter = 80) {
  const int N = logL.size();
  if (N == 0L) return 0.0;
  const double deltaMax = std::max(0.0, 1.0 - betaOld);
  if (deltaMax <= 0.0) return 0.0;

  // Try the max step first. If ESS at deltaMax >= target, return deltaMax.
  NumericVector lw(N);
  auto ess_at = [&](double d) {
    for (int i = 0; i < N; ++i) lw[i] = logwPrev[i] + d * logL[i];
    return smc_ess(lw);
  };
  const double essMax = ess_at(deltaMax);
  if (essMax >= targetESS) return deltaMax;

  // Bisect on [lo, hi] with f = ESS - targetESS, monotone decreasing in d.
  double lo = 0.0, hi = deltaMax;
  for (int it = 0; it < maxIter; ++it) {
    const double mid = 0.5 * (lo + hi);
    const double e = ess_at(mid);
    if (std::fabs(e - targetESS) < tol) return mid;
    if (e > targetESS) lo = mid; else hi = mid;
    if (hi - lo < 1e-12) return mid;
  }
  return 0.5 * (lo + hi);
}
