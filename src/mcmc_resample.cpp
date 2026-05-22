// Stratified, residual, and multinomial resampling (companions to
// smcSystematicResample in src/smc_kernel.cpp). All routines take
// normalised weights and return 1-based indices for direct R use.

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;


// Stratified resampling. Divides [0, 1) into N equal strata and draws one
// uniform from each stratum, then maps via the inverse CDF (cumulative
// weights). Lower variance than multinomial, slightly higher than
// systematic. O(N).
//
// [[Rcpp::export(name = "mcmcStratifiedResample")]]
IntegerVector mcmc_stratified_resample(const NumericVector& weights) {
  const int N = weights.size();
  IntegerVector idx(N);
  if (N == 0L) return idx;
  double c = weights[0];
  int j = 0;
  for (int i = 0; i < N; ++i) {
    const double u = (i + R::unif_rand()) / static_cast<double>(N);
    while (c < u && j < N - 1) {
      ++j;
      c += weights[j];
    }
    idx[i] = j + 1;
  }
  return idx;
}


// Residual resampling. Deterministic part: each particle gets
// floor(N * w_i) copies. Remainder of weights (length R) is resampled
// multinomially to fill the rest.
//
// [[Rcpp::export(name = "mcmcResidualResample")]]
IntegerVector mcmc_residual_resample(const NumericVector& weights) {
  const int N = weights.size();
  IntegerVector idx(N);
  if (N == 0L) return idx;

  std::vector<int> floor_counts(N);
  std::vector<double> resid(N);
  int filled = 0;
  for (int i = 0; i < N; ++i) {
    const double nw = N * weights[i];
    floor_counts[i] = static_cast<int>(std::floor(nw));
    resid[i]        = nw - floor_counts[i];
    for (int c = 0; c < floor_counts[i] && filled < N; ++c) {
      idx[filled++] = i + 1;
    }
  }

  if (filled < N) {
    // Renormalise residuals (sum = N - filled) for multinomial sampling.
    double resid_sum = 0.0;
    for (int i = 0; i < N; ++i) resid_sum += resid[i];
    if (resid_sum <= 0.0) {
      // Numerical edge: fill remainder uniformly.
      for (; filled < N; ++filled)
        idx[filled] = static_cast<int>(N * R::unif_rand()) + 1;
      return idx;
    }
    std::vector<double> cdf(N);
    cdf[0] = resid[0] / resid_sum;
    for (int i = 1; i < N; ++i) cdf[i] = cdf[i - 1] + resid[i] / resid_sum;

    while (filled < N) {
      const double u = R::unif_rand();
      // Binary search in cdf
      int lo = 0, hi = N - 1;
      while (lo < hi) {
        const int mid = (lo + hi) / 2;
        if (cdf[mid] < u) lo = mid + 1; else hi = mid;
      }
      idx[filled++] = lo + 1;
    }
  }
  return idx;
}


// Plain multinomial resampling via the inverse CDF. O(N log N) due to
// binary search; useful as a reference / when variance is not a concern.
//
// [[Rcpp::export(name = "mcmcMultinomialResample")]]
IntegerVector mcmc_multinomial_resample(const NumericVector& weights) {
  const int N = weights.size();
  IntegerVector idx(N);
  if (N == 0L) return idx;
  std::vector<double> cdf(N);
  cdf[0] = weights[0];
  for (int i = 1; i < N; ++i) cdf[i] = cdf[i - 1] + weights[i];
  const double total = cdf[N - 1];
  for (int i = 0; i < N; ++i) {
    const double u = total * R::unif_rand();
    int lo = 0, hi = N - 1;
    while (lo < hi) {
      const int mid = (lo + hi) / 2;
      if (cdf[mid] < u) lo = mid + 1; else hi = mid;
    }
    idx[i] = lo + 1;
  }
  return idx;
}
