// Sparse-grid Gauss-Hermite quadrature for adaptive NLME marginal integration.
//
// Smolyak combination of 1D physicists' Gauss-Hermite rules (weight exp(-z^2),
// 1D total mass sqrt(pi)). 1D nodes/weights via Golub-Welsch on the Hermite
// Jacobi matrix (diagonal 0, off-diagonal sqrt(i/2)) using LAPACK dstev.
//
// Conventions:
//   - Output `nodes` is z-space, batch-first [B, K] (column-major
//     NumericMatrix; matches dMod's batch-first invariant for sensitivities).
//   - 1D level l -> m(l) = 2l - 1 nodes.
//   - Smolyak: A(L, K) = sum_{i: |i| in [L-K+1, L]} (-1)^(L-|i|) C(K-1, L-|i|)
//     * (U^{i_1} (x) ... (x) U^{i_K}). For non-nested GH rules this leaves
//     signed weights after dedup; downstream signed-LSE handles aggregation.
//   - The R wrapper folds the z -> eta change of variable
//     (eta = eta_hat + sqrt(2) * L^{-T} z, L = chol(H_i)) into log_weights.

#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>

#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

#ifndef FCONE
#define FCONE
#endif


/* --------------------------------------------------------------------------
 * 1D physicists' Gauss-Hermite via Golub-Welsch.
 *
 * Jacobi matrix J_n is n x n symmetric tridiagonal: diagonal 0, off-diagonal
 * sqrt(i/2). dstev returns eigenvalues (= nodes) and eigenvectors. The weight
 * for node i is sqrt(pi) * (first eigenvector component)^2.
 * -------------------------------------------------------------------------- */
static void gh_1d(int n, std::vector<double>& nodes, std::vector<double>& weights) {
  nodes.resize(n);
  weights.resize(n);
  if (n == 1) {
    nodes[0]   = 0.0;
    weights[0] = std::sqrt(M_PI);
    return;
  }

  std::vector<double> d(n, 0.0);
  std::vector<double> e(n - 1);
  for (int i = 0; i < n - 1; ++i) {
    e[i] = std::sqrt(static_cast<double>(i + 1) / 2.0);
  }
  std::vector<double> Z(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
  std::vector<double> work(static_cast<std::size_t>(std::max(1, 2 * n - 2)), 0.0);

  int info = 0;
  int n_   = n;
  int ldz  = n;
  F77_CALL(dstev)("V", &n_, d.data(), e.data(), Z.data(), &ldz, work.data(), &info FCONE);
  if (info != 0) Rcpp::stop("gh_1d: LAPACK dstev failed (info = %d).", info);

  const double sqrt_pi = std::sqrt(M_PI);
  for (int i = 0; i < n; ++i) {
    nodes[i] = d[i];
    double v0 = Z[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + 0];
    weights[i] = sqrt_pi * v0 * v0;
  }
}


/* --------------------------------------------------------------------------
 * Smolyak coefficient: (-1)^(L - q) * C(K - 1, L - q) where q = sum(i_j).
 * -------------------------------------------------------------------------- */
static double smolyak_coef(int L, int K, int q) {
  int d = L - q;
  if (d < 0 || d > K - 1) return 0.0;
  double c = 1.0;
  for (int k = 1; k <= d; ++k) {
    c *= static_cast<double>(K - k) / static_cast<double>(k);
  }
  return ((d % 2 == 0) ? 1.0 : -1.0) * c;
}


/* --------------------------------------------------------------------------
 * Enumerate length-K positive-integer multi-indices with sum equal to q.
 * Result appended to `out`. Used inside the Smolyak combination loop.
 * -------------------------------------------------------------------------- */
static void enumerate_compositions(int K, int q,
                                   std::vector<int>& cur,
                                   int idx,
                                   int remaining,
                                   std::vector<std::vector<int>>& out) {
  if (idx == K - 1) {
    if (remaining >= 1) {
      cur[idx] = remaining;
      out.push_back(cur);
    }
    return;
  }
  int hi = remaining - (K - 1 - idx);
  for (int v = 1; v <= hi; ++v) {
    cur[idx] = v;
    enumerate_compositions(K, q, cur, idx + 1, remaining - v, out);
  }
}


// Build the K-D Smolyak sparse grid for physicists' Gauss-Hermite at depth
// `level`. Returns nodes [B, K] in z-space and signed weights (length B).
//
// For non-nested 1D rules (GH is non-nested), Smolyak weights are mixed-sign.
// Dedup merges weights at coincident nodes; remaining sign mixture is handled
// downstream via signed log-sum-exp in the per-subject evaluator.
//' @name sparse_grid_gh
//' @title Sparse-grid Gauss-Hermite quadrature nodes (Smolyak rule)
//' @description Builds the K-dimensional Smolyak sparse grid for physicists'
//'   Gauss-Hermite at depth `level`. Returns nodes `[B, K]` in z-space and
//'   signed weights (length `B`).
//' @param K Integer >= 1, problem dimension.
//' @param level Integer >= K, Smolyak depth (K+1..K+3 is the useful range).
//' @param deriv_mode Reserved for future Genz-Keister / adaptive refinement;
//'   currently ignored.
//' @return A list with `nodes` (B x K, batch-first), `weights` (length B,
//'   signed), and `K`, `level`.
//' @export
// [[Rcpp::export]]
List sparse_grid_gh(int K, int level, int deriv_mode = 0) {
  (void) deriv_mode;
  if (K < 1)     Rcpp::stop("sparse_grid_gh: K must be >= 1.");
  if (level < 1) Rcpp::stop("sparse_grid_gh: level must be >= 1.");

  /* K = 1 short-circuit: the only contributing multi-index is (level),
     so the grid is just the level-L 1D rule. */
  if (K == 1) {
    std::vector<double> n1, w1;
    gh_1d(2 * level - 1, n1, w1);
    int B = static_cast<int>(n1.size());
    NumericMatrix nodes(B, 1);
    NumericVector weights(B);
    for (int b = 0; b < B; ++b) {
      nodes(b, 0) = n1[b];
      weights[b]  = w1[b];
    }
    return List::create(_["nodes"]   = nodes,
                        _["weights"] = weights,
                        _["K"]       = K,
                        _["level"]   = level);
  }

  /* Precompute 1D rules at all levels 1..level. */
  std::vector<std::vector<double>> rule_nodes(level + 1);
  std::vector<std::vector<double>> rule_weights(level + 1);
  for (int l = 1; l <= level; ++l) {
    gh_1d(2 * l - 1, rule_nodes[l], rule_weights[l]);
  }

  /* Dedup map: key is the rounded z-vector (to suppress floating-point
     equality noise from Golub-Welsch; 1D nodes from dstev are reproducible
     to ~1e-15, so any rounding tolerance >= 1e-12 works). */
  auto key_of = [](const std::vector<double>& z) {
    std::vector<long long> key(z.size());
    for (std::size_t i = 0; i < z.size(); ++i) {
      key[i] = static_cast<long long>(std::llround(z[i] * 1e12));
    }
    return key;
  };
  std::map<std::vector<long long>, std::pair<std::vector<double>, double>> grid;

  /* Valid multi-index sum range: max(K, L-K+1) <= |i| <= L, since each
     i_j >= 1 forces |i| >= K. If level < K the grid is empty. */
  int lo = std::max(K, level - K + 1);
  int hi = level;
  if (lo > hi) {
    Rcpp::stop("sparse_grid_gh: level (%d) is below dimension K (%d); the "
               "minimum useful level is K (single-node grid). Increase `level`.",
               level, K);
  }

  std::vector<int> cur(K, 0);
  std::vector<std::vector<int>> midx_q;
  for (int q = lo; q <= hi; ++q) {
    midx_q.clear();
    enumerate_compositions(K, q, cur, 0, q, midx_q);
    double c = smolyak_coef(level, K, q);
    if (c == 0.0) continue;

    for (const auto& i : midx_q) {
      /* Tensor-product cartesian iteration over 1D nodes at levels i_1..i_K. */
      std::vector<int> sizes(K);
      for (int d = 0; d < K; ++d) sizes[d] = 2 * i[d] - 1;
      std::vector<int> idx(K, 0);

      while (true) {
        std::vector<double> z(K);
        double w = c;
        for (int d = 0; d < K; ++d) {
          z[d] = rule_nodes[i[d]][idx[d]];
          w   *= rule_weights[i[d]][idx[d]];
        }
        auto key = key_of(z);
        auto it  = grid.find(key);
        if (it == grid.end()) grid[key] = std::make_pair(z, w);
        else                  it->second.second += w;

        /* Increment the K-D index. */
        int d = K - 1;
        while (d >= 0) {
          ++idx[d];
          if (idx[d] < sizes[d]) break;
          idx[d] = 0;
          --d;
        }
        if (d < 0) break;
      }
    }
  }

  /* Materialize: drop zero-weight cancellation residuals (rare but possible
     at small floating-point noise levels). */
  int B = 0;
  for (auto& kv : grid) if (kv.second.second != 0.0) ++B;
  NumericMatrix nodes(B, K);
  NumericVector weights(B);
  int b = 0;
  for (auto& kv : grid) {
    if (kv.second.second == 0.0) continue;
    for (int d = 0; d < K; ++d) nodes(b, d) = kv.second.first[d];
    weights[b] = kv.second.second;
    ++b;
  }

  return List::create(_["nodes"]   = nodes,
                      _["weights"] = weights,
                      _["K"]       = K,
                      _["level"]   = level);
}
