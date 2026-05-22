// Shared subproblem helpers for the trust-region optimisers in
// trust_kernel.cpp and trustL1_kernel.cpp.
//
// `trust_sub` is a faithful port of the Moré-Sorensen subproblem solver in
// the historical R `trust()` (R/trust.R). The algorithm parametrises the
// Lagrange multiplier as `beep = sigma + lam_min`, so the smallest shifted
// eigenvalue `beta_j = vals_j - lam_min` is always non-negative and the
// degenerate index set is exactly `{j : beta_j == 0}` rather than a
// tolerance bucket. Branch selection then follows R's classical Moré-Sorensen
// formulation:
//   Newton              all(vals > 0) and ||H^{-1} g|| <= r
//   easy (incl hard-easy)   C2 > 0  OR  C1 > r^2     (sigma found via root)
//   hard-hard           C2 == 0 AND C1 <= r^2        (step lands on the
//                                                     min-eigenspace boundary)
// where C1 = sum_{j !in imin}(q_j / beta_j)^2, C2 = sum_{j in imin} q_j^2.
//
// BLAS: `q = vecs^T g` and the back-projection `p = +/- vecs * y` go through
// dgemv. Eigen decomposition uses LAPACK dsyevr.

#ifndef DMOD_TRUST_SUBPROBLEM_H
#define DMOD_TRUST_SUBPROBLEM_H

#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>

#ifndef FCONE
#define FCONE
#endif

extern "C" {
  void dsyevr_(const char* jobz, const char* range, const char* uplo,
               const int* n, double* a, const int* lda,
               const double* vl, const double* vu,
               const int* il, const int* iu,
               const double* abstol, int* m, double* w, double* z,
               const int* ldz, int* isuppz,
               double* work, const int* lwork,
               int* iwork, const int* liwork, int* info);
}

namespace dmod { namespace trust_internal {

// Symmetric K x K eigen via LAPACK dsyevr. `vals` (K, ascending),
// `vecs` (K*K, column-major). A is read-only (copied internally).
inline void eigen_sym_local(const double* A, int K, double* vals, double* vecs) {
  std::vector<double> A_copy(A, A + (std::size_t) K * K);
  std::vector<int>    isuppz(2 * K);
  int info = 0, m_out = 0;
  double abstol = 0.0;
  double wkopt;
  int    iwkopt;
  int    lwork  = -1;
  int    liwork = -1;
  dsyevr_("V", "A", "U", &K, A_copy.data(), &K,
          NULL, NULL, NULL, NULL, &abstol, &m_out,
          vals, vecs, &K, isuppz.data(),
          &wkopt, &lwork, &iwkopt, &liwork, &info);
  if (info != 0) throw std::runtime_error("dsyevr workspace query failed");
  lwork  = static_cast<int>(wkopt);
  liwork = iwkopt;
  std::vector<double> work(lwork);
  std::vector<int>    iwork(liwork);
  std::copy(A, A + (std::size_t) K * K, A_copy.begin());
  dsyevr_("V", "A", "U", &K, A_copy.data(), &K,
          NULL, NULL, NULL, NULL, &abstol, &m_out,
          vals, vecs, &K, isuppz.data(),
          work.data(), &lwork, iwork.data(), &liwork, &info);
  if (info != 0) throw std::runtime_error("dsyevr failed");
}

// Moré-Sorensen trust-region subproblem.
// Minimize g^T p + 1/2 p^T H p subject to ||p|| <= r, where H is given
// implicitly by its eigendecomposition (vals ascending, vecs column-major).
inline void trust_sub(int K, const double* g,
                      const double* vals, const double* vecs,
                      double r, double* p_out, double* predicted_red,
                      bool* is_newton, bool* is_hard, bool* is_easy) {
  *is_newton = false;
  *is_hard   = false;
  *is_easy   = false;

  const int    i_one  = 1;
  const double d_one  = 1.0;
  const double d_zero = 0.0;
  const double d_mone = -1.0;

  // q = vecs^T * g     (BLAS dgemv, T mode)
  std::vector<double> q(K);
  F77_CALL(dgemv)("T", &K, &K, &d_one, vecs, &K,
                  g, &i_one, &d_zero, q.data(), &i_one FCONE);

  // Newton: take if all eigenvalues positive AND ||H^{-1} g|| <= r.
  bool all_pos = true;
  for (int j = 0; j < K; ++j) {
    if (!(vals[j] > 0.0)) { all_pos = false; break; }
  }
  if (all_pos) {
    std::vector<double> y(K);
    double pn2 = 0.0;
    for (int j = 0; j < K; ++j) {
      y[j] = -q[j] / vals[j];
      pn2 += y[j] * y[j];
    }
    if (std::sqrt(pn2) <= r) {
      F77_CALL(dgemv)("N", &K, &K, &d_one, vecs, &K,
                      y.data(), &i_one, &d_zero, p_out, &i_one FCONE);
      double pred = 0.0;
      for (int j = 0; j < K; ++j) pred += 0.5 * q[j] * q[j] / vals[j];
      *predicted_red = pred;
      *is_newton = true;
      return;
    }
  }

  // Non-Newton path: shift eigenvalues so the smallest becomes 0.
  double lam_min = vals[0];
  for (int j = 1; j < K; ++j) if (vals[j] < lam_min) lam_min = vals[j];

  std::vector<double> beta(K);
  std::vector<unsigned char> imin(K);
  for (int j = 0; j < K; ++j) {
    beta[j] = vals[j] - lam_min;
    imin[j] = (beta[j] == 0.0) ? 1u : 0u;
  }

  // C1: contribution from non-degenerate eigenvectors at beep -> 0
  // C2: gradient mass in the min-eigenspace
  // C3: ||q||^2 = ||g||^2 (orthonormal eigenbasis)
  double C1 = 0.0, C2 = 0.0, C3 = 0.0;
  for (int j = 0; j < K; ++j) {
    double qj2 = q[j] * q[j];
    C3 += qj2;
    if (imin[j]) {
      C2 += qj2;
    } else {
      double t = q[j] / beta[j];
      C1 += t * t;
    }
  }

  std::vector<double> w(K, 0.0);

  if (C2 > 0.0 || C1 > r * r) {
    // Easy / hard-easy: solve for beep on (0, infty) such that ||p|| == r.
    *is_easy = true;
    *is_hard = (C2 == 0.0);

    // fred(beep) = sqrt(1 / sum_j (q_j / (beta_j + beep))^2) - 1/r
    // monotonically increasing in beep on (0, infty).
    auto fred = [&](double beep) -> double {
      if (beep == 0.0) {
        if (C2 > 0.0) return -1.0 / r;
        return std::sqrt(1.0 / C1) - 1.0 / r;
      }
      double s = 0.0;
      for (int j = 0; j < K; ++j) {
        double d = beta[j] + beep;
        double t = q[j] / d;
        s += t * t;
      }
      return std::sqrt(1.0 / s) - 1.0 / r;
    };

    // Bracket [beta_dn, beta_up] from R's trust(): conservative outer bounds
    // derived from C2 and C3 = ||g||^2.
    double beta_dn = std::sqrt(C2) / r;
    double beta_up = std::sqrt(C3) / r;

    double root;
    double f_up = fred(beta_up);
    double f_dn = fred(beta_dn);
    if (f_up <= 0.0) {
      root = beta_up;
    } else if (f_dn >= 0.0) {
      root = beta_dn;
    } else {
      // Bisection. fred is monotone increasing, so this converges robustly
      // to within ~50 iterations even when beta_dn / beta_up span many
      // orders of magnitude.
      double a = beta_dn, b = beta_up;
      for (int it = 0; it < 100; ++it) {
        double mid = 0.5 * (a + b);
        if (!(mid > a && mid < b)) { a = mid; b = mid; break; }
        double fm = fred(mid);
        if (fm > 0.0) b = mid; else a = mid;
        if ((b - a) <= 1e-14 * (std::fabs(b) + std::fabs(a) + 1.0)) break;
      }
      root = 0.5 * (a + b);
    }

    for (int j = 0; j < K; ++j) w[j] = q[j] / (beta[j] + root);

    // p_out = -vecs * w
    F77_CALL(dgemv)("N", &K, &K, &d_mone, vecs, &K,
                    w.data(), &i_one, &d_zero, p_out, &i_one FCONE);

    double m_val = 0.0;
    for (int j = 0; j < K; ++j) {
      double y_j = -w[j];
      m_val += q[j] * y_j + 0.5 * vals[j] * y_j * y_j;
    }
    *predicted_red = -m_val;
    return;
  }

  // Hard-hard: gradient orthogonal to min-eigenspace AND the off-min
  // pseudo-inverse step already fits inside the trust region. Take that
  // step and extend along the min-eigenspace direction to the boundary.
  *is_hard = true;
  *is_easy = false;

  for (int j = 0; j < K; ++j) w[j] = imin[j] ? 0.0 : (q[j] / beta[j]);

  // p_out = -vecs * w
  F77_CALL(dgemv)("N", &K, &K, &d_mone, vecs, &K,
                  w.data(), &i_one, &d_zero, p_out, &i_one FCONE);

  double pn2 = 0.0;
  for (int i = 0; i < K; ++i) pn2 += p_out[i] * p_out[i];
  double utry = std::sqrt(std::max(0.0, r * r - pn2));
  if (utry > 0.0) {
    int jmin = -1;
    for (int j = 0; j < K; ++j) if (imin[j]) { jmin = j; break; }
    if (jmin >= 0) {
      for (int i = 0; i < K; ++i) p_out[i] += utry * vecs[i + (std::size_t) jmin * K];
    }
  }

  double m_val = 0.0;
  for (int j = 0; j < K; ++j) {
    // y_j eigen-frame coords: -w[j] for non-min, +utry on the chosen min idx.
    double y_j = -w[j];
    m_val += q[j] * y_j + 0.5 * vals[j] * y_j * y_j;
  }
  // The added min-eigen contribution to m_val is zero in exact arithmetic
  // (lam_min ~= 0 and q[j] for imin is ~= 0), so we leave it out to avoid
  // amplifying floating-point noise.
  *predicted_red = -m_val;
}

}}  // namespace dmod::trust_internal

#endif  // DMOD_TRUST_SUBPROBLEM_H
