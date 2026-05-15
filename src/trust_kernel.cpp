// dMod trust-region optimizer (full C++ implementation).
//
// User-facing entry point is `trust()` (exported via Rcpp::export and
// auto-roxygen below). Replaces the previous R-level `trust()` wrapper
// in R/trust.R. Closed-form Moré-Sorensen subproblem on a working
// reduced space that drops coordinates pinned at a bound where the
// gradient pushes further into the bound (active-set treatment matching
// the legacy dMod::trust() semantics).
//
// Features supported in C++:
//   - rinit / rmax adaptive radius
//   - fterm / mterm termination
//   - minimize = TRUE/FALSE
//   - parscale rescaling (g and H divided by parscale and outer(ps, ps))
//   - parupper / parlower box bounds with name-or-scalar broadcast
//   - on_step(rho, accepted, iter, r) R callback per step decision
//   - printIter console logging
//   - traceFile CSV log (iter,value,p1,p2,...)
//   - blather per-iteration trace (argpath, argtry, steptype, accept,
//     r, rho, valpath, valtry, preddiff, stepnorm)
//
// objfun is an R closure: trust() does not accept `...`; bind any
// extra arguments in the closure on the R side.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>

using namespace Rcpp;

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

namespace {

// Symmetric K×K eigen via LAPACK dsyevr. `vals` (K, ascending),
// `vecs` (K*K, column-major). A is read-only (copied internally).
void eigen_sym_local(const double* A, int K, double* vals, double* vecs) {
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
  if (info != 0) throw std::runtime_error("trust: dsyevr workspace query failed");
  lwork  = static_cast<int>(wkopt);
  liwork = iwkopt;
  std::vector<double> work(lwork);
  std::vector<int>    iwork(liwork);
  std::copy(A, A + (std::size_t) K * K, A_copy.begin());
  dsyevr_("V", "A", "U", &K, A_copy.data(), &K,
          NULL, NULL, NULL, NULL, &abstol, &m_out,
          vals, vecs, &K, isuppz.data(),
          work.data(), &lwork, iwork.data(), &liwork, &info);
  if (info != 0) throw std::runtime_error("trust: dsyevr failed");
}

// Moré-Sorensen trust-region subproblem.
// Minimize g^T p + 1/2 p^T H p subject to ||p|| <= r, where H is given
// implicitly by its eigendecomposition (vals ascending, vecs column-major).
// Outputs:
//   p_out          : step (length K)
//   predicted_red  : -m(p) where m(p) = g^T p + 0.5 p^T H p  (positive
//                    for descent in the working frame)
//   is_newton      : true if the unconstrained Newton step lies in TR
//   is_hard        : true if the "hard case" branch was taken
//   is_easy        : true if the "easy" sub-branch of hard case was taken
//                    (only meaningful when is_hard is also true)
void trust_sub(int K, const double* g,
               const double* vals, const double* vecs,
               double r, double* p_out, double* predicted_red,
               bool* is_newton, bool* is_hard, bool* is_easy) {
  *is_newton = false;
  *is_hard   = false;
  *is_easy   = false;

  // q = V^T g  (decompose gradient along eigenbasis)
  std::vector<double> q(K, 0.0);
  for (int j = 0; j < K; ++j) {
    double s = 0.0;
    for (int i = 0; i < K; ++i) s += vecs[i + j * K] * g[i];
    q[j] = s;
  }

  // Try the Newton step if H positive definite.
  double lam_min = vals[0];
  if (lam_min > 1e-12) {
    std::vector<double> y(K);
    double pn2 = 0.0;
    for (int j = 0; j < K; ++j) {
      y[j] = -q[j] / vals[j];
      pn2 += y[j] * y[j];
    }
    if (std::sqrt(pn2) <= r) {
      for (int i = 0; i < K; ++i) {
        double s = 0.0;
        for (int j = 0; j < K; ++j) s += vecs[i + j * K] * y[j];
        p_out[i] = s;
      }
      double pred = 0.0;
      for (int j = 0; j < K; ++j) pred += 0.5 * q[j] * q[j] / vals[j];
      *predicted_red = pred;
      *is_newton = true;
      return;
    }
  }

  // Constrained: find sigma >= max(-lam_min, 0) such that ||p(sigma)|| = r.
  // Diagnose hard-vs-easy: lam_min <= 0 AND the eigenvector of lam_min is
  // (nearly) orthogonal to g (q_min ≈ 0) → hard case.
  // Otherwise the standard "easy" Lagrange-multiplier branch.
  //
  // Index set of eigenvalues at the minimum (within fuzz):
  std::vector<int> imin;
  for (int j = 0; j < K; ++j) {
    if (std::fabs(vals[j] - lam_min) < 1e-12 * (std::fabs(lam_min) + 1.0)) imin.push_back(j);
  }
  double C2 = 0.0;
  for (int j : imin) C2 += q[j] * q[j];

  if (lam_min <= 0.0 && C2 < 1e-24) {
    // True "hard-hard" case: lam_min <= 0 and q_min ≈ 0.
    *is_hard = true;
    *is_easy = false;

    // Set p as the constrained minimizer in the lam > lam_min subspace, then
    // shift by sqrt(r^2 - ||p||^2) along the lam_min eigenvector.
    std::vector<double> w(K, 0.0);
    for (int j = 0; j < K; ++j) {
      double beta = vals[j] - lam_min;
      if (beta > 1e-12) w[j] = -q[j] / beta;
    }
    double pn2 = 0.0;
    for (int j = 0; j < K; ++j) pn2 += w[j] * w[j];
    double extra2 = r * r - pn2;
    if (extra2 > 0.0 && !imin.empty()) {
      w[imin[0]] += std::sqrt(extra2);
    }
    for (int i = 0; i < K; ++i) {
      double s = 0.0;
      for (int j = 0; j < K; ++j) s += vecs[i + j * K] * w[j];
      p_out[i] = s;
    }
    double m_val = 0.0;
    for (int j = 0; j < K; ++j) m_val += q[j] * w[j] + 0.5 * vals[j] * w[j] * w[j];
    *predicted_red = -m_val;
    return;
  }

  // "Hard-easy" branch is the standard bracketed root-find on sigma.
  // We can't easily distinguish hard-easy from easy-easy without a
  // separate analytic check; convention from upstream trust:
  //   - hard case (lam_min <= 0) but C2 > 0 → "hard-easy"
  //   - lam_min > 0 → "easy-easy"
  if (lam_min <= 0.0) { *is_hard = true; *is_easy = true; }
  else                { *is_hard = false; *is_easy = false; }

  double sigma_lo = std::max(-lam_min, 0.0) + 1e-12;
  double sigma_hi = sigma_lo + 1.0;
  for (int it = 0; it < 50; ++it) {
    double pn2 = 0.0;
    for (int j = 0; j < K; ++j) {
      double d = vals[j] + sigma_hi;
      pn2 += (q[j] / d) * (q[j] / d);
    }
    if (std::sqrt(pn2) <= r) break;
    sigma_hi *= 2.0;
  }
  double sigma = 0.0;
  for (int it = 0; it < 60; ++it) {
    sigma = 0.5 * (sigma_lo + sigma_hi);
    double pn2 = 0.0;
    for (int j = 0; j < K; ++j) {
      double d = vals[j] + sigma;
      pn2 += (q[j] / d) * (q[j] / d);
    }
    double pn = std::sqrt(pn2);
    if (std::fabs(pn - r) < 1e-10 * r) break;
    if (pn > r) sigma_lo = sigma; else sigma_hi = sigma;
  }
  std::vector<double> y(K);
  for (int j = 0; j < K; ++j) y[j] = -q[j] / (vals[j] + sigma);
  for (int i = 0; i < K; ++i) {
    double s = 0.0;
    for (int j = 0; j < K; ++j) s += vecs[i + j * K] * y[j];
    p_out[i] = s;
  }
  double m_val = 0.0;
  for (int j = 0; j < K; ++j) m_val += q[j] * y[j] + 0.5 * vals[j] * y[j] * y[j];
  *predicted_red = -m_val;
}

}  // namespace


// Internal kernel. User-facing `trust()` lives in R/trust.R and calls
// this entry point after capturing `...` into a closure around objfun.
// [[Rcpp::export]]
List trust_impl(Function objfun,
           NumericVector parinit,
           double rinit,
           double rmax,
           Nullable<NumericVector> parscale = R_NilValue,
           int    iterlim = 100,
           double fterm   = 1e-6,
           double mterm   = 1e-6,
           bool   minimize = true,
           bool   blather  = false,
           Nullable<NumericVector>  parupper  = R_NilValue,
           Nullable<NumericVector>  parlower  = R_NilValue,
           bool   printIter = false,
           Nullable<CharacterVector> traceFile = R_NilValue,
           Nullable<Function>       on_step    = R_NilValue) {

  const int K = parinit.size();
  if (K == 0) stop("trust: parinit must be non-empty");
  if (!parinit.hasAttribute("names"))
    stop("trust: parinit must be a named numeric vector");
  CharacterVector parnames = parinit.names();
  for (int i = 0; i < K; ++i)
    if (!std::isfinite(parinit[i])) stop("trust: parinit not all finite");

  // ---- Build length-K parlower / parupper ----
  std::vector<double> pl(K, -std::numeric_limits<double>::infinity());
  std::vector<double> pu(K, +std::numeric_limits<double>::infinity());
  auto fill_bound = [&](Nullable<NumericVector> nv, std::vector<double>& out) {
    if (nv.isNull()) return;
    NumericVector v(nv);
    if (v.size() == 0) return;
    if (v.hasAttribute("names")) {
      CharacterVector vn = v.names();
      for (int j = 0; j < v.size(); ++j) {
        std::string nm = as<std::string>(vn[j]);
        for (int i = 0; i < K; ++i) {
          if (as<std::string>(parnames[i]) == nm) { out[i] = v[j]; break; }
        }
      }
    } else {
      std::fill(out.begin(), out.end(), v[0]);
    }
  };
  fill_bound(parupper, pu);
  fill_bound(parlower, pl);

  // ---- Validate parscale ----
  bool rescale = parscale.isNotNull();
  std::vector<double> ps(K, 1.0);
  if (rescale) {
    NumericVector v(parscale);
    if (v.size() != K) stop("trust: parscale and parinit not same length");
    for (int i = 0; i < K; ++i) {
      if (!(v[i] > 0))             stop("trust: parscale not all positive");
      if (!std::isfinite(v[i]) || !std::isfinite(1.0 / v[i]))
        stop("trust: parscale or 1/parscale not all finite");
      ps[i] = v[i];
    }
  }

  // ---- Clip parinit to bounds with warning ----
  std::vector<double> theta(K);
  bool any_above = false, any_below = false;
  for (int i = 0; i < K; ++i) {
    if (parinit[i] > pu[i])      { any_above = true; theta[i] = pu[i]; }
    else if (parinit[i] < pl[i]) { any_below = true; theta[i] = pl[i]; }
    else                          theta[i] = parinit[i];
  }
  if (any_above) Rf_warning("init above range");
  if (any_below) Rf_warning("init below range");

  // ---- Boundary indicator vectors (carry across iterations) ----
  // Initial values reflect parinit after clipping.
  std::vector<unsigned char> at_upper(K, 0), at_lower(K, 0);
  for (int i = 0; i < K; ++i) {
    at_upper[i] = (theta[i] >= pu[i]) ? 1 : 0;
    at_lower[i] = (theta[i] <= pl[i]) ? 1 : 0;
  }

  // ---- Initial evaluation ----
  NumericVector x_named(theta.begin(), theta.end());
  x_named.names() = parnames;
  List out_init = as<List>(objfun(x_named));
  double val = as<double>(out_init["value"]);
  NumericVector grad0   = as<NumericVector>(out_init["gradient"]);
  NumericMatrix Hmat0   = as<NumericMatrix>(out_init["hessian"]);
  if (!std::isfinite(val)) stop("parinit not feasible: value is not finite");

  std::vector<double> grad_full(grad0.begin(), grad0.end());
  std::vector<double> H_full((std::size_t) K * K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      H_full[i + (std::size_t) j * K] = Hmat0(i, j);

  int neval = 1;  // evaluation counter (incl. initial)

  // ---- printIter / traceFile setup ----
  int iter_width = static_cast<int>(std::to_string(iterlim).size());
  std::ofstream trace_ofs;
  std::string trace_path;
  if (traceFile.isNotNull()) {
    CharacterVector tf(traceFile);
    if (tf.size() > 0) trace_path = as<std::string>(tf[0]);
  }
  auto trace_write = [&](int it, double v, const NumericVector& x, bool head) {
    if (trace_path.empty()) return;
    if (head) {
      trace_ofs.open(trace_path.c_str());
      trace_ofs << "Iteration,Obj";
      for (int i = 0; i < K; ++i)
        trace_ofs << "," << as<std::string>(parnames[i]);
      trace_ofs << "\n";
    }
    trace_ofs << std::setprecision(15) << it << "," << v;
    for (int i = 0; i < K; ++i) trace_ofs << "," << x[i];
    trace_ofs << "\n";
  };
  if (printIter) {
    Rcpp::Rcout << "Iteration: " << std::setw(iter_width) << neval
                << "      Objective value: " << val << "\n";
  }
  if (!trace_path.empty()) trace_write(neval, val, x_named, /*head=*/true);

  // ---- Blather buffers ----
  std::vector<double> argpath_flat, argtry_flat;
  std::vector<std::string> steptype_v;
  std::vector<int>    accept_v;
  std::vector<double> r_v, rho_v, valpath_v, valtry_v, preddiff_v, stepnorm_v;

  // ---- Working reduced subproblem state (rebuilt on each accepted iter) ----
  std::vector<int>    active;
  std::vector<double> g_red, H_red, eigvals_red, eigvecs_red;
  double f_used = val;

  bool accept = true;
  double r = rinit;
  bool   converged = false;
  bool   is_terminate = false;
  int    iter = 0;
  int    n_fail = 0;

  std::vector<double> theta_try(K), p_full(K), p_red;

  for (iter = 1; iter <= iterlim; ++iter) {

    if (blather) {
      argpath_flat.insert(argpath_flat.end(), theta.begin(), theta.end());
      r_v.push_back(r);
      valpath_v.push_back(val);
    }

    // ---- Active-set / reduced subproblem build (only when accept) ----
    if (accept) {
      active.clear();
      for (int i = 0; i < K; ++i) {
        double gi = grad_full[i];
        bool drop;
        if (minimize) drop = (at_upper[i] && gi < 0.0) || (at_lower[i] && gi > 0.0);
        else          drop = (at_upper[i] && gi > 0.0) || (at_lower[i] && gi < 0.0);
        if (!drop) active.push_back(i);
      }
      int Kred = static_cast<int>(active.size());
      g_red.assign(Kred, 0.0);
      H_red.assign((std::size_t) Kred * Kred, 0.0);
      for (int ii = 0; ii < Kred; ++ii) {
        int i = active[ii];
        g_red[ii] = grad_full[i];
        for (int jj = 0; jj < Kred; ++jj) {
          int j = active[jj];
          H_red[ii + (std::size_t) jj * Kred] = H_full[i + (std::size_t) j * K];
        }
      }
      if (rescale) {
        for (int ii = 0; ii < Kred; ++ii) g_red[ii] /= ps[active[ii]];
        for (int jj = 0; jj < Kred; ++jj)
          for (int ii = 0; ii < Kred; ++ii)
            H_red[ii + (std::size_t) jj * Kred] /= ps[active[ii]] * ps[active[jj]];
      }
      f_used = val;
      if (!minimize) {
        for (auto& x : g_red) x = -x;
        for (auto& x : H_red) x = -x;
        f_used = -val;
      }
      if (Kred > 0) {
        eigvals_red.assign(Kred, 0.0);
        eigvecs_red.assign((std::size_t) Kred * Kred, 0.0);
        eigen_sym_local(H_red.data(), Kred, eigvals_red.data(), eigvecs_red.data());
      }
    }

    int Kred = static_cast<int>(active.size());

    // ---- Solve subproblem on the reduced space ----
    p_red.assign(Kred, 0.0);
    bool is_newton = false, is_hard = false, is_easy = false;
    double m_value = 0.0;        // g^T p + 0.5 p^T B p   (working frame)
    double pred_pos = 0.0;       // -m_value (positive for descent)
    if (Kred > 0) {
      trust_sub(Kred, g_red.data(), eigvals_red.data(), eigvecs_red.data(), r,
                p_red.data(), &pred_pos, &is_newton, &is_hard, &is_easy);
      // Recompute m_value directly from p_red and B_red for consistency:
      double gp = 0.0;
      for (int ii = 0; ii < Kred; ++ii) gp += g_red[ii] * p_red[ii];
      double pBp = 0.0;
      for (int ii = 0; ii < Kred; ++ii) {
        double s = 0.0;
        for (int jj = 0; jj < Kred; ++jj)
          s += H_red[ii + (std::size_t) jj * Kred] * p_red[jj];
        pBp += p_red[ii] * s;
      }
      m_value  = gp + 0.5 * pBp;
      pred_pos = -m_value;
    }

    // ---- Scatter step back to full K and apply parscale inverse ----
    std::fill(p_full.begin(), p_full.end(), 0.0);
    for (int ii = 0; ii < Kred; ++ii) {
      int i = active[ii];
      double s = p_red[ii];
      if (rescale) s /= ps[i];
      p_full[i] = s;
    }
    double stepnorm = 0.0;
    for (int i = 0; i < K; ++i) stepnorm += p_full[i] * p_full[i];
    stepnorm = std::sqrt(stepnorm);

    // ---- Propose theta_try, refresh at_upper/at_lower, clip ----
    for (int i = 0; i < K; ++i) theta_try[i] = theta[i] + p_full[i];
    for (int i = 0; i < K; ++i) {
      at_upper[i] = !(theta_try[i] < pu[i]) ? 1 : 0;
      at_lower[i] = !(theta_try[i] > pl[i]) ? 1 : 0;
      if (at_upper[i]) theta_try[i] = pu[i];
      if (at_lower[i]) theta_try[i] = pl[i];
    }

    // ---- Evaluate at theta_try ----
    NumericVector x_try(theta_try.begin(), theta_try.end());
    x_try.names() = parnames;
    List out_try;
    bool eval_ok = true;
    try {
      out_try = as<List>(objfun(x_try));
    } catch (...) {
      eval_ok = false;
    }
    double val_try = std::numeric_limits<double>::infinity();
    NumericVector grad_try;
    NumericMatrix Htry_mat;
    if (eval_ok) {
      val_try  = as<double>(out_try["value"]);
      grad_try = as<NumericVector>(out_try["gradient"]);
      Htry_mat = as<NumericMatrix>(out_try["hessian"]);
      if (!std::isfinite(val_try)) eval_ok = false;
    }
    neval++;
    if (printIter) {
      Rcpp::Rcout << "Iteration: " << std::setw(iter_width) << neval
                  << "      Objective value: " << val_try << "\n";
    }
    if (!trace_path.empty()) trace_write(neval, val_try, x_try, /*head=*/false);

    // ---- Acceptance / radius update ----
    double ftry_used = minimize ? val_try : -val_try;
    double rho;
    if (eval_ok && pred_pos > 0.0) {
      double actual_red = f_used - ftry_used;
      rho = actual_red / pred_pos;
    } else {
      rho = -std::numeric_limits<double>::infinity();
    }

    is_terminate = eval_ok &&
                   (std::fabs(ftry_used - f_used) < fterm ||
                    std::fabs(m_value)            < mterm);

    if (!eval_ok) {
      n_fail++;
      accept = false;
      r *= 0.25;
      if (n_fail >= 3) {
        Rf_warning("trust: objfun evaluation failed 3 times in a row");
        break;
      }
    } else {
      n_fail = 0;
      if (is_terminate) {
        accept = (ftry_used < f_used);
      } else if (rho < 0.25) {
        accept = false;
        r *= 0.25;
      } else {
        accept = true;
        if (rho > 0.75 && !is_newton) r = std::min(2.0 * r, rmax);
      }
    }

    // ---- Apply acceptance ----
    if (accept && eval_ok) {
      for (int i = 0; i < K; ++i) theta[i] = theta_try[i];
      val = val_try;
      grad_full.assign(grad_try.begin(), grad_try.end());
      for (int j = 0; j < K; ++j)
        for (int i = 0; i < K; ++i)
          H_full[i + (std::size_t) j * K] = Htry_mat(i, j);
    }

    // ---- on_step callback ----
    if (on_step.isNotNull()) {
      Function cb(on_step);
      cb(Named("rho")      = rho,
         Named("accepted") = accept,
         Named("iter")     = iter,
         Named("r")        = r);
    }

    // ---- Blather record (post-step) ----
    if (blather) {
      argtry_flat.insert(argtry_flat.end(), theta_try.begin(), theta_try.end());
      valtry_v.push_back(val_try);
      accept_v.push_back(accept ? 1 : 0);
      preddiff_v.push_back(m_value);   // sign-flipped below if !minimize
      stepnorm_v.push_back(stepnorm);
      rho_v.push_back(rho);
      std::string type_str;
      if (is_newton)               type_str = "Newton";
      else if (is_hard && is_easy) type_str = "hard-easy";
      else if (is_hard)            type_str = "hard-hard";
      else                          type_str = "easy-easy";
      steptype_v.push_back(type_str);
    }

    if (is_terminate) { converged = true; break; }
  }

  if (trace_ofs.is_open()) trace_ofs.close();

  // Loop counter post-exit: iter is in {break_point, iterlim+1}.
  int final_iter = (iter <= iterlim) ? iter : iterlim;
  if (!converged && final_iter == iterlim) {
    Rf_warning("Maximum number of iterations exceeded. Fit is not converged.");
  }

  // ---- Build result ----
  NumericVector arg_out(theta.begin(), theta.end());
  arg_out.names() = parnames;
  NumericVector grad_out(grad_full.begin(), grad_full.end());
  grad_out.names() = parnames;
  NumericMatrix Hess_out(K, K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      Hess_out(i, j) = H_full[i + (std::size_t) j * K];
  Hess_out.attr("dimnames") = List::create(parnames, parnames);

  List result = List::create(
      Named("argument")   = arg_out,
      Named("value")      = val,
      Named("gradient")   = grad_out,
      Named("hessian")    = Hess_out,
      Named("iterations") = final_iter,
      Named("converged")  = converged);

  if (blather) {
    int n_iters = final_iter;
    NumericMatrix argpath_M(n_iters, K);
    NumericMatrix argtry_M(n_iters, K);
    for (int it = 0; it < n_iters; ++it) {
      for (int c = 0; c < K; ++c) {
        argpath_M(it, c) = argpath_flat[(std::size_t) it * K + c];
        argtry_M (it, c) = argtry_flat [(std::size_t) it * K + c];
      }
    }
    argpath_M.attr("dimnames") = List::create(R_NilValue, parnames);
    argtry_M.attr("dimnames")  = List::create(R_NilValue, parnames);

    CharacterVector steptype_out(steptype_v.size());
    for (std::size_t i = 0; i < steptype_v.size(); ++i)
      steptype_out[i] = steptype_v[i];

    LogicalVector accept_out(accept_v.size());
    for (std::size_t i = 0; i < accept_v.size(); ++i)
      accept_out[i] = (accept_v[i] != 0);

    NumericVector preddiff_out(preddiff_v.size());
    for (std::size_t i = 0; i < preddiff_v.size(); ++i)
      preddiff_out[i] = minimize ? preddiff_v[i] : -preddiff_v[i];

    result["argpath"]  = argpath_M;
    result["argtry"]   = argtry_M;
    result["steptype"] = steptype_out;
    result["accept"]   = accept_out;
    result["r"]        = NumericVector(r_v.begin(), r_v.end());
    result["rho"]      = NumericVector(rho_v.begin(), rho_v.end());
    result["valpath"]  = NumericVector(valpath_v.begin(), valpath_v.end());
    result["valtry"]   = NumericVector(valtry_v.begin(), valtry_v.end());
    result["preddiff"] = preddiff_out;
    result["stepnorm"] = NumericVector(stepnorm_v.begin(), stepnorm_v.end());
  }

  return result;
}
