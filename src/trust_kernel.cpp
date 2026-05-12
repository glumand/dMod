// dMod trust C++ kernel.
// A minimal outer trust-region loop in C++ that calls back into R for the
// objective evaluation. Replaces the R-level trust loop's eigen + uniroot
// + tryCatch tower; uses the same Moré-Sorensen subproblem solver as the
// FOCEI inner trust (see src/focei_kernel.cpp).
//
// Supported features (subset of R `trust()`):
//   - basic outer trust loop (rinit / rmax)
//   - fterm / mterm convergence criteria
//   - minimize == TRUE (the dMod default)
//
// NOT supported here (falls back to R when the caller uses any of these):
//   - parscale rescaling, parupper / parlower bounds
//   - on_step, printIter, traceFile, blather

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>

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
  if (info != 0) throw std::runtime_error("trust_kernel: dsyevr workspace query failed");
  lwork  = static_cast<int>(wkopt);
  liwork = iwkopt;
  std::vector<double> work(lwork);
  std::vector<int>    iwork(liwork);
  std::copy(A, A + (std::size_t) K * K, A_copy.begin());
  dsyevr_("V", "A", "U", &K, A_copy.data(), &K,
          NULL, NULL, NULL, NULL, &abstol, &m_out,
          vals, vecs, &K, isuppz.data(),
          work.data(), &lwork, iwork.data(), &liwork, &info);
  if (info != 0) throw std::runtime_error("trust_kernel: dsyevr failed");
}

// Moré-Sorensen trust-region subproblem (cf. src/focei_kernel.cpp:126).
// Minimize g^T p + 1/2 p^T H p subject to ||p|| <= r.
void trust_sub(int K, const double* g,
               const double* vals, const double* vecs,
               double r, double* p_out, double* predicted_red,
               bool* is_newton) {
  std::vector<double> q(K, 0.0);
  for (int j = 0; j < K; ++j) {
    double s = 0.0;
    for (int i = 0; i < K; ++i) s += vecs[i + j * K] * g[i];
    q[j] = s;
  }
  double lam_min = vals[0];
  *is_newton = false;
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
  double pred = 0.0;
  for (int j = 0; j < K; ++j) {
    double d = vals[j] + sigma;
    pred += q[j] * q[j] / d - 0.5 * vals[j] * q[j] * q[j] / (d * d);
  }
  *predicted_red = pred;
}

}  // namespace


// [[Rcpp::export]]
List trust_kernel(Function objfun,
                  NumericVector parinit,
                  double rinit,
                  double rmax,
                  int    iterlim,
                  double fterm,
                  double mterm) {

  const int K = parinit.size();
  CharacterVector parnames = parinit.names();

  // Evaluate objfun at parinit.
  auto eval = [&](const NumericVector& x) -> List {
    return as<List>(objfun(x));
  };

  List out = eval(parinit);
  double val = as<double>(out["value"]);
  NumericVector grad = as<NumericVector>(out["gradient"]);
  NumericMatrix Hmat = as<NumericMatrix>(out["hessian"]);

  std::vector<double> theta(parinit.begin(), parinit.end());
  std::vector<double> g(grad.begin(), grad.end());
  std::vector<double> H((std::size_t) K * K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      H[i + (std::size_t) j * K] = Hmat(i, j);

  double r = rinit;
  bool   converged = false;
  int    iter = 0;

  std::vector<double> eigvals(K), eigvecs((std::size_t) K * K);
  std::vector<double> p_step(K), theta_try(K);

  for (iter = 1; iter <= iterlim; ++iter) {
    eigen_sym_local(H.data(), K, eigvals.data(), eigvecs.data());

    bool is_newton = false;
    double predicted_red = 0.0;
    trust_sub(K, g.data(), eigvals.data(), eigvecs.data(), r,
              p_step.data(), &predicted_red, &is_newton);

    for (int i = 0; i < K; ++i) theta_try[i] = theta[i] + p_step[i];
    NumericVector x_try(theta_try.begin(), theta_try.end());
    x_try.names() = parnames;

    List out_try;
    bool eval_ok = true;
    try {
      out_try = eval(x_try);
    } catch (...) {
      eval_ok = false;
    }
    double val_try = std::numeric_limits<double>::infinity();
    NumericVector grad_try;
    NumericMatrix Htry;
    if (eval_ok) {
      val_try = as<double>(out_try["value"]);
      grad_try = as<NumericVector>(out_try["gradient"]);
      Htry     = as<NumericMatrix>(out_try["hessian"]);
      if (!std::isfinite(val_try)) eval_ok = false;
    }

    double actual_red = val - val_try;
    double rho = (predicted_red > 0.0 && eval_ok)
                   ? actual_red / predicted_red
                   : -std::numeric_limits<double>::infinity();

    bool is_terminate = eval_ok &&
                        (std::fabs(actual_red)    < fterm ||
                         std::fabs(predicted_red) < mterm);

    bool accept;
    if (is_terminate) {
      accept = (val_try < val);
    } else if (rho < 0.25) {
      accept = false;
      r = r * 0.25;
    } else {
      accept = true;
      if (rho > 0.75 && !is_newton) r = std::min(2.0 * r, rmax);
    }

    if (accept && eval_ok) {
      for (int i = 0; i < K; ++i) theta[i] = theta_try[i];
      val  = val_try;
      g.assign(grad_try.begin(), grad_try.end());
      for (int j = 0; j < K; ++j)
        for (int i = 0; i < K; ++i)
          H[i + (std::size_t) j * K] = Htry(i, j);
    }
    if (is_terminate) { converged = true; break; }
  }

  NumericVector arg_out(theta.begin(), theta.end());
  arg_out.names() = parnames;
  NumericVector grad_out(g.begin(), g.end());
  grad_out.names() = parnames;
  NumericMatrix Hess_out(K, K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      Hess_out(i, j) = H[i + (std::size_t) j * K];
  Hess_out.attr("dimnames") = List::create(parnames, parnames);

  return List::create(
      Named("argument")   = arg_out,
      Named("value")      = val,
      Named("gradient")   = grad_out,
      Named("hessian")    = Hess_out,
      Named("iterations") = iter,
      Named("converged")  = converged,
      Named("r")          = r);
}
