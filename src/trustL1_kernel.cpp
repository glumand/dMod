// trust-region optimiser with L1 penalty
//
// Minimises  objfun(theta) + sum_i lambda_i * |theta_i - mu_i|
// (or one-sided  lambda_i * max(0, mu_i - theta_i)) over a subset of
// parameters named in mu, on top of the same Moré-Sorensen subproblem
// used by trust_impl. Two L1-specific eingriffe sitzen im inneren Loop:
//
//  1. L1 active set: a penalised coordinate sitting exactly on its kink
//     theta_i = mu_i is dropped from the reduced subproblem whenever
//     |grad_obj_i| <= lambda_i (two-sided) resp. -grad_obj_i <= lambda_i
//     (one-sided). Below that threshold the L1 force dominates and the
//     coordinate stays pinned.
//
//  2. Kink clamping: after the trust step is added to theta, any
//     penalised coordinate that crossed its kink is snapped back to
//     mu_i, so the next iteration can re-examine the active set.

#include <Rcpp.h>
#include "trust_subproblem.h"
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>

using namespace Rcpp;
using dmod::trust_internal::eigen_sym_local;
using dmod::trust_internal::trust_sub;


// [[Rcpp::export]]
List trustL1_impl(Function objfun,
                  NumericVector parinit,
                  NumericVector mu,
                  NumericVector lambda,
                  bool   one_sided,
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
                  Nullable<CharacterVector> traceFile = R_NilValue) {

  const int K = parinit.size();
  if (K == 0) stop("trustL1: parinit must be non-empty");
  if (!parinit.hasAttribute("names"))
    stop("trustL1: parinit must be a named numeric vector");
  CharacterVector parnames = parinit.names();
  for (int i = 0; i < K; ++i)
    if (!std::isfinite(parinit[i])) stop("trustL1: parinit not all finite");

  // ---- Build per-parameter L1 metadata from (mu, lambda) ----
  std::vector<unsigned char> has_l1(K, 0);
  std::vector<double>        mu_full(K, 0.0);
  std::vector<double>        lambda_full(K, 0.0);
  if (mu.size() > 0) {
    if (!mu.hasAttribute("names"))
      stop("trustL1: mu must be a named numeric vector");
    if (lambda.size() != mu.size())
      stop("trustL1: lambda must have the same length as mu");
    CharacterVector mnames = mu.names();
    for (int j = 0; j < mu.size(); ++j) {
      std::string nm = as<std::string>(mnames[j]);
      for (int i = 0; i < K; ++i) {
        if (as<std::string>(parnames[i]) == nm) {
          has_l1[i] = 1;
          mu_full[i] = mu[j];
          lambda_full[i] = lambda[j];
          break;
        }
      }
    }
  }

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
    if (v.size() != K) stop("trustL1: parscale and parinit not same length");
    for (int i = 0; i < K; ++i) {
      if (!(v[i] > 0))             stop("trustL1: parscale not all positive");
      if (!std::isfinite(v[i]) || !std::isfinite(1.0 / v[i]))
        stop("trustL1: parscale or 1/parscale not all finite");
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

  std::vector<unsigned char> at_upper(K, 0), at_lower(K, 0);
  for (int i = 0; i < K; ++i) {
    at_upper[i] = (theta[i] >= pu[i]) ? 1 : 0;
    at_lower[i] = (theta[i] <= pl[i]) ? 1 : 0;
  }

  // ---- L1 helpers ----
  auto l1_value = [&](const std::vector<double>& th) {
    double s = 0.0;
    for (int i = 0; i < K; ++i) {
      if (!has_l1[i]) continue;
      double d = th[i] - mu_full[i];
      if (one_sided) {
        if (d < 0.0) s += lambda_full[i] * (-d);
      } else {
        s += lambda_full[i] * std::fabs(d);
      }
    }
    return s;
  };
  // Subgradient at theta, used for the smooth part of the reduced
  // subproblem (away from the kink it is the unique gradient; at the
  // kink the coordinate is pinned by the active set and dropped).
  auto l1_grad = [&](int i, const std::vector<double>& th) -> double {
    if (!has_l1[i]) return 0.0;
    double d = th[i] - mu_full[i];
    if (one_sided) {
      if (d < 0.0) return -lambda_full[i];
      return 0.0;
    }
    if (d > 0.0) return  lambda_full[i];
    if (d < 0.0) return -lambda_full[i];
    return 0.0;
  };
  auto l1_pinned = [&](int i, double th_i, double grad_obj_i) {
    if (!has_l1[i] || th_i != mu_full[i]) return false;
    if (one_sided) return (-grad_obj_i) <= lambda_full[i];
    return std::fabs(grad_obj_i) <= lambda_full[i];
  };

  // ---- Initial evaluation (objfun only; L1 part added in C++) ----
  NumericVector x_named(theta.begin(), theta.end());
  x_named.names() = parnames;
  List out_init = as<List>(objfun(x_named));
  double val_obj = as<double>(out_init["value"]);
  NumericVector grad0   = as<NumericVector>(out_init["gradient"]);
  NumericMatrix Hmat0   = as<NumericMatrix>(out_init["hessian"]);
  if (!std::isfinite(val_obj)) stop("parinit not feasible: value is not finite");

  std::vector<double> grad_obj(grad0.begin(), grad0.end());
  std::vector<double> H_full((std::size_t) K * K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      H_full[i + (std::size_t) j * K] = Hmat0(i, j);

  double val = val_obj + l1_value(theta);

  int neval = 1;

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

  // ---- Working reduced subproblem state ----
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
        double gi = grad_obj[i];
        bool drop_bound;
        if (minimize) drop_bound = (at_upper[i] && gi < 0.0) || (at_lower[i] && gi > 0.0);
        else          drop_bound = (at_upper[i] && gi > 0.0) || (at_lower[i] && gi < 0.0);
        bool drop_l1 = l1_pinned(i, theta[i], gi);
        if (!drop_bound && !drop_l1) active.push_back(i);
      }
      int Kred = static_cast<int>(active.size());
      g_red.assign(Kred, 0.0);
      H_red.assign((std::size_t) Kred * Kred, 0.0);
      for (int ii = 0; ii < Kred; ++ii) {
        int i = active[ii];
        // Combined gradient: smooth objective + L1 subgradient at theta.
        g_red[ii] = grad_obj[i] + l1_grad(i, theta);
        for (int jj = 0; jj < Kred; ++jj) {
          int j = active[jj];
          // L1 Hessian contribution is zero — only the smooth part.
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
    double m_value = 0.0;
    double pred_pos = 0.0;
    if (Kred > 0) {
      trust_sub(Kred, g_red.data(), eigvals_red.data(), eigvecs_red.data(), r,
                p_red.data(), &pred_pos, &is_newton, &is_hard, &is_easy);
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

    // ---- Propose theta_try ----
    for (int i = 0; i < K; ++i) theta_try[i] = theta[i] + p_full[i];

    // ---- L1 kink clamping ----
    for (int i = 0; i < K; ++i) {
      if (!has_l1[i]) continue;
      double dnow = theta[i]     - mu_full[i];
      double dtry = theta_try[i] - mu_full[i];
      if (dnow * dtry < 0.0) theta_try[i] = mu_full[i];
      if (one_sided && theta_try[i] < mu_full[i]) theta_try[i] = mu_full[i];
    }

    // ---- Box-bounds clipping ----
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
    double val_obj_try = std::numeric_limits<double>::infinity();
    double val_try     = std::numeric_limits<double>::infinity();
    NumericVector grad_try;
    NumericMatrix Htry_mat;
    if (eval_ok) {
      val_obj_try = as<double>(out_try["value"]);
      grad_try    = as<NumericVector>(out_try["gradient"]);
      Htry_mat    = as<NumericMatrix>(out_try["hessian"]);
      if (!std::isfinite(val_obj_try)) eval_ok = false;
      else val_try = val_obj_try + l1_value(theta_try);
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
        Rf_warning("trustL1: objfun evaluation failed 3 times in a row");
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
      grad_obj.assign(grad_try.begin(), grad_try.end());
      for (int j = 0; j < K; ++j)
        for (int i = 0; i < K; ++i)
          H_full[i + (std::size_t) j * K] = Htry_mat(i, j);
    }

    if (blather) {
      argtry_flat.insert(argtry_flat.end(), theta_try.begin(), theta_try.end());
      valtry_v.push_back(val_try);
      accept_v.push_back(accept ? 1 : 0);
      preddiff_v.push_back(m_value);
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

  int final_iter = (iter <= iterlim) ? iter : iterlim;
  if (!converged && final_iter == iterlim) {
    Rf_warning("Maximum number of iterations exceeded. Fit is not converged.");
  }

  // ---- Build result ----
  NumericVector arg_out(theta.begin(), theta.end());
  arg_out.names() = parnames;
  // Return the combined (penalised) gradient at the final iterate so
  // the result is self-consistent with `value`.
  NumericVector grad_out(K);
  for (int i = 0; i < K; ++i) grad_out[i] = grad_obj[i] + l1_grad(i, theta);
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
