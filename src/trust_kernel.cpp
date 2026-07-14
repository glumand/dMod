// trust-region optimiser (Moré-Sorensen subproblem on a working
// reduced space with active-set treatment for box bounds).

#include <Rcpp.h>
#include "trust_subproblem.h"
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>

using namespace Rcpp;
using dmod::trust_internal::eigen_sym_local;
using dmod::trust_internal::trust_sub;

namespace {

// objfun (an R objfn closure such as normL2()'s) is free to return its
// gradient/hessian with names in whatever order its own internals settled
// on (e.g. normL2's par_names_global, built from deriv dimnames -- which
// need not match the order the caller's parinit happened to be typed in).
// trust_impl otherwise treats grad_full[i]/H_full[i,j] as belonging to
// parnames[i]/parnames[j] purely by position, so a silent name/position
// mismatch here scrambles the whole trust-region model (wrong descent
// direction, persistent rho mismatch) while still taking real steps in the
// correctly-labelled full parameter vector. Align by name defensively.

std::unordered_map<std::string, int> name_index_map(const CharacterVector& names) {
  std::unordered_map<std::string, int> m;
  m.reserve(names.size());
  for (int i = 0; i < names.size(); ++i) m[as<std::string>(names[i])] = i;
  return m;
}

// Reorder a gradient vector to match target_names order. Falls back to
// positional use only if the vector carries no names at all.
std::vector<double> align_grad(const NumericVector& v, const CharacterVector& target_names) {
  const int K = target_names.size();
  std::vector<double> out(K);
  if (!v.hasAttribute("names")) {
    if (v.size() != K)
      throw std::runtime_error("trust: objfun gradient length does not match parinit and has no names to align by.");
    std::copy(v.begin(), v.end(), out.begin());
    return out;
  }
  CharacterVector vn = v.names();
  auto idx = name_index_map(vn);
  for (int i = 0; i < K; ++i) {
    std::string nm = as<std::string>(target_names[i]);
    auto it = idx.find(nm);
    if (it == idx.end())
      throw std::runtime_error("trust: objfun gradient is missing parameter '" + nm + "'.");
    out[i] = v[it->second];
  }
  return out;
}

// Reorder a Hessian matrix (both rows and columns) to match target_names
// order. Falls back to positional use only if it carries no dimnames.
std::vector<double> align_hess(const NumericMatrix& H, const CharacterVector& target_names) {
  const int K = target_names.size();
  std::vector<double> out((std::size_t) K * K);
  RObject dn = H.attr("dimnames");
  if (dn.isNULL()) {
    if (H.nrow() != K || H.ncol() != K)
      throw std::runtime_error("trust: objfun hessian dim does not match parinit and has no dimnames to align by.");
    for (int j = 0; j < K; ++j)
      for (int i = 0; i < K; ++i)
        out[i + (std::size_t) j * K] = H(i, j);
    return out;
  }
  List dnl(dn);
  CharacterVector rn = dnl[0];
  CharacterVector cn = dnl[1];
  auto ridx = name_index_map(rn);
  auto cidx = name_index_map(cn);
  std::vector<int> ri(K), ci(K);
  for (int i = 0; i < K; ++i) {
    std::string nm = as<std::string>(target_names[i]);
    auto itr = ridx.find(nm);
    auto itc = cidx.find(nm);
    if (itr == ridx.end() || itc == cidx.end())
      throw std::runtime_error("trust: objfun hessian is missing parameter '" + nm + "'.");
    ri[i] = itr->second;
    ci[i] = itc->second;
  }
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      out[i + (std::size_t) j * K] = H(ri[i], ci[j]);
  return out;
}

}  // namespace


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
           Nullable<CharacterVector> traceFile = R_NilValue) {

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

  // Align to parnames order by name -- objfun's own internal parameter
  // order need not match parinit's (see comment on align_grad/align_hess).
  std::vector<double> grad_full = align_grad(grad0, parnames);
  std::vector<double> H_full    = align_hess(Hmat0, parnames);

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
    double val_try = std::numeric_limits<double>::infinity();
    std::vector<double> grad_try_aligned, H_try_aligned;
    try {
      out_try = as<List>(objfun(x_try));
      val_try = as<double>(out_try["value"]);
      if (!std::isfinite(val_try)) throw std::runtime_error("value not finite");
      NumericVector grad_try_raw = as<NumericVector>(out_try["gradient"]);
      NumericMatrix Htry_mat_raw = as<NumericMatrix>(out_try["hessian"]);
      grad_try_aligned = align_grad(grad_try_raw, parnames);
      H_try_aligned    = align_hess(Htry_mat_raw, parnames);
    } catch (...) {
      eval_ok = false;
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
      grad_full = grad_try_aligned;
      H_full    = H_try_aligned;
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
