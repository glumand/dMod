// dMod normL2 C++ kernel.
//
// Consumes a prdlist (one prdframe per condition) plus the per-condition
// errmodel output, plus a pre-built metadata table mapping data rows to
// prediction/errmodel cells. Iterates conditions (optionally OpenMP-parallel,
// `threads > 1`), partitions ALOQ vs BLOQ rows, calls accumulate_aloq /
// accumulate_bloq from residual_kernel for each partition, and scatters the
// per-condition local-parameter gradient/Hessian into a global gradient/Hessian
// vector indexed by `par_names_global`.
//
// Parallel safety: all R API access (extracting matrices, attributes, ints)
// happens in the serial phase 1 (gather). Phase 2 (accumulate) only touches
// plain std::vector data; each thread has its own grad/hess local buffer
// that is merged in a critical region.

#include "residual_kernel.h"

#include <Rcpp.h>
#include <vector>
#include <string>
#include <cstddef>
#include <cstring>
#include <stdexcept>

#ifdef _OPENMP
#  include <omp.h>
#endif

using namespace Rcpp;

namespace {

// Per-condition packed data, built outside the OpenMP region.
struct CondInputs {
  int n_data        = 0;
  int n_aloq        = 0;
  int n_bloq        = 0;
  int n_par_local   = 0;
  bool has_dsigma   = false;
  bool has_d2pred   = false;
  bool has_d2sigma  = false;
  // local par -> global par index (1-based; 0 means not in global, skipped)
  std::vector<int> par_idx_global;
  // ALOQ-first ordering for all arrays; the first n_aloq entries are ALOQ rows,
  // followed by n_bloq BLOQ rows. n_data == n_aloq + n_bloq.
  std::vector<double> pred;
  std::vector<double> y_data;
  std::vector<double> sigma;
  std::vector<double> lloq;
  std::vector<double> dpred;    // row-major [n_data, n_par_local]
  std::vector<double> dsigma;   // row-major [n_data, n_par_local], or empty
  std::vector<double> d2pred;   // row-major [n_data, n_par_local^2], or empty
  std::vector<double> d2sigma;  // row-major [n_data, n_par_local^2], or empty
};


inline int match_name(const CharacterVector& haystack, const std::string& s) {
  for (int i = 0; i < haystack.size(); ++i) {
    if (as<std::string>(haystack[i]) == s) return i;  // 0-based
  }
  return -1;
}

CondInputs gather_one_condition(
    NumericMatrix prdf,
    Nullable<NumericMatrix> err_mat_opt,
    const List& meta,
    const CharacterVector& par_names_global,
    bool want_d2pred) {
  CondInputs C;
  IntegerVector t_idx       = meta["t_idx_in_pred"];        // 1-based
  IntegerVector o_idx       = meta["o_idx_in_pred"];
  IntegerVector o_idx_d     = meta["o_idx_in_deriv"];
  IntegerVector bloq_mask   = meta["bloq_mask"];             // 0/1 per data row
  IntegerVector t_idx_err   = meta["t_idx_in_err"];
  IntegerVector o_idx_err   = meta["o_idx_in_err"];
  IntegerVector o_idx_err_d = meta["o_idx_in_err_deriv"];
  IntegerVector sigma_is_na = meta["sigma_is_na"];           // 0/1 per data row
  NumericVector y_in        = meta["y_data"];                // already LOQ-substituted
  NumericVector lloq_in     = meta["lloq"];
  NumericVector sigma_fixed = meta["sigma_fixed"];           // per-row, 0 where errmodel supplies

  const int n_data = t_idx.size();
  C.n_data = n_data;

  // Resolve par_local via the prdframe's deriv dimnames.
  RObject deriv_attr_sexp = prdf.attr("deriv");
  if (deriv_attr_sexp.isNULL())
    throw std::runtime_error("normL2_kernel: prdframe has no deriv attribute.");
  NumericVector dpred_flat(deriv_attr_sexp);
  IntegerVector deriv_dim = dpred_flat.attr("dim");
  List deriv_dimnames     = dpred_flat.attr("dimnames");
  CharacterVector par_local_names = deriv_dimnames[2];
  const int Dp0 = deriv_dim[0];
  const int Dp1 = deriv_dim[1];
  const int n_par_local = par_local_names.size();
  C.n_par_local = n_par_local;

  C.par_idx_global.resize(n_par_local);
  for (int p = 0; p < n_par_local; ++p) {
    std::string nm = as<std::string>(par_local_names[p]);
    int g = match_name(par_names_global, nm);
    C.par_idx_global[p] = (g < 0) ? 0 : (g + 1);  // 1-based
  }

  // ALOQ-first permutation: stable sort indices so ALOQ rows come first.
  std::vector<int> perm(n_data);
  int aw = 0, bw = 0;
  for (int i = 0; i < n_data; ++i) {
    if (bloq_mask[i] == 0) ++aw; else ++bw;
  }
  C.n_aloq = aw;
  C.n_bloq = bw;
  int a_pos = 0, b_pos = C.n_aloq;
  for (int i = 0; i < n_data; ++i) {
    perm[bloq_mask[i] == 0 ? a_pos++ : b_pos++] = i;
  }

  // Gather scalar arrays.
  C.pred.resize(n_data);
  C.y_data.resize(n_data);
  C.sigma.resize(n_data);
  C.lloq.resize(n_data);
  for (int j = 0; j < n_data; ++j) {
    const int i = perm[j];
    const int ti = t_idx[i] - 1;
    const int oi = o_idx[i] - 1;
    C.pred[j]   = prdf(ti, oi);
    C.y_data[j] = y_in[i];
    C.lloq[j]   = lloq_in[i];
    if (sigma_is_na[i]) {
      if (err_mat_opt.isNull())
        throw std::runtime_error("normL2_kernel: NA sigma but no errmodel for this condition.");
      NumericMatrix em(err_mat_opt.get());
      const int tei = t_idx_err[i] - 1;
      const int oei = o_idx_err[i] - 1;
      C.sigma[j] = em(tei, oei);
    } else {
      C.sigma[j] = sigma_fixed[i];
    }
  }

  // Gather dpred.
  C.dpred.assign((std::size_t) n_data * n_par_local, 0.0);
  for (int j = 0; j < n_data; ++j) {
    const int i = perm[j];
    const int ti = t_idx[i] - 1;
    const int od = o_idx_d[i] - 1;
    double* row = C.dpred.data() + (std::size_t) j * n_par_local;
    for (int p = 0; p < n_par_local; ++p) {
      row[p] = dpred_flat[ti + od * Dp0 + p * Dp0 * Dp1];
    }
  }

  // Gather d2pred if requested + present.
  if (want_d2pred) {
    RObject d2_sexp = prdf.attr("deriv2");
    if (!d2_sexp.isNULL()) {
      NumericVector d2_attr(d2_sexp);
      IntegerVector d2_dim = d2_attr.attr("dim");
      if (d2_dim.size() == 4) {
        const int D0 = d2_dim[0], D1 = d2_dim[1], D2 = d2_dim[2];
        const std::size_t N2 = (std::size_t) n_par_local * (std::size_t) n_par_local;
        C.d2pred.assign((std::size_t) n_data * N2, 0.0);
        // R d2 attr col-major: [time, obs, par1, par2]
        // We assume dimnames(d2)[[3]] == dimnames(d2)[[4]] == par_local_names.
        for (int j = 0; j < n_data; ++j) {
          const int i = perm[j];
          const int ti = t_idx[i] - 1;
          const int od = o_idx_d[i] - 1;
          double* dst = C.d2pred.data() + (std::size_t) j * N2;
          for (int p = 0; p < n_par_local; ++p) {
            for (int q = 0; q < n_par_local; ++q) {
              dst[(std::size_t) p * n_par_local + q] =
                  d2_attr[ti + od * D0 + p * D0 * D1 + q * D0 * D1 * D2];
            }
          }
        }
        C.has_d2pred = true;
      }
    }
  }

  // Gather dsigma / d2sigma if errmodel has deriv (mapped from err-par names
  // to local pars).
  if (err_mat_opt.isNotNull()) {
    NumericMatrix em(err_mat_opt.get());
    RObject ed_sexp = em.attr("deriv");
    if (!ed_sexp.isNULL()) {
      NumericVector ed_attr(ed_sexp);
      IntegerVector ed_dim = ed_attr.attr("dim");
      const int Ed0 = ed_dim[0];
      const int Ed1 = ed_dim[1];
      List ed_dimnames = ed_attr.attr("dimnames");
      CharacterVector err_par_names = ed_dimnames[2];
      const int n_err_par = err_par_names.size();
      std::vector<int> err_to_local(n_err_par, -1);
      for (int q = 0; q < n_err_par; ++q) {
        std::string nm = as<std::string>(err_par_names[q]);
        for (int p = 0; p < n_par_local; ++p) {
          if (as<std::string>(par_local_names[p]) == nm) {
            err_to_local[q] = p; break;
          }
        }
      }
      C.dsigma.assign((std::size_t) n_data * n_par_local, 0.0);
      bool any = false;
      for (int j = 0; j < n_data; ++j) {
        const int i = perm[j];
        if (!sigma_is_na[i]) continue;  // fixed sigma row: no dsigma
        const int tei = t_idx_err[i] - 1;
        const int oed_one = o_idx_err_d[i];
        if (oed_one <= 0) continue;
        const int oed = oed_one - 1;
        double* row = C.dsigma.data() + (std::size_t) j * n_par_local;
        for (int q = 0; q < n_err_par; ++q) {
          const int lp = err_to_local[q];
          if (lp < 0) continue;
          row[lp] = ed_attr[tei + oed * Ed0 + q * Ed0 * Ed1];
          any = true;
        }
      }
      C.has_dsigma = any;

      // d2sigma: same err-par -> local-par mapping; gather only when caller
      // requested deriv2 (= want_d2pred) and when errmodel exposes deriv2.
      if (want_d2pred) {
        RObject ed2_sexp = em.attr("deriv2");
        if (!ed2_sexp.isNULL()) {
          NumericVector ed2_attr(ed2_sexp);
          IntegerVector ed2_dim = ed2_attr.attr("dim");
          if (ed2_dim.size() == 4) {
            const int E0 = ed2_dim[0], E1 = ed2_dim[1], E2 = ed2_dim[2];
            // Assume dimnames(deriv2)[[3]] == dimnames(deriv2)[[4]] ==
            // err_par_names (the err deriv parameter order). We reuse the
            // err_to_local mapping computed above.
            const std::size_t N2 = (std::size_t) n_par_local
                                   * (std::size_t) n_par_local;
            C.d2sigma.assign((std::size_t) n_data * N2, 0.0);
            bool any2 = false;
            for (int j = 0; j < n_data; ++j) {
              const int i = perm[j];
              if (!sigma_is_na[i]) continue;
              const int tei = t_idx_err[i] - 1;
              const int oed_one = o_idx_err_d[i];
              if (oed_one <= 0) continue;
              const int oed = oed_one - 1;
              double* dst = C.d2sigma.data() + (std::size_t) j * N2;
              for (int p = 0; p < n_err_par; ++p) {
                const int lp = err_to_local[p];
                if (lp < 0) continue;
                for (int q = 0; q < n_err_par; ++q) {
                  const int lq = err_to_local[q];
                  if (lq < 0) continue;
                  const double v = ed2_attr[tei + oed * E0
                                            + p * E0 * E1
                                            + q * E0 * E1 * E2];
                  dst[(std::size_t) lp * n_par_local + lq] = v;
                  if (v != 0.0) any2 = true;
                }
              }
            }
            C.has_d2sigma = any2;
          }
        }
      }
    }
  }

  return C;
}

}  // namespace


namespace {

dmod::BloqMode parse_normL2_bloq_mode(const std::string& s) {
  if (s == "M1")     return dmod::BloqMode::M1;
  if (s == "M3")     return dmod::BloqMode::M3;
  if (s == "M4NM")   return dmod::BloqMode::M4NM;
  if (s == "M4BEAL") return dmod::BloqMode::M4BEAL;
  if (s == "NONE")   return dmod::BloqMode::NONE;
  throw std::runtime_error("normL2_kernel: unknown bloq_mode '" + s +
                           "' (expected one of M1, M3, M4NM, M4BEAL).");
}

}  // namespace


// [[Rcpp::export]]
List normL2_kernel(
    List prediction,
    Nullable<List> err_list_opt,
    List meta_list,
    CharacterVector par_names_global,
    double bessel,
    bool deriv2_requested,
    int threads,
    std::string bloq_mode = "M3") {

  const int n_cond = prediction.size();
  if ((int) meta_list.size() != n_cond)
    throw std::runtime_error("normL2_kernel: meta_list size mismatch.");
  const int n_par_global = par_names_global.size();

  const dmod::BloqMode bmode = parse_normL2_bloq_mode(bloq_mode);

  List err_list;
  if (err_list_opt.isNotNull()) err_list = err_list_opt.get();

  // Phase 1: gather (serial, R API allowed)
  std::vector<CondInputs> conds;
  conds.reserve(n_cond);
  for (int c = 0; c < n_cond; ++c) {
    NumericMatrix prdf = as<NumericMatrix>(prediction[c]);
    Nullable<NumericMatrix> err_opt = R_NilValue;
    if (err_list_opt.isNotNull()) {
      SEXP em_sexp = err_list[c];
      if (em_sexp != R_NilValue) {
        err_opt = Nullable<NumericMatrix>(NumericMatrix(em_sexp));
      }
    }
    conds.push_back(gather_one_condition(
        prdf, err_opt, as<List>(meta_list[c]),
        par_names_global, deriv2_requested));
  }

  // Phase 2: accumulate (parallel-over-conditions)
  std::vector<double> grad_global(n_par_global, 0.0);
  std::vector<double> hess_global((std::size_t) n_par_global * n_par_global, 0.0);
  double value_global = 0.0;

  int n_threads = std::max(1, threads);
#ifndef _OPENMP
  n_threads = 1;
#endif

  dmod::AccumOpts base_opts;
  base_opts.use_deriv2_exact = deriv2_requested;
  base_opts.bloq_mode        = bmode;
  base_opts.bessel           = bessel;

#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads)
#endif
  {
    std::vector<double> grad_local(n_par_global, 0.0);
    std::vector<double> hess_local((std::size_t) n_par_global * n_par_global, 0.0);
    double value_local = 0.0;

#ifdef _OPENMP
    #pragma omp for schedule(static) nowait
#endif
    for (int c = 0; c < n_cond; ++c) {
      const CondInputs& C = conds[c];
      const int npl = C.n_par_local;

      // Per-condition gradient + Hessian (LOCAL parameter space).
      std::vector<double> grad_cond(npl, 0.0);
      std::vector<double> hess_cond((std::size_t) npl * npl, 0.0);
      double value_cond = 0.0;

      dmod::AccumOpts opts = base_opts;
      opts.sigma_depends_on_par = C.has_dsigma;
      opts.d2sigma_present      = C.has_d2sigma;

      // ALOQ rows
      if (C.n_aloq > 0) {
        dmod::accumulate_aloq_residual(
            C.n_aloq, npl,
            C.pred.data(), C.dpred.data(),
            C.has_d2pred ? C.d2pred.data() : nullptr,
            C.y_data.data(), C.sigma.data(),
            C.has_dsigma  ? C.dsigma.data()  : nullptr,
            C.has_d2sigma ? C.d2sigma.data() : nullptr,
            C.lloq.data(),
            opts,
            value_cond, grad_cond.data(), hess_cond.data());
      }
      // BLOQ rows (offset = n_aloq)
      if (C.n_bloq > 0) {
        const std::size_t off_n  = (std::size_t) C.n_aloq;
        const std::size_t off_n2 = off_n * npl;
        const std::size_t off_n3 = off_n * (std::size_t) npl * (std::size_t) npl;
        dmod::accumulate_bloq_residual(
            C.n_bloq, npl,
            C.pred.data()   + off_n,
            C.dpred.data()  + off_n2,
            C.has_d2pred ? (C.d2pred.data() + off_n3) : nullptr,
            C.y_data.data() + off_n,
            C.sigma.data()  + off_n,
            C.has_dsigma  ? (C.dsigma.data()  + off_n2) : nullptr,
            C.has_d2sigma ? (C.d2sigma.data() + off_n3) : nullptr,
            C.lloq.data()   + off_n,
            opts,
            value_cond, grad_cond.data(), hess_cond.data());
      }

      value_local += value_cond;
      // Scatter local -> global gradient/Hessian.
      for (int p = 0; p < npl; ++p) {
        const int g = C.par_idx_global[p] - 1;
        if (g < 0) continue;
        grad_local[g] += grad_cond[p];
      }
      for (int p2 = 0; p2 < npl; ++p2) {
        const int g2 = C.par_idx_global[p2] - 1;
        if (g2 < 0) continue;
        for (int p1 = 0; p1 < npl; ++p1) {
          const int g1 = C.par_idx_global[p1] - 1;
          if (g1 < 0) continue;
          hess_local[(std::size_t) g1 + (std::size_t) g2 * n_par_global] +=
              hess_cond[(std::size_t) p1 + (std::size_t) p2 * npl];
        }
      }
    }

    // Merge thread-local accumulators into global.
#ifdef _OPENMP
    #pragma omp critical
#endif
    {
      value_global += value_local;
      for (int p = 0; p < n_par_global; ++p) grad_global[p] += grad_local[p];
      const std::size_t H = (std::size_t) n_par_global * n_par_global;
      for (std::size_t k = 0; k < H; ++k) hess_global[k] += hess_local[k];
    }
  }

  // Package result as objlist-shaped list.
  NumericVector grad_R(grad_global.begin(), grad_global.end());
  grad_R.names() = par_names_global;
  NumericMatrix hess_R(n_par_global, n_par_global);
  std::memcpy(&hess_R(0, 0), hess_global.data(),
              sizeof(double) * (std::size_t) n_par_global * n_par_global);
  hess_R.attr("dimnames") = List::create(par_names_global, par_names_global);

  return List::create(
      Named("value")    = value_global,
      Named("gradient") = grad_R,
      Named("hessian")  = hess_R);
}
