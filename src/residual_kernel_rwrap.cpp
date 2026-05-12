// Rcpp wrappers around the dMod shared residual kernel.
// Test-facing entry points (not part of the production hot path).

#include "residual_kernel.h"

#include <Rcpp.h>
#include <vector>
#include <stdexcept>
#include <cstring>

using namespace Rcpp;

namespace {

dmod::BloqMode parse_bloq_mode(const std::string& s) {
  if (s == "NONE")    return dmod::BloqMode::NONE;
  if (s == "M1")      return dmod::BloqMode::M1;
  if (s == "M3")      return dmod::BloqMode::M3;
  if (s == "M4NM")    return dmod::BloqMode::M4NM;
  if (s == "M4BEAL")  return dmod::BloqMode::M4BEAL;
  throw std::runtime_error("residual_kernel: unknown bloq_mode '" + s + "'");
}

// Repacks R's column-major [n_obs, n_par] matrix into row-major storage.
std::vector<double> repack_dpred(const NumericMatrix& M, int n_obs, int n_par) {
  if (M.nrow() != n_obs || M.ncol() != n_par) {
    throw std::runtime_error(
        "residual_kernel: dpred shape mismatch (need [n_obs, n_par]).");
  }
  std::vector<double> out((std::size_t) n_obs * (std::size_t) n_par, 0.0);
  for (int i = 0; i < n_obs; ++i) {
    for (int k = 0; k < n_par; ++k) {
      out[(std::size_t) i * n_par + k] = M(i, k);
    }
  }
  return out;
}

// Repacks R's column-major 3D array [n_obs, n_par, n_par] into row-major
// (n_obs blocks of n_par x n_par row-major).
std::vector<double> repack_d2pred(const NumericVector& arr, int n_obs, int n_par) {
  IntegerVector dim = arr.attr("dim");
  if (dim.size() != 3 || dim[0] != n_obs || dim[1] != n_par || dim[2] != n_par) {
    throw std::runtime_error(
        "residual_kernel: d2pred shape mismatch (need [n_obs, n_par, n_par]).");
  }
  const std::size_t N2 = (std::size_t) n_par * (std::size_t) n_par;
  std::vector<double> out((std::size_t) n_obs * N2, 0.0);
  // Source col-major: arr[i + p*n_obs + q*n_obs*n_par]
  // Dest row-major:   out[i*N2 + p*n_par + q]
  for (int i = 0; i < n_obs; ++i) {
    for (int p = 0; p < n_par; ++p) {
      for (int q = 0; q < n_par; ++q) {
        out[(std::size_t) i * N2 + (std::size_t) p * n_par + q] =
            arr[i + p * n_obs + q * n_obs * n_par];
      }
    }
  }
  return out;
}

dmod::AccumOpts parse_opts(const List& opts_list) {
  dmod::AccumOpts opts;
  opts.use_deriv2_exact     = as<bool>(opts_list["use_deriv2_exact"]);
  opts.bloq_mode            = parse_bloq_mode(as<std::string>(opts_list["bloq_mode"]));
  opts.sigma_depends_on_par = as<bool>(opts_list["sigma_depends_on_par"]);
  opts.d2sigma_present      = as<bool>(opts_list["d2sigma_present"]);
  opts.bessel               = as<double>(opts_list["bessel"]);
  opts.aloq_part1           = as<bool>(opts_list["aloq_part1"]);
  opts.aloq_part2           = as<bool>(opts_list["aloq_part2"]);
  opts.aloq_part3           = as<bool>(opts_list["aloq_part3"]);
  opts.bloq_part1           = as<bool>(opts_list["bloq_part1"]);
  opts.bloq_part2           = as<bool>(opts_list["bloq_part2"]);
  opts.bloq_part3           = as<bool>(opts_list["bloq_part3"]);
  return opts;
}

List run_kernel(bool aloq,
                NumericVector pred, NumericMatrix dpred,
                Nullable<NumericVector> d2pred_in,
                NumericVector y_data, NumericVector sigma,
                Nullable<NumericMatrix> dsigma_in,
                Nullable<NumericVector> lloq_in,
                List opts_list) {
  const int n_obs = pred.size();
  const int n_par = dpred.ncol();

  if (y_data.size() != n_obs || sigma.size() != n_obs) {
    throw std::runtime_error(
        "residual_kernel: pred, y_data, sigma must have length n_obs.");
  }

  dmod::AccumOpts opts = parse_opts(opts_list);

  std::vector<double> dpred_rm  = repack_dpred(dpred, n_obs, n_par);
  std::vector<double> d2pred_rm;
  const double* d2pred_ptr = nullptr;
  if (d2pred_in.isNotNull()) {
    NumericVector d2_in(d2pred_in.get());
    d2pred_rm = repack_d2pred(d2_in, n_obs, n_par);
    d2pred_ptr = d2pred_rm.data();
  }

  std::vector<double> dsigma_rm;
  const double* dsigma_ptr = nullptr;
  if (dsigma_in.isNotNull()) {
    NumericMatrix dsigma_mat(dsigma_in.get());
    dsigma_rm = repack_dpred(dsigma_mat, n_obs, n_par);
    dsigma_ptr = dsigma_rm.data();
  }

  std::vector<double> lloq;
  const double* lloq_ptr = nullptr;
  if (lloq_in.isNotNull()) {
    NumericVector lloq_v(lloq_in.get());
    if (lloq_v.size() != n_obs) {
      throw std::runtime_error("residual_kernel: lloq must have length n_obs.");
    }
    lloq.assign(lloq_v.begin(), lloq_v.end());
    lloq_ptr = lloq.data();
  }

  double value = 0.0;
  std::vector<double> grad(n_par, 0.0);
  std::vector<double> hess((std::size_t) n_par * (std::size_t) n_par, 0.0);

  if (aloq) {
    dmod::accumulate_aloq_residual(
        n_obs, n_par,
        pred.begin(), dpred_rm.data(), d2pred_ptr,
        y_data.begin(), sigma.begin(),
        dsigma_ptr, /*d2sigma=*/ nullptr, lloq_ptr,
        opts, value, grad.data(), hess.data());
  } else {
    dmod::accumulate_bloq_residual(
        n_obs, n_par,
        pred.begin(), dpred_rm.data(), d2pred_ptr,
        y_data.begin(), sigma.begin(),
        dsigma_ptr, /*d2sigma=*/ nullptr, lloq_ptr,
        opts, value, grad.data(), hess.data());
  }

  // Return Hessian as a col-major n_par x n_par matrix.
  NumericMatrix H(n_par, n_par);
  std::memcpy(&H(0, 0), hess.data(), sizeof(double) * (std::size_t) n_par * n_par);

  return List::create(
      Named("value")    = value,
      Named("gradient") = NumericVector(grad.begin(), grad.end()),
      Named("hessian")  = H);
}

}  // anonymous namespace


// [[Rcpp::export]]
List residual_kernel_aloq(NumericVector pred, NumericMatrix dpred,
                          Nullable<NumericVector> d2pred,
                          NumericVector y_data, NumericVector sigma,
                          Nullable<NumericMatrix> dsigma,
                          Nullable<NumericVector> lloq,
                          List opts) {
  return run_kernel(/*aloq=*/ true,
                    pred, dpred, d2pred, y_data, sigma,
                    dsigma, lloq, opts);
}

// [[Rcpp::export]]
List residual_kernel_bloq(NumericVector pred, NumericMatrix dpred,
                          Nullable<NumericVector> d2pred,
                          NumericVector y_data, NumericVector sigma,
                          Nullable<NumericMatrix> dsigma,
                          Nullable<NumericVector> lloq,
                          List opts) {
  return run_kernel(/*aloq=*/ false,
                    pred, dpred, d2pred, y_data, sigma,
                    dsigma, lloq, opts);
}
