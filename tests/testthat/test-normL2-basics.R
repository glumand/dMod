# Behavioral tests for normL2() against mathematical ground truth.
#
# Every value / gradient claim is checked against a closed-form analytic
# formula (Gaussian log-likelihood, decay sensitivities composed by chain
# rule, etc.). No finite-difference reference is used.
#
# Hessian and second-order chain-rule semantics are covered by
# test-deriv2-normL2.R; this file deliberately does not duplicate that.
#
# Every block runs under both objfn backends (R reference and C++ kernel)
# via for_each_backend(); the info= tag distinguishes failures.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


## ---- Value properties ---------------------------------------------------

test_that("normL2 value equals sum(wr^2) + sum(log(2*pi*sigma^2)) at a known point", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  prd_at <- bench$prd_id(times = data$C1$time, pars = bench$outerpars_id)
  closed <- truth_nll_aloq(prd_at$C1[, "y"], data$C1$value, data$C1$sigma)

  for_each_backend(function(cpp) {
    obj <- normL2(data, bench$prd_id)
    o <- obj(bench$outerpars_id)
    expect_equal(o$value, closed, tolerance = 1e-10,
                 info = paste0("cpp=", cpp))
  })
})


test_that("normL2 value scales as 1/sigma^2 when sigma is rescaled", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  # Identical synthetic data at two sigma scales. The chi-square part scales
  # as 1/c^2; the log-sigma part shifts by 2 * n * log(c).
  c1 <- 0.1; c2 <- 0.2
  data1 <- fx_decay_data(sigma = c1)
  data2 <- data1
  data2$C1$value <- data1$C1$value
  data2$C1$sigma <- c2

  for_each_backend(function(cpp) {
    o1 <- normL2(data1, bench$prd_id)(bench$outerpars_id)
    o2 <- normL2(data2, bench$prd_id)(bench$outerpars_id)
    chi1 <- o1$value - sum(log(2 * pi * data1$C1$sigma^2))
    chi2 <- o2$value - sum(log(2 * pi * data2$C1$sigma^2))
    expect_equal(chi2 * (c2 / c1)^2, chi1, tolerance = 1e-10,
                 info = paste0("cpp=", cpp))
  })
})


## ---- Gradient properties ------------------------------------------------

test_that("normL2 gradient equals 2 * Jt * (pred - y) / sigma^2 (analytic decay sens)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)
  # Perturb away from the (statistical) minimum so the gradient is non-trivial.
  pars <- bench$outerpars_id + c(A = 0.15, k = -0.1)

  # Analytic reference: for y = A0 * exp(-k * t),
  #   dy/dA0 = exp(-k*t),  dy/dk = -t * A0 * exp(-k*t).
  # normL2 has value sum(wr^2) + log_term with wr = (pred - obs) / sigma,
  # so dvalue/dtheta = 2 * sum_i wr_i * d(pred_i)/dtheta / sigma_i.
  t <- data$C1$time
  obs <- data$C1$value
  sigma_vec <- data$C1$sigma
  pred <- pars[["A"]] * exp(-pars[["k"]] * t)
  wr <- (pred - obs) / sigma_vec
  J_A <- exp(-pars[["k"]] * t)
  J_k <- -t * pars[["A"]] * exp(-pars[["k"]] * t)
  g_ref <- c(A = 2 * sum(wr * J_A / sigma_vec),
             k = 2 * sum(wr * J_k / sigma_vec))

  for_each_backend(function(cpp) {
    obj <- normL2(data, bench$prd_id)
    g_ana <- obj(pars)$gradient
    expect_equal(unname(g_ana[names(g_ref)]), unname(g_ref),
                 tolerance = 1e-4, info = paste0("cpp=", cpp))
  })
})


## ---- Sigma source equivalence -------------------------------------------

# Build a "prediction chain with sigma pass-through" plus a matching
# constant-sigma errmodel. The errmodel parameter (sigma_y) must appear in
# the inner-parameter set of the prediction chain so normL2's call into
# errmodel sees it. Mirrors the standard dMod pattern.
.build_errmodel_chain <- function(bench, mn_suffix) {
  .dmod_with_fx_workdir({
    e_const <- Y(c(y = "sigma_y"), f = bench$gfn, attach.input = FALSE,
                 condition = "C1",
                 modelname = paste0("fx_decay_err_", mn_suffix),
                 compile = TRUE)
    pfn_sig <- P(eqnvec(A = "A", k = "k", sigma_y = "sigma_y"),
                 condition = "C1",
                 modelname = paste0("fx_decay_p_sig_", mn_suffix),
                 compile = TRUE)
    prd_sig <- bench$gfn * bench$xfn * pfn_sig
  })
  list(prd = prd_sig, e = e_const)
}

test_that("sigma from data column == sigma from errmodel (constant case, no bessel)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  sigma_const <- 0.1

  data_col <- fx_decay_data(sigma = sigma_const)
  data_em  <- data_col
  data_em$C1$sigma <- NA_real_

  ec <- .build_errmodel_chain(bench, "src")
  pars_em <- c(bench$outerpars_id, sigma_y = sigma_const)

  for_each_backend(function(cpp) {
    o_col <- normL2(data_col, bench$prd_id, use.bessel = FALSE)(bench$outerpars_id)
    o_em  <- normL2(data_em,  ec$prd, errmodel = ec$e,
                    use.bessel = FALSE)(pars_em)
    expect_equal(o_em$value, o_col$value, tolerance = 1e-10,
                 info = paste0("cpp=", cpp))
  })
})


## ---- Bessel correction --------------------------------------------------

test_that("bessel correction inflates the chi-square term by exactly factor^2", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  sigma_const <- 0.1
  data_em <- fx_decay_data(sigma = sigma_const)
  data_em$C1$sigma <- NA_real_

  ec <- .build_errmodel_chain(bench, "bes")
  pars_em <- c(bench$outerpars_id, sigma_y = sigma_const)

  # Reproduce normL2's bessel factor so we can predict the exact inflation.
  n     <- nrow(data_em$C1)
  p_all <- union(getParameters(ec$prd), getParameters(ec$e))
  p_err <- setdiff(cOde::getSymbols(unlist(getEquations(ec$e))),
                   names(unlist(getEquations(ec$e))))
  bessel <- sqrt(n / (n - length(p_all) + length(p_err)))

  for_each_backend(function(cpp) {
    o_no <- normL2(data_em, ec$prd, errmodel = ec$e,
                   use.bessel = FALSE)(pars_em)
    o_bs <- normL2(data_em, ec$prd, errmodel = ec$e,
                   use.bessel = TRUE)(pars_em)
    # log_sigma_term is per-data-point; sigma is constant across n rows.
    log_sigma_term <- n * log(2 * pi * sigma_const^2)
    chi_no <- o_no$value - log_sigma_term
    chi_bs <- o_bs$value - log_sigma_term
    expect_equal(chi_bs, chi_no * bessel^2, tolerance = 1e-9,
                 info = paste0("cpp=", cpp))
  })
})


## ---- Multi-condition aggregation ----------------------------------------

test_that("normL2 sums per-condition contributions across two conditions", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  # Build a second-condition trafo so the prd chain knows about C2 (same
  # identity mapping; conditions only differ in the data).
  .dmod_with_fx_workdir({
    pfn_C2 <- P(eqnvec(A = "A", k = "k"), condition = "C2",
                modelname = "fx_decay_p_id_C2", compile = TRUE)
  })
  prd_multi <- bench$gfn * bench$xfn * (bench$pfn_id + pfn_C2)
  data_multi <- fx_decay_data_multi(
    parslist = list(C1 = c(A = 1.0, k = 0.5), C2 = c(A = 1.0, k = 1.0)),
    sigma = 0.1)

  pars <- c(A = 1.0, k = 0.7)  # away from both truths

  for_each_backend(function(cpp) {
    o_joint <- normL2(data_multi, prd_multi)(pars)
    o_C1 <- normL2(data_multi["C1"], prd_multi)(pars)
    o_C2 <- normL2(data_multi["C2"], prd_multi)(pars)
    expect_equal(o_joint$value, o_C1$value + o_C2$value, tolerance = 1e-10,
                 info = paste0("cpp=", cpp))
  })
})
