# ============================================================================
# Behavioral tests for normL2() against mathematical ground truth.
#
# Every value / gradient claim is checked against a closed-form analytic
# formula (Gaussian log-likelihood, decay sensitivities composed by chain
# rule, errmodel propagation). No finite-difference reference is used.
#
# Hessian and second-order chain-rule semantics live in test-deriv2.R.
#
# All blocks run under both objfn backends (R reference and C++ kernel)
# via for_each_backend(); the info= tag distinguishes failures.
# ============================================================================

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


# ---- Basics: value / gradient -------------------------------------------

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


test_that("normL2 gradient equals 2 * Jt * (pred - y) / sigma^2 (analytic decay sens)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)
  pars <- bench$outerpars_id + c(A = 0.15, k = -0.1)

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


# ---- Sigma source equivalence -------------------------------------------

# Build a "prediction chain with sigma pass-through" plus a matching
# constant-sigma errmodel. The errmodel parameter (sigma_y) must appear in
# the inner-parameter set of the prediction chain so normL2's call into
# errmodel sees it.
.build_const_errmodel_chain <- function(bench, mn_suffix) {
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

  ec <- .build_const_errmodel_chain(bench, "src")
  pars_em <- c(bench$outerpars_id, sigma_y = sigma_const)

  for_each_backend(function(cpp) {
    o_col <- normL2(data_col, ec$prd, use.bessel = FALSE)(pars_em)
    o_em  <- normL2(data_em,  ec$prd, errmodel = ec$e,
                    use.bessel = FALSE)(pars_em)
    expect_equal(o_em$value, o_col$value, tolerance = 1e-10,
                 info = paste0("cpp=", cpp))
  })
})


# ---- Bessel correction --------------------------------------------------

test_that("bessel correction inflates the chi-square term by exactly factor^2", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  sigma_const <- 0.1
  data_em <- fx_decay_data(sigma = sigma_const)
  data_em$C1$sigma <- NA_real_

  ec <- .build_const_errmodel_chain(bench, "bes")
  pars_em <- c(bench$outerpars_id, sigma_y = sigma_const)

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
    log_sigma_term <- n * log(2 * pi * sigma_const^2)
    chi_no <- o_no$value - log_sigma_term
    chi_bs <- o_bs$value - log_sigma_term
    expect_equal(chi_bs, chi_no * bessel^2, tolerance = 1e-9,
                 info = paste0("cpp=", cpp))
  })
})


# ---- Multi-condition aggregation ----------------------------------------

test_that("normL2 sums per-condition contributions across two conditions", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  .dmod_with_fx_workdir({
    pfn_C2 <- P(eqnvec(A = "A", k = "k"), condition = "C2",
                modelname = "fx_decay_p_id_C2", compile = TRUE)
  })
  prd_multi <- bench$gfn * bench$xfn * (bench$pfn_id + pfn_C2)
  data_multi <- fx_decay_data_multi(
    parslist = list(C1 = c(A = 1.0, k = 0.5), C2 = c(A = 1.0, k = 1.0)),
    sigma = 0.1)

  pars <- c(A = 1.0, k = 0.7)

  for_each_backend(function(cpp) {
    o_joint <- normL2(data_multi, prd_multi)(pars)
    o_C1 <- normL2(data_multi["C1"], prd_multi)(pars)
    o_C2 <- normL2(data_multi["C2"], prd_multi)(pars)
    expect_equal(o_joint$value, o_C1$value + o_C2$value, tolerance = 1e-10,
                 info = paste0("cpp=", cpp))
  })
})


# ---- BLOQ: closed-form value (M3) ---------------------------------------

test_that("normL2 with BLOQ rows adds -2 * sum(log Phi(-wr_bloq)) (M3) over ALOQ value", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  data <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                             times = seq(0, 10, by = 1))
  pars <- bench$outerpars_id

  prd_at <- bench$prd_id(times = data$C1$time, pars = pars,
                         deriv = FALSE)$C1
  pred <- prd_at[, "y"]
  sigma_vec <- data$C1$sigma
  lloq_vec  <- data$C1$lloq
  val_post  <- pmax(data$C1$value, lloq_vec)
  is_bloq   <- val_post <= lloq_vec

  closed_aloq <- truth_nll_aloq(pred[!is_bloq], data$C1$value[!is_bloq],
                                sigma_vec[!is_bloq])
  closed_bloq <- truth_nll_bloq_m3(pred[is_bloq], lloq_vec[is_bloq],
                                   sigma_vec[is_bloq])
  expected <- closed_aloq + closed_bloq

  for_each_backend(function(cpp) {
    obj <- normL2(data, bench$prd_id)
    o   <- obj(pars)
    expect_equal(o$value, expected, tolerance = 1e-4,
                 info = paste0("cpp=", cpp))
  })
})


test_that("adding a BLOQ row to data strictly increases the normL2 value", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  raw_pars <- bench$outerpars_id
  data_no <- fx_decay_data(pars = raw_pars, sigma = 0.05,
                           times = seq(0, 10, by = 1), seed = 7L)
  data_yes <- data_no
  data_yes$C1$lloq <- 0.1

  for_each_backend(function(cpp) {
    o_no  <- normL2(data_no,  bench$prd_id)(raw_pars)
    o_yes <- normL2(data_yes, bench$prd_id)(raw_pars)
    expect_gt(o_yes$value, o_no$value,
              label = paste0("cpp=", cpp, " bloq monotonicity"))
  })
})


test_that("normL2 gradient on a BLOQ dataset equals analytic ALOQ + M3 BLOQ contributions", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  pars <- bench$outerpars_id + c(A = 0.2, k = -0.1)

  t <- data$C1$time; obs <- data$C1$value; sigma_vec <- data$C1$sigma
  lloq_vec <- data$C1$lloq
  pred <- pars[["A"]] * exp(-pars[["k"]] * t)
  J_A  <- exp(-pars[["k"]] * t)
  J_k  <- -t * pars[["A"]] * exp(-pars[["k"]] * t)
  val_post <- pmax(obs, lloq_vec)
  is_bloq  <- val_post <= lloq_vec
  wr_aloq <- (pred[!is_bloq] - obs[!is_bloq]) / sigma_vec[!is_bloq]
  s_aloq  <- sigma_vec[!is_bloq]
  grad_aloq_A <- 2 * sum(wr_aloq * J_A[!is_bloq] / s_aloq)
  grad_aloq_k <- 2 * sum(wr_aloq * J_k[!is_bloq] / s_aloq)
  if (any(is_bloq)) {
    wr_b <- (pred[is_bloq] - lloq_vec[is_bloq]) / sigma_vec[is_bloq]
    G    <- exp(stats::dnorm(-wr_b, log = TRUE) -
                  stats::pnorm(-wr_b, log.p = TRUE))
    s_b  <- sigma_vec[is_bloq]
    grad_bloq_A <- 2 * sum(G * J_A[is_bloq] / s_b)
    grad_bloq_k <- 2 * sum(G * J_k[is_bloq] / s_b)
  } else {
    grad_bloq_A <- 0; grad_bloq_k <- 0
  }
  g_ref <- c(A = grad_aloq_A + grad_bloq_A, k = grad_aloq_k + grad_bloq_k)

  for_each_backend(function(cpp) {
    obj <- normL2(data, bench$prd_id)
    g_ana <- obj(pars)$gradient
    expect_equal(unname(g_ana[names(g_ref)]), unname(g_ref),
                 tolerance = 1e-3, info = paste0("cpp=", cpp))
  })
})


test_that("R and C++ backends give the same value on a BLOQ dataset", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  pars <- bench$outerpars_id

  with_cpp_backend(FALSE, {
    v_R <- normL2(data, bench$prd_id)(pars)$value
  })
  with_cpp_backend(TRUE, {
    v_C <- normL2(data, bench$prd_id)(pars)$value
  })
  expect_equal(v_C, v_R, tolerance = 1e-9)
})


test_that("normL2(opt.BLOQ = ...) selects the BLOQ method on both backends", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  pars  <- bench$outerpars_id

  prd_at    <- bench$prd_id(times = data$C1$time, pars = pars,
                            deriv = FALSE)$C1
  pred      <- prd_at[, "y"]
  sigma_vec <- data$C1$sigma
  lloq_vec  <- data$C1$lloq
  val_post  <- pmax(data$C1$value, lloq_vec)
  is_bloq   <- val_post <= lloq_vec

  aloq_value <- truth_nll_aloq(pred[!is_bloq], data$C1$value[!is_bloq],
                               sigma_vec[!is_bloq])
  bloq_m3    <- truth_nll_bloq_m3(pred[is_bloq], lloq_vec[is_bloq],
                                  sigma_vec[is_bloq])
  bloq_m4    <- truth_nll_bloq_m4(pred[is_bloq], lloq_vec[is_bloq],
                                  sigma_vec[is_bloq])
  w0_aloq <- pred[!is_bloq] / sigma_vec[!is_bloq]
  m4beal_aloq_correction <- 2 * sum(stats::pnorm(w0_aloq, log.p = TRUE))

  expected <- c(
    M1     = aloq_value,
    M3     = aloq_value + bloq_m3,
    M4NM   = aloq_value + bloq_m4,
    M4BEAL = aloq_value + m4beal_aloq_correction + bloq_m4
  )

  for_each_backend(function(cpp) {
    for (mode in names(expected)) {
      obj <- normL2(data, bench$prd_id, opt.BLOQ = mode)
      o   <- obj(pars)
      expect_equal(o$value, expected[[mode]], tolerance = 1e-3,
                   info = paste0("cpp=", cpp, " mode=", mode))
    }
  })

  for (mode in c("M1", "M3", "M4NM", "M4BEAL")) {
    with_cpp_backend(FALSE, {
      g_R <- normL2(data, bench$prd_id, opt.BLOQ = mode)(pars)$gradient
    })
    with_cpp_backend(TRUE, {
      g_C <- normL2(data, bench$prd_id, opt.BLOQ = mode)(pars)$gradient
    })
    expect_equal(g_C, g_R, tolerance = 1e-9,
                 info = paste0("gradient parity, mode=", mode))
  }
})


test_that("normL2 rejects unknown opt.BLOQ values", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  expect_error(normL2(data, bench$prd_id, opt.BLOQ = "M2"),
               "should be one of")
})


# ---- errmodel: proportional sigma ---------------------------------------

# Build a prediction chain plus a proportional-error errmodel:
#   sigma(y) = srel * y    where y = A (the observable).
.build_prop_errmodel_chain <- function(mn_suffix) {
  bench <- fx_decay_compiled()
  .dmod_with_fx_workdir({
    e_prop <- Y(c(y = "srel * y"), f = bench$gfn, attach.input = FALSE,
                condition = "C1",
                modelname = paste0("fx_decay_err_prop_", mn_suffix),
                compile = TRUE)
    pfn_prop <- P(eqnvec(A = "A", k = "k", srel = "srel"),
                  condition = "C1",
                  modelname = paste0("fx_decay_p_prop_", mn_suffix),
                  compile = TRUE)
    prd_prop <- bench$gfn * bench$xfn * pfn_prop
  })
  list(prd = prd_prop, e = e_prop)
}


test_that("normL2 with sigma = srel*y matches the proportional-error log-likelihood", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .build_prop_errmodel_chain("val")

  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data(pars = pars[c("A", "k")], sigma = 0.05)
  data$C1$sigma <- NA_real_

  prd_at <- ec$prd(times = data$C1$time, pars = pars, deriv = FALSE)$C1
  pred <- prd_at[, "y"]
  sigma_pred <- pars[["srel"]] * pred
  expected <- truth_nll_aloq(pred, data$C1$value, sigma_pred)

  for_each_backend(function(cpp) {
    obj <- normL2(data, ec$prd, errmodel = ec$e, use.bessel = FALSE)
    o   <- obj(pars)
    expect_equal(o$value, expected, tolerance = 1e-4,
                 info = paste0("cpp=", cpp))
  })
})


test_that("normL2 gradient with proportional errmodel follows the analytic closed form", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .build_prop_errmodel_chain("grad")

  data <- fx_decay_data(pars = c(A = 1.0, k = 0.5), sigma = 0.05)
  data$C1$sigma <- NA_real_
  pars <- c(A = 1.2, k = 0.45, srel = 0.08)

  t <- data$C1$time; obs <- data$C1$value
  A <- pars[["A"]]; k <- pars[["k"]]; srel <- pars[["srel"]]
  pred <- A * exp(-k * t)
  sigma_vec <- srel * pred
  wr <- (pred - obs) / sigma_vec
  dpred <- cbind(A = exp(-k * t), k = -t * A * exp(-k * t), srel = 0)
  dsigma <- srel * dpred
  dsigma[, "srel"] <- pred
  dlog_sigma <- dsigma / sigma_vec
  dwr <- (dpred * sigma_vec - (pred - obs) * dsigma) / sigma_vec^2
  g_ref <- colSums(2 * wr * dwr) + colSums(2 * dlog_sigma)

  for_each_backend(function(cpp) {
    obj <- normL2(data, ec$prd, errmodel = ec$e, use.bessel = FALSE)
    g_ana <- obj(pars)$gradient
    expect_equal(unname(g_ana[names(g_ref)]), unname(g_ref),
                 tolerance = 1e-3, info = paste0("cpp=", cpp))
  })
})


test_that("rows with explicit sigma keep it; NA rows fall through to errmodel", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .build_prop_errmodel_chain("mix")

  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data(pars = pars[c("A", "k")], sigma = 0.05)
  n <- nrow(data$C1)
  data$C1$sigma <- ifelse(seq_len(n) <= n %/% 2, NA_real_, 0.05)

  for_each_backend(function(cpp) {
    obj <- normL2(data, ec$prd, errmodel = ec$e, use.bessel = FALSE)
    o   <- obj(pars)
    expect_true(is.finite(o$value),
                label = paste0("cpp=", cpp, " mixed sigma finite"))
  })

  data_na <- data; data_na$C1 <- data_na$C1[is.na(data$C1$sigma), ]
  data_ex <- data; data_ex$C1 <- data_ex$C1[!is.na(data$C1$sigma), ]
  # Use a common time grid across all three so the adaptive integrator
  # produces bit-identical predictions at the shared data times.
  all_times <- sort(unique(data$C1$time))
  with_cpp_backend(FALSE, {
    v_full <- normL2(data,    ec$prd, errmodel = ec$e, times = all_times,
                     use.bessel = FALSE)(pars)$value
    v_na   <- normL2(data_na, ec$prd, errmodel = ec$e, times = all_times,
                     use.bessel = FALSE)(pars)$value
    v_ex   <- normL2(data_ex, ec$prd, errmodel = ec$e, times = all_times,
                     use.bessel = FALSE)(pars)$value
  })
  expect_equal(v_full, v_na + v_ex, tolerance = 1e-9)
})


test_that("getParameters(normL2(..., errmodel = ec$e)) includes errmodel pars", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .build_prop_errmodel_chain("pars")
  data <- fx_decay_data()
  data$C1$sigma <- NA_real_
  obj <- normL2(data, ec$prd, errmodel = ec$e)
  expect_true("srel" %in% attr(obj, "parameters"))
})
