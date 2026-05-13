# Behavioral tests for normL2() on datasets with BLOQ (below-LLOQ) rows.
#
# normL2() forwards `opt.BLOQ` to both backends (R nll() and the C++
# normL2_kernel). The first three blocks below cover the default M3 path;
# the trailing block exercises all four modes (M1, M3, M4NM, M4BEAL) via
# the normL2() public API on both backends.
#
# Hessian / second-order semantics live in test-deriv2-normL2.R.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


## ---- Closed-form BLOQ contribution (M3) --------------------------------

test_that("normL2 with BLOQ rows adds -2 * sum(log Phi(-wr_bloq)) (M3) over ALOQ value", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  # Build a dataset where late time points fall below LLOQ. With the truth
  # k = 0.5, A(t) = exp(-0.5 t) reaches 0.05 around t = 6, so an LLOQ of 0.1
  # censors the tail.
  data <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                             times = seq(0, 10, by = 1))
  pars <- bench$outerpars_id

  # Hand-compute the closed-form M3 contribution at the truth predictions.
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
    # Closed-form comparison; tolerance reflects ODE integrator rtol.
    expect_equal(o$value, expected, tolerance = 1e-4,
                 info = paste0("cpp=", cpp))
  })
})


## ---- BLOQ rows monotonicity --------------------------------------------

test_that("adding a BLOQ row to data strictly increases the normL2 value", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  # Two datasets sharing the same noise sample: one without LLOQ, one with
  # an LLOQ that censors the late tail.
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


## ---- Gradient on BLOQ dataset: closed form (ALOQ + M3 BLOQ) -----------

test_that("normL2 gradient on a BLOQ dataset equals analytic ALOQ + M3 BLOQ contributions", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  # Evaluate away from the (statistical) minimum.
  pars <- bench$outerpars_id + c(A = 0.2, k = -0.1)

  # Closed-form gradient pieces:
  # Predictions and analytic sensitivities for A(t) = A0 * exp(-k*t).
  t <- data$C1$time; obs <- data$C1$value; sigma_vec <- data$C1$sigma
  lloq_vec <- data$C1$lloq
  pred <- pars[["A"]] * exp(-pars[["k"]] * t)
  J_A  <- exp(-pars[["k"]] * t)
  J_k  <- -t * pars[["A"]] * exp(-pars[["k"]] * t)
  val_post <- pmax(obs, lloq_vec)
  is_bloq  <- val_post <= lloq_vec
  # ALOQ part: 2 * sum_i wr_i * J_i / sigma_i
  wr_aloq <- (pred[!is_bloq] - obs[!is_bloq]) / sigma_vec[!is_bloq]
  s_aloq  <- sigma_vec[!is_bloq]
  grad_aloq_A <- 2 * sum(wr_aloq * J_A[!is_bloq] / s_aloq)
  grad_aloq_k <- 2 * sum(wr_aloq * J_k[!is_bloq] / s_aloq)
  # BLOQ M3 part: d/dtheta [-2 sum log Phi(-wr)] = sum 2 * (phi(wr)/Phi(-wr)) * (J/sigma)
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


## ---- Backend agreement on BLOQ dataset --------------------------------

test_that("R and C++ backends give the same value on a BLOQ dataset", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  pars <- bench$outerpars_id

  # Sanity cross-check: even though every other test already runs both
  # backends, this one explicitly compares R against C++ on a BLOQ dataset.
  with_cpp_backend(FALSE, {
    v_R <- normL2(data, bench$prd_id)(pars)$value
  })
  with_cpp_backend(TRUE, {
    v_C <- normL2(data, bench$prd_id)(pars)$value
  })
  expect_equal(v_C, v_R, tolerance = 1e-9)
})


## ---- opt.BLOQ wiring: M1 / M3 / M4NM / M4BEAL via normL2() -------------

test_that("normL2(opt.BLOQ = ...) selects the BLOQ method on both backends", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  pars  <- bench$outerpars_id

  # Closed-form per-mode reference at these pars.
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
  # M4BEAL adds 2 * sum(log Phi(w0)) on the ALOQ rows.
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

  # Cross-backend agreement on the gradient (the new bloq_mode path has to
  # match R for every method we propagate).
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
