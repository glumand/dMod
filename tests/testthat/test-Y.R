# Behavioral tests for Y() (observation function).
#
# Verifies:
#   * value: observable g(states) evaluates correctly
#   * composition: (Y * Xs)(...) equals Y applied to Xs output
#   * derivMode: "symbolic" and "dual" backends agree numerically
#   * attach.input: pass-through of inputs alongside outputs
#   * gradient: analytic chain rule on y = A^2 (no numDeriv)
#
# Second-order chain rule is covered by test-deriv2-Y.R.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


## ---- Value: linear observable ------------------------------------------

test_that("Y(y = A) on linear decay matches A(t) directly", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  prd_obs <- bench$gfn * bench$xfn * bench$pfn_id

  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.7, k = 0.42)
  out <- prd_obs(times = times, pars = pars, deriv = FALSE)

  closed <- pars[["A"]] * exp(-pars[["k"]] * times)
  expect_equal(out$C1[, "y"], closed, tolerance = 1e-5)
})


## ---- Value: nonlinear observable ---------------------------------------

test_that("Y(y = A^2) evaluates the closed-form (A(t))^2", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)
  bench <- fx_decay_compiled()

  g_sq <- Y(c(y = "A^2"), f = bench$xfn, condition = NULL,
            attach.input = FALSE,
            modelname = "test_Y_sq", compile = TRUE)
  prd_sq <- g_sq * bench$xfn * bench$pfn_id

  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.4, k = 0.3)
  out <- prd_sq(times = times, pars = pars, deriv = FALSE)

  closed <- (pars[["A"]] * exp(-pars[["k"]] * times))^2
  expect_equal(out$C1[, "y"], closed, tolerance = 1e-5)
})


## ---- derivMode parity --------------------------------------------------

test_that("Y derivMode 'symbolic' and 'dual' agree on a nonlinear observable", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)
  bench <- fx_decay_compiled()

  g_sym <- Y(c(y = "A^2"), f = bench$xfn, condition = NULL,
             attach.input = FALSE, derivMode = "symbolic",
             modelname = "test_Y_dm_sym", compile = TRUE)
  g_dual <- Y(c(y = "A^2"), f = bench$xfn, condition = NULL,
              attach.input = FALSE, derivMode = "dual",
              modelname = "test_Y_dm_dual", compile = TRUE)

  prd_sym  <- g_sym  * bench$xfn * bench$pfn_id
  prd_dual <- g_dual * bench$xfn * bench$pfn_id

  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.4, k = 0.3)
  o_sym  <- prd_sym (times = times, pars = pars, deriv = TRUE)
  o_dual <- prd_dual(times = times, pars = pars, deriv = TRUE)

  expect_equal(o_sym$C1[, "y"], o_dual$C1[, "y"], tolerance = 1e-8)
  expect_equal(attr(o_sym$C1, "deriv"), attr(o_dual$C1, "deriv"),
               tolerance = 1e-8)
})


## ---- attach.input ------------------------------------------------------

test_that("Y with attach.input = TRUE returns inputs and outputs", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)
  bench <- fx_decay_compiled()

  g_with_input <- Y(c(y = "A"), f = bench$xfn, condition = NULL,
                    attach.input = TRUE,
                    modelname = "test_Y_attach", compile = TRUE)
  prd_in <- g_with_input * bench$xfn * bench$pfn_id

  times <- c(0, 1, 2)
  pars  <- c(A = 1.0, k = 0.5)
  out <- prd_in(times = times, pars = pars, deriv = FALSE)

  cn <- colnames(out$C1)
  expect_true("y" %in% cn)
  expect_true("A" %in% cn)
  # When attach.input is TRUE, the observable column should equal the input
  # column for the identity observable.
  expect_equal(out$C1[, "y"], out$C1[, "A"], tolerance = 1e-12)
})


## ---- Gradient: analytic chain rule on y = A^2 --------------------------

test_that("Y gradient on y = A^2 follows the analytic chain rule dy/dtheta = 2 A * dA/dtheta", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)
  bench <- fx_decay_compiled()

  g_sq <- Y(c(y = "A^2"), f = bench$xfn, condition = NULL,
            attach.input = FALSE,
            modelname = "test_Y_grad", compile = TRUE)
  prd_sq <- g_sq * bench$xfn * bench$pfn_id

  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.0, k = 0.5)
  out <- prd_sq(times = times, pars = pars, deriv = TRUE)
  d <- attr(out$C1, "deriv")  # [time, var, par]

  # Closed form: A(t) = A0 * exp(-k * t), y(t) = A(t)^2.
  #   dy/dA0 = 2 * A * (dA/dA0) = 2 * A0 * exp(-k*t) * exp(-k*t) = 2 * A0 * exp(-2 k t)
  #   dy/dk  = 2 * A * (dA/dk)  = 2 * A0 * exp(-k*t) * (-t * A0 * exp(-k*t))
  #                              = -2 * t * A0^2 * exp(-2 k t)
  A0 <- pars[["A"]]; k <- pars[["k"]]
  ref_dA <- 2 * A0 * exp(-2 * k * times)
  ref_dk <- -2 * times * A0^2 * exp(-2 * k * times)
  expect_equal(d[, "y", "A"], ref_dA, tolerance = 1e-5)
  expect_equal(d[, "y", "k"], ref_dk, tolerance = 1e-5)
})


# ============================================================================
# Edge case: Y with pure-numeric observable (no outer parameters)
# ============================================================================

test_that("Y with pure-numeric observable composes with an Xs prediction", {
  withr::local_dir(tempdir())
  f <- as.eqnvec(c(A = "-k*A"))
  m <- odemodel(f, modelname = "noparam_y_ode", compile = TRUE,
                solver = "CppODE")
  x <- Xs(m)

  g <- Y(c(y1 = "1.0"), f = NULL, states = c("A"),
         parameters = character(0),
         derivMode = "symbolic", compile = FALSE,
         modelname = "noparam_y_obs")

  out <- (g * x)(seq(0, 5, length.out = 3), c(A = 1.0, k = 0.1))
  pred <- out[[1]]
  expect_true(all(pred[, "y1"] == 1.0))
  expect_equal(pred[, "time"], c(0, 2.5, 5))
})
