# Behavioral tests for Xs() (ODE-solver-backed prediction function).
#
# Verifies state trajectories and first-order sensitivities against
# closed-form analytical solutions of simple ODE systems:
#   * linear decay      A -> 0        rate k * A
#   * linear cascade    A -> B -> C   rates k1 * A, k2 * B
#
# Hessian and second-order chain-rule semantics are covered by
# test-deriv2-Xs.R.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


## ---- State trajectories: closed-form -----------------------------------

test_that("Xs on linear decay matches A0 * exp(-k * t)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  # bench$xfn alone is condition-agnostic. To get a prdlist we route the
  # ODE solver via the identity trafo (yields a single-condition prdlist).
  prd_states <- bench$xfn * bench$pfn_id

  # Include t = 0 explicitly; Xs prepends it to the time grid otherwise.
  times <- c(0, 0.5, 1, 2, 5, 10)
  pars  <- c(A = 1.7, k = 0.42)
  out <- prd_states(times = times, pars = pars, deriv = FALSE)

  closed <- pars[["A"]] * exp(-pars[["k"]] * times)
  expect_equal(out$C1[, "A"], closed, tolerance = 1e-5)
})


test_that("Xs on two-step cascade matches the closed-form A(t), B(t), C(t)", {
  skip_if_no_compile()
  testthat::skip_if_not_installed("CppODE")
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  reactions <- eqnlist() |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "C", "k2 * B")
  m  <- odemodel(reactions, modelname = "test_xs_cascade", compile = TRUE)
  xf <- Xs(m)
  pf <- P(eqnvec(A = "A", B = "0", C = "0", k1 = "k1", k2 = "k2"),
          condition = "C1",
          modelname = "test_xs_cascade_p", compile = TRUE)
  prd <- xf * pf

  times <- c(0, 0.5, 1, 2, 4)
  pars <- c(A = 1.0, k1 = 0.7, k2 = 0.3)
  out <- prd(times = times, pars = pars, deriv = FALSE)

  closed <- truth_two_step(times, pars["A"], pars["k1"], pars["k2"])
  expect_equal(out$C1[, "A"], closed$A, tolerance = 1e-5)
  expect_equal(out$C1[, "B"], closed$B, tolerance = 1e-5)
  expect_equal(out$C1[, "C"], closed$C, tolerance = 1e-5)
})


## ---- First-order sensitivities -----------------------------------------

test_that("Xs sensitivities on linear decay match analytical d/d(A0,k)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  prd_states <- bench$xfn * bench$pfn_id

  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.0, k = 0.5)
  out <- prd_states(times = times, pars = pars, deriv = TRUE)
  d <- attr(out$C1, "deriv")  # shape [time, var, par]

  # dA/dA0 = exp(-k * t),   dA/dk = -t * A0 * exp(-k * t)
  expect_equal(d[, "A", "A"], exp(-pars[["k"]] * times),     tolerance = 1e-5)
  expect_equal(d[, "A", "k"], -times * pars[["A"]] * exp(-pars[["k"]] * times),
               tolerance = 1e-5)
})


test_that("Xs predictions across conditions are independent and parameter-local", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  .dmod_with_fx_workdir({
    pfn_C2 <- P(eqnvec(A = "A_C2", k = "k_C2"), condition = "C2",
                modelname = "test_xs_p_C2", compile = TRUE)
  })
  prd_multi <- bench$xfn * (bench$pfn_id + pfn_C2)
  times <- c(0, 1, 2, 5)
  pars <- c(A = 1.0, k = 0.5, A_C2 = 2.0, k_C2 = 1.0)
  out <- prd_multi(times = times, pars = pars, deriv = TRUE)

  expect_equal(out$C1[, "A"], 1.0 * exp(-0.5 * times), tolerance = 1e-5)
  expect_equal(out$C2[, "A"], 2.0 * exp(-1.0 * times), tolerance = 1e-5)

  # The deriv tensor is parameter-local per condition: C1's deriv lists
  # only (A, k), C2's only (A_C2, k_C2). The composition + chain rule
  # propagates this sparsity (no cross-condition coupling).
  dC1 <- attr(out$C1, "deriv")
  dC2 <- attr(out$C2, "deriv")
  expect_setequal(dimnames(dC1)[[3]], c("A", "k"))
  expect_setequal(dimnames(dC2)[[3]], c("A_C2", "k_C2"))
})


## ---- Sensitivities at non-integer times: analytic closed form ---------

test_that("Xs sensitivities match the analytical decay formulas at non-integer times", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  prd_states <- bench$xfn * bench$pfn_id
  times <- c(0, 0.7, 1.4, 3.5)
  pars  <- c(A = 1.3, k = 0.65)

  out <- prd_states(times = times, pars = pars, deriv = TRUE)
  d <- attr(out$C1, "deriv")  # [time, var, par], rows aligned to `times`

  # Analytical: dA/dA0 = exp(-k*t), dA/dk = -t * A0 * exp(-k*t).
  expect_equal(d[, "A", "A"], exp(-pars[["k"]] * times), tolerance = 1e-5)
  expect_equal(d[, "A", "k"],
               -times * pars[["A"]] * exp(-pars[["k"]] * times),
               tolerance = 1e-5)
})
