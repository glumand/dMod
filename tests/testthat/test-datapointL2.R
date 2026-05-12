# Behavioral tests for datapointL2() (validation-point L2 constraint).
#
# datapointL2 penalises the deviation between a predicted observable at a
# given time and a target parameter value:
#
#   val = ((prediction[name, t] - target_par) / sigma)^2
#
# It reads prediction from the shared `env` set by an upstream objective
# (typically normL2). Tests use the normL2 + datapointL2 composition so
# env$prediction is populated naturally.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


## ---- Closed-form value -------------------------------------------------

test_that("datapointL2 value equals ((pred - target) / sigma)^2 at a known point", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  obj_main <- normL2(data, bench$prd_id)
  obj_val  <- datapointL2(name = "y", time = 2.0, value = "newpoint",
                          sigma = 0.05, condition = "C1")
  obj <- obj_main + obj_val

  # Target value `newpoint` deliberately != model prediction at t = 2 so the
  # constraint has a non-zero contribution.
  pars <- c(bench$outerpars_id, newpoint = 0.5)
  o    <- obj(pars)
  o_main <- obj_main(bench$outerpars_id)

  pred_at_t <- unname(bench$prd_id(times = c(0, 2), pars = bench$outerpars_id,
                                   deriv = FALSE)$C1[2, "y"])
  expected_val <- ((pred_at_t - pars[["newpoint"]]) / 0.05)^2
  # Tolerance reflects the ODE integrator's rtol; the prediction at t=2 in
  # the obj_val call comes from the (denser) timesD grid used inside normL2,
  # not a clean 2-point integration.
  expect_equal(unname(o$value - o_main$value), expected_val, tolerance = 1e-3)
})


## ---- Additive composition with normL2 ----------------------------------

test_that("normL2 + datapointL2 value equals the sum of the parts", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  obj_main <- normL2(data, bench$prd_id)
  obj_val  <- datapointL2(name = "y", time = 5.0, value = "vtarget",
                          sigma = 0.1, condition = "C1")
  obj <- obj_main + obj_val
  pars <- c(bench$outerpars_id, vtarget = 0.08)

  # +.objfn evaluates main first, then datapointL2 reuses env$prediction.
  # We reproduce the sum by hand: main on its own pars, datapointL2 on a
  # tail call with env from main.
  v_main <- obj_main(bench$outerpars_id)
  env_main <- attr(v_main, "env")
  v_pt   <- obj_val(pars = pars, env = env_main)

  v_total <- obj(pars)
  expect_equal(v_total$value, v_main$value + v_pt$value, tolerance = 1e-10)
})


## ---- Gradient: closed form on the validation-point contribution -------

test_that("datapointL2 contribution to the combined gradient follows the closed form", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  obj_main <- normL2(data, bench$prd_id)
  obj_val  <- datapointL2(name = "y", time = 3.0, value = "target",
                          sigma = 0.05, condition = "C1")
  obj <- obj_main + obj_val
  pars <- c(A = 1.0, k = 0.5, target = 0.2)

  # Closed form of the datapointL2 contribution at (A, k, target):
  # val_pt = ((pred(3) - target) / sigma_pt)^2 with pred(t) = A * exp(-k*t).
  # dpt/dA       = 2 (pred - target) * exp(-k*t) / sigma_pt^2
  # dpt/dk       = 2 (pred - target) * (-t * A * exp(-k*t)) / sigma_pt^2
  # dpt/dtarget  = -2 (pred - target) / sigma_pt^2
  t_pt <- 3.0; sigma_pt <- 0.05
  pred_pt <- pars[["A"]] * exp(-pars[["k"]] * t_pt)
  r <- pred_pt - pars[["target"]]
  contrib <- c(A      = 2 * r * exp(-pars[["k"]] * t_pt) / sigma_pt^2,
               k      = 2 * r * (-t_pt * pars[["A"]] *
                                   exp(-pars[["k"]] * t_pt)) / sigma_pt^2,
               target = -2 * r / sigma_pt^2)

  g_combined <- obj(pars)$gradient
  g_main     <- obj_main(pars[c("A", "k")])$gradient

  # The combined gradient at (A, k) equals main + datapoint contribution
  # at A and k. The `target` entry is purely from datapointL2.
  expect_equal(unname(g_combined[["A"]]),
               unname(g_main[["A"]] + contrib[["A"]]), tolerance = 1e-3)
  expect_equal(unname(g_combined[["k"]]),
               unname(g_main[["k"]] + contrib[["k"]]), tolerance = 1e-3)
  expect_equal(unname(g_combined[["target"]]), unname(contrib[["target"]]),
               tolerance = 1e-3)
})
