# Behavioral tests for mstrust() (multi-start trust) and profile()
# (parameter profile likelihood).
#
# These tests use synthetic quadratic / Gaussian objectives so the expected
# behaviour is unambiguous. Real-fit scenarios are covered by the FOCEI /
# ECM suite which exercises both functions as part of larger pipelines.


# Both mstrust() and profile() require an objfn (class "objfn") rather than
# a plain R function. We use constraintL2 (a convex quadratic centred at
# `mu`) as our synthetic objective; min value = 0 at p = mu, Hessian = 2 I.


## ---- mstrust returns the global minimum on a convex problem -----------

test_that("mstrust returns the global minimum across random starts on a convex quadratic", {
  target <- c(a = 1.0, b = -0.5)
  obj <- constraintL2(mu = target, sigma = 1)

  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)
  fits <- mstrust(objfun = obj, center = c(a = 0, b = 0),
                  studyname = "test_mstrust_convex",
                  rinit = 1, rmax = 10, iterlim = 100,
                  fits = 5, sd = 1, cores = 1, output = FALSE)

  pf <- as.parframe(fits)
  best <- pf[which.min(pf$value), ]
  expect_lt(abs(best$a - target[["a"]]), 1e-6)
  expect_lt(abs(best$b - target[["b"]]), 1e-6)
  expect_equal(min(pf$value), 0, tolerance = 1e-8)
})


## ---- profile increases monotonically away from the minimum ------------

test_that("profile on a 1D quadratic increases monotonically on both sides", {
  # 2D objective with a single parameter to profile.
  obj <- constraintL2(mu = c(theta = 0.0, nuisance = 0.0), sigma = 1)

  prof <- profile(obj = obj, pars = c(theta = 0, nuisance = 0),
                  whichPar = "theta",
                  limits = c(lower = -2, upper = 2),
                  method = "integrate",
                  verbose = FALSE, cores = 1)
  pf <- as.data.frame(prof)
  # The integrate method writes parframe rows; the profiled coordinate is
  # stored in the `theta` column. `value` should be smallest at theta ~ 0
  # and increase monotonically on either side.
  centered <- pf[order(pf$theta), ]
  left  <- centered[centered$theta <= 0, ]
  right <- centered[centered$theta >= 0, ]
  expect_true(all(diff(left$value)  <= 1e-6))
  expect_true(all(diff(right$value) >= -1e-6))
})


## ---- profile chi^2 crossing at theta = +/- sigma * z(0.95) -----------

test_that("profile on a 1D Gaussian crosses chi^2 = 3.84 at +/- z(0.95)*sigma", {
  # constraintL2(mu = 0, sigma = sigma) value at theta is (theta/sigma)^2.
  # For sigma = 0.4, the 95% chi-square threshold qchisq(0.95, 1) ~ 3.84 is
  # crossed at theta = +/- sqrt(3.84) * sigma = +/- 1.96 * sigma = +/- 0.784.
  sigma <- 0.4
  obj <- constraintL2(mu = c(theta = 0, nuisance = 0), sigma = sigma)
  prof <- profile(obj, pars = c(theta = 0, nuisance = 0),
                  whichPar = "theta",
                  limits = c(lower = -3, upper = 3),
                  method = "integrate", verbose = FALSE, cores = 1)
  ci <- confint(prof, level = 0.95, val.column = "value")
  expected_half_width <- qnorm(0.975) * sigma
  half_width <- (ci$upper - ci$lower) / 2
  expect_lt(abs(half_width - expected_half_width) / expected_half_width,
            0.10)
})
