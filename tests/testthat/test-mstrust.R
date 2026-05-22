# ============================================================================
# mstrust() (multi-start trust), profile() (profile likelihood),
# vcov() (variance-covariance), and confint.parframe() (profile CIs).
#
# All four interact closely. We use convex-quadratic / Gaussian objectives
# (constraintL2) so the expected behaviour is closed-form. Real-fit
# scenarios are covered by the FOCEI / ECM suite in test-focei.R and
# test-nlmefit.R.
# ============================================================================


# ---- mstrust ------------------------------------------------------------

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


# ---- profile ------------------------------------------------------------

test_that("profile on a 1D quadratic increases monotonically on both sides", {
  obj <- constraintL2(mu = c(theta = 0.0, nuisance = 0.0), sigma = 1)

  prof <- profile(obj = obj, pars = c(theta = 0, nuisance = 0),
                  whichPar = "theta",
                  limits = c(lower = -2, upper = 2),
                  method = "integrate",
                  verbose = FALSE, cores = 1)
  pf <- as.data.frame(prof)
  centered <- pf[order(pf$theta), ]
  left  <- centered[centered$theta <= 0, ]
  right <- centered[centered$theta >= 0, ]
  expect_true(all(diff(left$value)  <= 1e-6))
  expect_true(all(diff(right$value) >= -1e-6))
})


test_that("profile on a 1D Gaussian crosses chi^2 = 3.84 at +/- z(0.95)*sigma", {
  # constraintL2(mu = 0, sigma = sigma) value at theta is (theta/sigma)^2.
  # 95% chi-square threshold qchisq(0.95, 1) ~ 3.84 crosses at
  # theta = +/- sqrt(3.84) * sigma = +/- 1.96 * sigma.
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


# ---- vcov ---------------------------------------------------------------

test_that("vcov(fit) equals solve(0.5 * H) for a quadratic with known Hessian", {
  # constraintL2 with sigma = 1: obj(p) = sum((p - mu)^2), Hessian = 2 I.
  # After trust converges to p = mu, vcov = (0.5 * 2 I)^-1 = I.
  mu <- c(a = 0.5, b = -0.3)
  obj <- constraintL2(mu = mu, sigma = 1)
  fit <- trust(obj, parinit = c(a = 0, b = 0),
               rinit = 1, rmax = 10, iterlim = 50, printIter = FALSE)
  V <- vcov(fit)
  expect_equal(unname(V), diag(2), tolerance = 1e-6)

  # sigma != 1: H = 2 / sigma^2, vcov = sigma^2 * I.
  mu2 <- c(a = 0.0, b = 0.0)
  obj2 <- constraintL2(mu = mu2, sigma = 0.5)
  fit2 <- trust(obj2, parinit = c(a = 0.2, b = -0.2),
                rinit = 1, rmax = 10, iterlim = 50, printIter = FALSE)
  V2 <- vcov(fit2)
  expect_equal(unname(V2), 0.5^2 * diag(2), tolerance = 1e-6)
})


# ---- confint ------------------------------------------------------------

test_that("confint.parframe yields half-width = z(0.95)*sigma on a 1D Gaussian profile", {
  testthat::skip_on_cran()
  # constraintL2's value = sum((p - mu)^2 / sigma^2); for one parameter and
  # sigma = 0.5, value(theta) = (theta / 0.5)^2 = 4 theta^2. The 95%
  # chi-square threshold delta = qchisq(0.95, 1) ~ 3.841 is crossed at
  # |theta| = sqrt(3.841)/2 ~ 0.98 ~ 1.96 * 0.5.
  sigma <- 0.5
  obj <- constraintL2(mu = c(theta = 0, nuisance = 0), sigma = sigma)
  prof <- profile(obj, pars = c(theta = 0, nuisance = 0),
                  whichPar = "theta",
                  limits = c(lower = -3, upper = 3),
                  method = "integrate", verbose = FALSE, cores = 1)
  # confint.parframe defaults val.column = "data"; our profile parframe
  # uses "value", so pass it explicitly.
  ci <- confint(prof, level = 0.95, val.column = "value")
  half_width <- (ci$upper - ci$lower) / 2
  expected <- qnorm(0.975) * sigma
  expect_lt(abs(half_width - expected) / expected, 0.10)
})
