# Behavioral tests for vcov() (variance-covariance matrix from a fit)
# and confint.parframe() (profile-likelihood confidence intervals).
#
# vcov: closed-form on a quadratic objective whose Hessian is known.
# confint: profile a 1D Gaussian and check that the 95 % CI half-width
# equals z(0.95) ~ 1.96 sigma, the analytic Wald interval.


## ---- vcov: closed form on a known Hessian ------------------------------

test_that("vcov(fit) equals solve(0.5 * H) for a quadratic with known Hessian", {
  # constraintL2 with sigma = 1 yields obj(p) = sum((p - mu)^2), Hessian
  # = 2 * I (diagonal of 2). After trust converges to p = mu, vcov should
  # be (0.5 * 2 I)^-1 = I = sigma^2 * I.
  mu <- c(a = 0.5, b = -0.3)
  obj <- constraintL2(mu = mu, sigma = 1)
  fit <- trust(obj, parinit = c(a = 0, b = 0),
               rinit = 1, rmax = 10, iterlim = 50, printIter = FALSE)
  V <- vcov(fit)
  expect_equal(unname(V), diag(2), tolerance = 1e-6)

  # Now use sigma != 1: H = 2 / sigma^2, vcov = sigma^2 * I.
  mu2 <- c(a = 0.0, b = 0.0)
  obj2 <- constraintL2(mu = mu2, sigma = 0.5)
  fit2 <- trust(obj2, parinit = c(a = 0.2, b = -0.2),
                rinit = 1, rmax = 10, iterlim = 50, printIter = FALSE)
  V2 <- vcov(fit2)
  expect_equal(unname(V2), 0.5^2 * diag(2), tolerance = 1e-6)
})


## ---- confint.parframe: 95% CI half-width is z(0.95) * sigma ----------

test_that("confint.parframe yields half-width = z(0.95)*sigma on a 1D Gaussian profile", {
  testthat::skip_on_cran()
  # 1D Gaussian: obj(theta) = ((theta - 0) / sigma)^2 with the convention
  # that 'sigma' here is the standard deviation. constraintL2's value =
  # sum((p - mu)^2 / sigma^2), so for one parameter and sigma = 0.5,
  # value(theta) = (theta / 0.5)^2 = 4 theta^2. The 95 % chi-square
  # threshold is delta = qchisq(0.95, 1) ~= 3.841. Solving 4 theta^2 = 3.841
  # gives |theta| = sqrt(3.841)/2 ~ 0.98 ~ 1.96 * 0.5.
  sigma <- 0.5
  obj <- constraintL2(mu = c(theta = 0, nuisance = 0), sigma = sigma)
  prof <- profile(obj, pars = c(theta = 0, nuisance = 0),
                  whichPar = "theta",
                  limits = c(lower = -3, upper = 3),
                  method = "integrate", verbose = FALSE, cores = 1)
  # NB: confint.parframe defaults val.column = "data"; our profile parframe
  # uses "value", so we pass it explicitly.
  ci <- confint(prof, level = 0.95, val.column = "value")
  half_width <- (ci$upper - ci$lower) / 2
  expected <- qnorm(0.975) * sigma
  expect_lt(abs(half_width - expected) / expected, 0.10)
})
