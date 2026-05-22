# Behavioral tests for trustL1() (trust optimizer with L1 penalty).
#
# trustL1 minimises  objfun(p) + lambda * sum(|p - mu_L1|)
# on top of a smooth objfun. For a 1D quadratic objfun (p - x0)^2 the
# closed-form optimum is the soft-thresholding solution:
#
#   |x0 - mu_L1| <= lambda / 2:  p* = mu_L1
#   else:                        p* = x0 - sign(x0 - mu_L1) * lambda / 2
#
# We verify trustL1 lands on these analytical points across both regimes.


# Quadratic objfn centred at `target` with unit "sigma". Returns an
# objlist so trustL1's `+.objlist` accumulation works correctly.
.quadratic_objfn <- function(target) {
  d <- length(target)
  function(p, ...) {
    gr <- 2 * (p - target); names(gr) <- names(p)
    hs <- 2 * diag(d); dimnames(hs) <- list(names(p), names(p))
    objlist(value = sum((p - target)^2), gradient = gr, hessian = hs)
  }
}


## ---- Soft-threshold above the kink: |x0 - mu| > lambda/2 ---------------

test_that("trustL1 lands on soft(x0, lambda/2) when above the L1 kink", {
  x0 <- 1.5
  lambda <- 0.4
  obj <- .quadratic_objfn(c(x = x0))

  fit <- trustL1(obj, parinit = c(x = 0), mu = c(x = 0),
                 lambda = lambda, rinit = 0.5, rmax = 5, iterlim = 100)
  # |x0 - mu_L1| = 1.5 > lambda/2 = 0.2, so the soft-threshold is active.
  expected <- x0 - sign(x0 - 0) * lambda / 2
  expect_equal(unname(fit$argument[["x"]]), expected, tolerance = 1e-4)
})


## ---- Multivariate: soft-threshold acts elementwise --------------------

test_that("trustL1 with diagonal quadratic applies soft-threshold elementwise", {
  x0 <- c(a = 1.0, b = 0.05, c = -0.8)
  lambda <- 0.3
  obj <- .quadratic_objfn(x0)
  mu_L1 <- c(a = 0, b = 0, c = 0)

  fit <- trustL1(obj, parinit = c(a = 0, b = 0, c = 0),
                 mu = mu_L1, lambda = lambda,
                 rinit = 0.5, rmax = 5, iterlim = 200)

  # Each coordinate: soft(x0_i, lambda/2) w.r.t. mu_i = 0.
  soft <- function(z, t) if (abs(z) <= t) 0 else z - sign(z) * t
  expected <- vapply(x0, soft, t = lambda / 2, FUN.VALUE = 0.0)
  expect_equal(unname(fit$argument[names(x0)]),
               unname(expected), tolerance = 5e-3)
})


## ---- Lambda = 0 reduces to plain trust ---------------------------------

test_that("trustL1 with lambda = 0 converges to the unpenalized minimum", {
  x0 <- c(a = 0.7, b = -1.2)
  obj <- .quadratic_objfn(x0)
  fit <- trustL1(obj, parinit = c(a = 0, b = 0),
                 mu = c(a = 0, b = 0), lambda = 0,
                 rinit = 0.5, rmax = 5, iterlim = 100)
  expect_equal(unname(fit$argument[names(x0)]), unname(x0),
               tolerance = 1e-5)
})


## ---- Sub-kink coords are pinned to mu exactly -------------------------

test_that("trustL1 pins sub-kink coords to mu exactly via the active set", {
  # x0 places 'b' and 'd' clearly inside the soft-threshold dead zone
  # (|x0_i| < lambda/2 = 0.25); 'a' and 'c' are firmly outside.
  x0     <- c(a = 1.5, b = 0.10, c = -0.8, d = -0.05)
  lambda <- 0.5
  obj    <- .quadratic_objfn(x0)
  mu_L1  <- c(a = 0, b = 0, c = 0, d = 0)

  fit <- trustL1(obj, parinit = c(a = 0, b = 0, c = 0, d = 0),
                 mu = mu_L1, lambda = lambda,
                 rinit = 0.5, rmax = 5, iterlim = 200)

  # b and d must land at mu = 0 exactly (active-set pinning, not just
  # numerically close).
  expect_identical(unname(fit$argument[["b"]]), 0)
  expect_identical(unname(fit$argument[["d"]]), 0)

  # a and c follow the analytic soft-threshold solution.
  expect_equal(unname(fit$argument[["a"]]),
               1.5 - lambda / 2, tolerance = 1e-6)
  expect_equal(unname(fit$argument[["c"]]),
               -0.8 + lambda / 2, tolerance = 1e-6)

  # Reported value equals smooth-part value at theta plus L1 contribution.
  smooth_at_theta <- sum((fit$argument - x0)^2)
  l1_at_theta    <- lambda * sum(abs(fit$argument - mu_L1))
  expect_equal(fit$value, smooth_at_theta + l1_at_theta, tolerance = 1e-8)
})


## ---- Sub-kink coords stay pinned when parinit starts off the kink -----

test_that("trustL1 pulls sub-kink coords to mu when started away from it", {
  # Same dead-zone setup, but initialise far from mu so the optimiser has
  # to *cross* the kink to pin the parameter.
  x0     <- c(a = 1.5, b = 0.10)
  lambda <- 0.5
  obj    <- .quadratic_objfn(x0)
  mu_L1  <- c(a = 0, b = 0)

  fit <- trustL1(obj, parinit = c(a = 2.0, b = 1.0),
                 mu = mu_L1, lambda = lambda,
                 rinit = 0.5, rmax = 5, iterlim = 200)

  expect_identical(unname(fit$argument[["b"]]), 0)
  expect_equal(unname(fit$argument[["a"]]),
               1.5 - lambda / 2, tolerance = 1e-6)
})


## ---- One-sided penalty acts as a lower wall at mu ---------------------

test_that("trustL1 one-sided penalty pins coords pushing below mu", {
  # x0 pulls both coords below mu = 0. With one.sided = TRUE the penalty
  # only fires for theta_i < mu_i, so the optimum is the projection of
  # the soft-threshold solution onto [mu, +Inf).
  x0     <- c(a = -1.0, b = 0.5)
  lambda <- 0.4
  obj    <- .quadratic_objfn(x0)
  mu_L1  <- c(a = 0, b = 0)

  fit <- trustL1(obj, parinit = c(a = 0.5, b = 0.5),
                 mu = mu_L1, one.sided = TRUE, lambda = lambda,
                 rinit = 0.5, rmax = 5, iterlim = 200)

  # a wants to go to -1.0 but is held at the lower wall mu = 0.
  expect_identical(unname(fit$argument[["a"]]), 0)
  # b is above mu, penalty inactive, so b lands at the unpenalised min.
  expect_equal(unname(fit$argument[["b"]]), 0.5, tolerance = 1e-6)
})
