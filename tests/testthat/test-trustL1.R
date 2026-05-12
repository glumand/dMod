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
