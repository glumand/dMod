# Behavioral tests for constraintL2() on the scalar / diagonal path.
#
# Soft Gaussian L2 constraint:
#
#   value      = sum_i ((p_i - mu_i) / sigma_i)^2  [+ 2 * sum log(sigma) if est]
#   gradient_i = 2 * (p_i - mu_i) / sigma_i^2
#   Hessian    = 2 * diag(1 / sigma_i^2)
#
# When `sigma` is character, the sigma parameters are themselves estimated
# (log-parametrised) so the gradient/Hessian pick up extra terms; this file
# focuses on the fixed-sigma path and the diagonal multi-element case. The
# MVN-Omega path lives in test-constraintL2-mvn.R and is unchanged.

# Note: constraintL2 *not* a typical fit objective (no prediction frame
# needed); we call it directly. Backend parametrisation still applies so
# both R and C++ paths are exercised.


## ---- Scalar sigma, single parameter ------------------------------------

test_that("constraintL2 with scalar sigma equals (p - mu)^2 / sigma^2", {
  mu    <- c(theta = 0.5)
  sigma <- 0.1
  obj   <- constraintL2(mu = mu, sigma = sigma)

  for_each_backend(function(cpp) {
    o <- obj(c(theta = 0.8))
    expected <- ((0.8 - 0.5) / sigma)^2
    expect_equal(unname(o$value), expected, tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    expect_equal(unname(o$gradient[["theta"]]),
                 2 * (0.8 - 0.5) / sigma^2, tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    expect_equal(unname(o$hessian[1, 1]), 2 / sigma^2,
                 tolerance = 1e-12, info = paste0("cpp=", cpp))
  })
})


## ---- Vector sigma, multiple parameters ---------------------------------

test_that("constraintL2 sums per-parameter contributions when mu has length > 1", {
  mu    <- c(a = 0.0, b = 1.0, c = -0.5)
  sigma <- c(a = 1.0, b = 0.5, c = 0.2)
  obj   <- constraintL2(mu = mu, sigma = sigma)
  p     <- c(a = 0.3, b = 1.4, c = -0.6)

  for_each_backend(function(cpp) {
    o <- obj(p)
    expected <- sum(((p - mu) / sigma)^2)
    expect_equal(unname(o$value), unname(expected), tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    expect_equal(unname(o$gradient[names(p)]),
                 unname(2 * (p - mu) / sigma^2), tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    expect_equal(unname(diag(o$hessian)),
                 unname(2 / sigma^2), tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    # Off-diagonal Hessian entries are zero for the diagonal-sigma case.
    H <- o$hessian
    expect_lt(max(abs(H[upper.tri(H)])), 1e-12,
              label = paste0("cpp=", cpp, " off-diag"))
  })
})


## ---- Value is zero at mu, gradient is zero at mu -----------------------

test_that("constraintL2 has value 0 and gradient 0 at the prior mean", {
  mu  <- c(a = 0.1, b = 0.2, c = 0.3)
  sg  <- c(a = 0.5, b = 0.5, c = 0.5)
  obj <- constraintL2(mu = mu, sigma = sg)

  for_each_backend(function(cpp) {
    o <- obj(mu)
    expect_equal(unname(o$value), 0, tolerance = 1e-14,
                 info = paste0("cpp=", cpp))
    expect_lt(max(abs(o$gradient)), 1e-14,
              label = paste0("cpp=", cpp))
  })
})


## ---- Subset of parameters: only matching names contribute --------------

test_that("constraintL2 only constrains parameters whose names appear in mu", {
  mu  <- c(a = 0.0, b = 0.0)
  sg  <- 1.0
  obj <- constraintL2(mu = mu, sigma = sg)
  # `c` is not in mu -> no constraint on it.
  p   <- c(a = 0.4, b = -0.3, c = 100)

  for_each_backend(function(cpp) {
    o <- obj(p)
    expect_equal(unname(o$value), 0.4^2 + (-0.3)^2, tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    expect_equal(unname(o$gradient[["c"]]), 0, tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
  })
})


## ---- Composition: + with another constraintL2 sums values -------------

test_that("(constraintL2(mu1) + constraintL2(mu2))(pars) sums per-element", {
  mu1 <- c(a = 0.0, b = 0.0); mu2 <- c(a = 1.0, b = 1.0)
  obj <- constraintL2(mu = mu1, sigma = 1, attr.name = "c1") +
         constraintL2(mu = mu2, sigma = 1, attr.name = "c2")
  p <- c(a = 0.3, b = 0.7)

  for_each_backend(function(cpp) {
    v_total <- obj(p)$value
    v_a <- (0.3 - 0)^2 + (0.7 - 0)^2 + (0.3 - 1)^2 + (0.7 - 1)^2
    expect_equal(unname(v_total), v_a, tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
  })
})
