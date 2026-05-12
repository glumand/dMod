# Behavioral tests for priorL2() (exponential L2 prior).
#
# priorL2 computes  value = exp(lambda) * ||p - mu||^2  with derivatives w.r.t.
# both p (the parameters constrained by mu) and lambda (the strength of the
# constraint, parametrised log-scale). All assertions are direct closed-form
# checks against the documented formula.


## ---- Value: zero at p == mu --------------------------------------------

test_that("priorL2 value is zero when all constrained parameters equal mu", {
  mu <- c(A = 0.5, B = -0.3)
  obj <- priorL2(mu = mu, lambda = "lambda")
  p   <- c(A = 0.5, B = -0.3, lambda = 0)
  expect_equal(obj(p)$value, 0, tolerance = 1e-14)
})


## ---- Value: closed form on a non-zero displacement ---------------------

test_that("priorL2 value equals exp(lambda) * sum((p - mu)^2)", {
  mu <- c(A = 0.0, B = 0.0, C = 0.0)
  lambda_val <- 0.7
  obj <- priorL2(mu = mu, lambda = "lambda")
  p   <- c(A = 0.4, B = -0.2, C = 1.1, lambda = lambda_val)
  expected <- exp(lambda_val) * sum((p[names(mu)] - mu)^2)
  expect_equal(unname(obj(p)$value), unname(expected), tolerance = 1e-12)
})


## ---- Gradient w.r.t. p -------------------------------------------------

test_that("priorL2 gradient w.r.t. p is 2 * exp(lambda) * (p - mu)", {
  mu <- c(A = 0.1, B = -0.4)
  lambda_val <- -0.3
  obj <- priorL2(mu = mu, lambda = "lambda")
  p   <- c(A = 0.7, B = 0.2, lambda = lambda_val)
  g <- obj(p)$gradient

  expect_equal(unname(g[["A"]]),
               2 * exp(lambda_val) * (p[["A"]] - mu[["A"]]),
               tolerance = 1e-12)
  expect_equal(unname(g[["B"]]),
               2 * exp(lambda_val) * (p[["B"]] - mu[["B"]]),
               tolerance = 1e-12)
})


## ---- Gradient w.r.t. lambda --------------------------------------------

test_that("priorL2 gradient w.r.t. lambda equals exp(lambda) * sum((p - mu)^2)", {
  mu <- c(A = 0.0, B = 0.0)
  lambda_val <- 0.5
  obj <- priorL2(mu = mu, lambda = "lambda")
  p   <- c(A = 1.0, B = -0.5, lambda = lambda_val)
  g <- obj(p)$gradient
  expected <- exp(lambda_val) * sum((p[names(mu)] - mu)^2)
  expect_equal(unname(g[["lambda"]]), unname(expected), tolerance = 1e-12)
})


## ---- Hessian: diagonal with closed-form entries ------------------------

test_that("priorL2 Hessian: diag(p) = 2 * exp(lambda); cross terms via product rule", {
  mu <- c(A = 0.0, B = 0.0)
  lambda_val <- 0.0
  obj <- priorL2(mu = mu, lambda = "lambda")
  p   <- c(A = 0.3, B = -0.1, lambda = lambda_val)

  H <- obj(p)$hessian
  # d^2 / dA^2 = 2 * exp(lambda)
  expect_equal(unname(H["A", "A"]), 2 * exp(lambda_val), tolerance = 1e-12)
  expect_equal(unname(H["B", "B"]), 2 * exp(lambda_val), tolerance = 1e-12)
  # d^2 / dA dB = 0
  expect_equal(unname(H["A", "B"]), 0, tolerance = 1e-12)
  # d^2 / dA dlambda = 2 * exp(lambda) * (A - mu_A)
  expect_equal(unname(H["A", "lambda"]),
               2 * exp(lambda_val) * (p[["A"]] - mu[["A"]]),
               tolerance = 1e-12)
  # d^2 / dlambda^2 = exp(lambda) * sum((p - mu)^2)
  expect_equal(unname(H["lambda", "lambda"]),
               exp(lambda_val) * sum((p[names(mu)] - mu)^2),
               tolerance = 1e-12)
})


## ---- Parameter outside mu does not contribute --------------------------

test_that("parameters not in mu are unconstrained (zero gradient entries)", {
  mu <- c(A = 0.0)
  obj <- priorL2(mu = mu, lambda = "lambda")
  p   <- c(A = 0.2, OTHER = 1e6, lambda = 0)
  g <- obj(p)$gradient
  expect_equal(unname(g[["OTHER"]]), 0, tolerance = 1e-12)
})
