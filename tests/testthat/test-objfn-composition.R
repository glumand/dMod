# Behavioral tests for objective and function composition operators.
#
# Operators tested:
#   * `+.objfn`   sums two objective functions (value, gradient, Hessian)
#   * `*.fn`      composes prediction functions  (f * g)(x) = f(g(x))
#   * `+.fn`      merges parameter trafos across conditions
#   * getParameters() consistency through composition

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


## ---- +.objfn: value sums elementwise -----------------------------------

test_that("(constraintL2(mu1) + constraintL2(mu2))(p)$value = sum of parts", {
  mu1 <- c(a = 0, b = 0)
  mu2 <- c(a = 1, b = -1)
  o1 <- constraintL2(mu1, sigma = 1, attr.name = "p1")
  o2 <- constraintL2(mu2, sigma = 1, attr.name = "p2")
  obj <- o1 + o2

  p <- c(a = 0.3, b = 0.4)
  v_total <- obj(p)$value
  v_parts <- o1(p)$value + o2(p)$value
  expect_equal(unname(v_total), unname(v_parts), tolerance = 1e-12)
})


test_that("(o1 + o2)$gradient sums per-parameter contributions", {
  mu1 <- c(a = 0, b = 0); mu2 <- c(a = 1, b = -1)
  o1 <- constraintL2(mu1, sigma = 1, attr.name = "p1")
  o2 <- constraintL2(mu2, sigma = 1, attr.name = "p2")
  obj <- o1 + o2
  p <- c(a = 0.3, b = 0.4)

  g_total <- obj(p)$gradient
  g_parts <- o1(p)$gradient + o2(p)$gradient
  expect_equal(unname(g_total[names(p)]),
               unname(g_parts[names(p)]), tolerance = 1e-12)
})


test_that("(o1 + o2)$hessian sums per-parameter contributions blockwise", {
  mu1 <- c(a = 0, b = 0); mu2 <- c(a = 1, b = -1)
  o1 <- constraintL2(mu1, sigma = 1)
  o2 <- constraintL2(mu2, sigma = 1)
  obj <- o1 + o2
  p <- c(a = 0.3, b = 0.4)

  H_total <- obj(p)$hessian
  H_parts <- o1(p)$hessian + o2(p)$hessian
  expect_equal(unname(H_total[names(p), names(p)]),
               unname(H_parts[names(p), names(p)]), tolerance = 1e-12)
})


## ---- *.fn chains Jacobians on prediction functions --------------------

test_that("(g * x)(times, pars) equals g(x(times, pars)) for the decay chain", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.0, k = 0.5)

  out_chain <- (bench$gfn * bench$xfn * bench$pfn_id)(times, pars,
                                                     deriv = FALSE)
  # Compute the same via the intermediate state directly. We can't easily
  # call gfn(prdframe, pars) outside the chain, so instead we check that
  # the chained y-value equals the closed-form expression for our observable
  # y = A and the analytical decay solution.
  expect_equal(out_chain$C1[, "y"], pars[["A"]] * exp(-pars[["k"]] * times),
               tolerance = 1e-5)
})


## ---- constraintL2 + normL2: prior penalty added exactly ---------------

test_that("normL2 + constraintL2 value equals normL2 + prior penalty (closed)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)
  mu    <- c(A = 0.0, k = 0.0)
  sigma <- 1.0
  prior <- constraintL2(mu = mu, sigma = sigma)

  obj_main <- normL2(data, bench$prd_id)
  obj      <- obj_main + prior

  p <- c(A = 1.2, k = 0.4)
  v_total <- obj(p)$value
  v_main  <- obj_main(p)$value
  v_prior <- sum((p - mu)^2)  # since sigma = 1
  expect_equal(unname(v_total), unname(v_main + v_prior), tolerance = 1e-9)
})


## ---- Parameter set union through composition --------------------------

test_that("attr(o1 + o2, 'parameters') = union(parameters(o1), parameters(o2))", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  obj_main <- normL2(data, bench$prd_id)
  prior_extra <- constraintL2(mu = c(A = 0.0, gamma = 0.0), sigma = 1.0)
  obj <- obj_main + prior_extra

  expect_setequal(attr(obj, "parameters"),
                  union(attr(obj_main, "parameters"),
                        attr(prior_extra, "parameters")))
})


## ---- Order invariance of +.objfn ---------------------------------------

test_that("(o1 + o2)(p) equals (o2 + o1)(p) at the value level", {
  o1 <- constraintL2(c(a = 0), sigma = 1)
  o2 <- constraintL2(c(a = 2), sigma = 1)
  v12 <- (o1 + o2)(c(a = 0.5))$value
  v21 <- (o2 + o1)(c(a = 0.5))$value
  expect_equal(unname(v12), unname(v21), tolerance = 1e-12)
})
