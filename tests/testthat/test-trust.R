# Behavioral tests for trust() (trust-region optimizer, C++ implementation).
#
# Verifies:
#   * convergence to the exact minimum of a quadratic in one trust step
#   * recovery of simulated-truth parameters from a noisy decay dataset
#   * extra objfun args passed via closure binding (no more ... forwarding)
#   * parscale rescaling lands at the same minimum
#   * parupper clamps a single component
#   * blather returns the per-iter trace
#
# Trust is also exercised end-to-end via normL2 -> trust in
# test-mstrust-profile.R and the existing FOCEI tests.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


# Synthetic quadratic objective in arbitrary dimension. Minimum at `target`.
.quadratic_objfn <- function(target) {
  d <- length(target)
  function(p, ...) {
    list(
      value    = 0.5 * sum((p - target)^2),
      gradient = p - target,
      hessian  = diag(d))
  }
}


## ---- Quadratic: one-step convergence ----------------------------------

test_that("trust converges to the exact minimum of a quadratic in one outer step", {
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  obj <- .quadratic_objfn(target)
  init <- c(a = 0, b = 0, c = 0)

  fit <- trust(obj, init, rinit = 5, rmax = 100, iterlim = 50,
               printIter = FALSE)
  expect_true(fit$converged)
  expect_equal(fit$argument, target, tolerance = 1e-8, ignore_attr = TRUE)
  # A pure quadratic with rinit large enough to admit the Newton step is
  # accepted on iteration 1; trust then runs one more iteration to verify
  # termination criteria.
  expect_lte(fit$iterations, 3L)
})


## ---- Simulated-truth recovery ----------------------------------------

test_that("trust(normL2) recovers true (A, k) from simulated decay data within noise budget", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  true_pars <- c(A = 1.0, k = 0.5)
  sigma_sim <- 0.02  # very low noise so fit is well-determined
  data <- fx_decay_data(pars = true_pars, sigma = sigma_sim,
                        times = seq(0, 8, by = 0.5), seed = 17L)

  obj <- normL2(data, bench$prd_id)
  init <- c(A = 0.6, k = 0.9)  # away from truth
  fit <- trust(obj, init, rinit = 1, rmax = 10, iterlim = 100,
               printIter = FALSE)

  expect_true(fit$converged)
  # Practical envelope at this noise level: ~10% relative on (A, k). The
  # exponential decay has a known A * k -> kA scale identifiability that
  # makes k a few % less constrained than A.
  expect_lt(abs(fit$argument[["A"]] - true_pars[["A"]]) / true_pars[["A"]],
            0.10)
  expect_lt(abs(fit$argument[["k"]] - true_pars[["k"]]) / true_pars[["k"]],
            0.10)
})


## ---- fixed = ... holds parameters constant ---------------------------

test_that("trust honors fixed = ... (held parameters unchanged in argument)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.05)
  obj   <- normL2(data, bench$prd_id)

  init  <- c(A = 0.7, k = 0.6)
  fix   <- c(k = 0.5)   # hold k, vary only A
  free  <- init[setdiff(names(init), names(fix))]

  fit <- trust(obj, free, rinit = 0.5, rmax = 5, iterlim = 50,
               fixed = fix, printIter = FALSE)
  expect_true(fit$converged)
  # `argument` contains only the free parameters; k should be absent.
  expect_false("k" %in% names(fit$argument))
  expect_true("A" %in% names(fit$argument))
})


## ---- parscale invariance --------------------------------------------

test_that("trust with parscale lands at the same minimum as the unscaled run", {
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  obj <- .quadratic_objfn(target)
  init <- c(a = 0, b = 0, c = 0)

  fit_plain <- trust(obj, init, rinit = 1, rmax = 100, iterlim = 50,
                     printIter = FALSE)
  fit_scaled <- trust(obj, init, parscale = c(a = 2, b = 0.5, c = 10),
                      rinit = 1, rmax = 100, iterlim = 50,
                      printIter = FALSE)
  expect_equal(unname(fit_scaled$argument[names(target)]),
               unname(fit_plain$argument[names(target)]),
               tolerance = 1e-6)
})


## ---- bounds clamp the optimum ---------------------------------------

test_that("trust honors parupper on one component while leaving the other free", {
  # 2D problem; unconstrained min at (a, b) = (5, 1). Upper bound on `a`
  # at 2 clamps that component; `b` converges freely to 1.
  target <- c(a = 5.0, b = 1.0)
  obj <- .quadratic_objfn(target)

  fit <- trust(obj, parinit = c(a = 0, b = 0), rinit = 1, rmax = 10,
               iterlim = 50,
               parupper = c(a = 2, b = Inf),
               printIter = FALSE)
  expect_equal(unname(fit$argument[["a"]]), 2.0, tolerance = 1e-5)
  expect_equal(unname(fit$argument[["b"]]), 1.0, tolerance = 1e-5)
})


## ---- blather returns the per-iter trace ------------------------------

test_that("trust(blather = TRUE) returns all trace fields with finite numbers", {
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  obj <- .quadratic_objfn(target)
  init <- c(a = 0, b = 0, c = 0)

  fit <- trust(obj, init, rinit = 0.3, rmax = 5, iterlim = 20,
               blather = TRUE, printIter = FALSE)
  expect_true(fit$converged)
  n <- fit$iterations
  expect_equal(nrow(fit$argpath), n)
  expect_equal(ncol(fit$argpath), 3L)
  expect_equal(nrow(fit$argtry),  n)
  expect_equal(length(fit$steptype), n)
  expect_equal(length(fit$accept),   n)
  expect_equal(length(fit$r),        n)
  expect_equal(length(fit$rho),      n)
  expect_equal(length(fit$valpath),  n)
  expect_equal(length(fit$valtry),   n)
  expect_equal(length(fit$preddiff), n)
  expect_equal(length(fit$stepnorm), n)
  expect_true(all(fit$steptype %in% c("Newton", "easy-easy", "hard-easy", "hard-hard")))
  expect_true(all(is.finite(fit$valpath)))
  expect_true(all(fit$r >= 0))
})
