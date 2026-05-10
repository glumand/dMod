test_that("Y deriv2 (AD) reproduces analytical observation Hessian", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  # Plain Y at the head of the chain: g = a*x^2 + b*x.
  # No upstream parameter transformation, so theta = c(states, params)
  # = c(x, a, b). Y returns d2y/dtheta^2 directly.
  gfn <- Y(c(y = "a*x^2 + b*x"), states = "x", parameters = c("a", "b"),
           modelname = paste0("y_d2_", as.integer(Sys.time())),
           compile = TRUE, deriv2 = TRUE, derivMode = "dual",
           attach.input = FALSE)

  # Build a synthetic prediction "out" with matching parameters attribute.
  times <- c(0.0, 0.5, 1.0)
  out <- cbind(time = times, x = c(1.0, 0.7, 0.4))
  class(out) <- c("prdframe", "matrix")
  pars <- as.parvec(c(a = 0.3, b = -0.5))

  res <- gfn(out, pars, deriv = TRUE, deriv2 = TRUE)[[1]]

  # Analytical references
  a <- pars["a"]; b <- pars["b"]
  xv <- out[, "x"]
  ref_dy_dx <- 2 * a * xv + b
  ref_dy_da <- xv^2
  ref_dy_db <- xv

  d1 <- attr(res, "deriv")
  d2 <- attr(res, "deriv2")
  # First-order
  expect_equal(d1[, "y", "x"], unname(ref_dy_dx), tolerance = 1e-10)
  expect_equal(d1[, "y", "a"], unname(ref_dy_da), tolerance = 1e-10)
  expect_equal(d1[, "y", "b"], unname(ref_dy_db), tolerance = 1e-10)

  # Hessian per timepoint
  for (i in seq_along(times)) {
    H_ref <- matrix(0, 3, 3, dimnames = list(c("x", "a", "b"), c("x", "a", "b")))
    H_ref["x", "x"] <- 2 * a
    H_ref["x", "a"] <- 2 * xv[i]; H_ref["a", "x"] <- 2 * xv[i]
    H_ref["x", "b"] <- 1;          H_ref["b", "x"] <- 1
    expect_equal(d2[i, "y", c("x", "a", "b"), c("x", "a", "b")], H_ref,
                 tolerance = 1e-10)
  }
})

test_that("Y deriv2 (symbolic) reproduces analytical observation Hessian", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  gfn <- Y(c(y = "a*x^2 + b*x"), states = "x", parameters = c("a", "b"),
           modelname = paste0("y_d2_sym_", as.integer(Sys.time())),
           compile = TRUE, deriv2 = TRUE, derivMode = "symbolic",
           attach.input = FALSE)

  times <- c(0.5, 1.0)
  out <- cbind(time = times, x = c(0.6, 0.3))
  class(out) <- c("prdframe", "matrix")
  pars <- as.parvec(c(a = -0.2, b = 0.4))

  res <- gfn(out, pars, deriv = TRUE, deriv2 = TRUE)[[1]]
  a <- pars["a"]; b <- pars["b"]
  xv <- out[, "x"]
  d2 <- attr(res, "deriv2")
  for (i in seq_along(times)) {
    H_ref <- matrix(0, 3, 3, dimnames = list(c("x", "a", "b"), c("x", "a", "b")))
    H_ref["x", "x"] <- 2 * a
    H_ref["x", "a"] <- 2 * xv[i]; H_ref["a", "x"] <- 2 * xv[i]
    H_ref["x", "b"] <- 1;          H_ref["b", "x"] <- 1
    expect_equal(d2[i, "y", c("x", "a", "b"), c("x", "a", "b")], H_ref,
                 tolerance = 1e-10)
  }
})
