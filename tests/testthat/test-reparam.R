context("Xs.CppODE theta-sensitivity path (heap vs stack AD slab)")

# Tests that the unified Xs.CppODE path (Phi'(theta) as sens1ini) produces
# identical sensitivities whether the underlying AD slab is heap-allocated
# (default nStack = Inf) or compile-time stack-allocated (finite nStack).
# Per-condition varying theta counts are covered by a two-condition setup.

test_that("Heap and stack AD slabs match on a single-condition linear model", {

  setwd(tempdir())

  f <- c(A = "-k1*A + k2*B",
         B =  "k1*A - k2*B")

  # Default heap slab vs explicit stack slab (nStack = 4 covers {A,B,log_k1,log_k2}).
  mod_v1 <- odemodel(f, modelname = "rep_v1", solver = "CppODE")
  mod_v2 <- odemodel(f, modelname = "rep_v2", solver = "CppODE", nStack = 4L)

  # Same parameter transformation for both.
  trafo <- c(A = "A", B = "B", k1 = "exp(log_k1)", k2 = "exp(log_k2)")
  p1 <- P(trafo, modelname = "rep_trafo_v1", compile = TRUE)
  p2 <- P(trafo, modelname = "rep_trafo_v2", compile = TRUE)

  tight <- list(atol = 1e-10, rtol = 1e-10)
  x1 <- Xs(mod_v1, optionsSens = tight) * p1
  x2 <- Xs(mod_v2, optionsSens = tight) * p2

  theta <- c(A = 1.0, B = 0.2, log_k1 = log(0.5), log_k2 = log(0.3))
  times <- seq(0, 3, length.out = 7)

  pred1 <- x1(times, theta,
              conditions = NULL,
              deriv = TRUE)
  pred2 <- x2(times, theta,
              conditions = NULL,
              deriv = TRUE)

  d1 <- getDerivs(pred1)[[1]]
  d2 <- getDerivs(pred2)[[1]]

  arr1 <- attr(d1, "deriv")
  arr2 <- attr(d2, "deriv")
  expect_equal(dim(arr1), dim(arr2))
  expect_equal(dimnames(arr1)[[2]], dimnames(arr2)[[2]])
  expect_equal(dimnames(arr1)[[3]], dimnames(arr2)[[3]])
  expect_equal(as.numeric(arr1), as.numeric(arr2), tolerance = 1e-6)
})


test_that("Heap/stack parity holds with per-condition varying theta subsets", {

  setwd(tempdir())

  f <- c(A = "-k1*A + k2*B",
         B =  "k1*A - k2*B")

  mod_v1 <- odemodel(f, modelname = "repmulti_v1", solver = "CppODE")
  # Stack upper bound: any condition may activate up to 4 thetas.
  mod_v2 <- odemodel(f, modelname = "repmulti_v2", solver = "CppODE", nStack = 4L)

  # Condition "closed" uses log_k1; condition "open" uses log_k_open instead.
  # Global theta set has 5 elements; each condition activates 4.
  trafo_closed <- c(A = "A", B = "B", k1 = "exp(log_k1)",      k2 = "exp(log_k2)")
  trafo_open   <- c(A = "A", B = "B", k1 = "exp(log_k_open)",  k2 = "exp(log_k2)")

  p1 <-
    P(trafo_closed, condition = "closed",
      modelname = "repmulti_trafo_cl_v1", compile = TRUE) +
    P(trafo_open,   condition = "open",
      modelname = "repmulti_trafo_op_v1", compile = TRUE)

  p2 <-
    P(trafo_closed, condition = "closed",
      modelname = "repmulti_trafo_cl_v2", compile = TRUE) +
    P(trafo_open,   condition = "open",
      modelname = "repmulti_trafo_op_v2", compile = TRUE)

  tight <- list(atol = 1e-10, rtol = 1e-10)
  x1 <- Xs(mod_v1, optionsSens = tight) * p1
  x2 <- Xs(mod_v2, optionsSens = tight) * p2

  theta <- c(A = 1.0, B = 0.2,
             log_k1 = log(0.5), log_k_open = log(0.8), log_k2 = log(0.3))
  times <- seq(0, 2, length.out = 5)

  pred1 <- x1(times, theta, deriv = TRUE)
  pred2 <- x2(times, theta, deriv = TRUE)

  for (cond in c("closed", "open")) {
    arr1 <- attr(pred1[[cond]], "deriv")
    arr2 <- attr(pred2[[cond]], "deriv")
    expect_equal(dim(arr1), dim(arr2),
                 info = paste("shape mismatch for condition", cond))
    expect_equal(dim(arr2)[3], 4L,
                 info = paste("ncol mismatch for condition", cond))
    expect_equal(dimnames(arr1)[[3]], dimnames(arr2)[[3]],
                 info = paste("theta names mismatch for condition", cond))
    expect_equal(as.numeric(arr1), as.numeric(arr2), tolerance = 1e-6,
                 info = paste("values mismatch for condition", cond))
  }
})
