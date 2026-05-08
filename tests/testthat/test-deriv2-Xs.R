test_that("Xs.CppODE deriv2 reproduces linear-decay analytical Hessian", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  setwd(tempdir())

  f <- c(x = "-k * x")  # x(t) = x0 * exp(-k*t)
  m <- odemodel(f, modelname = paste0("decay_d2_", as.integer(Sys.time())),
                solver = "CppODE", deriv2 = TRUE, verbose = FALSE)
  xfn <- Xs(m, optionsOde = list(atol = 1e-10, rtol = 1e-10),
            optionsSens = list(atol = 1e-10, rtol = 1e-10))

  times <- c(0.0, 0.5, 1.0, 1.5)
  pars <- c(x = 2.0, k = 0.7)
  pred <- xfn(times, pars, deriv = TRUE, deriv2 = TRUE)[[1]]

  x0 <- pars["x"]; k <- pars["k"]
  ref_J  <- cbind(x  = exp(-k * times),
                  k  = -x0 * times * exp(-k * times))
  ref_H_kk <-  x0 * times^2 * exp(-k * times)
  ref_H_xk <- -times * exp(-k * times)
  ref_H_xx <- rep(0, length(times))

  d1 <- attr(pred, "deriv")
  d2 <- attr(pred, "deriv2")

  expect_equal(d1[, "x", "x"], unname(ref_J[, "x"]), tolerance = 1e-7)
  expect_equal(d1[, "x", "k"], unname(ref_J[, "k"]), tolerance = 1e-7)

  expect_equal(d2[, "x", "x", "x"], unname(ref_H_xx), tolerance = 1e-7)
  expect_equal(d2[, "x", "k", "k"], unname(ref_H_kk), tolerance = 1e-7)
  expect_equal(d2[, "x", "x", "k"], unname(ref_H_xk), tolerance = 1e-7)
  expect_equal(d2[, "x", "k", "x"], unname(ref_H_xk), tolerance = 1e-7)
})

test_that("odemodel(deriv2 = TRUE) emits <m>, <m>_s, <m>_s2 and Xs.CppODE dispatches", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  td <- file.path(tempdir(), paste0("triple_", as.integer(Sys.time())))
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  setwd(td)

  f <- c(x = "-k * x")
  modelname <- "triple"
  m <- odemodel(f, modelname = modelname, solver = "CppODE",
                deriv2 = TRUE, verbose = FALSE)
  expect_true(file.exists(paste0(modelname,    ".cpp")))
  expect_true(file.exists(paste0(modelname, "_s.cpp")))
  expect_true(file.exists(paste0(modelname, "_s2.cpp")))
  expect_false(isTRUE(attr(m$extended,  "deriv2")))
  expect_true(isTRUE(attr(m$extended2, "deriv2")))

  # 1st-order call should run via the cheaper _s extension and not return d2y;
  # 2nd-order call uses _s2 and returns the 4D Hessian array.
  xfn <- Xs(m, optionsOde = list(atol = 1e-12, rtol = 1e-12),
            optionsSens = list(atol = 1e-12, rtol = 1e-12))
  times <- c(0.0, 0.5, 1.0)
  pars <- c(x = 2.0, k = 0.7)
  r1 <- xfn(times, pars, deriv = TRUE, deriv2 = FALSE)[[1]]
  r2 <- xfn(times, pars, deriv = TRUE, deriv2 = TRUE )[[1]]

  expect_null(attr(r1, "deriv2"))
  expect_equal(dim(attr(r2, "deriv2")), c(length(times), 1L, 2L, 2L))

  # Both paths must agree on the 1st-order sensitivities up to integrator
  # tolerance (we use 1e-12 atol/rtol so a strict 1e-9 here is realistic).
  expect_equal(attr(r1, "deriv"), attr(r2, "deriv"), tolerance = 1e-9)
})

test_that("Xs.deSolve refuses deriv2 = TRUE", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  setwd(tempdir())

  f <- c(x = "-k*x")
  m <- odemodel(f, modelname = paste0("decay_des_", as.integer(Sys.time())),
                solver = "deSolve", verbose = FALSE)
  xfn <- Xs(m)
  expect_error(xfn(c(0, 1), c(x = 1, k = 0.5), deriv2 = TRUE),
               "Xs.deSolve")
})
