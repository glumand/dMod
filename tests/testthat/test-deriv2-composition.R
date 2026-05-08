test_that("(Y * Xs)(times, pars, deriv2 = TRUE) chain-rule matches analytical", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  setwd(tempdir())

  f <- c(x = "-k * x")
  m <- odemodel(f,
                modelname = paste0("compose_decay_", as.integer(Sys.time())),
                solver = "CppODE", deriv2 = TRUE, verbose = FALSE)
  xfn <- Xs(m, optionsOde = list(atol = 1e-10, rtol = 1e-10),
            optionsSens = list(atol = 1e-10, rtol = 1e-10))
  gfn <- Y(c(y = "a*x^2 + b*x"), f = f, parameters = c("a", "b"),
           modelname = paste0("compose_obs_", as.integer(Sys.time())),
           compile = TRUE, deriv2 = TRUE, attach.input = FALSE)
  prd <- gfn * xfn

  times <- c(0.0, 0.5, 1.0)
  pars <- c(x = 2.0, k = 0.7, a = 0.3, b = -0.5)
  res <- prd(times, pars, deriv = TRUE, deriv2 = TRUE)[[1]]
  d2 <- attr(res, "deriv2")

  # y = a*x0^2*exp(-2kt) + b*x0*exp(-kt). Closed-form Hessian wrt (x0, k):
  x0 <- pars["x"]; k <- pars["k"]; a <- pars["a"]; b <- pars["b"]
  for (i in seq_along(times)) {
    t <- times[i]
    e1 <- exp(-k * t); e2 <- exp(-2 * k * t)
    H_x0_x0 <- 2 * a * e2
    H_k_k   <- 4 * t^2 * a * x0^2 * e2 + t^2 * b * x0 * e1
    H_x0_k  <- -4 * t * a * x0 * e2 - t * b * e1
    expect_equal(d2[i, "y", "x", "x"], unname(H_x0_x0), tolerance = 1e-7)
    expect_equal(d2[i, "y", "k", "k"], unname(H_k_k),   tolerance = 1e-7)
    expect_equal(d2[i, "y", "x", "k"], unname(H_x0_k),  tolerance = 1e-7)
    expect_equal(d2[i, "y", "k", "x"], unname(H_x0_k),  tolerance = 1e-7)
  }
})
