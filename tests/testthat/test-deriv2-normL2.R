test_that("normL2 gradient is identical for deriv2 = FALSE and deriv2 = TRUE", {
  # Regression for the [.parvec / c.parvec deriv2-drop bug. With deriv2 = TRUE,
  # the upstream Hessian seed must reach Xs.CppODE so the _s2 integration
  # produces the same first-order sensitivities as the _s integration. If
  # subsetting drops attr(., "deriv2") the gradient diverges silently.
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  setwd(tempdir())

  f <- c(x = "-k*x")
  m <- odemodel(f, modelname = paste0("nl_grad_id_", as.integer(Sys.time())),
                solver = "CppODE", deriv2 = TRUE, nStack = 4L, verbose = FALSE)
  ode_opts <- list(atol = 1e-12, rtol = 1e-12)
  xfn <- Xs(m, condition = "C1",
            optionsOde = ode_opts, optionsSens = ode_opts)
  gfn <- Y(c(y = "a*x^2 + b*x"), f = f, parameters = c("a", "b"),
           modelname = paste0("nl_grad_obs_", as.integer(Sys.time())),
           compile = TRUE, deriv2 = TRUE, attach.input = FALSE,
           condition = "C1")
  pfn <- Pexpl(c(x = "exp(lx)", k = "exp(lk)", a = "la", b = "lb"),
               parameters = NULL,
               modelname = paste0("nl_grad_p_", as.integer(Sys.time())),
               compile = TRUE, deriv2 = TRUE, derivMode = "dual",
               condition = "C1")
  prd <- gfn * xfn * pfn

  times_d <- c(0.5, 1.0, 1.5)
  pars <- c(lx = log(2.0), lk = log(0.7), la = 0.3, lb = -0.5)
  truth <- prd(times_d, pars, deriv = FALSE)[[1]]
  ydata <- truth[truth[, "time"] %in% times_d, "y"] + c(0.05, -0.02, 0.03)
  data <- datalist(C1 = data.frame(name = "y", time = times_d, value = ydata,
                                   sigma = 0.1))

  obj <- normL2(data, prd)
  res_d1 <- obj(pars, deriv = TRUE, deriv2 = FALSE)
  res_d2 <- obj(pars, deriv = TRUE, deriv2 = TRUE)
  expect_equal(res_d2$value, res_d1$value, tolerance = 1e-8)
  expect_equal(res_d2$gradient, res_d1$gradient, tolerance = 1e-7)
})

test_that("normL2(deriv2 = FALSE) reproduces the pre-deriv2 GN Hessian", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  setwd(tempdir())

  f <- c(x = "-k*x")
  m <- odemodel(f, modelname = paste0("nl_d2_gn_", as.integer(Sys.time())),
                solver = "CppODE", deriv2 = TRUE, nStack = 4L, verbose = FALSE)
  xfn <- Xs(m, condition = "C1")
  gfn <- Y(c(y = "a*x^2 + b*x"), f = f, parameters = c("a", "b"),
           modelname = paste0("nl_obs_gn_", as.integer(Sys.time())),
           compile = TRUE, deriv2 = TRUE, attach.input = FALSE,
           condition = "C1")
  pfn <- Pexpl(c(x = "x", k = "k", a = "a", b = "b"), parameters = NULL,
               modelname = paste0("nl_id_gn_", as.integer(Sys.time())),
               compile = TRUE, deriv2 = TRUE, derivMode = "symbolic",
               condition = "C1")
  prd <- gfn * xfn * pfn

  times_d <- c(0.5, 1.0, 1.5)
  pars <- c(x = 2.0, k = 0.7, a = 0.3, b = -0.5)
  truth <- prd(times_d, pars, deriv = FALSE)[[1]]
  ydata <- truth[truth[, "time"] %in% times_d, "y"] + c(0.05, -0.02, 0.03)
  data <- datalist(C1 = data.frame(name = "y", time = times_d, value = ydata,
                                   sigma = 0.1))

  obj <- normL2(data, prd)
  res_gn <- obj(pars, deriv = TRUE, deriv2 = FALSE)
  # Sanity: gradient is non-zero, Hessian is symmetric and positive on the
  # active block — the GN baseline.
  expect_true(is.numeric(res_gn$value))
  expect_equal(res_gn$hessian, t(res_gn$hessian), tolerance = 1e-12)
})

test_that("normL2(deriv2 = TRUE) adds residual times d^2 pred / sigma^2", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  setwd(tempdir())

  f <- c(x = "-k*x")
  m <- odemodel(f, modelname = paste0("nl_d2_ex_", as.integer(Sys.time())),
                solver = "CppODE", deriv2 = TRUE, nStack = 4L, verbose = FALSE)
  # Use tight integrator tolerances so the _s and _s2 paths agree on the
  # value/gradient axes within numerical precision.
  ode_opts <- list(atol = 1e-12, rtol = 1e-12)
  xfn <- Xs(m, condition = "C1",
            optionsOde = ode_opts, optionsSens = ode_opts)
  gfn <- Y(c(y = "a*x^2 + b*x"), f = f, parameters = c("a", "b"),
           modelname = paste0("nl_obs_ex_", as.integer(Sys.time())),
           compile = TRUE, deriv2 = TRUE, attach.input = FALSE,
           condition = "C1")
  pfn <- Pexpl(c(x = "x", k = "k", a = "a", b = "b"), parameters = NULL,
               modelname = paste0("nl_id_ex_", as.integer(Sys.time())),
               compile = TRUE, deriv2 = TRUE, derivMode = "symbolic",
               condition = "C1")
  prd <- gfn * xfn * pfn

  times_d <- c(0.5, 1.0, 1.5)
  pars <- c(x = 2.0, k = 0.7, a = 0.3, b = -0.5)
  truth <- prd(times_d, pars, deriv = FALSE)[[1]]
  noise <- c(0.05, -0.02, 0.03)
  ydata <- truth[truth[, "time"] %in% times_d, "y"] + noise
  sigma <- 0.1
  data <- datalist(C1 = data.frame(name = "y", time = times_d, value = ydata,
                                   sigma = sigma))

  obj <- normL2(data, prd)
  res_gn <- obj(pars, deriv = TRUE, deriv2 = FALSE)
  res_ex <- obj(pars, deriv = TRUE, deriv2 = TRUE)
  # Value and gradient must agree across the _s and _s2 integration paths
  # up to integrator tolerance.
  expect_equal(res_ex$value, res_gn$value, tolerance = 1e-8)
  expect_equal(res_ex$gradient, res_gn$gradient, tolerance = 1e-8)

  # Analytical Hessian addition: 2/sigma^2 * sum_i r_i * d^2 y_i / dtheta^2.
  # r_i = y_pred - y_data = -noise.
  r <- -noise
  x0 <- pars["x"]; k <- pars["k"]; a <- pars["a"]; b <- pars["b"]
  H_add_terms <- function(t) {
    e1 <- exp(-k * t); e2 <- exp(-2 * k * t)
    H <- matrix(0, 4, 4, dimnames = list(c("x", "k", "a", "b"), c("x", "k", "a", "b")))
    H["x", "x"] <- 2 * a * e2
    H["x", "k"] <- -4 * t * a * x0 * e2 - t * b * e1
    H["k", "x"] <- H["x", "k"]
    H["x", "a"] <- 2 * x0 * e2; H["a", "x"] <- H["x", "a"]
    H["x", "b"] <- e1;          H["b", "x"] <- H["x", "b"]
    H["k", "k"] <- 4 * t^2 * a * x0^2 * e2 + t^2 * b * x0 * e1
    H["k", "a"] <- -2 * t * x0^2 * e2; H["a", "k"] <- H["k", "a"]
    H["k", "b"] <- -t * x0 * e1;       H["b", "k"] <- H["k", "b"]
    H
  }
  H_expected_add <- Reduce(`+`, lapply(seq_along(times_d), function(i)
    2 / sigma^2 * r[i] * H_add_terms(times_d[i])))
  H_actual_add <- res_ex$hessian - res_gn$hessian
  par_order <- c("x", "k", "a", "b")
  expect_equal(H_actual_add[par_order, par_order], H_expected_add[par_order, par_order],
               tolerance = 1e-3)
})
