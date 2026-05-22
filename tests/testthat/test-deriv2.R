# ============================================================================
# deriv2 (Hessian) propagation through the dMod function stack.
#
# Section order follows the composition chain:
#   Xs -> Y -> Pexpl -> Pequil -> (Y * Xs) -> res -> normL2 -> constraintL2
# ============================================================================


# ---- Xs -------------------------------------------------------------------

test_that("Xs.CppODE deriv2 reproduces linear-decay analytical Hessian", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
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

  xfn <- Xs(m, optionsOde = list(atol = 1e-12, rtol = 1e-12),
            optionsSens = list(atol = 1e-12, rtol = 1e-12))
  times <- c(0.0, 0.5, 1.0)
  pars <- c(x = 2.0, k = 0.7)
  r1 <- xfn(times, pars, deriv = TRUE, deriv2 = FALSE)[[1]]
  r2 <- xfn(times, pars, deriv = TRUE, deriv2 = TRUE )[[1]]

  expect_null(attr(r1, "deriv2"))
  expect_equal(dim(attr(r2, "deriv2")), c(length(times), 1L, 2L, 2L))
  expect_equal(attr(r1, "deriv"), attr(r2, "deriv"), tolerance = 1e-9)
})

test_that("Xs.deSolve refuses deriv2 = TRUE", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  f <- c(x = "-k*x")
  m <- odemodel(f, modelname = paste0("decay_des_", as.integer(Sys.time())),
                solver = "deSolve", verbose = FALSE, compile = FALSE)
  xfn <- Xs(m); compile(xfn)
  expect_error(xfn(c(0, 1), c(x = 1, k = 0.5), deriv2 = TRUE),
               "Xs.deSolve")
})


# ---- Y --------------------------------------------------------------------

test_that("Y deriv2 (AD) reproduces analytical observation Hessian", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  gfn <- Y(c(y = "a*x^2 + b*x"), states = "x", parameters = c("a", "b"),
           modelname = paste0("y_d2_", as.integer(Sys.time())),
           compile = TRUE, deriv2 = TRUE, derivMode = "dual",
           attach.input = FALSE)

  times <- c(0.0, 0.5, 1.0)
  out <- cbind(time = times, x = c(1.0, 0.7, 0.4))
  class(out) <- c("prdframe", "matrix")
  pars <- as.parvec(c(a = 0.3, b = -0.5))

  res <- gfn(out, pars, deriv = TRUE, deriv2 = TRUE)[[1]]

  a <- pars["a"]; b <- pars["b"]
  xv <- out[, "x"]
  ref_dy_dx <- 2 * a * xv + b
  ref_dy_da <- xv^2
  ref_dy_db <- xv

  d1 <- attr(res, "deriv")
  d2 <- attr(res, "deriv2")
  expect_equal(d1[, "y", "x"], unname(ref_dy_dx), tolerance = 1e-10)
  expect_equal(d1[, "y", "a"], unname(ref_dy_da), tolerance = 1e-10)
  expect_equal(d1[, "y", "b"], unname(ref_dy_db), tolerance = 1e-10)

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


# ---- Pexpl ----------------------------------------------------------------

test_that("Pexpl deriv2 (AD) reproduces analytical Hessian", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  trafo <- c(a = "exp(la)", b = "la^2 + lb", c = "la*lb")
  p <- Pexpl(trafo, parameters = NULL,
             modelname = paste0("ad_pexpl_d2_", as.integer(Sys.time())),
             compile = TRUE, deriv2 = TRUE, derivMode = "dual")

  pars <- c(la = 0.3, lb = 0.5)
  pinner <- p(pars, deriv = TRUE, deriv2 = TRUE)[[1]]

  la <- pars["la"]; lb <- pars["lb"]
  J_ref <- matrix(c(exp(la), 2 * la, lb,
                    0,       1,      la),
                  3, 2,
                  dimnames = list(c("a", "b", "c"), c("la", "lb")))
  H_ref <- array(0, c(3, 2, 2),
                 dimnames = list(c("a", "b", "c"), c("la", "lb"), c("la", "lb")))
  H_ref["a", "la", "la"] <- exp(la)
  H_ref["b", "la", "la"] <- 2
  H_ref["c", "la", "lb"] <- 1
  H_ref["c", "lb", "la"] <- 1

  expect_equal(attr(pinner, "deriv"),
               J_ref[rownames(attr(pinner, "deriv")), , drop = FALSE],
               tolerance = 1e-10)
  expect_equal(attr(pinner, "deriv2")[c("a", "b", "c"), , ], H_ref,
               tolerance = 1e-10)
})

test_that("Pexpl deriv2 (symbolic) reproduces analytical Hessian", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  trafo <- c(a = "exp(la)", b = "la^2 + lb", c = "la*lb")
  p <- Pexpl(trafo, parameters = NULL,
             modelname = paste0("sym_pexpl_d2_", as.integer(Sys.time())),
             compile = TRUE, deriv2 = TRUE, derivMode = "symbolic")

  pars <- c(la = -0.4, lb = 0.7)
  pinner <- p(pars, deriv = TRUE, deriv2 = TRUE)[[1]]

  la <- pars["la"]; lb <- pars["lb"]
  H_ref <- array(0, c(3, 2, 2),
                 dimnames = list(c("a", "b", "c"), c("la", "lb"), c("la", "lb")))
  H_ref["a", "la", "la"] <- exp(la)
  H_ref["b", "la", "la"] <- 2
  H_ref["c", "la", "lb"] <- 1
  H_ref["c", "lb", "la"] <- 1

  expect_equal(attr(pinner, "deriv2")[c("a", "b", "c"), , ], H_ref,
               tolerance = 1e-10)
})

test_that("Pexpl deriv2 (AD) handles identity pass-through entries", {
  # Regression: with `parameters` supplied, identity entries (e.g. la = "la")
  # are appended to the trafo; the compiled AD2 entry must propagate dual2nd
  # values through them without corrupting derivatives.
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  trafo <- c(a = "exp(la)", b = "la^2 + lb")
  p <- Pexpl(trafo, parameters = c("la", "lb"),
             modelname = paste0("id_pexpl_d2_", as.integer(Sys.time())),
             compile = TRUE, deriv2 = TRUE, derivMode = "dual")

  pars <- c(la = 0.3, lb = 0.5)
  pinner <- p(pars, deriv = TRUE, deriv2 = TRUE)[[1]]

  la <- pars["la"]
  d2 <- attr(pinner, "deriv2")
  expect_equal(unname(d2["a", "la", "la"]), unname(exp(la)), tolerance = 1e-10)
  expect_equal(unname(d2["b", "la", "la"]), 2,                tolerance = 1e-10)
  expect_equal(unname(d2["la", "la", "la"]), 0,               tolerance = 1e-12)
  expect_equal(unname(d2["lb", "lb", "lb"]), 0,               tolerance = 1e-12)
  expect_equal(unname(d2["a", "lb", "lb"]), 0,                tolerance = 1e-12)
})

test_that("[.parvec and c.parvec propagate deriv2 attributes", {
  vals <- c(a = 1.0, b = 2.0, c = 3.0)
  J <- matrix(c(1, 0, 2, 1, 0, 3), 3, 2,
              dimnames = list(c("a", "b", "c"), c("p1", "p2")))
  H <- array(seq_len(3 * 2 * 2), c(3, 2, 2),
             dimnames = list(c("a", "b", "c"), c("p1", "p2"), c("p1", "p2")))
  pv <- as.parvec(vals, deriv = J, deriv2 = H)

  sub <- pv[c("a", "c")]
  expect_equal(rownames(attr(sub, "deriv")), c("a", "c"))
  expect_equal(dimnames(attr(sub, "deriv2"))[[1]], c("a", "c"))
  expect_equal(dim(attr(sub, "deriv2")), c(2, 2, 2))
  expect_equal(attr(sub, "deriv2")[1, , ], H["a", , ], tolerance = 1e-12)
  expect_equal(attr(sub, "deriv2")[2, , ], H["c", , ], tolerance = 1e-12)

  cat <- c(pv[c("a")], pv[c("b", "c")])
  expect_equal(dim(attr(cat, "deriv2")), c(3, 2, 2))
  expect_equal(attr(cat, "deriv2")[c("a", "b", "c"), , ], H,
               tolerance = 1e-12)
})

test_that("Pexpl(deriv2 = FALSE) refuses deriv2 = TRUE at call time", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  trafo <- c(a = "exp(la)")
  p <- Pexpl(trafo, parameters = NULL,
             modelname = paste0("nod2_pexpl_", as.integer(Sys.time())),
             compile = TRUE, deriv2 = FALSE, derivMode = "dual")
  expect_error(p(c(la = 0.1), deriv2 = TRUE),
               "deriv2 = FALSE")
})


# ---- Pequil ---------------------------------------------------------------

test_that("Pequil deriv2 reproduces analytical equilibrium Hessian", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  trafo <- c(x = "-k * x + s")  # x* = s/k
  p <- Pequil(trafo, parameters = c("k", "s"),
              modelname = paste0("equil_lin_d2_", as.integer(Sys.time())),
              compile = TRUE, deriv2 = TRUE, attach.input = FALSE,
              verbose = FALSE)

  pars <- c(k = 0.5, s = 2.0, x = 1.0)
  pinner <- p(pars, deriv = TRUE, deriv2 = TRUE)[[1]]

  k <- pars["k"]; s <- pars["s"]
  expect_equal(as.numeric(pinner), as.numeric(s / k), tolerance = 1e-6)

  J_ref <- matrix(c(-s / k^2, 1 / k), nrow = 1,
                  dimnames = list("x", c("k", "s")))
  expect_equal(attr(pinner, "deriv")[, c("k", "s"), drop = FALSE],
               J_ref, tolerance = 1e-6)

  H_ref <- array(0, c(1, 2, 2),
                 dimnames = list("x", c("k", "s"), c("k", "s")))
  H_ref["x", "k", "k"] <- 2 * s / k^3
  H_ref["x", "k", "s"] <- -1 / k^2
  H_ref["x", "s", "k"] <- -1 / k^2
  H_ref["x", "s", "s"] <- 0
  expect_equal(attr(pinner, "deriv2")[, c("k", "s"), c("k", "s"), drop = FALSE],
               H_ref, tolerance = 1e-6)
})

test_that("Pequil(deriv2 = TRUE) emits <m>, <m>_s, <m>_s2 and dispatches", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  td <- file.path(tempdir(), paste0("equil_triple_", as.integer(Sys.time())))
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  setwd(td)

  trafo <- c(x = "-k*x + s")
  base <- "equil_triple"
  p <- Pequil(trafo, parameters = c("k", "s"),
              modelname = base,
              compile = TRUE, deriv2 = TRUE, attach.input = FALSE,
              verbose = FALSE)
  expect_true(file.exists(paste0(base,     ".cpp")))
  expect_true(file.exists(paste0(base, "_s.cpp")))
  expect_true(file.exists(paste0(base, "_s2.cpp")))

  pars <- c(k = 0.5, s = 2.0, x = 1.0)
  r1 <- p(pars, deriv = TRUE, deriv2 = FALSE)[[1]]
  r2 <- p(pars, deriv = TRUE, deriv2 = TRUE)[[1]]

  expect_null(attr(r1, "deriv2"))
  expect_false(is.null(attr(r2, "deriv2")))
  expect_equal(attr(r1, "deriv"), attr(r2, "deriv"), tolerance = 1e-6)
})

test_that("Pequil(deriv2 = FALSE) refuses deriv2 = TRUE at call time", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  trafo <- c(x = "-k*x + s")
  p <- Pequil(trafo, parameters = c("k", "s"),
              modelname = paste0("equil_nod2_", as.integer(Sys.time())),
              compile = TRUE, deriv2 = FALSE, attach.input = FALSE,
              verbose = FALSE)
  expect_error(p(c(k = 0.5, s = 2, x = 1), deriv2 = TRUE),
               "deriv2 = TRUE")
})


# ---- Composition ----------------------------------------------------------

test_that("(Y * Xs)(times, pars, deriv2 = TRUE) chain-rule matches analytical", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
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


# ---- res ------------------------------------------------------------------

test_that("res() reduces 4D deriv2 attribute to 3D [n_residuals, p, p]", {
  times <- c(0.0, 0.5, 1.0)
  obs   <- c("y1", "y2")
  pars  <- c("p1", "p2")

  out <- cbind(time = times,
               y1 = c(1.0, 2.0, 3.0),
               y2 = c(0.1, 0.2, 0.3))
  d2 <- array(NA_real_, c(length(times), length(obs), length(pars), length(pars)),
              dimnames = list(NULL, obs, pars, pars))
  for (i in seq_along(times))
    for (j in seq_along(obs))
      for (k1 in seq_along(pars))
        for (k2 in seq_along(pars))
          d2[i, j, k1, k2] <- 1000 * i + 100 * j + 10 * k1 + k2
  attr(out, "deriv2") <- d2
  attr(out, "deriv") <- array(0, c(length(times), length(obs), length(pars)),
                              dimnames = list(NULL, obs, pars))

  data_df <- data.frame(name = c("y1", "y2", "y1"),
                        time = c(0.0, 0.5, 1.0),
                        value = c(1.05, 0.18, 2.95),
                        sigma = c(0.1, 0.1, 0.1),
                        lloq = -Inf)

  ro <- res(data_df, out)
  rd2 <- attr(ro, "deriv2")
  expect_equal(dim(rd2), c(3, 2, 2))
  expect_equal(rd2[1, , ], matrix(c(1111, 1121, 1112, 1122), 2, 2,
                                  dimnames = list(pars, pars)))
  expect_equal(rd2[2, , ], matrix(c(2211, 2221, 2212, 2222), 2, 2,
                                  dimnames = list(pars, pars)))
  expect_equal(rd2[3, , ], matrix(c(3111, 3121, 3112, 3122), 2, 2,
                                  dimnames = list(pars, pars)))
})


# ---- normL2 ---------------------------------------------------------------

test_that("normL2 gradient is identical for deriv2 = FALSE and deriv2 = TRUE", {
  # Regression for the [.parvec / c.parvec deriv2-drop bug. With deriv2 = TRUE,
  # the upstream Hessian seed must reach Xs.CppODE so the _s2 integration
  # produces the same first-order sensitivities as the _s integration. If
  # subsetting drops attr(., "deriv2") the gradient diverges silently.
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
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
  withr::local_dir(tempdir())
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
  expect_true(is.numeric(res_gn$value))
  expect_equal(res_gn$hessian, t(res_gn$hessian), tolerance = 1e-12)
})

test_that("normL2(deriv2 = TRUE) adds residual times d^2 pred / sigma^2", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  f <- c(x = "-k*x")
  m <- odemodel(f, modelname = paste0("nl_d2_ex_", as.integer(Sys.time())),
                solver = "CppODE", deriv2 = TRUE, nStack = 4L, verbose = FALSE)
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
  expect_equal(res_ex$value, res_gn$value, tolerance = 1e-8)
  expect_equal(res_ex$gradient, res_gn$gradient, tolerance = 1e-8)

  # Analytical Hessian addition: 2/sigma^2 * sum_i r_i * d^2 y_i / dtheta^2.
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


# ---- constraintL2 ---------------------------------------------------------

test_that("constraintL2 deriv2 adds gi . dP2 chain term after Pexpl", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  # Pexpl: a = exp(la). constraintL2(mu = mu_a, sigma = s) on `a`,
  # composed via attr(p, "deriv") and attr(p, "deriv2") from Pexpl.
  pfn <- Pexpl(c(a = "exp(la)"), parameters = NULL,
               modelname = paste0("c2_pexpl_", as.integer(Sys.time())),
               compile = TRUE, deriv2 = TRUE, derivMode = "dual")

  mu <- c(a = 1.0); sg <- 0.5
  cfn <- constraintL2(mu = mu, sigma = sg)

  pars <- c(la = 0.3)
  pinner <- pfn(pars, deriv = TRUE, deriv2 = TRUE)[[1]]

  res_gn <- cfn(pinner, deriv = TRUE, deriv2 = FALSE)
  res_ex <- cfn(pinner, deriv = TRUE, deriv2 = TRUE)

  # L = ((exp(la) - mu) / sg)^2 with closed-form gradient/Hessian.
  la <- pars["la"]
  ref_grad <- unname(2 * (exp(la) - mu) * exp(la) / sg^2)
  ref_hess_GN <- unname(2 * exp(la)^2 / sg^2)
  ref_hess_exact <- ref_hess_GN + unname(2 * (exp(la) - mu) * exp(la) / sg^2)

  expect_equal(unname(res_gn$value), unname(((exp(la) - mu) / sg)^2),
               tolerance = 1e-10)
  expect_equal(unname(res_gn$gradient["la"]), ref_grad, tolerance = 1e-10)
  expect_equal(unname(res_gn$hessian["la", "la"]), ref_hess_GN,
               tolerance = 1e-10)
  expect_equal(unname(res_ex$hessian["la", "la"]), ref_hess_exact,
               tolerance = 1e-10)
})
