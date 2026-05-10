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

  # 1st-order must agree between the cheap and the expensive path.
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
