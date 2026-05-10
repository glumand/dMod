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
  # Regression for the CppODE dual2nd identity-output bug. When the caller
  # supplies `parameters` to Pexpl, identity entries (e.g. la = "la") are
  # appended to the trafo. The compiled AD2 entry must propagate dual2nd
  # values through these without corrupting derivatives.
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
  # Closed-form: only d^2 a/dla^2 = exp(la) and d^2 b/dla^2 = 2 are non-zero.
  # Identity entries la, lb contribute zero throughout.
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

  # c.parvec preserves deriv2 along axis 1.
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
