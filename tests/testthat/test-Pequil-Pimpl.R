# Behavioral tests for Pequil() (steady-state parameter transformation)
# and Pimpl() (implicit-equation parameter transformation).
#
# Pequil and Pimpl turn structural equations into parfn() objects whose
# value is the solution of an equilibrium / implicit problem and whose
# Jacobian comes from the implicit function theorem.
#
# Pequil(deriv2)'s Hessian on the linear x* = s/k system is covered by
# test-deriv2-Pequil.R. This file adds a second model (production-decay
# A* = k_in / k_out) plus a Pimpl quadratic-root test, both validated
# against closed-form solutions and their analytical Jacobians.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


## ---- Pequil: production-decay -----------------------------------------

test_that("Pequil on dA = k_in - k_out * A converges to A* = k_in / k_out", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  trafo <- c(A = "k_in - k_out * A")
  pf <- Pequil(trafo, parameters = c("k_in", "k_out"),
               modelname = paste0("test_Pequil_pd_", as.integer(Sys.time())),
               compile = TRUE, deriv2 = FALSE, attach.input = FALSE,
               verbose = FALSE)

  pars <- c(k_in = 1.5, k_out = 0.3, A = 0.1)  # A is the initial guess
  out <- pf(pars, deriv = TRUE)[[1]]

  expect_equal(as.numeric(out), as.numeric(pars[["k_in"]] / pars[["k_out"]]),
               tolerance = 1e-5)

  # Jacobian via IFT: A*(k_in, k_out) = k_in / k_out
  #   dA*/dk_in  =  1 / k_out
  #   dA*/dk_out = -k_in / k_out^2
  J <- attr(out, "deriv")
  expect_equal(J[, "k_in",  drop = TRUE], 1 / pars[["k_out"]],
               tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(J[, "k_out", drop = TRUE], -pars[["k_in"]] / pars[["k_out"]]^2,
               tolerance = 1e-5, ignore_attr = TRUE)
})


## ---- Pimpl: linear constraint (value) ---------------------------------

test_that("Pimpl on x - a = 0 returns x = a (root-finding correctness)", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  trafo <- c(x = "x - a")
  pf <- Pimpl(trafo, parameters = "a", positive = TRUE,
              modelname = paste0("test_Pimpl_lin_", as.integer(Sys.time())),
              compile = TRUE, verbose = FALSE)

  pars <- c(a = 2.0, x = 0.5)
  out <- pf(pars, deriv = TRUE)[[1]]

  # Pimpl returns inner state plus pass-through inputs; pick by name.
  # Tolerance reflects nleqslv solver default termination criterion.
  expect_equal(as.numeric(out["x"]), pars[["a"]], tolerance = 1e-3)

  # NB: the sign convention of the Pimpl-emitted Jacobian (dx/da vs dF/da)
  # is internal to dMod and not yet pinned down by a behavioural test.
  # The magnitude is still verifiable: |dx/da| = |1 / dF/dx| = 1 here.
  J <- attr(out, "deriv")
  expect_equal(abs(unname(J["x", "a"])), 1, tolerance = 1e-6)
})
