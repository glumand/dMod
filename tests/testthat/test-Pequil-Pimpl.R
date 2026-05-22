# Behavioral tests for Pequil() (steady-state parameter transformation)
# and Pimpl() (implicit-equation parameter transformation).
#
# Pequil and Pimpl turn structural equations into parfn() objects whose
# value is the solution of an equilibrium / implicit problem and whose
# Jacobian comes from the implicit function theorem.
#
# Pequil(deriv2)'s Hessian on the linear x* = s/k system is covered by
# test-deriv2.R. This file adds a second model (production-decay
# A* = k_in / k_out) plus a Pimpl quadratic-root test, both validated
# against closed-form solutions and their analytical Jacobians.
#
# The resetWarmStarts() cache-management tests live at the bottom (also
# specific to Pequil/Pimpl, which are the only parfns that warm-start).

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

  # IFT: x*(a) = a, so dx*/da = +1. (The legacy fn used to recycle a single
  # Jacobian entry across all columns, producing the wrong sign here; the
  # current implementation pins the sign down to the standard IFT result.)
  J <- attr(out, "deriv")
  expect_equal(unname(J["x", "a"]), 1, tolerance = 1e-6)
})


## ---- Pimpl: IFT second-order on a 1D quadratic SS ---------------------

test_that("Pimpl(deriv2) matches FD on a 1D mass-action steady state", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # A + B <-> C, with totA, totB substituted in: solve k1*(totA-C)*(totB-C) - km*C = 0
  trafo <- c(C = "k1*(totA - C)*(totB - C) - km*C")
  pf <- Pimpl(trafo, parameters = c("k1","km","totA","totB"),
              positive = TRUE, deriv2 = TRUE,
              modelname = paste0("test_Pimpl_d2_", as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsNleqslv = list(nstarts = 5, ftol = 1e-12, xtol = 1e-12))

  p0 <- c(k1 = 1, km = 0.5, totA = 2, totB = 3, C = 0.5)
  inputs <- c("k1","km","totA","totB")
  v <- pf(p0, deriv2 = TRUE)[[1]]
  J_an <- attr(v, "deriv")
  H_an <- attr(v, "deriv2")

  # Closed form at C* = 1.5: dxdp = -dfdp/dfdx
  expect_equal(unname(J_an["C", inputs]),
               c(0.3, -0.6, 0.6, 0.2), tolerance = 1e-6)

  # FD-check the Hessian against the analytical Jacobian
  jac_at <- function(par_in) {
    pp <- p0; pp[names(par_in)] <- par_in
    attr(pf(pp, deriv2 = FALSE)[[1]], "deriv")["C", inputs, drop = FALSE]
  }
  H_fd <- matrix(0, length(inputs), length(inputs), dimnames = list(inputs, inputs))
  h <- 1e-5
  for (j in seq_along(inputs)) {
    pj <- p0[inputs]; pj[j] <- pj[j] + h
    pn <- p0[inputs]; pn[j] <- pn[j] - h
    H_fd[, j] <- (jac_at(pj) - jac_at(pn)) / (2*h)
  }
  H_fd <- 0.5 * (H_fd + t(H_fd))
  expect_equal(unname(H_an["C", inputs, inputs]), unname(H_fd),
               tolerance = 1e-5)
})


## ---- Pimpl: IFT second-order with coupled n_dep = 2 -------------------

test_that("Pimpl(deriv2) matches FD on a 2-state coupled SS", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  trafo <- c(x1 = "a*x1 - b*x2 - 1",
             x2 = "x1*x2 - c")
  pf <- Pimpl(trafo, parameters = c("a","b","c"),
              positive = TRUE, deriv2 = TRUE,
              modelname = paste0("test_Pimpl_d2_2_", as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsNleqslv = list(nstarts = 5, ftol = 1e-12, xtol = 1e-12))

  p0 <- c(a = 2, b = 0.5, c = 1, x1 = 0.5, x2 = 0.5)
  inputs <- c("a","b","c")
  v <- pf(p0, deriv2 = TRUE)[[1]]
  J_an <- attr(v, "deriv")
  H_an <- attr(v, "deriv2")

  jac_at <- function(par_in) {
    pp <- p0; pp[names(par_in)] <- par_in
    attr(pf(pp, deriv2 = FALSE)[[1]], "deriv")[c("x1","x2"), inputs, drop = FALSE]
  }
  H_fd <- array(0, c(2, length(inputs), length(inputs)),
                dimnames = list(c("x1","x2"), inputs, inputs))
  h <- 1e-5
  for (j in seq_along(inputs)) {
    pj <- p0[inputs]; pj[j] <- pj[j] + h
    pn <- p0[inputs]; pn[j] <- pn[j] - h
    H_fd[, , j] <- (jac_at(pj) - jac_at(pn)) / (2*h)
  }
  H_fd <- 0.5 * (H_fd + aperm(H_fd, c(1,3,2)))
  expect_equal(unname(H_an[c("x1","x2"), inputs, inputs]),
               unname(H_fd), tolerance = 1e-4)
})


## ---- Pimpl: IFT second-order through CQ-elim reconstruction -----------

test_that("Pimpl(deriv2) propagates Hessian through CQ-eliminated species", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # A <-> B; CQ A + B = total_1; A gets eliminated, B stays dependent
  el <- eqnlist()
  el <- addReaction(el, "A", "B", "k*A")
  el <- addReaction(el, "B", "A", "km*B")

  pf <- Pimpl(el, parameters = c("k","km"),
              positive = FALSE, deriv2 = TRUE,
              modelname = paste0("test_Pimpl_d2_cq_", as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsNleqslv = list(nstarts = 1, ftol = 1e-12, xtol = 1e-12))

  p0 <- c(k = 2, km = 0.5, total_1 = 3, A = 0.5, B = 0.5)
  inputs <- c("k","km","total_1")
  v <- pf(p0, deriv2 = TRUE)[[1]]
  J_an <- attr(v, "deriv")
  H_an <- attr(v, "deriv2")

  # Closed form: B* = k*total_1 / (k+km); A* = total_1 - B*
  Bs <- p0[["k"]] * p0[["total_1"]] / (p0[["k"]] + p0[["km"]])
  expect_equal(as.numeric(v["B"]), Bs, tolerance = 1e-6)
  expect_equal(as.numeric(v["A"]), p0[["total_1"]] - Bs, tolerance = 1e-6)

  # A and B's Hessians are linked by A = total_1 - B  =>  d^2 A = -d^2 B
  expect_equal(unname(H_an["A", inputs, inputs]),
               -unname(H_an["B", inputs, inputs]), tolerance = 1e-10)

  # FD-check on B's analytical Hessian
  jac_at <- function(par_in) {
    pp <- p0; pp[names(par_in)] <- par_in
    attr(pf(pp, deriv2 = FALSE)[[1]], "deriv")["B", inputs, drop = FALSE]
  }
  H_fd <- matrix(0, length(inputs), length(inputs), dimnames = list(inputs, inputs))
  h <- 1e-5
  for (j in seq_along(inputs)) {
    pj <- p0[inputs]; pj[j] <- pj[j] + h
    pn <- p0[inputs]; pn[j] <- pn[j] - h
    H_fd[, j] <- (jac_at(pj) - jac_at(pn)) / (2*h)
  }
  H_fd <- 0.5 * (H_fd + t(H_fd))
  expect_equal(unname(H_an["B", inputs, inputs]), unname(H_fd),
               tolerance = 1e-5)
})


## ---- Pimpl: singular df/dx triggers a diagnostic error ----------------

test_that("Pimpl uses the pseudoinverse when df/dx is rank-deficient", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Two-state system with a hidden conservation law that .detect_and_substitute_cq
  # cannot pick up because there is no `eqnlist` smatrix:
  #   f1 = x1 + x2 - s
  #   f2 = x1 + x2 - s    (same as f1 -> df/dx has rank 1)
  # The constraint manifold is the line x1 + x2 = s. The IFT does not apply
  # uniquely, but the Moore-Penrose pseudoinverse gives the minimum-norm
  # sensitivity (movement on the manifold, perpendicular to the null space)
  # which downstream callers can still use; Pimpl warns and proceeds.
  trafo <- c(x1 = "x1 + x2 - s",
             x2 = "x1 + x2 - s")
  pf <- Pimpl(trafo, parameters = "s",
              positive = FALSE, deriv2 = FALSE,
              modelname = paste0("test_Pimpl_sing_", as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsNleqslv = list(nstarts = 1, ftol = 1e-12, xtol = 1e-12))

  p0 <- c(s = 1, x1 = 0.5, x2 = 0.5)
  expect_warning(out <- pf(p0, deriv = TRUE)[[1]], "rank-deficient")
  J <- attr(out, "deriv")
  expect_true(is.matrix(J))
  # Pseudoinverse projects equal-weighted sensitivity onto x1 + x2: the
  # constraint says dx1/ds + dx2/ds = 1, with the minimum-norm split.
  expect_equal(sum(J[c("x1", "x2"), "s"]), 1, tolerance = 1e-6)
})


## ---- resetWarmStarts: clear Pequil / Pimpl cache ----------------------

# Both Pequil and Pimpl cache the previous root (+ sensitivities) for
# warm-starting the next call. After a structural change or a basin
# jump, that cache becomes a liability -- resetWarmStarts() drops it.

test_that("resetWarmStarts clears Pequil's cache by name", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  pf <- Pequil(c(A = "k_in - k_out * A"),
               parameters = c("k_in", "k_out"),
               modelname = paste0("test_rws_pequil_", as.integer(Sys.time())),
               compile = TRUE, verbose = FALSE, attach.input = FALSE)
  pars <- c(k_in = 1.5, k_out = 0.3, A = 0.1)
  pf(pars)

  reset_env <- environment(attr(pf, "resetWarmStart"))$cache_ref
  expect_false(is.null(reset_env$yini))

  labels <- resetWarmStarts(pf, verbose = FALSE)
  expect_length(labels, 1L)
  expect_match(labels, "^Pequil\\(")
  expect_null(reset_env$yini)
  expect_null(reset_env$last_hash)
})


test_that("resetWarmStarts clears Pimpl's cache by name", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  trafo <- c(x = "x - a")
  pf <- Pimpl(trafo, parameters = "a",
              modelname = paste0("test_rws_pimpl_", as.integer(Sys.time())),
              compile = TRUE, verbose = FALSE)
  pf(c(a = 1, x = 0.5))

  reset_env <- environment(attr(pf, "resetWarmStart"))$cache_ref
  expect_false(is.null(reset_env$guess))

  labels <- resetWarmStarts(pf, verbose = FALSE)
  expect_length(labels, 1L)
  expect_match(labels, "^Pimpl\\(")
  expect_null(reset_env$guess)
  expect_null(reset_env$failed_pv_hash)
})


test_that("resetWarmStarts walks into composed functions", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  p1 <- Pequil(c(A = "k1 - kdg1 * A"),
               parameters = c("k1", "kdg1"),
               modelname = paste0("test_rws_p1_", as.integer(Sys.time())),
               compile = TRUE, verbose = FALSE, attach.input = FALSE)
  p2 <- Pequil(c(B = "k2 - kdg2 * B"),
               parameters = c("k2", "kdg2"),
               modelname = paste0("test_rws_p2_", as.integer(Sys.time())),
               compile = TRUE, verbose = FALSE, attach.input = FALSE)

  p1(c(k1 = 1, kdg1 = 0.5, A = 0.1))
  p2(c(k2 = 2, kdg2 = 0.4, B = 0.2))

  ref1 <- environment(attr(p1, "resetWarmStart"))$cache_ref
  ref2 <- environment(attr(p2, "resetWarmStart"))$cache_ref
  expect_false(is.null(ref1$yini))
  expect_false(is.null(ref2$yini))

  wrap <- function(pars) list(p1 = p1(pars), p2 = p2(pars))
  labels <- resetWarmStarts(wrap, verbose = FALSE)
  expect_setequal(sub("\\(.*", "", labels), c("Pequil", "Pequil"))
  expect_length(labels, 2L)
  expect_null(ref1$yini)
  expect_null(ref2$yini)
})


test_that("resetWarmStarts is idempotent and silent when nothing to reset", {
  f <- function(x) x + 1
  labels <- resetWarmStarts(f, verbose = FALSE)
  expect_length(labels, 0L)
})


test_that("resetWarmStarts rejects non-function inputs", {
  expect_error(resetWarmStarts(1L), "must be a function")
  expect_error(resetWarmStarts(NULL), "must be a function")
})
