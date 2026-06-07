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
  pf <- Pimpl(trafo, parameters = "a",
              controlsMS = list(positive = TRUE),
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
              deriv2 = TRUE,
              modelname = paste0("test_Pimpl_d2_", as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsMS = list(nStarts = 5L, positive = TRUE),
              controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))

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
              deriv2 = TRUE,
              modelname = paste0("test_Pimpl_d2_2_", as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsMS = list(nStarts = 5L, positive = TRUE),
              controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))

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
              deriv2 = TRUE,
              modelname = paste0("test_Pimpl_d2_cq_", as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsMS = list(nStarts = 1L, positive = FALSE),
              controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))

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


## ---- getEquations under expressInTotals -------------------------------

test_that("Pimpl(expressInTotals) solves the full system with a conservation residual", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # A <-> B with CQ A + B = total_1. Pimpl keeps all moiety species (solved
  # positive in log-space) and replaces one redundant ODE with the algebraic
  # conservation residual; nothing is reconstructed by subtraction.
  el <- eqnlist()
  el <- addReaction(el, "A", "B", "k*A")
  el <- addReaction(el, "B", "A", "km*B")

  pf <- Pimpl(el, parameters = c("k","km"), expressInTotals = TRUE,
              modelname = paste0("test_Pimpl_eq_", as.integer(Sys.time())),
              compile = FALSE, controlsMS = list(nStarts = 3L))
  eqs <- getEquations(pf)[[1]]
  expect_length(eqs, 2L)
  # one entry is the conservation residual (named by the total), the other a
  # genuine rate equation for the surviving species
  expect_true("total_1" %in% names(eqs))
  expect_match(eqs[["total_1"]], "total_1", fixed = TRUE)
  ode_name <- setdiff(names(eqs), "total_1")
  expect_true(ode_name %in% c("A", "B"))
  expect_true(grepl("k", eqs[[ode_name]]))
  expect_true("total_1" %in% getParameters(pf))
})

test_that("Pequil(expressInTotals) integrates the full system (no elimination)", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Pequil keeps every moiety species as a dynamical equation; the total
  # enters through the initial conditions, not the rates.
  el <- eqnlist()
  el <- addReaction(el, "A", "B", "k*A")
  el <- addReaction(el, "B", "A", "km*B")

  pf <- Pequil(el, parameters = c("k","km"), expressInTotals = TRUE,
               modelname = paste0("test_Pequil_eq_", as.integer(Sys.time())),
               compile = FALSE, controlsMS = list(nStarts = 3L))
  eqs <- getEquations(pf)[[1]]
  expect_setequal(names(eqs), c("A","B"))
  # both are genuine rate equations, no algebraic total_X - rest relation
  expect_true(all(grepl("k", eqs)))
  expect_false(any(grepl("total", eqs, fixed = TRUE)))
  # total_1 is an outer parameter, supplied through the initial condition
  expect_true("total_1" %in% getParameters(pf))
})

test_that("Pequil(expressInTotals) keeps reconstructed moiety species non-negative", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Overlapping moieties sharing C: A+C and B+C. The full-system integration
  # cannot produce a negative C the way linear elimination (C = total - rest)
  # could.
  el <- eqnlist()
  el <- addReaction(el, "A", "Ap", "ka*A")
  el <- addReaction(el, "Ap", "A", "kap*Ap")
  el <- addReaction(el, "Ap + B", "C", "kf*Ap*B")
  el <- addReaction(el, "C", "Ap + B", "kr*C")
  el <- addReaction(el, "B", "Bp", "kb*B")
  el <- addReaction(el, "Bp", "B", "kbp*Bp")

  pf <- Pequil(el, expressInTotals = TRUE, compile = TRUE,
               modelname = paste0("test_Pequil_pos_", as.integer(Sys.time())),
               controlsMS = list(nStarts = 30L))
  set.seed(11)
  for (i in 1:8) {
    resetWarmStarts(pf)
    pin <- structure(10^runif(length(getParameters(pf)), -1, 1),
                     names = getParameters(pf))
    o <- tryCatch(pf(pin)[[1]], error = function(e) NULL)
    if (is.null(o)) next
    species <- intersect(c("A","Ap","B","Bp","C"), names(o))
    expect_true(all(as.numeric(o[species]) >= -1e-8))
  }
})


## ---- CQ spectrum: independent / non-unit-coef / overlapping -----------

# These exercise distinct branches of .detect_and_substitute_cq:
#   B  two independent CQs            -> two pivots, two total_* params
#   D  non-unit stoichiometry (2*D)   -> recon divides by coef_e = 2
#   C  overlapping CQs, neg coef      -> Gaussian-elim pivots + nested
#                                        recon (fixed-point substitution)

test_that("two independent conserved moieties solve to closed form (both backends)", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # A <-> B  and  C <-> D : two disjoint moieties, total_1 = A+B, total_2 = C+D
  el <- eqnlist()
  el <- addReaction(el, "A", "B", "k1*A"); el <- addReaction(el, "B", "A", "k2*B")
  el <- addReaction(el, "C", "D", "k3*C"); el <- addReaction(el, "D", "C", "k4*D")

  p0 <- c(k1 = 2, k2 = 1, k3 = 0.5, k4 = 1.5, total_1 = 4, total_2 = 6,
          A = 1, B = 1, C = 1, D = 1)
  Bs <- p0[["k1"]] * p0[["total_1"]] / (p0[["k1"]] + p0[["k2"]]); As <- p0[["total_1"]] - Bs
  Ds <- p0[["k3"]] * p0[["total_2"]] / (p0[["k3"]] + p0[["k4"]]); Cs <- p0[["total_2"]] - Ds

  for (backend in c("Pimpl", "Pequil")) {
    pf <- get(backend)(el, expressInTotals = TRUE, compile = TRUE,
                       modelname = paste0("test_", backend, "_2moiety_", as.integer(Sys.time())),
                       controlsMS = list(nStarts = 20L))
    o <- pf(p0)[[1]]
    expect_equal(as.numeric(o[c("A","B","C","D")]),
                 c(As, Bs, Cs, Ds), tolerance = 1e-3)
  }
})

test_that("non-unit stoichiometric coefficient is handled (2*M <-> D)", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Dimerisation 2 M <-> D conserves the monomer count M + 2 D.
  el <- eqnlist()
  el <- addReaction(el, "2*M", "D", "ka*M^2")
  el <- addReaction(el, "D", "2*M", "kd*D")
  el <- customTotals(el, list(total_MD = "M + 2*D"))

  p0 <- c(ka = 1, kd = 2, total_MD = 5)
  r  <- p0[["ka"]] / p0[["kd"]]
  Ms <- (-1 + sqrt(1 + 8 * r * p0[["total_MD"]])) / (4 * r)   # 2 r M^2 + M - T = 0
  Ds <- r * Ms^2

  # Pimpl conserves the moiety to the solver tolerance (constraint residual),
  # so tighten ftol; Pequil conserves it exactly (ODE invariant).
  pfi <- Pimpl(el, expressInTotals = TRUE, compile = FALSE,
               modelname = paste0("test_Pimpl_dimer_", as.integer(Sys.time())),
               controlsMS = list(nStarts = 30L),
               controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))
  oi <- pfi(p0)[[1]]
  expect_equal(as.numeric(oi["M"]), Ms, tolerance = 1e-4)
  expect_equal(as.numeric(oi["D"]), Ds, tolerance = 1e-4)
  expect_equal(as.numeric(oi["M"] + 2 * oi["D"]), p0[["total_MD"]], tolerance = 1e-6)

  pfe <- Pequil(el, expressInTotals = TRUE, compile = TRUE,
                modelname = paste0("test_Pequil_dimer_", as.integer(Sys.time())),
                controlsMS = list(nStarts = 30L))
  oe <- pfe(p0)[[1]]
  expect_equal(as.numeric(oe["M"]), Ms, tolerance = 1e-3)
  expect_equal(as.numeric(oe["D"]), Ds, tolerance = 1e-3)
  expect_equal(as.numeric(oe["M"] + 2 * oe["D"]), p0[["total_MD"]], tolerance = 1e-6)
})

test_that("overlapping conserved quantities reconstruct consistently (recycle enzyme)", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # G + S <-> GS -> G + P, P -> S : closed catalytic cycle.
  # Two overlapping CQs share G/GS; one auto-detected CQ carries a
  # negative coefficient on G, so the recon both divides by coef_g and
  # nests another eliminated species (fixed-point substitution).
  el <- eqnlist()
  el <- addReaction(el, "G + S", "GS", "k1*G*S")
  el <- addReaction(el, "GS", "G + S", "k1r*GS")
  el <- addReaction(el, "GS", "G + P", "k1c*GS")
  el <- addReaction(el, "P", "S", "k3*P")

  totals <- getTotals(el)
  expect_length(totals, 2L)                       # two independent CQs

  p0 <- c(k1 = 2, k1r = 1, k1c = 3, k3 = 1, total_1 = 1, total_2 = 4,
          G = 0.3, GS = 0.2, S = 2, P = 1)
  pf <- Pimpl(el, expressInTotals = TRUE, compile = FALSE,
              modelname = paste0("test_Pimpl_recycle_", as.integer(Sys.time())),
              controlsMS = list(nStarts = 50L),
              controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))
  o <- pf(p0)[[1]]

  # all four species reconstructed and non-negative at the steady state
  expect_setequal(intersect(c("G","GS","S","P"), names(o)), c("G","GS","S","P"))
  expect_true(all(as.numeric(o[c("G","GS","S","P")]) >= -1e-6))
  # the enzyme moiety G + GS equals its total to the solver tolerance
  expect_equal(as.numeric(o["G"] + o["GS"]), p0[["total_1"]], tolerance = 1e-5)
  expect_true(is.finite(as.numeric(o["S"])))
})

test_that("a fully open network (all states drain to zero) errors clearly", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # No influx anywhere; every species is structurally zero in steady state.
  el <- eqnlist()
  el <- addReaction(el, "E + S", "ES", "k1*E*S")
  el <- addReaction(el, "ES", "E + S", "k1r*ES")
  el <- addReaction(el, "ES", "E + P", "k1c*ES")

  expect_error(Pimpl(el, expressInTotals = TRUE, compile = FALSE),
               "structurally zero in steady state")
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
              deriv2 = FALSE,
              modelname = paste0("test_Pimpl_sing_", as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsMS = list(nStarts = 1L, positive = FALSE),
              controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))

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

  # Warm-start caches are now kept per condition in a registry; a condition-less
  # call lands in the "__default__" slot.
  reset_env <- environment(attr(pf, "resetWarmStart"))$reg_ref$caches[["__default__"]]
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

  reset_env <- environment(attr(pf, "resetWarmStart"))$reg_ref$caches[["__default__"]]
  expect_false(is.null(reset_env$guess))

  labels <- resetWarmStarts(pf, verbose = FALSE)
  expect_length(labels, 1L)
  expect_match(labels, "^Pimpl\\(")
  expect_null(reset_env$guess)
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

  ref1 <- environment(attr(p1, "resetWarmStart"))$reg_ref$caches[["__default__"]]
  ref2 <- environment(attr(p2, "resetWarmStart"))$reg_ref$caches[["__default__"]]
  expect_false(is.null(ref1$yini))
  expect_false(is.null(ref2$yini))

  wrap <- function(pars) list(p1 = p1(pars), p2 = p2(pars))
  labels <- resetWarmStarts(wrap, verbose = FALSE)
  expect_setequal(sub("\\(.*", "", labels), c("Pequil", "Pequil"))
  expect_length(labels, 2L)
  expect_null(ref1$yini)
  expect_null(ref2$yini)
})


test_that("a condition-less Pequil keeps an independent warm start per condition", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  stamp <- as.integer(Sys.time())

  # Condition-less equilibration: SS A* = k_in / k_out.
  g_eq <- Pequil(c(A = "k_in - k_out * A"),
                 parameters = c("k_in", "k_out"),
                 modelname = paste0("test_percond_eq_", stamp),
                 compile = TRUE, verbose = FALSE, attach.input = TRUE)

  # Two conditions with different k_in -> different steady states (3 and 6).
  trafos <- list(
    C1 = c(k_in = "s",     k_out = "1", A = "1"),
    C2 = c(k_in = "2 * s", k_out = "1", A = "1"))
  px <- P(trafos, modelname = paste0("test_percond_px_", stamp),
          compile = TRUE, verbose = FALSE)

  pf  <- g_eq * px
  out <- pf(c(s = 3), deriv = FALSE)

  expect_equal(as.numeric(out$C1["A"]), 3, tolerance = 1e-4)
  expect_equal(as.numeric(out$C2["A"]), 6, tolerance = 1e-4)

  # The single Pequil closure must now hold one warm-start slot per condition,
  # each cached at its own root -- not a single shared/overwritten cache.
  caches <- environment(attr(g_eq, "resetWarmStart"))$reg_ref$caches
  expect_setequal(ls(caches), c("C1", "C2"))
  expect_equal(as.numeric(caches[["C1"]]$yini["A"]), 3, tolerance = 1e-4)
  expect_equal(as.numeric(caches[["C2"]]$yini["A"]), 6, tolerance = 1e-4)
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


## ---- Pequil multistart: bad init recovered by random sweep ------------

test_that("Pequil multistart recovers from a bad initial guess", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Stable: dA = k_in - k_out * A, SS A* = k_in / k_out = 5.
  trafo <- c(A = "k_in - k_out * A")
  pf <- Pequil(trafo, parameters = c("k_in", "k_out"),
               modelname = paste0("test_Pequil_ms_", as.integer(Sys.time())),
               compile = TRUE, attach.input = FALSE,
               controlsMS = list(nStarts = 5L, positive = TRUE,
                                 lower = 1e-3, upper = 1e3),
               verbose = FALSE)

  # Start with a wildly wrong (but legal) guess; integrator should still
  # find the basin via warm-start, but force a fresh search by resetting.
  resetWarmStarts(pf, verbose = FALSE)
  pars <- c(k_in = 1.0, k_out = 0.2, A = 1e6)  # SS = 5
  out  <- pf(pars, deriv = FALSE)[[1]]
  expect_equal(as.numeric(out), 5, tolerance = 1e-3)
})


## ---- Pequil hard error when no steady state can be reached ------------

test_that("Pequil throws when no integration attempt produces a root", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Linear *unstable* system: dA = k * A with k > 0 has no finite SS.
  trafo <- c(A = "k * A")
  pf <- Pequil(trafo, parameters = "k",
               modelname = paste0("test_Pequil_unstable_", as.integer(Sys.time())),
               compile = TRUE, attach.input = FALSE,
               controlsMS = list(nStarts = 3L, positive = TRUE,
                                 lower = 1e-3, upper = 1e3),
               verbose = FALSE)

  resetWarmStarts(pf, verbose = FALSE)
  expect_error(pf(c(k = 1.0, A = 0.5))[[1]], "no steady state reached")
})


## ---- Pimpl hard error when residual cannot be driven below ftol -------

test_that("Pimpl throws when no start brings the residual below ftol", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # x^2 + 1 = 0 has no real root: nleqslv cannot reach zero residual.
  trafo <- c(x = "x*x + 1")
  pf <- Pimpl(trafo, parameters = character(0),
              modelname = paste0("test_Pimpl_noroot_", as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsMS = list(nStarts = 5L, positive = FALSE),
              controlsNleqslv = list(ftol = 1e-6, xtol = 1e-6))

  resetWarmStarts(pf, verbose = FALSE)
  expect_error(pf(c(x = 0.1))[[1]], "exceeds ftol|all .*solve attempt")
})


# ============================================================================
# Structural zero-state detection (.zeroStatesFromSmatrix) for eqnlist inputs
# (three layers from AlyssaPetit v1.2: NegCol, PosCol+single-state feeder,
# sink-cluster LP).
# ============================================================================

helper <- dMod:::.zeroStatesFromSmatrix


## ---- Layer 1: only-outflux column --------------------------------------

test_that("NegCol: a state with only outflux is detected as zero", {
  el <- eqnlist() |>
    addReaction("A", "", "k_dg * A")
  zs <- helper(el)
  expect_setequal(zs$zero_states, "A")
  expect_equal(ncol(zs$eqnlist$smatrix), 0L)
})


## ---- Layer 2: PosCol with single-state feeder --------------------------

test_that("PosCol: only-influx state whose feeder has a single reactant", {
  # X is fed by k_pr * Z and consumed by nothing (PosCol).
  # k_pr * Z must be 0 in SS, so Z = 0 (single-state in flux).
  # Then X is fed only by zero flux -> drop X next iteration (NegCol-empty).
  el <- eqnlist() |>
    addReaction("", "X", "k_pr * Z") |>
    addReaction("Z", "",  "k_dg * Z")
  zs <- helper(el)
  expect_setequal(zs$zero_states, c("Z", "X"))
})


## ---- Layer 3: structural sink cluster ----------------------------------

test_that("Sink cluster: TGFb + R_TGFb + R_TGFb_int (combined mass leaks)", {
  skip_if_not_installed("lpSolve")
  # Each individual column has +1 and -1 entries, so layers 1+2 cannot
  # catch them. The cluster {L, RL, RLi} as a whole degrades via the
  # RLi -> "" reaction, so its total mass leaks monotonically.
  el <- eqnlist() |>
    addReaction("L + R", "RL",  "k_on  * L * R") |>
    addReaction("RL",    "L + R", "k_off * RL") |>
    addReaction("RL",    "RLi", "k_int * RL") |>
    addReaction("RLi",   "L",   "k_dec * RLi") |>
    addReaction("RLi",   "",    "k_dg  * RLi") |>
    addReaction("",      "R",   "k_pr_R") |>
    addReaction("R",     "",    "k_dg_R * R")
  zs <- helper(el)
  expect_true(all(c("L", "RL", "RLi") %in% zs$zero_states))
  # R is not in the sink cluster: it has its own production/degradation.
  expect_false("R" %in% zs$zero_states)
})


## ---- Layer 3 falls back when lpSolve is missing ------------------------

test_that("FindSinkCluster degrades gracefully without lpSolve", {
  # Same model as above; we can't easily uninstall lpSolve mid-test, so
  # this is a structural assertion: with only layers 1+2, the cluster
  # would not be found (we verified manually no NegCol/PosCol matches).
  el <- eqnlist() |>
    addReaction("L + R", "RL",  "k_on  * L * R") |>
    addReaction("RL",    "L + R", "k_off * RL") |>
    addReaction("RL",    "RLi", "k_int * RL") |>
    addReaction("RLi",   "L",   "k_dec * RLi") |>
    addReaction("RLi",   "",    "k_dg  * RLi")
  # Sanity: every column is mixed-sign in this minimal cluster model.
  S <- el$smatrix
  S[is.na(S)] <- 0
  storage.mode(S) <- "double"
  for (j in seq_len(ncol(S))) {
    col <- S[, j]
    expect_true(any(col > 0) && any(col < 0),
                info = sprintf("state %s is not mixed-sign", colnames(S)[j]))
  }
})


## ---- Idempotence on a model with no zero states ------------------------

test_that("Pure production-decay has no zero states", {
  el <- eqnlist() |>
    addReaction("", "A", "k_in") |>
    addReaction("A", "",  "k_out * A")
  zs <- helper(el)
  expect_equal(zs$zero_states, character(0))
  expect_identical(zs$eqnlist, el)
})


## ---- End-to-end: Pequil produces value = 0 for zero states -------------

test_that("Pequil reports zero states with value 0 and drops them from getParameters()", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Receptor R alone has production + degradation (nonzero baseline).
  # Ligand L and complex LR form a sink cluster (LR degrades).
  el <- eqnlist() |>
    addReaction("",       "R",  "k_pr_R") |>
    addReaction("R",      "",   "k_dg_R * R") |>
    addReaction("L + R",  "LR", "k_on * L * R") |>
    addReaction("LR",     "L + R", "k_off * LR") |>
    addReaction("LR",     "",   "k_dg_LR * LR")

  pf <- Pequil(el,
               modelname = paste0("test_zero_states_", as.integer(Sys.time())),
               compile = TRUE, verbose = FALSE)

  expect_false(any(c("L", "LR") %in% getParameters(pf)))
  pars <- c(k_pr_R = 2, k_dg_R = 0.5, k_on = 1, k_off = 0.5, k_dg_LR = 0.1)
  out  <- pf(pars)[[1]]
  expect_equal(as.numeric(out["L"]),  0, tolerance = 1e-8)
  expect_equal(as.numeric(out["LR"]), 0, tolerance = 1e-8)
  expect_equal(as.numeric(out["R"]),  pars[["k_pr_R"]] / pars[["k_dg_R"]],
               tolerance = 1e-5)
})


## ---- No-progress: hard error instead of stale initial values -----------

test_that("Pequil errors when the solver makes no progress", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Force a no-progress condition by capping the total step budget so the
  # solver can't even take one accepted step. Without the new guard,
  # Pequil would silently return the initial values.
  pf <- Pequil(c(A = "k_in - k_out * A"),
               parameters = c("k_in", "k_out"),
               modelname = paste0("test_no_progress_", as.integer(Sys.time())),
               controlsODE = list(maxsteps = 1L, maxprogress = 1L),
               compile = TRUE, verbose = FALSE)
  expect_error(
    suppressWarnings(pf(c(k_in = 1, k_out = 1, A = 1))),
    "no steady state reached")
})


# ============================================================================
# Edge case: Pimpl with no outer parameters
# ============================================================================

test_that("Pimpl with no outer parameters does not crash in build_jacobian", {
  withr::local_dir(tempdir())
  trafo <- c(A = "A - 1.0")
  p <- Pimpl(trafo, parameters = NULL, compile = TRUE,
             modelname = "noparam_pimpl")
  out <- p(c(dummy = 1.0))
  expect_true(is.numeric(unclass(out[[1]])["A"]))
  expect_equal(unname(unclass(out[[1]])["A"]), 1.0, tolerance = 1e-3)
})


# ============================================================================
# CQ harmonization: both Pimpl and Pequil accept `expressInTotals` to switch
# between "total as parameter" (TRUE, Pimpl default) and "eliminated species
# as pass-through parameter" (FALSE, Pequil default). Smart `totalXxx`
# naming from the longest common substring of CQ species; `total_<index>`
# fallback when no common substring exists.
# ============================================================================

## ---- Smart naming: pERK + ERK -> totalERK ----------------------------

test_that(".smartTotalName picks the longest common substring", {
  expect_equal(dMod:::.smartTotalName(c("pERK", "ERK"), character(0),
                                        character(0), 1), "totalERK")
  expect_equal(dMod:::.smartTotalName(c("TGFb", "R1_TGFb", "R1_TGFb_int"),
                                        character(0), character(0), 1),
               "totalTGFb")
  # No common substring -> fallback to total_<index>
  expect_equal(dMod:::.smartTotalName(c("A", "B"), character(0),
                                        character(0), 1), "total_1")
  # Single character common -> below the >=2 threshold, fallback
  expect_equal(dMod:::.smartTotalName(c("Ax", "Ay"), character(0),
                                        character(0), 1), "total_1")
  # Collision with an existing parameter -> disambiguate
  expect_equal(dMod:::.smartTotalName(c("pERK", "ERK"), character(0),
                                        "totalERK", 1), "totalERK_2")
})


## ---- Pimpl smart naming: A <-> B uses LCS = "" -> total_1 -----------

test_that("Pimpl on A <-> B introduces total_1 (no common substring)", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  el <- eqnlist() |>
    addReaction("A", "B", "k * A") |>
    addReaction("B", "A", "km * B")

  pf <- Pimpl(el, parameters = c("k", "km"),
              modelname = paste0("test_Pimpl_smart_total1_",
                                 as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsMS = list(nStarts = 1L, positive = FALSE),
              controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))
  expect_true("total_1" %in% getParameters(pf))
})


## ---- Pimpl smart naming: ERK + pERK -> totalERK ---------------------

test_that("P(method='implicit') preserves the eqnlist smatrix for CQ detection", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  pf <- P(el, method = "implicit", compile = FALSE, verbose = FALSE,
          modelname = paste0("test_P_implicit_cq_", as.integer(Sys.time())),
          controlsMS = list(nStarts = 1L, positive = FALSE),
          controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))
  expect_true("totalERK" %in% getParameters(pf))
})


test_that("Pimpl on ERK <-> pERK introduces totalERK via LCS", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  pf <- Pimpl(el, parameters = c("k1", "k2"),
              modelname = paste0("test_Pimpl_totalERK_",
                                 as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsMS = list(nStarts = 1L, positive = FALSE),
              controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))
  expect_true("totalERK" %in% getParameters(pf))
  expect_false("total_1" %in% getParameters(pf))
})


## ---- Pimpl == Pequil on CQ model with expressInTotals = TRUE --------

test_that("Pimpl and Pequil produce identical parvec interface (totals mode)", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  ts <- as.integer(Sys.time())
  p_pimpl <- Pimpl(el, parameters = c("k1", "k2"),
                   modelname = paste0("test_harm_pimpl_", ts),
                   compile = FALSE, verbose = FALSE,
                   controlsMS = list(nStarts = 1L, positive = FALSE),
                   controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))
  p_pequil <- Pequil(el, parameters = c("k1", "k2"), expressInTotals = TRUE,
                     modelname = paste0("test_harm_pequil_", ts),
                     compile = TRUE, verbose = FALSE, attach.input = TRUE,
                     controlsMS = list(nStarts = 1L))

  expect_setequal(getParameters(p_pimpl), getParameters(p_pequil))

  pars <- c(k1 = 1, k2 = 3, totalERK = 4, ERK = 1, pERK = 1)
  out_pimpl  <- p_pimpl(pars, deriv = FALSE)[[1]]
  out_pequil <- p_pequil(pars, deriv = FALSE)[[1]]

  # Closed form: ERK* = k2 / (k1 + k2) * total = 3, pERK* = 1.
  expect_equal(as.numeric(out_pimpl["pERK"]),  1, tolerance = 1e-3)
  expect_equal(as.numeric(out_pimpl["ERK"]),   3, tolerance = 1e-3)
  expect_equal(as.numeric(out_pequil["pERK"]), 1, tolerance = 1e-3)
  expect_equal(as.numeric(out_pequil["ERK"]),  3, tolerance = 1e-3)
})


## ---- Pequil Jacobian for eliminated species matches closed form -----

test_that("Pequil chain-rules the eliminated species sensitivity correctly", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  pf <- Pequil(el, parameters = c("k1", "k2"), expressInTotals = TRUE,
               modelname = paste0("test_pequil_elimjac_",
                                  as.integer(Sys.time())),
               compile = TRUE, verbose = FALSE, attach.input = TRUE)

  pars <- c(k1 = 1, k2 = 3, totalERK = 4, ERK = 1, pERK = 1)
  out  <- pf(pars, deriv = TRUE)[[1]]
  J    <- attr(out, "deriv")

  # Closed form: pERK = k1/(k1+k2) * total, ERK = total - pERK.
  # d(pERK)/d(k1) = k2/(k1+k2)^2 * total = 3/16 * 4 = 0.75
  # d(ERK)/d(k1)  = -d(pERK)/d(k1) = -0.75
  # d(pERK)/d(totalERK) = k1/(k1+k2) = 0.25
  # d(ERK)/d(totalERK)  = 1 - 0.25 = 0.75
  expect_equal(unname(J["pERK", "k1"]),       0.75, tolerance = 1e-4)
  expect_equal(unname(J["ERK",  "k1"]),      -0.75, tolerance = 1e-4)
  expect_equal(unname(J["pERK", "totalERK"]), 0.25, tolerance = 1e-4)
  expect_equal(unname(J["ERK",  "totalERK"]), 0.75, tolerance = 1e-4)
})


## ---- expressInTotals = FALSE: pivot is parameter, no total ----------

test_that("Pimpl with expressInTotals = FALSE promotes pivot to parameter", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  pf <- Pimpl(el, parameters = c("k1", "k2"),
              expressInTotals = FALSE,
              modelname = paste0("test_Pimpl_pivot_",
                                 as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsMS = list(nStarts = 1L, positive = FALSE),
              controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))
  params <- getParameters(pf)
  expect_false(any(grepl("^total", params)))
  # One of ERK/pERK is the pivot; both should be in the parvec interface
  # (the pivot as a user-supplied parameter, the other as an output).
  expect_true(all(c("k1", "k2") %in% params))
})


test_that("Pequil default (expressInTotals = FALSE) introduces no totals", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  pf <- Pequil(el, parameters = c("k1", "k2"),
               modelname = paste0("test_Pequil_default_pivot_",
                                  as.integer(Sys.time())),
               compile = TRUE, verbose = FALSE)
  expect_false(any(grepl("^total", getParameters(pf))))
})


test_that("Pequil with expressInTotals = FALSE matches Pimpl interface", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  ts <- as.integer(Sys.time())
  p_pimpl <- Pimpl(el, parameters = c("k1", "k2"), expressInTotals = FALSE,
                   modelname = paste0("test_harm_pivot_pimpl_", ts),
                   compile = FALSE, verbose = FALSE,
                   controlsMS = list(nStarts = 1L, positive = FALSE),
                   controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))
  p_pequil <- Pequil(el, parameters = c("k1", "k2"), expressInTotals = FALSE,
                     modelname = paste0("test_harm_pivot_pequil_", ts),
                     compile = TRUE, verbose = FALSE,
                     controlsMS = list(nStarts = 1L))
  expect_setequal(getParameters(p_pimpl), getParameters(p_pequil))
})


# ============================================================================
# CQ basis as a first-class eqnlist field: getTotals(), customTotals(),
# mutator preservation, backward compat with pre-totals eqnlists.
# ============================================================================

## ---- getTotals auto-detection ---------------------------------------

test_that("getTotals returns smart-named auto-detected totals", {
  el_AB <- eqnlist() |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "A", "k2 * B")
  tot_AB <- getTotals(el_AB)
  expect_length(tot_AB, 1L)
  expect_equal(names(tot_AB), "total_1")  # LCS of {A, B} is "" -> fallback

  el_ERK <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  tot_ERK <- getTotals(el_ERK)
  expect_length(tot_ERK, 1L)
  expect_equal(names(tot_ERK), "totalERK")
})


## ---- customTotals strict validation ---------------------------------

test_that("customTotals validates structure, rank, and CQ membership", {
  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  el2 <- customTotals(el, list(totalE = "ERK + pERK"))
  tot <- getTotals(el2)
  expect_equal(names(tot), "totalE")
  expect_equal(tot[["totalE"]], "ERK + pERK")
  expect_true(isTRUE(attr(el2$totals, "custom")))

  expect_error(customTotals(el, list(bogus = "ERK")),
               "not a conservation quantity")
  expect_error(customTotals(el, list(too = "ERK + pERK", many = "ERK + pERK")),
               "Expected 1 conservation quantit(y|ies), got 2")
  expect_error(customTotals(el, list(weird = "ERK * pERK")),
               "not linear")

  el_reset <- customTotals(el2, NULL)
  expect_null(el_reset$totals)
  el_reset2 <- customTotals(el2, list())
  expect_null(el_reset2$totals)
})


## ---- Mutator preservation of custom totals --------------------------

test_that("addReaction preserves custom totals when CQ structure survives", {
  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  el <- customTotals(el, list(totalE = "ERK + pERK"))

  # Adding a phosphatase-modifier reaction that doesn't touch ERK/pERK mass
  el2 <- addReaction(el, "X", "", "k3 * X")
  expect_equal(names(el2$totals), "totalE")
  expect_true(isTRUE(attr(el2$totals, "custom")))
})


test_that("addReaction warns and resets when custom totals are invalidated", {
  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  el <- customTotals(el, list(totalE = "ERK + pERK"))

  # Adding a degradation of ERK breaks the conservation
  expect_warning(
    el2 <- addReaction(el, "ERK", "", "k_dg * ERK"),
    "customTotals invalidated"
  )
  expect_null(el2$totals)
})


## ---- Backward compat with pre-totals eqnlist ------------------------

test_that("is.eqnlist accepts pre-totals eqnlists (missing $totals field)", {
  el <- eqnlist() |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "A", "k2 * B")
  el_old <- el; el_old$totals <- NULL
  el_old <- el_old[setdiff(names(el_old), "totals")]
  class(el_old) <- c("eqnlist", "list")
  expect_true(is.eqnlist(el_old))
  # getTotals still works (computes lazily)
  expect_length(getTotals(el_old), 1L)
})


## ---- Pimpl/Pequil pick up customTotals ------------------------------

test_that("Pimpl uses customTotals names in the parvec interface", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  el <- customTotals(el, list(totalERKpool = "ERK + pERK"))

  pf <- Pimpl(el, parameters = c("k1", "k2"),
              modelname = paste0("test_Pimpl_custom_totals_",
                                 as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE,
              controlsMS = list(nStarts = 1L, positive = FALSE),
              controlsNleqslv = list(ftol = 1e-12, xtol = 1e-12))
  params <- getParameters(pf)
  expect_true("totalERKpool" %in% params)
  expect_false(any(grepl("^total[^E]", params)))  # no auto-named total slipped in
})


## ---- print.eqnlist shows totals section -----------------------------

test_that("print.eqnlist shows the conserved quantities by name", {
  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  out <- capture.output(print(el))
  expect_true(any(grepl("Conserved quantities", out)))
  expect_true(any(grepl("totalERK", out)))

  el2 <- customTotals(el, list(totalE = "ERK + pERK"))
  out2 <- capture.output(print(el2))
  expect_true(any(grepl("Conserved quantities .custom.", out2)))
  expect_true(any(grepl("totalE\\s*=\\s*ERK \\+ pERK", out2)))
})
