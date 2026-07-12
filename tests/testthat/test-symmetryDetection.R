# Behavioral tests for symmetryDetection() (structural non-identifiabilities).
#
# symmetryDetection drives the self-contained Python module
# symmetryDetectionVersion2 via reticulate. Engines:
#   method = "observability"  -> observability-identifiability matrix (rank)
#   method = "polynomial"      -> polynomial Lie-symmetry ansatz (generators)
#   method = "scaling"        -> integer-kernel scaling symmetries
# The API is purely symbolic: f / g / trafo are coercible with as.eqnvec.
# Reference: Merkt et al. 2015 (PRE 92, 012920).


.sympy_works <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  isTRUE(tryCatch({
    sympy <- reticulate::import("sympy", convert = TRUE)
    nzchar(sympy[["__version__"]])
  }, error = function(e) FALSE))
}


.sd_module <- function() {
  code_dir <- system.file("code", package = "dMod")
  sysmod <- reticulate::import("sys", convert = TRUE)
  if (!(code_dir %in% sysmod$path)) sysmod$path <- c(code_dir, sysmod$path)
  reticulate::import("symmetryDetectionVersion2", convert = TRUE)
}


.canonical <- function() {
  eqnlist() |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "A", "k2 * B")
}


# exact symbolic equality of two expression strings via sympy
.sym_expr_equal <- function(a, b) {
  spy <- reticulate::import("sympy", convert = TRUE)
  as.character(spy$simplify(spy$sympify(paste0("(", a, ") - (", b, ")")))) == "0"
}


# numeric tangent of a reported direction at point `pt` over `coords`. $generator
# always holds the tangent components xi_i directly (a scaling's integer weight
# w_i is expanded to xi_i = w_i * z_i at the finalisation boundary), so evaluating
# each component at the point gives the tangent.
.sym_tangent <- function(d, pt, coords) {
  v <- setNames(numeric(length(coords)), coords)
  for (nm in names(d$generator))
    v[nm] <- eval(parse(text = d$generator[[nm]]), pt)
  v
}


test_that("liesym finds and verifies the canonical scaling symmetry", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  res <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "polynomial", reduceCQ = FALSE,
                           polynomial = polynomialControl(ansatz = "uni", pMax = 1L))

  gens <- res$symmetries
  expect_gte(length(gens), 1L)
  scaling <- gens[grepl("scaling", vapply(gens, function(r) r$type, character(1)))]
  expect_gte(length(scaling), 1L)
  expect_true(all(c("A", "B", "alpha") %in% names(scaling[[1]]$generator)))
  expect_true(isTRUE(scaling[[1]]$verified))
})


test_that("liesym exact and legacy float solvers agree on the canonical case", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  ex  <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "polynomial", reduceCQ = FALSE,
                           polynomial = polynomialControl(ansatz = "uni", pMax = 1L, exact = TRUE))
  leg <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "polynomial", reduceCQ = FALSE,
                           polynomial = polynomialControl(ansatz = "uni", pMax = 1L, exact = FALSE))
  expect_gte(length(ex$symmetries), 1L)          # non-vacuous: both actually find a symmetry
  expect_equal(length(ex$symmetries), length(leg$symmetries))
})


test_that("liesym handles a log10 observable with an offset parameter", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  res <- symmetryDetection(.canonical(), eqnvec(Aobs = "log10(A) + c"),
                           method = "polynomial", reduceCQ = FALSE,
                           polynomial = polynomialControl(ansatz = "uni", pMax = 1L))

  expect_gte(length(res$symmetries), 1L)
  gen <- res$symmetries[[1]]
  expect_true(all(c("A", "B", "c") %in% names(gen$generator)))
  expect_true(isTRUE(gen$verified))
})


test_that("a trafo that fixes a parameter removes its symmetry", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # substituting alpha = 1 through the trafo removes the calibration scaling
  res <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           trafo = eqnvec(alpha = "1"), method = "polynomial",
                           reduceCQ = FALSE, polynomial = polynomialControl(ansatz = "uni", pMax = 1L))
  expect_equal(length(res$symmetries), 0L)
})


test_that("a trafo renames a rate and is substituted into the model", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # k1 -> kf renames the forward rate; the canonical scaling is unaffected
  res <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           trafo = eqnvec(k1 = "kf"), method = "observability",
                           reduceCQ = FALSE)
  expect_false(res$identifiable)
  supp <- unlist(lapply(res$symmetries, function(d) d$support))
  expect_true(all(c("A", "B", "alpha") %in% supp))
  expect_false("k1" %in% supp)
})


test_that("observability flags the canonical non-identifiability", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  res <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "observability", reduceCQ = FALSE)

  expect_false(res$identifiable)
  expect_lt(res$rank, res$dim)
  supports <- lapply(res$symmetries, function(d) d$support)
  expect_true(any(vapply(supports,
                         function(s) all(c("A", "B", "alpha") %in% s), logical(1))))
})


test_that("observability reports a fully observed model as identifiable", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  res <- symmetryDetection(.canonical(), eqnvec(o1 = "A", o2 = "B"),
                           method = "observability", reduceCQ = FALSE)
  expect_true(res$identifiable)
  expect_equal(res$rank, res$dim)
})


test_that("scaling finds the canonical scaling exactly via the integer kernel", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  res <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "scaling", reduceCQ = FALSE)
  expect_equal(length(res$symmetries), 1L)
  v <- res$symmetries[[1]]$weights
  expect_true(all(c("A", "B", "alpha") %in% names(v)))
  expect_equal(as.integer(v[["A"]]), as.integer(v[["B"]]))
  expect_equal(as.integer(v[["alpha"]]), -as.integer(v[["A"]]))
})


test_that("observability scales to a deep chain via the exact modular engine", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  eq <- eqnlist() |>
    addReaction("",  "A", "kin")   |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "C", "k2 * B") |>
    addReaction("C", "D", "k3 * C") |>
    addReaction("D", "",  "k4 * D")
  res <- symmetryDetection(eq, eqnvec(yD = "scale * D"), reduceCQ = FALSE)

  expect_false(res$identifiable)
  expect_lt(res$rank, res$dim)
  expect_gte(length(res$symmetries), 1L)
})


test_that("observability rejects a non-rational observable", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  expect_error(
    symmetryDetection(.canonical(), eqnvec(Aobs = "sqrt(A)"),
                      method = "observability", reduceCQ = FALSE),
    "rational")
})


test_that("reconstruct observability returns exact rational directions", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  res <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "observability", reconstruct = TRUE,
                           reduceCQ = FALSE)

  expect_false(res$identifiable)
  expect_equal(res$info$engine, "modular")
  expect_true(all(vapply(res$symmetries,
                         function(d) is.logical(d$explicit) && length(d$support) > 0,
                         logical(1))))
  # the calibration scaling, peeled exactly as an integer toric direction
  scal <- Filter(function(d) all(c("A", "B", "alpha") %in% names(d$generator)),
                 res$symmetries)
  expect_gte(length(scal), 1L)
  d <- scal[[1]]
  expect_identical(d$type, "scaling")
  # A and B carry equal weight, opposite to the readout gain alpha
  w <- vapply(d$weights, as.integer, integer(1))
  expect_equal(unname(w[["A"]]), unname(w[["B"]]))
  expect_equal(unname(w[["alpha"]]), -unname(w[["A"]]))
})


test_that("reconstruct observability reconstructs a non-monomial direction", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  res <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "observability", reconstruct = TRUE,
                           reduceCQ = FALSE)
  # the conserved-quantity direction, reported as its canonical polynomial-
  # primitive (affine) generator: dB = A + B, dk1 = k2, dk2 = -k2 -- the same
  # direction as the raw dk2 = 1, dk1 = -1, dB = -(A + B)/k2, cleared of its 1/k2 gauge
  d <- Filter(function(d) all(c("B", "k1", "k2") %in% names(d$generator)) &&
                isTRUE(d$explicit), res$symmetries)
  expect_gte(length(d), 1L)
  d <- d[[1]]
  expect_equal(d$type, "affine")
  expect_true(.sym_expr_equal(d$generator[["k2"]], "-k2"))
  expect_true(.sym_expr_equal(d$generator[["k1"]], "k2"))
  expect_true(.sym_expr_equal(d$generator[["B"]], "A + B"))

  # the scaling engine alone finds only the scaling symmetry, not this direction
  sc <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                          method = "scaling", reduceCQ = FALSE)
  expect_equal(length(sc$symmetries), 1L)
})


test_that("the closed-form and support-only verdicts agree", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  ana <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "observability", reconstruct = TRUE,
                           reduceCQ = FALSE)
  sup <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "observability", reconstruct = FALSE,
                           reduceCQ = FALSE)
  expect_equal(ana$rank, sup$rank)
  expect_equal(ana$dim, sup$dim)
  # scalings are peeled exactly and reported in closed form in both modes;
  # only the residual (non-scaling) directions stay support-only without reconstruct
  general <- Filter(function(d) !isTRUE(d$type == "scaling"), sup$symmetries)
  expect_false(any(vapply(general, function(d) isTRUE(d$explicit), logical(1))))
})


test_that("reconstruct observability handles a deep chain", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  eq <- eqnlist() |>
    addReaction("",  "A", "kin")   |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "C", "k2 * B") |>
    addReaction("C", "D", "k3 * C") |>
    addReaction("D", "",  "k4 * D")
  res <- symmetryDetection(eq, eqnvec(yD = "scale * D"),
                           method = "observability", reconstruct = TRUE,
                           reduceCQ = FALSE)

  expect_false(res$identifiable)
  expect_gte(length(res$symmetries), 1L)
  # the scaling directions of the chain are peeled and returned in closed form
  expect_true(any(vapply(res$symmetries,
                         function(d) isTRUE(d$explicit), logical(1))))
})


test_that("a fixed parameter is excluded from z and removes its symmetry", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  free  <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                             method = "observability", reconstruct = TRUE,
                             reduceCQ = FALSE)
  fixed <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                             method = "observability", reconstruct = TRUE,
                             fixed = "alpha", reduceCQ = FALSE)

  expect_equal(fixed$dim, free$dim - 1L)
  expect_lt(length(fixed$symmetries), length(free$symmetries))
  expect_false("alpha" %in% unlist(lapply(fixed$symmetries,
                                          function(d) d$support)))
})


test_that("observability, liesym and scaling agree on the scaling symmetry", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  abAlpha <- function(s) all(c("A", "B", "alpha") %in% s)

  obs <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "observability", reduceCQ = FALSE)
  lie <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "polynomial", reduceCQ = FALSE,
                           polynomial = polynomialControl(ansatz = "uni", pMax = 1L))
  scl <- symmetryDetection(.canonical(), eqnvec(Aobs = "alpha * A"),
                           method = "scaling", reduceCQ = FALSE)

  expect_true(any(vapply(obs$symmetries,
                         function(d) abAlpha(d$support), logical(1))))
  expect_true(any(vapply(lie$symmetries,
                         function(g) abAlpha(names(g$generator)), logical(1))))
  expect_true(any(vapply(scl$symmetries,
                         function(d) abAlpha(d$support), logical(1))))
})


test_that("a steady-state expression initial condition is seeded with its duals", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # x' = b - a x rests at x* = b/a; with that initial condition (and no dose) the
  # trajectory is constant, so only the product s*b/a is observed (rank 1 of 3)
  res <- symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
                           trafo = eqnvec(x = "b/a"), method = "observability",
                           reconstruct = TRUE, reduceCQ = FALSE)
  expect_false(res$identifiable)
  expect_equal(res$rank, 1L)
  expect_equal(res$dim, 3L)
})


test_that("a pre-equilibrated model with a known dose is identifiable", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  ident <- function(m) {
    ev <- addEvent(eventlist(), var = "x", time = 0, value = "dose", method = m)
    symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
                      trafo = eqnvec(x = "b/a"), events = ev,
                      conditions = data.frame(dose = 2, row.names = "stim"),
                      method = "observability", reduceCQ = FALSE)$identifiable
  }
  # an absolute dose (replace or add) fixes the scale; a relative one does not
  expect_true(ident("replace"))
  expect_true(ident("add"))
  expect_false(ident("multiply"))
})


test_that("equilibrate reproduces the explicit steady state", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # no symbolic initial: the modular solver must find the interior point
  # x* = b/a of f = 0 and seed it with its IFT parameter-duals, giving the same
  # rank 1/3 and non-identifiable directions as the explicit trafo x = b/a
  f <- eqnvec(x = "b - a*x")
  g <- eqnvec(y = "s*x")
  res <- symmetryDetection(f, g, method = "observability",
                           equilibrate = TRUE, reconstruct = TRUE,
                           reduceCQ = FALSE)
  expect_false(res$identifiable)
  expect_equal(res$rank, 1L)
  expect_equal(res$dim, 3L)
  # every reported direction is an exact null direction of the lone observed
  # gradient grad(s*b/a) = (-s*b/a^2, s/a, b/a) in coordinates (a, b, s)
  vals <- list(a = 2, b = 3, s = 5)
  grad <- c(a = -vals$s * vals$b / vals$a^2, b = vals$s / vals$a, s = vals$b / vals$a)
  for (d in res$symmetries) {
    v <- .sym_tangent(d, vals, c("a", "b", "s"))
    expect_equal(sum(grad * v), 0, tolerance = 1e-9)
  }
})


test_that("equilibrate handles a steady state with no rational form", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # x* = sqrt(b/a) is irrational, so steadyStates() cannot return it, but the
  # saturated variety has an interior GF(p) point whenever b/a is a quadratic
  # residue; the modular solver finds it and the verdict is still reached
  res <- symmetryDetection(eqnvec(x = "b - a*x^2"), eqnvec(y = "s*x"),
                           method = "observability", equilibrate = TRUE,
                           reconstruct = TRUE, reduceCQ = FALSE)
  expect_false(res$identifiable)
  expect_equal(res$rank, 1L)
  expect_equal(res$dim, 3L)
})


test_that("equilibrate ignores a state initial condition given in trafo", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # a state-named trafo entry is an initial condition; equilibrate solves the
  # steady state from f = 0 instead, so it is dropped with a warning
  expect_warning(
    res <- symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
                             method = "observability", equilibrate = TRUE,
                             trafo = eqnvec(x = "x0"), reduceCQ = FALSE),
    "ignored")
  expect_false(res$identifiable)
  expect_equal(res$rank, 1L)
  expect_equal(res$dim, 3L)
  # a params-only trafo is a substitution and the constraint runs through it
  res <- symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
                           equilibrate = TRUE, trafo = eqnvec(a = "k1*k2"),
                           method = "observability", reduceCQ = FALSE)
  expect_false(res$identifiable)
  expect_equal(res$rank, 1L)
  expect_equal(res$dim, 4L)
})


test_that("equilibrate applies a t0 dose on top of the resting state", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # the steady state x* = b/a is solved event-free; a replace dose at t0 then
  # perturbs the start point, so the relaxation is observed and (a, b, s) become
  # identifiable. Without the dose only the resting product s*b/a is seen.
  dose <- addEvent(eventlist(), var = "x", time = 0, value = "dose",
                   method = "replace")
  res <- symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
                           method = "observability", equilibrate = TRUE,
                           events = dose,
                           conditions = data.frame(dose = 2, row.names = "stim"))
  expect_true(res$identifiable)
  expect_equal(res$rank, 3L)

  res0 <- symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
                            method = "observability", equilibrate = TRUE)
  expect_false(res0$identifiable)
  expect_equal(res0$rank, 1L)
})


test_that("the compiled steady-state seed matches the symbolic solve", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # The compiled linear-elimination seed must reproduce the symbolic solve
  # bit-for-bit mod p, including the implicit-function parameter duals and t0
  # event composition. setSteadyStateForceSympy toggles the symbolic path.
  sd <- .sd_module()
  primes <- c(2147483647, 2147483629, 2147483587)
  cases <- list(
    list(model = "x = b - a*x", states = "x", params = c("a", "b"), ev = NULL),
    list(model = c("x = b - a*x", "y = c*x - d*y"), states = c("x", "y"),
         params = c("a", "b", "c", "d"), ev = NULL),
    list(model = "x = b - a*x", states = "x", params = c("a", "b", "dose"),
         ev = list(list(var = "x", value = "dose*2", method = "replace"))))

  for (cs in cases) for (p in primes) for (seed in 4:8) {
    pv <- as.list(setNames((seed * 7 + seq_along(cs$params) * 131 + 11) %% p,
                           cs$params))
    sd$setSteadyStateForceSympy(FALSE)
    fast <- sd$solveSteadyStateModular(cs$model, cs$states, cs$params, pv, p,
                                       t0events = cs$ev)
    sd$setSteadyStateForceSympy(TRUE)
    symb <- sd$solveSteadyStateModular(cs$model, cs$states, cs$params, pv, p,
                                       t0events = cs$ev)
    sd$setSteadyStateForceSympy(FALSE)
    expect_identical(fast$ok, symb$ok)
    expect_identical(fast$xstar, symb$xstar)
    expect_identical(fast$dx, symb$dx)
  }
})


test_that("the fast numeric solve matches the symbolic solve on a coupled steady state", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # A steady state that is NOT generically linear (a bilinear A*B core, so the
  # genericLinear/linTerms fast path does not fire) must still solve bit-for-bit mod p
  # via the fast dict-polynomial numeric elimination (_solve_states_fast): linear
  # states eliminated by modular arithmetic, the small coupled residual by Groebner.
  # setSteadyStateForceSympy toggles the reference symbolic path (_solve_states_modular).
  sd <- .sd_module()
  model <- c("A = kp  - kd*A  - kf*A*B + kr*C",
             "B = kpB - kdB*B - kf*A*B + kr*C",
             "C = kf*A*B - kr*C - kdC*C")
  states <- c("A", "B", "C")
  params <- c("kp", "kd", "kf", "kr", "kpB", "kdB", "kdC")
  on.exit(sd$setSteadyStateForceSympy(FALSE))

  primes <- c(2147483647, 2147483629, 2147483587)
  nfeasible <- 0L
  for (p in primes) for (seed in 3:12) {
    pv <- as.list(setNames((seed * (seq_along(params) + 3L) + 7L) %% p, params))
    sd$setSteadyStateForceSympy(FALSE)
    fast <- sd$solveSteadyStateModular(model, states, params, pv, p, jointMode = TRUE)
    sd$setSteadyStateForceSympy(TRUE)
    symb <- sd$solveSteadyStateModular(model, states, params, pv, p, jointMode = TRUE)
    expect_identical(fast$ok, symb$ok)
    if (isTRUE(fast$ok)) {
      expect_identical(fast$valBy, symb$valBy)
      expect_identical(fast$dfJx, symb$dfJx)
      expect_identical(fast$dfJt, symb$dfJt)
      nfeasible <- nfeasible + 1L
    }
  }
  expect_gt(nfeasible, 0L)   # at least some points exercised the coupled solve
})


test_that("the C++ steady-state seed reproduces the symbolic verdict", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # the compiled C++ seed (the steady-state IC tape for a plain linear resting
  # state, the symSteadyStateSeed kernel for a free-Hill-exponent recast one) must
  # reproduce the symbolic per-point solve.
  sd <- .sd_module()
  dirset <- function(r) sort(vapply(r$symmetries,
    function(d) paste(sort(d$support), collapse = "+"), character(1)))
  agree <- function(...) {
    cpp <- symmetryDetection(...)
    old <- options(dMod.symSeedCpp = FALSE)
    sd$setForceConstraintSeed(TRUE)
    on.exit({ options(old); sd$setForceConstraintSeed(FALSE) })
    sym <- symmetryDetection(...)
    expect_identical(cpp$rank, sym$rank)
    expect_identical(cpp$identifiable, sym$identifiable)
    expect_identical(dirset(cpp), dirset(sym))
  }

  # plain linear resting state with a t0 dose
  agree(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"), method = "observability",
        equilibrate = TRUE,
        events = addEvent(eventlist(), var = "x", time = 0, value = "dose",
                          method = "replace"),
        conditions = data.frame(dose = 2, row.names = "stim"))

  # free Hill exponent (power recast): the seed solves with E = FB^nh generic
  fb <- eqnvec(R = "k_pr/(1 + k_fb*FB^nh) - k_dg*R", FB = "k_pf*R - k_df*FB",
               u = "0")
  agree(fb, eqnvec(R_obs = "s1*R", FB_obs = "s2*FB"), method = "observability",
        equilibrate = TRUE, forcings = "u",
        events = addEvent(eventlist(), var = "u", time = -1, value = "var_u",
                          method = "replace"),
        conditions = data.frame(var_u = c(0, 1), row.names = c("c1", "c2")),
        reconstruct = TRUE, reduceCQ = FALSE)
})


test_that("the analytic-segment kernel agrees across thread counts", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # the per-condition gap-propagation build runs over `cores` OpenMP threads; the
  # stacked nullspace must not depend on the thread count
  f <- as.eqnvec(c(A = "-(k1 + u*k2)*A + b", u = "0"))
  g <- as.eqnvec(c(y = "s*A"))
  ev <- eventlist() |>
    addEvent(var = "A", time = 0, value = "dose", method = "add") |>
    addEvent(var = "u", time = -1, value = "var_u", method = "replace")
  grid <- data.frame(var_u = c(0, 1, 0, 1), dose = c(1, 1, 2, 2),
                     row.names = c("c1", "c2", "c3", "c4"))
  run <- function(cores)
    symmetryDetection(f, g, method = "observability", events = ev,
                      conditions = grid, forcings = "u",
                      reconstruct = TRUE, reduceCQ = FALSE, cores = cores)
  serial <- run(1L)
  parallel <- run(4L)
  expect_identical(serial$rank, parallel$rank)
  expect_identical(serial$identifiable, parallel$identifiable)
  expect_equal(length(serial$symmetries), length(parallel$symmetries))
})

test_that("the batched chain kernel matches the serial path (joint + gap)", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # equilibrate + a later-event gap + reconstruct routes the reconstruction through
  # the batched per-(point, condition) chain kernel (kchunk). DMOD_SYM_NOCHUNK forces
  # the serial per-point fallback; the two must be byte-identical.
  f  <- as.eqnvec(c(R = "kpr - kdg*R + kon*u*R", u = "0"))
  g  <- as.eqnvec(c(y = "scale*R"))
  ev <- eventlist() |>
    addEvent(var = "u", time = 0,  value = "init_u", method = "replace") |>
    addEvent(var = "u", time = 60, value = "0",      method = "replace")
  cg <- data.frame(init_u = 1, row.names = "Ctrl")
  run <- function(nochunk) {
    old <- Sys.getenv("DMOD_SYM_NOCHUNK")
    Sys.setenv(DMOD_SYM_NOCHUNK = if (nochunk) "1" else "")
    on.exit(Sys.setenv(DMOD_SYM_NOCHUNK = old))
    symmetryDetection(f, g, method = "observability", equilibrate = TRUE,
                      events = ev, conditions = cg, forcings = "u",
                      reconstruct = TRUE, reduceCQ = FALSE)
  }
  batched <- run(FALSE)
  serial  <- run(TRUE)
  lineOf <- function(o) sort(vapply(o$symmetries, dMod:::.sym_direction_line,
                                    character(1)))
  expect_identical(batched$rank, serial$rank)
  expect_identical(batched$identifiable, serial$identifiable)
  expect_identical(lineOf(batched), lineOf(serial))
})


test_that("observability recovers the Michaelis-Menten enzyme symmetries", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # the substrate-depletion assay S' = -kcat*Etot*S/(Km+S) read out on an
  # arbitrary scale y = s*S has two structural non-identifiabilities that are not
  # obvious from the equations: the turnover and the enzyme amount enter only as
  # the product Vmax = kcat*Etot, and the molar units of S, Km, Etot trade against
  # the readout scale s. Both are derived in the vignette; here we check the two
  # closed-form directions annihilate exactly the observable invariants
  # {s*S, s*Km, s*kcat*Etot}.
  res <- symmetryDetection(eqnvec(S = "-kcat*Etot*S/(Km + S)"), eqnvec(y = "s*S"),
                           method = "observability", reconstruct = TRUE,
                           reduceCQ = FALSE)
  expect_false(res$identifiable)
  expect_equal(res$rank, 3L)
  expect_equal(res$dim, 5L)
  expect_length(res$symmetries, 2L)
  expect_true(all(vapply(res$symmetries,
                         function(d) isTRUE(d$explicit), logical(1))))

  pt <- list(S = 2, Etot = 3, Km = 5, s = 7, kcat = 11)
  coords <- c("S", "Etot", "Km", "s", "kcat")
  # directional derivative of each observable invariant must vanish on every
  # reported direction (the invariants span the rank-3 identifiable subspace)
  dInv <- function(v) c(
    sS    = pt$s * v["S"]  + pt$S * v["s"],
    sKm   = pt$s * v["Km"] + pt$Km * v["s"],
    sVmax = pt$kcat * pt$Etot * v["s"] + pt$s * pt$Etot * v["kcat"] +
            pt$s * pt$kcat * v["Etot"])
  for (d in res$symmetries)
    expect_equal(unname(dInv(.sym_tangent(d, pt, coords))), c(0, 0, 0),
                 tolerance = 1e-9)
})


test_that("observability finds the transcription-translation rate curve", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # gene expression observed only through the protein: transcription ktx and
  # translation ktl trade on the nonlinear curve ktx*ktl = const (the hidden mRNA
  # scales along it), giving one structural non-identifiability whose closed-form
  # direction must annihilate the observable invariants {ktx*ktl, ktl*m}
  res <- symmetryDetection(eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p"),
                           eqnvec(y = "p"), method = "observability",
                           reconstruct = TRUE, reduceCQ = FALSE)
  expect_false(res$identifiable)
  expect_equal(res$rank, 5L)
  expect_equal(res$dim, 6L)
  expect_length(res$symmetries, 1L)
  d <- res$symmetries[[1]]
  expect_true(isTRUE(d$explicit))
  # the curve lives in the two rates and the hidden mRNA only
  expect_setequal(names(d$generator), c("m", "ktx", "ktl"))

  pt <- list(m = 2, p = 3, ktx = 5, ktl = 7, dm = 11, dp = 13)
  v <- .sym_tangent(d, pt, c("m", "p", "ktx", "ktl", "dm", "dp"))
  # d(ktx*ktl) and d(ktl*m) vanish along the direction
  expect_equal(unname(pt$ktl * v["ktx"] + pt$ktx * v["ktl"]), 0, tolerance = 1e-9)
  expect_equal(unname(pt$m * v["ktl"] + pt$ktl * v["m"]), 0, tolerance = 1e-9)
})


test_that("scaling directions are peeled exactly past the interpolation cap", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # a wide product readout y = s*a*b*c*d*e*x of a decaying state: the decay rate k
  # and the product s*a*b*c*d*e*x0 are identifiable (rank 2 of 8), leaving six
  # scaling non-identifiabilities. Each couples several coordinates, so the dense
  # interpolation would fall back to support-only; the integer-kernel peel returns
  # them all exactly, in both closed-form and support mode.
  f <- eqnvec(x = "-k * x")
  g <- eqnvec(y = "s*a*b*c*d*e*x")
  res <- symmetryDetection(f, g, method = "observability", reconstruct = TRUE,
                           reduceCQ = FALSE)
  expect_false(res$identifiable)
  expect_equal(res$dim, 8L)
  expect_equal(res$rank, 2L)
  expect_length(res$symmetries, 6L)
  expect_true(all(vapply(res$symmetries, function(d)
    isTRUE(d$type == "scaling") && isTRUE(d$explicit), logical(1))))

  # every scaling leaves the product invariant: the signed weights of its factors
  # (the free initial value x0 is the coordinate "x") sum to zero
  factors <- c(s = 1, a = 1, b = 1, c = 1, d = 1, e = 1, x = 1)
  inv <- function(d) {
    w <- setNames(numeric(length(factors)), names(factors))
    for (nm in names(d$weights)) if (nm %in% names(factors))
      w[nm] <- as.numeric(d$weights[[nm]])
    sum(w * factors)
  }
  expect_true(all(vapply(res$symmetries,
                         function(d) abs(inv(d)) < 1e-9, logical(1))))

  # the default (support) mode reports the same scalings in closed form
  sup <- symmetryDetection(f, g, method = "observability", reduceCQ = FALSE)
  expect_true(all(vapply(sup$symmetries, function(d)
    isTRUE(d$type == "scaling") && !is.null(d$generator), logical(1))))
  expect_length(sup$symmetries, 6L)
})


test_that("peeled scalings and per-entry reconstruction combine on two moieties", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # two independent conserved moieties read out through a shared scale: the
  # readout scaling is peeled exactly, and each moiety's conserved direction is a
  # non-scaling rational reconstructed over only its own (narrow) variables
  eq <- eqnlist() |>
    addReaction("A1", "B1", "k1 * A1") |> addReaction("B1", "A1", "k2 * B1") |>
    addReaction("A2", "B2", "k3 * A2") |> addReaction("B2", "A2", "k4 * B2")
  res <- symmetryDetection(eq, eqnvec(y1 = "s * A1", y2 = "s * A2"),
                           method = "observability", reconstruct = TRUE,
                           reduceCQ = FALSE)
  expect_equal(res$rank, 6L)
  expect_equal(res$dim, 9L)
  expect_length(res$symmetries, 3L)
  types <- vapply(res$symmetries, function(d) d$type, character(1))
  expect_equal(sum(types == "scaling"), 1L)
  expect_equal(sum(types == "affine"), 2L)   # each moiety direction is affine once canonicalised
  expect_true(all(vapply(res$symmetries,
                         function(d) isTRUE(d$explicit), logical(1))))

  # each conserved direction satisfies its moiety relation: the canonical affine
  # generator shifts B_i along the conserved total A_i + B_i
  gen <- Filter(function(d) d$type == "affine", res$symmetries)
  for (d in gen) {
    expect_true(any(grepl("^B[12]$", names(d$generator))))
    expect_true(.sym_expr_equal(d$generator[[grep("^B[12]$", names(d$generator),
                                                value = TRUE)]],
                                if ("B1" %in% names(d$generator)) "A1 + B1"
                                else "A2 + B2"))
  }
})


test_that("equilibrate supports a free Hill/power exponent", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  ev <- addEvent(eventlist(), var = "u", time = -1, value = "1", method = "replace")
  cond <- data.frame(var = 1, row.names = "stim")
  # cooperative self-inhibition with a free Hill coefficient nhill, pre-equilibrated
  r <- symmetryDetection(
    eqnvec(p = "kpr/(1+kinh*p^nhill) - dp*p + kin*u", u = "0"), eqnvec(y = "s*p"),
    method = "observability", equilibrate = TRUE, forcings = "u", events = ev,
    conditions = cond, reconstruct = FALSE, reduceCQ = FALSE)
  supp <- unlist(lapply(r$symmetries, function(d) d$support))
  expect_false("nhill" %in% supp)            # the Hill coefficient is identifiable
  expect_equal(r$dim, 6L)
  # a base without a linear turnover term has a steady state that is a fractional
  # root; the inverted recast solves E = x^q and keeps x, log(x) generic. The
  # symmetry is a parameter-weighted scaling recovered exactly over Q(q) by the
  # toric peel (no rational fit): weight 1 - q on dp, integers elsewhere.
  r2 <- symmetryDetection(
    eqnvec(x = "kpr - dp*x^q + kin*u", u = "0"), eqnvec(y = "s*x"),
    method = "observability", equilibrate = TRUE, forcings = "u", events = ev,
    conditions = cond, reconstruct = TRUE, reduceCQ = FALSE)
  expect_length(r2$symmetries, 1L)
  d2 <- r2$symmetries[[1]]
  expect_true(isTRUE(d2$explicit))
  expect_equal(d2$type, "scaling")                   # a q-weighted toric scaling
  expect_false("q" %in% d2$support)                  # the exponent is identifiable
  # the free-exponent recast makes dp's scaling weight 1 - q; the readout scale s and
  # the synthesis rates carry integer weights, and the input rate kin scales along
  expect_true(.sym_expr_equal(d2$weights[["dp"]], "1 - q"))
  expect_true(.sym_expr_equal(d2$weights[["kpr"]], "1"))
  expect_true(.sym_expr_equal(d2$weights[["s"]], "-1"))
  exprs2 <- as.character(unlist(d2$weights))
  expect_false(any(grepl("log\\(", exprs2)))         # the direction is rational
  expect_false(any(grepl("_E_|_L_", exprs2)))         # no internal recast symbol leaks
  # a free exponent WITHOUT equilibrate is handled by the transient recast: E = p^nhill
  # and L = log(p) become free-initial-value coordinates tied to (p, nhill) by the
  # recast relation, so the observability system stays rational. The transient model
  # (free initial value p, no steady state) is MORE identifiable than the equilibrated
  # one: the single symmetry is the parameter-weighted scaling p->L*p, kpr->L*kpr,
  # kinh->kinh*L^-nhill, s->s/L, recovered in closed form with kinh weighted by nhill.
  r3 <- symmetryDetection(eqnvec(p = "kpr/(1+kinh*p^nhill) - dp*p"),
                          eqnvec(y = "s*p"), method = "observability",
                          equilibrate = FALSE, reduceCQ = FALSE, reconstruct = TRUE)
  expect_equal(r3$dim, 6L)                            # p(0), kpr, kinh, dp, nhill, s
  expect_equal(r3$rank, 5L)
  expect_length(r3$symmetries, 1L)
  d3 <- r3$symmetries[[1]]
  expect_true(isTRUE(d3$explicit))
  expect_false("nhill" %in% d3$support)              # the Hill exponent is identifiable
  expect_setequal(d3$support, c("kinh", "kpr", "p", "s"))
  gen3 <- as.character(unlist(d3$generator))
  expect_true(any(grepl("nhill", gen3)))             # kinh carries the nhill-weight
  expect_false(any(grepl("_E_|_L_", gen3)))          # no internal recast symbol leaks
  expect_false(any(grepl("log\\(", gen3)))           # this direction is rational in nhill
})


test_that("equilibrate reconstructs Hill directions in closed form via the recast", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  ev <- addEvent(eventlist(), var = "u", time = -1, value = "var_u", method = "replace")
  cond <- data.frame(var_u = c(0, 1), row.names = c("ctrl", "stim"))
  # two switch values make the Hill self-inhibition non-identifiable; the recast
  # coordinates E = p^nhill, L = log(p) are sampled and back-substituted so the
  # directions, including the transcendental ones, are returned in closed form
  r <- symmetryDetection(
    eqnvec(p = "kpr/(1+kinh*p^nhill) - dp*p + kin*u", u = "0"), eqnvec(y = "s*p"),
    method = "observability", equilibrate = TRUE, events = ev,
    conditions = cond, reconstruct = TRUE, reduceCQ = FALSE)
  expect_true(all(vapply(r$symmetries,
                         function(d) isTRUE(d$explicit), logical(1))))
  exprs <- unlist(lapply(r$symmetries, function(d) as.character(unlist(d$generator))))
  expect_true(any(grepl("p\\*\\*nhill", exprs)))   # the power E -> base^exp
  expect_true(any(grepl("log\\(p\\)", exprs)))      # the transcendental L -> log(base)
  expect_false(any(grepl("_E_|_L_", exprs)))        # no internal recast symbol leaks
})


# The Hill term above has a STATE base (p^nhill). A Michaelis constant K^n is a
# PARAMETER base -- it appears only under the exponent -- which the recast handles
# differently (the "Michaelis co-scaling" relation, xi_base = -exp/base and -1/base).
# That case had no coverage, which is why the parameter-base regression slipped
# through. These two tests guard the verdict and the ground-truth direction.
test_that("a parameter-base (Michaelis) Hill exponent is non-identifiable (dim 7)", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  hill <- eqnlist() |>
    addReaction("0",  "FB", "k_pr_FB")                     |>
    addReaction("FB", "0",  "d_FB * FB")                   |>
    addReaction("0",  "x",  "k_pr_x * K^n / (K^n + FB^n)") |>
    addReaction("x",  "0",  "d_x * x")
  obs <- eqnvec(xobs = "scale * x")
  # verdict: K (a Michaelis constant under the exponent) is a coordinate and is
  # non-identifiable, so dim is 7. A regression that drops K to a non-coordinate
  # (the pre-Michaelis-fix bug) gives dim 6 and misses the K non-identifiability.
  r <- symmetryDetection(hill, obs, method = "observability",
                         equilibrate = TRUE, reduceCQ = FALSE)
  expect_false(r$identifiable)
  expect_equal(r$dim, 7L)
  supp <- unlist(lapply(r$symmetries, function(d) d$support))
  expect_true("K" %in% supp)
  expect_true("n" %in% supp)

  # ground truth from the symbolic engine (base^n stays an atom, no recast, so it is
  # independent of the modular recast): n trades against K with the log(base) factor,
  # xi_n = 1, xi_K = (K/n) log((k_pr_FB/d_FB)/K).
  ss <- eqnvec(FB = "k_pr_FB/d_FB",
               x  = "k_pr_x*K^n/(d_x*(K^n + (k_pr_FB/d_FB)^n))")
  sres <- symmetryDetection(hill, obs, method = "observability",
                            symEngine = "symbolic", trafo = ss)
  nd <- Filter(function(d) "n" %in% d$support, sres$symmetries)
  expect_length(nd, 1L)
  v <- nd[[1]]$generator
  expect_true(.sym_expr_equal(v[["n"]], "1"))
  expect_true(.sym_expr_equal(v[["K"]], "(-K*log(K) + K*log(k_pr_FB/d_FB))/n"))
})


test_that("modular reconstruct recovers the parameter-base (Michaelis) Hill direction", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # The forward reconstruction closes the leaf-dependent residual recast directions
  # (each on a distinct physical anchor); the n-direction comes back in closed form
  # with its log(base) factor. The reported gauge differs from the symbolic engine's
  # (a valid alternative basis of the same residual space), so this only asserts that
  # every direction closes and a log factor appears -- not the exact representative.
  hill <- eqnlist() |>
    addReaction("0",  "FB", "k_pr_FB")                     |>
    addReaction("FB", "0",  "d_FB * FB")                   |>
    addReaction("0",  "x",  "k_pr_x * K^n / (K^n + FB^n)") |>
    addReaction("x",  "0",  "d_x * x")
  r <- symmetryDetection(hill, eqnvec(xobs = "scale * x"), method = "observability",
                         equilibrate = TRUE, reconstruct = TRUE, reduceCQ = FALSE)
  expect_true(all(vapply(r$symmetries,
                         function(d) isTRUE(d$explicit), logical(1))))
  exprs <- unlist(lapply(r$symmetries, function(d) as.character(unlist(d$generator))))
  expect_true(any(grepl("log\\(", exprs)))          # the n-direction's log(base) factor
})


test_that("equilibrate without reduceCQ uses the held-variable moiety parameterisation", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # A conserved moiety A + B under equilibrate makes f = 0 rank-deficient. reduceCQ
  # eliminates the pivot and reports the moiety freedom on a `total` parameter;
  # reduceCQ = FALSE instead holds one pivot's resting value under the pivot's own
  # name (dMod/deSolve convention: a parameter named like a state is its initial
  # value), keeps every species a coordinate, and reports the SAME identifiability
  # (rank/dim/direction count) with the pivot's initial value in place of the total.
  # u modulates the A->B rate (conservation preserved, transient informative). The
  # pool is a shared parameter, so a second condition does NOT double the direction.
  eq <- eqnlist() |>
    addReaction("A", "B", "k1 * A * (1 + kin*u)") |>
    addReaction("B", "A", "k2 * B") |>
    addReaction("0", "u", "0")
  obs <- eqnvec(y = "scale * A")
  ev  <- addEvent(eventlist(), var = "u", time = 0, value = "dose", method = "replace")
  cond <- data.frame(dose = c(1, 2), row.names = c("d1", "d2"))
  args <- list(eq, obs, method = "observability", equilibrate = TRUE,
               forcings = "u", events = ev, conditions = cond, reconstruct = TRUE)

  rT <- do.call(symmetryDetection, c(args, list(reduceCQ = TRUE)))
  rF <- do.call(symmetryDetection, c(args, list(reduceCQ = FALSE)))

  # same identifiability structure, just a different parameterisation of the moiety
  expect_equal(rF$rank, rT$rank)
  expect_equal(rF$dim, rT$dim)
  expect_equal(length(rF$symmetries), length(rT$symmetries))

  suppF <- unlist(lapply(rF$symmetries, function(d) d$support))
  suppT <- unlist(lapply(rT$symmetries, function(d) d$support))
  expect_true("A" %in% suppF)                      # held: pivot initial value appears
  expect_false(any(grepl("^total", suppF)))        # held: no `total` coordinate
  expect_true(any(grepl("^total", suppT)))         # reduceCQ: the total does appear
  expect_false("A" %in% suppT)                     # reduceCQ: pivot A is eliminated

  # every direction closes in the held parameterisation too: the pool scaling is
  # peeled onto the pivot's initial-value parameter (weight matching the total's), not
  # left as an unreconstructed support-only direction
  expect_true(all(vapply(rF$symmetries,
                         function(d) isTRUE(d$explicit), logical(1))))
  pool <- Filter(function(d) "A" %in% d$support, rF$symmetries)
  expect_true(length(pool) >= 1L)
  expect_true(.sym_expr_equal(pool[[1]]$generator[["A"]], "A"))     # weight +1

  # a single condition works too, and does not error on the rank-deficient f = 0
  cond1 <- data.frame(dose = 1, row.names = "d1")
  r1 <- do.call(symmetryDetection,
                c(list(eq, obs, method = "observability", equilibrate = TRUE,
                       forcings = "u", events = ev, conditions = cond1,
                       reconstruct = TRUE, reduceCQ = FALSE)))
  expect_true("A" %in% unlist(lapply(r1$symmetries, function(d) d$support)))
})


test_that("the verification gate accepts correct and rejects wrong closed forms", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  sd <- .sd_module()

  # exact rational evaluation over GF(q), checked against an independent R inverse
  # at a small prime where the arithmetic stays exact in doubles
  qs <- 100003
  modinv <- function(a, m) {
    a <- a %% m; t0 <- 0; t1 <- 1; r0 <- m; r1 <- a
    while (r1 != 0) { qq <- r0 %/% r1
      h <- t0 - qq * t1; t0 <- t1; t1 <- h
      h <- r0 - qq * r1; r0 <- r1; r1 <- h }
    ((t0 %% m) + m) %% m
  }
  ref <- ((-18 %% qs) * modinv(13, qs)) %% qs
  expect_equal(dMod:::.sym_eval_modq("-(A + B)/k2", list(A = 7, B = 11, k2 = 13),
                                     qs, sd), as.integer(ref))
  expect_equal(dMod:::.sym_eval_modq("-1", list(A = 7), qs, sd), as.integer(qs - 1))

  # a fake reduction whose free column 1 has null vector (a : -3, b : 1); the
  # matching closed form certifies, a wrong one is rejected
  fakeKcall <- function(point, p, Nt)
    list(ok = TRUE, R = matrix(c(1L, 3L), 1, 2), pivots = 0L, rank = 1L, dim = 2L)
  znames <- c("a", "b"); leafNames <- c("a", "b"); point0 <- c(7, 11)
  good <- list(vector = list(a = "-3", b = "1"), type = "general", reconstruct = TRUE)
  bad  <- list(vector = list(a = "-5", b = "1"), type = "general", reconstruct = TRUE)
  expect_true(dMod:::.sym_verify_direction(good, 1L, znames, leafNames, point0, 1L,
                                           fakeKcall, sd))
  expect_false(dMod:::.sym_verify_direction(bad, 1L, znames, leafNames, point0, 1L,
                                            fakeKcall, sd))
})


test_that("symSparsePoly recovers sparse polynomials by Ben-Or-Tiwari", {
  # a small prime keeps the sequence construction exact in R doubles (products of
  # two ~2^31 residues would exceed 2^53); the kernel itself is u64-exact for any
  # prime, so this validates the algorithm without R precision loss
  p <- 100003
  powmodR <- function(a, e, m) {
    r <- 1; a <- a %% m
    while (e > 0) { if (e %% 2 == 1) r <- (r * a) %% m; a <- (a * a) %% m; e <- e %/% 2 }
    r
  }
  # exact recovery of an arbitrary sparse polynomial in `nv` variables
  recover <- function(terms, bases, len, degree) {
    nv <- length(bases)
    mval <- function(e) { v <- 1; for (j in seq_len(nv)) v <- (v * powmodR(bases[j], e[j], p)) %% p; v }
    seqv <- vapply(0:(len - 1L), function(k) {
      acc <- 0
      for (tt in terms) acc <- (acc + (tt$c %% p) * powmodR(mval(tt$e), k, p)) %% p
      acc
    }, numeric(1))
    g <- as.matrix(expand.grid(rep(list(0:degree), nv)))
    g <- g[rowSums(g) <= degree, , drop = FALSE]
    monoRes <- apply(g, 1, mval)
    res <- dMod:::symSparsePoly(as.integer(seqv), matrix(as.integer(g), ncol = nv),
                                as.integer(monoRes), p)
    res$got <- if (res$status == "ok" && res$nterms > 0)
      data.frame(res$exps, c = res$coeffs) else NULL
    res
  }
  key <- function(df) { o <- do.call(order, df); df[o, , drop = FALSE] }

  # 3*x^2*y + 5*y^3 + 7
  r1 <- recover(list(list(e = c(2, 1), c = 3), list(e = c(0, 3), c = 5),
                     list(e = c(0, 0), c = 7)), c(2, 3), 8L, 3L)
  expect_identical(r1$status, "ok"); expect_equal(r1$nterms, 3L)
  expect_equal(unname(as.matrix(key(r1$got))),
               unname(as.matrix(key(data.frame(X1 = c(0, 0, 2), X2 = c(0, 3, 1),
                                               c = c(7, 5, 3))))))

  # a single constant term
  r2 <- recover(list(list(e = c(0, 0), c = 9)), c(2, 3), 6L, 2L)
  expect_identical(r2$status, "ok"); expect_equal(r2$nterms, 1L)
  expect_equal(r2$coeffs, 9L)

  # five terms in three variables, recovered exactly
  r3 <- recover(list(list(e = c(1, 0, 0), c = 2), list(e = c(0, 2, 1), c = -3),
                     list(e = c(3, 0, 0), c = 4), list(e = c(0, 0, 2), c = 5),
                     list(e = c(1, 1, 1), c = -6)), c(2, 3, 5), 12L, 4L)
  expect_identical(r3$status, "ok"); expect_equal(r3$nterms, 5L)

  # too short a sequence for the recovered order -> needmore, never a wrong answer
  nv <- 2; bases <- c(2, 3)
  mval <- function(e) (powmodR(bases[1], e[1], p) * powmodR(bases[2], e[2], p)) %% p
  short <- vapply(0:2, function(k) (3 * powmodR(mval(c(2, 1)), k, p) +
                                      5 * powmodR(mval(c(0, 3)), k, p) +
                                      7 * powmodR(mval(c(0, 0)), k, p)) %% p, numeric(1))
  g <- as.matrix(expand.grid(0:3, 0:3)); g <- g[rowSums(g) <= 3, , drop = FALSE]
  r4 <- dMod:::symSparsePoly(as.integer(short), matrix(as.integer(g), ncol = 2),
                             as.integer(apply(g, 1, mval)), p)
  expect_identical(r4$status, "needmore")
})


test_that("symMonoResidues computes modular monomial values including Laurent", {
  p <- 100003; b <- c(2, 3, 5)
  modinv <- function(a, m) {
    a <- a %% m; t0 <- 0; t1 <- 1; r0 <- m; r1 <- a
    while (r1 != 0) { qq <- r0 %/% r1
      h <- t0 - qq * t1; t0 <- t1; t1 <- h
      h <- r0 - qq * r1; r0 <- r1; r1 <- h }
    ((t0 %% m) + m) %% m
  }
  E <- matrix(c(1, 0, 0,  0, 1, 0,  2, 1, 0,  1, 0, -1,  0, 0, -1),
              ncol = 3, byrow = TRUE)
  got <- dMod:::symMonoResidues(matrix(as.integer(E), ncol = 3), as.integer(b), p)
  i5 <- modinv(5, p)
  expect_equal(got, as.integer(c(2, 3, (4 * 3) %% p, (2 * i5) %% p, i5)))
})


test_that(".sym_sparse_entry recovers a wide Laurent entry exactly", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  sd <- .sd_module()
  # a degree-1 numerator over a single-monomial denominator, 8 relevant variables
  # (beyond the dense cap): -(v1 + ... + v7)/v8
  vars <- paste0("v", 1:8)
  target <- "-(v1 + v2 + v3 + v4 + v5 + v6 + v7)/v8"
  nz <- 9L; f <- 8L; sc <- 0L
  # a fake reduction whose free column f carries `target` at support column sc
  fakeKcall <- function(point, p, Nt) {
    tv <- sd$evalRationalMod(target, as.list(vars), as.list(as.numeric(point[1:8])),
                             as.integer(p))
    if (is.null(tv)) return(list(ok = FALSE))
    R <- matrix(0L, 1L, nz)
    R[1, f + 1L] <- as.integer((p - (as.numeric(tv) %% p)) %% p)
    list(ok = TRUE, R = R, pivots = 0L, rank = 1L, dim = nz)
  }
  e <- dMod:::.sym_sparse_entry(1:8, sc, f, rep(1, 8), vars, 1L, fakeKcall, 0L)
  expect_false(is.null(e))
  expect_true(.sym_expr_equal(e, target))
})


test_that(".sym_general_rational_entry recovers a multi-term denominator", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  sd <- .sd_module()
  # a multi-term denominator with no constant term, 8 relevant variables:
  # (v1 + ... + v6 + 1) / (v7 + 2*v8)
  vars <- paste0("v", 1:8)
  target <- "(v1 + v2 + v3 + v4 + v5 + v6 + 1)/(v7 + 2*v8)"
  nz <- 9L; f <- 8L; sc <- 0L
  fakeKcall <- function(point, p, Nt) {
    tv <- sd$evalRationalMod(target, as.list(vars), as.list(as.numeric(point[1:8])),
                             as.integer(p))
    if (is.null(tv)) return(list(ok = FALSE))
    R <- matrix(0L, 1L, nz)
    R[1, f + 1L] <- as.integer((p - (as.numeric(tv) %% p)) %% p)
    list(ok = TRUE, R = R, pivots = 0L, rank = 1L, dim = nz)
  }
  point0 <- c(7, 11, 13, 17, 19, 23, 29, 31)
  e <- dMod:::.sym_general_rational_entry(1:8, sc, f, point0, vars, 1L, fakeKcall, 0L)
  expect_false(is.null(e))
  expect_true(.sym_expr_equal(e, target))
})


test_that("log-coordinate gauge reconstructs a parameter-weighted scaling", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  sd  <- .sd_module()
  spy <- reticulate::import("sympy", convert = TRUE)

  # one weighted scaling xi = (a: -n*a, b: b) whose weight is the model parameter
  # n. The free-column gauge carries the rational n*a/b; the log gauge carries the
  # constant weight, so back-substitution recovers the sparse c*z form.
  znames <- c("a", "b"); nz <- 2L; zSlots <- c(0L, 1L)
  leafNames <- c("a", "b", "n")
  P <- dMod:::.symPrimes[1]
  # free column b (=1) carries v_a = -n*a/b at the single pivot row (column a = 0)
  fakeKcall <- function(point, p, Nt) {
    val <- sd$evalRationalMod("n*a/b", as.list(leafNames),
                              as.list(as.numeric(point[1:3])), as.integer(p))
    if (is.null(val)) return(list(ok = FALSE))
    R <- matrix(0L, 1L, nz)
    R[1, 2L] <- as.integer(as.numeric(val) %% p)
    list(ok = TRUE, R = R, pivots = 0L, rank = 1L, dim = nz)
  }
  point0 <- c(a = 7, b = 11, n = 13)
  ref <- fakeKcall(point0, P, 1L)
  zvals0 <- as.numeric(point0[zSlots + 1L])

  lg <- dMod:::.sym_logcoord_gauge(1L, matrix(0L, 0L, nz), P, nz, list(ref = ref),
                                   zvals0)
  expect_equal(length(lg$anchors), 1L)

  pool <- dMod:::.sym_pool()
  dir <- dMod:::.sym_interpolate_direction(
    lg$anchors[1], ref, ref$pivots, znames, zSlots, leafNames, 3L,
    as.numeric(point0), pool, 100L, 1L, fakeKcall, spy, NULL,
    lg$residueFns[[1]], reconstControl())
  e <- dir$entry
  expect_true(isTRUE(e$closedForm))
  e$vector <- dMod:::.sym_logcoord_backsub(e$vector, spy)

  # back-transformed entries are the sparse c*z form (the b entry couples only n)
  expect_true(.sym_expr_equal(e$vector[["a"]], "a"))
  expect_true(.sym_expr_equal(e$vector[["b"]], "-b/n"))

  # the direction lies in the nullspace; a corrupted back-transform is rejected
  expect_true(dMod:::.sym_verify_in_nullspace(
    e, lg$anchors[1], znames, leafNames, as.numeric(point0), 1L, fakeKcall, pool,
    200L, nz, sd))
  bad <- e; bad$vector[["a"]] <- paste0("2*(", e$vector[["a"]], ")")
  expect_false(dMod:::.sym_verify_in_nullspace(
    bad, lg$anchors[1], znames, leafNames, as.numeric(point0), 1L, fakeKcall, pool,
    300L, nz, sd))
})


test_that("sparse reconstruction matches the dense path when forced", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # route a genuine rational entry (the conserved direction -(A + B)/k2) through
  # the sparse path by dropping the dense cap; it must reproduce the dense result
  # (both report the same canonical affine generator: dB = A + B, dk1 = k2)
  eq <- eqnlist() |> addReaction("A", "B", "k1 * A") |> addReaction("B", "A", "k2 * B")
  res <- symmetryDetection(eq, eqnvec(Aobs = "alpha * A"), method = "observability",
                           reconstruct = TRUE, reduceCQ = FALSE,
                           control = reconstControl(relevanceCap = 0L))
  d <- Filter(function(x) all(c("B", "k1", "k2") %in% names(x$generator)) &&
                isTRUE(x$explicit), res$symmetries)
  expect_gte(length(d), 1L)
  expect_true(.sym_expr_equal(d[[1]]$generator[["B"]], "A + B"))
  expect_true(.sym_expr_equal(d[[1]]$generator[["k1"]], "k2"))
})


test_that("multi-condition observability identifies a switch-gated rate", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # u gates the rate k1 + u*k2; one switch value cannot separate k1 and k2, but
  # two values (set per condition by the switch event) identify both
  f  <- eqnvec(A = "-(k1 + u*k2)*A", u = "0")
  g  <- eqnvec(y = "A")
  ev <- addEvent(eventlist(), var = "u", time = -1, value = "var_u", method = "replace")
  cg <- data.frame(var_u = c(0, 1), row.names = c("ctrl", "stim"))

  multi <- symmetryDetection(f, g, method = "observability",
                             events = ev, conditions = cg, reduceCQ = FALSE)
  expect_true(multi$identifiable)
  expect_equal(multi$info$conditions, 2L)
  expect_false(any(vapply(multi$symmetries,
                          function(d) "u" %in% d$support, logical(1))))
})


test_that("multi-condition observability spans a union parameter space", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # a condition-specific rate name (a knockdown reparametrisation): control uses
  # k, the perturbed condition k_knd; both rates and the shared free initial
  # value are identifiable across the union of conditions
  f  <- eqnvec(A = "-k*A")
  g  <- eqnvec(y = "A")
  cg <- data.frame(k = c("k", "k_knd"), row.names = c("ctrl", "knd"),
                   stringsAsFactors = FALSE)

  res <- symmetryDetection(f, g, method = "observability", conditions = cg,
                           reduceCQ = FALSE)
  expect_equal(res$dim, 3L)
  expect_true(res$identifiable)
})


test_that("multi-condition observability keeps a confounder common to all conditions", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # the calibration scale s and the free initial value A(0) confound in every
  # condition (only s*A is seen); the inert switch u never enters a direction
  f  <- eqnvec(A = "-k*A", u = "0")
  g  <- eqnvec(y = "s*A")
  ev <- addEvent(eventlist(), var = "u", time = -1, value = "var_u", method = "replace")
  cg <- data.frame(var_u = c(0, 1), row.names = c("a", "b"))

  res <- symmetryDetection(f, g, method = "observability", events = ev,
                           conditions = cg, reduceCQ = FALSE)
  expect_false(res$identifiable)
  supp <- unlist(lapply(res$symmetries, function(d) d$support))
  expect_true(all(c("A", "s") %in% supp))
  expect_false("u" %in% supp)
})


test_that("a post-t0 event opens a second segment, propagated exactly", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # decay observed on an unknown scale; a dose at t0 and a second add-dose at t = 1
  # open two segments. The state is propagated exactly across the gap (no free carry
  # coordinates), so the coordinate space is just A(0), s, k. The second dose's jump
  # in the readout is s * d1, which reveals the scale s; with both doses known the
  # model is then fully identifiable. (The second-segment dose must reach the kernel
  # for this: it is applied to the propagated state at the segment boundary.)
  f  <- eqnvec(A = "-k*A")
  g  <- eqnvec(y = "s*A")
  ev <- eventlist() |>
    addEvent(var = "A", time = 0, value = "d0", method = "add") |>
    addEvent(var = "A", time = 1, value = "d1", method = "add")
  cg <- data.frame(d0 = 1, d1 = 1, row.names = "c1")

  r <- symmetryDetection(f, g, method = "observability", events = ev,
                         conditions = cg, reconstruct = TRUE)
  expect_equal(r$info$conditions, 1L)
  expect_equal(r$info$segments, 2L)
  expect_equal(r$dim, 3L)             # A(0), s, k -- no free carry coordinates
  expect_true(r$identifiable)         # the second dose's jump reveals the scale s
})


test_that("with no post-t0 events the analysis collapses to a single segment", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # only a t0 dose: exactly one segment, identical to the pre-segment behaviour
  # (the scale s and the free initial value A(0) confound)
  f  <- eqnvec(A = "-k*A")
  g  <- eqnvec(y = "s*A")
  ev <- addEvent(eventlist(), var = "A", time = 0, value = "d0", method = "add")
  cg <- data.frame(d0 = 1, row.names = "c1")

  r <- symmetryDetection(f, g, method = "observability", events = ev,
                         conditions = cg)
  expect_equal(r$info$segments, 1L)
  expect_equal(r$dim, 3L)             # A(0), s, k
  expect_false(r$identifiable)
})


test_that("exact propagation identifies a transient-channel parameter", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # An inhibitor pre-incubation [-30, 0) feeds the stimulus phase. The inhibition
  # strength kinh affects only the inhibitor-relaxed boundary state at t = 0, not
  # any directly observed quantity. The expansion starts at the earliest event
  # (-30) and the state is propagated exactly across the gap to the stimulus, so
  # kinh becomes identifiable through that transient channel.
  f <- eqnvec(x = "kpr/(1 + kinh*inh) - kdeg*x + kstim*stim", inh = "0", stim = "0")
  g <- eqnvec(y = "s*x")
  ev <- eventlist() |>
    addEvent(var = "inh",  time = -30, value = "1", method = "replace") |>
    addEvent(var = "stim", time = 0,   value = "1", method = "replace")
  cg <- data.frame(row.names = "c1")

  r30 <- symmetryDetection(f, g, method = "observability", equilibrate = TRUE,
                           events = ev, conditions = cg,
                           forcings = c("inh", "stim"))
  supp30 <- unlist(lapply(r30$symmetries, function(d) d$support))
  expect_equal(r30$info$segments, 2L)
  expect_true(r30$info$gapOrderUsed >= 1L)       # the gap series carries kinh
  expect_false("kinh" %in% supp30)          # propagated: kinh identifiable
})


test_that("equilibrate seeds the first segment and propagates to a later one", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # S0 is seeded from the resting steady state; S1 is seeded by propagating S0
  # exactly across the gap and applying the t = 1 add-dose. A known t0 dose pins
  # the scale, then the relaxation pins the rate and the production.
  f  <- eqnvec(x = "b - a*x")
  g  <- eqnvec(y = "s*x")
  ev <- eventlist() |>
    addEvent(var = "x", time = 0, value = "dose", method = "replace") |>
    addEvent(var = "x", time = 1, value = "d1",   method = "add")
  cg <- data.frame(dose = 2, d1 = 1, row.names = "stim")

  r <- symmetryDetection(f, g, method = "observability", equilibrate = TRUE,
                         events = ev, conditions = cg)
  expect_equal(r$info$conditions, 1L)
  expect_equal(r$info$segments, 2L)
  expect_true(r$identifiable)         # known dose pins s, then b; a from relaxation
})


test_that("a free-Hill-exponent feedback closes as a weighted scaling", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # An unobserved feedback species FB inhibits the observed node R through a Hill
  # term with a FREE exponent nh, and is itself produced from R (a loop). Rescaling
  # FB's synthesis k_pf is compensated by k_fb with weight -nh (the invariant is
  # k_fb * k_pf^nh), a scaling whose weight is a parameter. It is not an integer
  # scaling and, under equilibrate, its free-column representative couples the whole
  # loop; it is recovered from its minimal support in log coordinates as the exact
  # monomial entry k_fb = -k_fb * nh / k_pf.
  f <- eqnvec(R  = "k_pr/(1 + k_fb*FB^nh) - k_dg*R + k_stim*u",
              FB = "k_pf*R - k_df*FB", u = "0")
  g <- eqnvec(R_obs = "s*R")
  r <- symmetryDetection(f, g, method = "observability", equilibrate = TRUE,
                         forcings = "u",
                         events = addEvent(eventlist(), var = "u", time = 0,
                                           value = "dose", method = "replace"),
                         conditions = data.frame(dose = 1, row.names = "c1"),
                         reconstruct = TRUE, reduceCQ = FALSE)
  # every non-identifiable direction is in closed form
  expect_true(all(vapply(r$symmetries,
                         function(d) isTRUE(d$explicit), logical(1))))
  # the weighted-scaling direction couples exactly the inhibition strength and the
  # feedback synthesis rate, with the Hill exponent as the weight
  hill <- Filter(function(d) setequal(d$support, c("k_fb", "k_pf")),
                 r$symmetries)
  expect_length(hill, 1L)
  expect_true(isTRUE(hill[[1]]$explicit))
  # the Hill exponent is the weight, so it appears in the direction (in whichever
  # entry is not the normalised anchor)
  expect_match(paste(unlist(hill[[1]]$weights), collapse = " "), "nh")
})

test_that("implicit steady state finds non-scaling multi-condition directions", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # two states with a non-monomial resting state x_i* = v_i*K/(k - v_i), observed
  # only through the SUM s*(x1+x2) so neither is individually pinned; two conditions
  # differ in k (enters as a DIFFERENCE), so a valid direction can move the resting
  # states non-proportionally across conditions. A single shared state column would
  # miss these and over-report identifiability; the per-condition joint coordinates
  # of the implicit determining system report the full rank-2/6 verdict (four
  # non-identifiable directions, two of them non-scaling).
  f <- eqnvec(x1 = "v1 - k*x1/(K + x1)", x2 = "v2 - k*x2/(K + x2)")
  g <- eqnvec(y = "s*(x1 + x2)")
  grid <- data.frame(k = c("k1", "k2"), row.names = c("c1", "c2"),
                     stringsAsFactors = FALSE)
  r <- symmetryDetection(f, g, method = "observability", equilibrate = TRUE,
                         conditions = grid, reduceCQ = FALSE, reconstruct = FALSE)
  expect_false(r$identifiable)
  expect_equal(r$dim, 6L)
  expect_equal(r$rank, 2L)
  expect_length(r$symmetries, 4L)
})

test_that("a per-condition trafo list matches the equivalent condition grid", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  f <- eqnvec(A = "-p*A"); g <- eqnvec(y = "s*A")
  # a symbol-rename per condition, once as a grid and once as a trafo list
  r_grid <- symmetryDetection(f, g, method = "observability",
    conditions = data.frame(p = c("pa", "pb"), row.names = c("c1", "c2")),
    reduceCQ = FALSE)
  r_traf <- symmetryDetection(f, g, method = "observability",
    trafo = list(c1 = eqnvec(p = "pa"), c2 = eqnvec(p = "pb")), reduceCQ = FALSE)
  expect_equal(r_traf$info$conditions, 2L)          # the list alone sets the conditions
  expect_equal(r_traf$rank, r_grid$rank)
  expect_equal(r_traf$dim, r_grid$dim)

  # a numeric bake per condition matches a numeric grid cell
  r_bake  <- symmetryDetection(f, g, method = "observability",
    trafo = list(eqnvec(p = "1/2"), eqnvec(p = "3/2")), reduceCQ = FALSE)
  r_bgrid <- symmetryDetection(f, g, method = "observability",
    conditions = data.frame(p = c(0.5, 1.5), row.names = c("c1", "c2")),
    reduceCQ = FALSE)
  expect_equal(r_bake$rank, r_bgrid$rank)
  expect_equal(r_bake$dim, r_bgrid$dim)
})

test_that("trafo = steadyStates() matches equilibrate (explicit vs implicit route)", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  withr::local_dir(withr::local_tempdir())

  # x' = b - a*x has the rational resting state x* = b/a. The explicit route
  # (steadyStates() fed as a trafo initial condition) and the implicit route
  # (equilibrate = TRUE, the f = 0 tangency constraint) must agree exactly.
  eq <- eqnlist() |> addReaction("", "x", "b") |> addReaction("x", "", "a*x")
  g  <- eqnvec(y = "s*x")
  ss <- steadyStates(eq, testSteady = "fast")
  skip_if_not(is.character(ss) && identical(unname(ss[["x"]]), "b/a"))

  r_trafo <- symmetryDetection(eq, g, method = "observability", trafo = ss,
                               reconstruct = TRUE, reduceCQ = FALSE)
  r_equil <- symmetryDetection(eq, g, method = "observability", equilibrate = TRUE,
                               reconstruct = TRUE, reduceCQ = FALSE)
  expect_equal(r_trafo$rank, r_equil$rank)
  expect_equal(r_trafo$dim, r_equil$dim)
  expect_equal(r_trafo$identifiable, r_equil$identifiable)
  # same support set on both routes
  supp <- function(r) sort(unique(unlist(lapply(r$symmetries, function(d) d$support))))
  expect_equal(supp(r_trafo), supp(r_equil))
})

test_that("scaling uses events (fixed weight) and conditions (lattice intersection)", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  eq <- .canonical(); g <- eqnvec(Aobs = "alpha * A")

  # baseline: the calibration scaling A, B, alpha
  base <- symmetryDetection(eq, g, method = "scaling", reduceCQ = FALSE)
  expect_equal(length(base$symmetries), 1L)

  # a known dose pins A's absolute value, so it can no longer scale (weight 0)
  dose <- addEvent(eventlist(), var = "A", time = 0, value = "2", method = "replace")
  pinned <- symmetryDetection(eq, g, method = "scaling", events = dose, reduceCQ = FALSE)
  expect_equal(length(pinned$symmetries), 0L)

  # multi-condition via a trafo list: a second condition that fixes alpha drops the
  # shared scaling (the intersection of the per-condition lattices)
  drop <- symmetryDetection(eq, g, method = "scaling", reduceCQ = FALSE,
                            trafo = list(c1 = eqnvec(k1 = "k1"), c2 = eqnvec(alpha = "1")))
  expect_equal(length(drop$symmetries), 0L)

  # multi-condition where both conditions keep it: the scaling survives
  keep <- symmetryDetection(eq, g, method = "scaling", reduceCQ = FALSE,
                            conditions = data.frame(k1 = c("ka", "kb"),
                                                    row.names = c("c1", "c2")))
  expect_equal(length(keep$symmetries), 1L)
})

test_that("the toric peel recovers a parameter-weighted (Hill) scaling over Q(nhill)", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # a free Hill exponent makes the inhibition k_inh * p^nhill invariant under
  # p -> lam*p, k_inh -> lam^(-nhill)*k_inh: a scaling whose WEIGHT is the exponent.
  # The toric peel imposes c_E = nhill*c_base and recovers it exactly over Q(nhill),
  # as a scaling with weight -nhill on k_inh -- no rational fit, no sampling.
  ev <- addEvent(eventlist(), var = "u", time = -1, value = "var_u", method = "replace")
  cond <- data.frame(var_u = c(0, 1), row.names = c("ctrl", "stim"))
  r <- symmetryDetection(
    eqnvec(p = "kpr/(1 + kinh*p^nhill) - dp*p + kin*u", u = "0"), eqnvec(y = "s*p"),
    method = "observability", equilibrate = TRUE, events = ev, conditions = cond,
    reconstruct = TRUE, reduceCQ = FALSE)
  hill <- Filter(function(d) isTRUE(d$type == "scaling") &&
                   "kinh" %in% names(d$generator), r$symmetries)
  expect_length(hill, 1L)
  expect_true(isTRUE(hill[[1]]$explicit))
  expect_true(.sym_expr_equal(hill[[1]]$weights[["kinh"]], "-nhill"))  # the exponent weight
})


test_that("both observability engines report the same canonical generator", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # gene expression observed at the protein level: transcription and translation
  # are confined to the hyperbola ktx*ktl = const. The modular engine peels this
  # as an exact integer scaling; the symbolic engine reconstructs it as a rational
  # ("disguised") direction. The classifier canonicalises both to the same scaling
  # generator, so the two engines must now agree exactly.
  gene <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
  mod <- symmetryDetection(gene, eqnvec(y = "p"), method = "observability",
                           reconstruct = TRUE)
  sym <- symmetryDetection(gene, eqnvec(y = "p"), method = "observability",
                           symEngine = "symbolic")
  expect_length(mod$symmetries, 1L)
  expect_length(sym$symmetries, 1L)
  dmod <- mod$symmetries[[1]]; dsym <- sym$symmetries[[1]]
  expect_equal(dmod$type, "scaling")
  expect_equal(dsym$type, "scaling")
  # identical canonical weights, hence an identical rendered generator line
  expect_equal(dmod$weights[order(names(dmod$weights))],
               dsym$weights[order(names(dsym$weights))])
  expect_identical(dMod:::.sym_direction_line(dmod), dMod:::.sym_direction_line(dsym))
})


test_that("a pure constant shift is classified as a translation", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # a conserved moiety A + B with only the total observed: shifting A up and B
  # down leaves the total unchanged -- a pure (constant) translation
  moi <- eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B")
  r <- symmetryDetection(moi, eqnvec(y = "A + B"), method = "observability",
                         reconstruct = TRUE, reduceCQ = FALSE)
  ab <- Filter(function(d) all(c("A", "B") %in% names(d$generator)), r$symmetries)
  expect_gte(length(ab), 1L)
  expect_equal(ab[[1]]$type, "translation")
  expect_true(.sym_expr_equal(ab[[1]]$generator[["A"]], "1"))
  expect_true(.sym_expr_equal(ab[[1]]$generator[["B"]], "-1"))
})


test_that("method = 'translation' peels the exact additive lattice", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # the conserved moiety A + B and both rates are additively non-identifiable
  # (the output A + B is constant): the translation engine peels all three
  # constant directions directly, reconstruction-free
  moi <- eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B")
  tr <- symmetryDetection(moi, eqnvec(y = "A + B"), method = "translation",
                          reduceCQ = FALSE)
  expect_equal(tr$method, "translation")
  expect_length(tr$symmetries, 3L)
  expect_true(all(vapply(tr$symmetries,
                         function(d) isTRUE(d$type == "translation"), logical(1))))
  supp <- lapply(tr$symmetries, function(d) d$support)
  expect_true(any(vapply(supp, function(s) setequal(s, c("A", "B")), logical(1))))

  # a pure scaling model has NO constant direction: its tangent w*z varies with the
  # sample point and drops out of the intersection, so the lattice is empty
  gene <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
  tg <- symmetryDetection(gene, eqnvec(y = "p"), method = "translation")
  expect_length(tg$symmetries, 0L)
})


test_that("certifyPoly flags a direction that is a polynomial Lie symmetry", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # harmonic oscillator with only the conserved energy observed: the rotation of
  # the initial state is a genuine (affine) Lie point symmetry -> certified
  rot <- eqnvec(x1 = "-x2", x2 = "x1")
  r <- symmetryDetection(rot, eqnvec(y = "x1^2 + x2^2"), method = "observability",
                         reconstruct = TRUE, reduceCQ = FALSE,
                         control = reconstControl(certifyPoly = TRUE))
  aff <- Filter(function(d) d$type == "affine", r$symmetries)
  expect_gte(length(aff), 1L)
  expect_true(isTRUE(aff[[1]]$certified))

  # default (certifyPoly = FALSE) attaches no certificate
  r0 <- symmetryDetection(rot, eqnvec(y = "x1^2 + x2^2"), method = "observability",
                          reconstruct = TRUE, reduceCQ = FALSE)
  aff0 <- Filter(function(d) d$type == "affine", r0$symmetries)
  expect_true(length(aff0) >= 1L && !isTRUE(aff0[[1]]$certified))
})


test_that("the symbolic engine handles multiple conditions and single-time events", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # switch-gated rate: two switch values identify k1 and k2; the symbolic engine
  # stacks the two conditions and must reach the same verdict as the modular one
  f  <- eqnvec(A = "-(k1 + u*k2)*A", u = "0")
  ev <- addEvent(eventlist(), var = "u", time = -1, value = "var_u", method = "replace")
  cg <- data.frame(var_u = c(0, 1), row.names = c("ctrl", "stim"))
  mod <- symmetryDetection(f, eqnvec(y = "A"), method = "observability",
                           events = ev, conditions = cg, reduceCQ = FALSE)
  sym <- symmetryDetection(f, eqnvec(y = "A"), method = "observability",
                           events = ev, conditions = cg, reduceCQ = FALSE,
                           symEngine = "symbolic")
  expect_equal(sym$rank, mod$rank)
  expect_equal(sym$dim, mod$dim)
  expect_true(sym$identifiable)
  expect_equal(sym$info$conditions, 2L)
  expect_false(any(vapply(sym$symmetries,
                          function(d) "u" %in% d$support, logical(1))))

  # a single switch value cannot separate k1 and k2 -> non-identifiable
  s1 <- symmetryDetection(f, eqnvec(y = "A"), method = "observability",
                          events = ev, conditions = data.frame(var_u = 1, row.names = "stim"),
                          reduceCQ = FALSE, symEngine = "symbolic")
  expect_false(s1$identifiable)

  # a condition-specific rate rename spans the union parameter space (k, k_knd, A)
  cgk <- data.frame(k = c("k", "k_knd"), row.names = c("ctrl", "knd"),
                    stringsAsFactors = FALSE)
  uk <- symmetryDetection(eqnvec(A = "-k*A"), eqnvec(y = "A"),
                          method = "observability", conditions = cgk,
                          reduceCQ = FALSE, symEngine = "symbolic")
  expect_equal(uk$dim, 3L)
  expect_true(uk$identifiable)
})


test_that("the symbolic engine rejects equilibrate and later-event gaps", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  gene <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
  # equilibrate is a modular feature; symbolic wants an explicit steady state via trafo
  expect_error(symmetryDetection(gene, eqnvec(y = "p"), method = "observability",
                                 equilibrate = TRUE, symEngine = "symbolic"),
               "equilibrate")
  # later events (gaps) need the modular kernel
  ev2 <- eventlist() |>
    addEvent(var = "u", time = -1, value = "var_u", method = "replace") |>
    addEvent(var = "A", time = 0, value = "dose", method = "add")
  cg2 <- data.frame(var_u = c(0, 1), dose = c(1, 1), row.names = c("c1", "c2"))
  expect_error(symmetryDetection(eqnvec(A = "-(k1 + u*k2)*A", u = "0"), eqnvec(y = "A"),
                                 method = "observability", events = ev2, conditions = cg2,
                                 forcings = "u", reduceCQ = FALSE, symEngine = "symbolic"),
               "single-segment")
})
