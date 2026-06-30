
#' Search for structural non-identifiabilities of a model
#'
#' @description Detects structural non-identifiabilities of a reaction network
#'   plus its observation map, driven through the Python module
#'   `symmetryDetectionVersion2` via `reticulate`. The model `f`, observables `g`
#'   and optional transformation `trafo` are given as equations (anything
#'   [as.eqnvec] accepts). The engine is chosen with `method`:
#'
#'   * `"observability"` (default): the rank of the observability-identifiability
#'     matrix over a finite field. Full column rank means structurally locally
#'     identifiable; otherwise the nullspace gives the non-identifiable directions
#'     and the symbols they involve. Exact and scalable, defined only for rational
#'     right-hand sides and observables. The sensible default.
#'   * `"liesym"`: the polynomial Lie-symmetry ansatz of Merkt et al. (2015),
#'     returning the explicit closed-form generator and transformation. More
#'     expensive (grows with `pMax`); use it when the explicit symmetry is wanted.
#'   * `"scaling"`: the scaling symmetries that multiply each state and parameter
#'     by a power of one common factor, from the integer kernel of the
#'     monomial-exponent conditions. Returns only the scaling part.
#'
#' @param f The model right-hand sides: an [eqnlist], an [eqnvec], or a named
#'   character vector keyed by state name (anything [as.eqnvec] accepts). With an
#'   [eqnlist], conserved quantities can be reduced first (`reduceCQ`).
#' @param g The observation functions, as an [eqnvec] or named character vector.
#' @param trafo Optional parameter transformation as an [eqnvec] or named
#'   character vector. A parameter-named entry (inner symbol = expression in the
#'   outer parameters) is substituted into `f` and `g`; a state-named entry is
#'   that state's initial condition (a parameter, number or rational expression).
#'   A solved steady state from [steadyStates] can be passed whole.
#' @param method One of `"observability"` (default), `"liesym"` or `"scaling"`;
#'   see Description.
#' @param parameters Character vector of extra symbols to treat as parameters.
#' @param forcings Character vector of externally driven (input) state names.
#'   For `"observability"` a forcing is an integrated state with default initial
#'   value 0 that is excluded from the `f = 0` steady-state constraint (so the
#'   constraint solves over the non-forcing states with every forcing held at
#'   0); for `"liesym"` and `"scaling"` the same states are the inputs whose
#'   infinitesimal ansatz is fixed. Listed explicitly here; an `RHS = 0` state
#'   that is not a forcing is instead treated as a constant and substituted by
#'   its value.
#' @param events Optional [eventlist]. For `"observability"` the event value is
#'   looked up in `conditions` when it names a grid column, so the boolean
#'   switches and the stimulus levels of a condition grid reach the analysis.
#'   Events at or before `t0` set the pre-stimulus regime and compose the first
#'   segment's initial condition; events after `t0` split the post-stimulus
#'   timeline into analytic segments, one local Taylor jet each, all stacked, with
#'   the state propagated exactly across each gap (Details).
#' @param conditions Optional data frame of experimental conditions (a dMod
#'   condition grid: one row per condition, columns named by model symbols or by
#'   event-value placeholders). When supplied, `"observability"` runs the
#'   multi-condition analysis: each condition is compiled with its own symbol
#'   substitutions (a numeric cell bakes a symbol to a constant, a symbol cell
#'   renames it, e.g. a knockdown-specific rate) and its event-defined initial
#'   values, and the observability matrices of all conditions are stacked. A
#'   direction is reported non-identifiable only when it is unobservable in every
#'   condition (the intersection nullspace), and known inputs are constants of
#'   each condition rather than free symbols, so they never appear in a
#'   reconstructed direction. The directions are returned in closed form.
#' @param fixed Character vector of symbols that are known and therefore not
#'   unknowns. For `"observability"` they are excluded from the coordinates `z`:
#'   a fixed parameter is a known constant, and a fixed state keeps its dynamics
#'   but carries no unknown initial value. For `"liesym"` their infinitesimal
#'   ansatz is set to zero.
#' @param t0 Optional numeric start time of the first segment's Taylor expansion
#'   (default: the earliest event time, else 0). Events at `t0` are the first
#'   stimulus (composed into the start-point initial condition); the distinct
#'   event times after `t0` are the boundaries of the later segments. A dynamic
#'   state's event perturbs its value by `replace` / `add` / `multiply`;
#'   constant-state (switch) events are regime substitutions.
#' @param equilibrate Logical. For `"observability"`, start the states at a steady
#'   state of `f` with forcings held at 0, solved exactly over a finite field
#'   without a symbolic solution (parameters are never solved for); t0 events such
#'   as a dose are applied on top. A free power or Hill exponent is supported,
#'   whether or not its base has a linear turnover term: the term `base^exp` is
#'   recast to a rational coordinate `E = base^exp` with a companion `L =
#'   log(base)`, so the verdict (which parameters, including the exponent, are
#'   identifiable) is exact. With `closedForm`, `E` and `L` are
#'   sampled as coordinates and back-substituted, so a non-identifiable direction
#'   is returned in closed form even when it involves `base^exp` or `log(base)`.
#'   State initial conditions given in `trafo` are ignored when equilibrating.
#' @param reduceCQ Logical. If the model is an [eqnlist] with conserved
#'   quantities, eliminate one pivot species per conserved quantity first.
#' @param closedForm Logical. For `"observability"`, reconstruct the non-scaling
#'   directions as exact rational functions of the parameters, not just their
#'   support. Scaling directions are peeled exactly from the integer kernel and
#'   reported in closed form regardless of this flag. A narrow entry is fit
#'   densely; an entry coupling many variables is recovered by sparse
#'   (Ben-Or-Tiwari) interpolation. Every reconstructed direction is certified
#'   against the nullspace at a fresh prime; one that cannot be reconstructed or
#'   certified is returned with its support and `closedForm = FALSE` instead.
#' @param cores Number of threads. For `"observability"` the per-condition Taylor
#'   build is parallelised over conditions and the reconstruction's sample bank
#'   over evaluation points; for `"liesym"` and `"scaling"` it is the worker count
#'   of the symbolic search.
#' @param control A [reconstControl()] list tuning the `"observability"` engine's
#'   saturation and closed-form reconstruction (relevance caps, fit degrees, term
#'   and gap-order caps). Raise the caps to recover wide or high-degree directions.
#' @param liesym A [liesymControl()] list tuning the `"liesym"` engine: the
#'   infinitesimal ansatz and degree, the extra Lie-derivative order, the symbolic
#'   backend and verification, and the evaluation point.
#' @param scaling A [scalingControl()] list tuning the `"scaling"` engine: the
#'   symbolic backend and the evaluation point.
#'
#' @return For `"observability"`: a list with `identifiable`, `rank`, `dim`,
#'   `lieOrderUsed` and `nonIdentifiable`. Each direction carries a `type`:
#'   `"scaling"` directions give the integer toric weights `c_i` (the symmetry is
#'   `z_i -> lambda^{c_i} z_i`); `"general"` directions give the reconstructed
#'   tangent components as rational expressions (or only the `support` when
#'   `closedForm = FALSE`, or with a `reason` naming the `reconstControl` cap that
#'   stopped a requested closed form). The multi-condition path adds `conditions`, the number
#'   of distinct conditions, `segments`, the number of stacked analytic segments
#'   (events after `t0` split each condition into several), and `gapOrderUsed`, the
#'   gap power-series order at which the rank saturated (0 with no gaps). For
#'   `"scaling"`: `count` and the same `nonIdentifiable` weights. For
#'   `"liesym"`: a list of generators, each with `infinitesimals`,
#'   `transformation`, `type` and `verified`. The result has class
#'   `symmetryDetection`; `print()` gives a one-line summary and `summary()`
#'   lists the directions.
#'
#' @details
#' **Segments and exact propagation.** With events after `t0`, the post-stimulus
#' timeline is split at the distinct event times into analytic segments, one local
#' Taylor jet per segment over a shared coordinate space, all stacked (across
#' conditions) like conditions. A later segment's start state is the previous
#' segment's state flowed across the inter-event gap and mapped by the event. That
#' finite-time nonlinear flow is transcendental, so it is carried not by plugging
#' in the numeric gap length but by keeping the gap length a *formal power-series
#' variable* \eqn{\Delta t}: the kernel propagates the state as its exact Taylor
#' series in \eqn{\Delta t} (the coefficients it already computes), and an
#' identifiability direction must annihilate every \eqn{\Delta t} coefficient (a
#' polynomial vanishing for generic timing). The gap order is raised alongside the
#' Lie order until the rank stabilises (`gapOrderUsed`), giving the **exact
#' identifiability for generic inter-event timing**, transient channels included.
#' Truncating the gap order early is sound (a conservative lower bound on the rank).
#'
#' This makes parameters identifiable that only enter through a transient between
#' events (e.g. an inhibitor pre-incubation whose relaxed state seeds the stimulus
#' phase): apply it at equilibrium (`t0` at the inhibitor time) so the gap is
#' propagated. A `replace`/`add`/`multiply` event on a dynamic state is applied
#' exactly to the propagated state. A dose on a conserved-moiety *pivot* species
#' eliminated by `reduceCQ` is not seen (it is no longer a state); keep that
#' species (`reduceCQ = FALSE`) to dose it.
#'
#' @note The interface, defaults and output structure may still change between
#'   releases.
#'
#' @references \[1\]
#' <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.012920>
#'
#' @examples
#' \dontrun{
#' # The canonical scaling symmetry, found by all three engines: a reversible
#' # reaction observed only through alpha * A leaves the absolute scale free.
#' eq <- eqnlist() |>
#'   addReaction("A", "B", "k1 * A") |>
#'   addReaction("B", "A", "k2 * B")
#'
#' summary(symmetryDetection(eq, eqnvec(Aobs = "alpha * A")))          # observability
#' summary(symmetryDetection(eq, eqnvec(Aobs = "alpha * A"), method = "liesym"))
#' summary(symmetryDetection(eq, eqnvec(Aobs = "alpha * A"), method = "scaling"))
#'
#' # A steady-state initial condition as an expression: a state-named trafo entry
#' # is the initial condition, seeded with its parameter sensitivities.
#' summary(symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
#'                           trafo = eqnvec(x = "b/a"), closedForm = TRUE))
#'
#' # Several experimental conditions: a switch held at one generic value makes
#' # k1 and k2 look unidentifiable, but two values (set per condition by an
#' # event, read from the grid) identify both.
#' fu <- eqnvec(A = "-(k1 + u * k2) * A", u = "0")
#' events <- addEvent(eventlist(), var = "u", time = -1, value = "var_u",
#'                    method = "replace")
#' grid <- data.frame(var_u = c(0, 1), row.names = c("ctrl", "stim"))
#' symmetryDetection(fu, eqnvec(y = "A"), method = "observability",
#'                   events = events, conditions = grid)$identifiable
#'
#' # Pre-equilibrated model with a dose event: the steady state b/a is the
#' # relaxation attractor; a known dose makes (a, b, s) identifiable.
#' dose <- addEvent(eventlist(), var = "x", time = 0, value = "dose",
#'                  method = "replace")
#' summary(symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
#'                           method = "observability", events = dose,
#'                           conditions = data.frame(dose = 2, row.names = "stim")))
#'
#' # Full worked script (scaling, closed-form directions, an enzyme assay, a
#' # transcription-translation rate curve, conditions and steady states):
#' file.edit(system.file("examples", "symmetryDetection.R", package = "dMod"))
#' }
#' @export
symmetryDetection <- function(f = NULL, g = NULL, trafo = NULL,
                              method = c("observability", "liesym", "scaling"),
                              parameters = NULL, fixed = NULL, forcings = NULL,
                              events = NULL, conditions = NULL,
                              t0 = NULL, equilibrate = FALSE,
                              reduceCQ = TRUE,
                              closedForm = FALSE, cores = 1,
                              control = reconstControl(),
                              liesym = liesymControl(),
                              scaling = scalingControl()) {

  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required for symmetryDetection().")
  method <- match.arg(method)
  equilibrate <- isTRUE(equilibrate)

  # warn about arguments that do not apply to the chosen engine instead of
  # silently ignoring them
  supplied <- setdiff(names(match.call())[-1], "")
  applies <- switch(method,
    observability = c("events", "conditions", "t0", "equilibrate", "control"),
    liesym  = "liesym",
    scaling = "scaling")
  methodSpecific <- c("events", "conditions", "t0", "equilibrate",
                      "control", "liesym", "scaling")
  ignored <- intersect(supplied, setdiff(methodSpecific, applies))
  if (length(ignored))
    warning("symmetryDetection(): argument(s) ", paste(ignored, collapse = ", "),
            " do not apply to method = \"", method, "\" and are ignored.",
            call. = FALSE)

  # model right-hand sides as an eqnvec; f is eqnlist, eqnvec or named character
  if (is.null(f))
    stop("Provide the model right-hand sides via `f` ",
         "(eqnlist, eqnvec or named character vector).")
  feqnlist <- if (inherits(f, "eqnlist")) f else NULL
  fdyn <- as.eqnvec(f)
  states <- names(fdyn)

  gobs <- if (is.null(g)) NULL else as.eqnvec(g)
  parameters <- parameters %||% character(0)

  # state-named trafo entries are initial conditions, parameter-named entries are
  # substitutions applied to f, g and the initial-condition expressions
  initial <- NULL
  if (!is.null(trafo)) {
    trafo <- as.eqnvec(trafo)
    icEntries  <- intersect(names(trafo), states)
    subEntries <- setdiff(names(trafo), states)
    subs <- trafo[subEntries]
    sub <- function(e) {
      if (is.null(e) || !length(e) || !length(subs)) return(e)
      e <- as.eqnvec(e)
      setNames(replaceSymbols(names(subs), subs, e), names(e))
    }
    fdyn <- sub(fdyn)
    gobs <- sub(gobs)
    if (length(icEntries)) initial <- sub(trafo[icEntries])
  }

  # equilibrate solves the initial conditions from f = 0, so any given in `trafo`
  # are dropped with a warning
  if (equilibrate && !is.null(initial)) {
    warning("symmetryDetection(): equilibrate solves the steady state from f = 0; ",
            "the initial condition(s) for ",
            paste(names(initial), collapse = ", "), " in `trafo` are ignored.",
            call. = FALSE)
    initial <- NULL
  }

  # conserved-quantity reduction (eqnlist input only)
  if (!is.null(feqnlist) && isTRUE(reduceCQ)) {
    totals <- getTotals(feqnlist)
    if (length(totals)) {
      cq <- .detect_and_substitute_cq(totals, TRUE, fdyn, names(fdyn),
                                      parameters, expressInTotals = TRUE)
      fdyn <- cq$f[setdiff(names(cq$f), cq$elim_states)]
      parameters <- cq$parameters
      states <- names(fdyn)
      if (length(cq$cq_info)) {
        keys <- vapply(cq$cq_info, function(ci) ci$elim_state, character(1))
        vals <- vapply(cq$cq_info, function(ci) unname(ci$recon_expr), character(1))
        recon <- function(e) {
          if (is.null(e) || !length(e)) return(e)
          e <- as.eqnvec(e)
          setNames(replaceSymbols(keys, vals, e), names(e))
        }
        gobs    <- recon(gobs)
        initial <- recon(initial)
      }
    }
  }

  toLines <- function(e) {
    if (is.null(e)) return(NULL)
    e <- as.eqnvec(e)
    if (!length(e)) return(NULL)
    as.character(vapply(seq_along(e),
                        function(i) paste(names(e)[i], "=", e[i]),
                        character(1)))
  }

  code_dir <- system.file("code", package = "dMod")
  sysmod <- reticulate::import("sys", convert = TRUE)
  if (!(code_dir %in% sysmod$path)) sysmod$path <- c(code_dir, sysmod$path)
  sd <- reticulate::import("symmetryDetectionVersion2", convert = TRUE)

  # event times must be numeric (an event at a parameter time has no place in the
  # local Taylor jet at t0)
  if (!is.null(events) && nrow(as.data.frame(events))) {
    et <- as.character(as.data.frame(events)$time)
    if (any(is.na(suppressWarnings(as.numeric(et)))))
      stop("event `time` must be numeric; a parameter time is not supported: ",
           paste(et[is.na(suppressWarnings(as.numeric(et)))], collapse = ", "),
           call. = FALSE)
  }

  if (method == "observability") {
    # grid-substitution targets include symbols that appear only in initial
    # values or event values, so the grid can fix such parameters too
    extraSyms <- c(if (!is.null(initial)) getSymbols(as.character(as.eqnvec(initial))),
                   if (!is.null(events)) getSymbols(as.character(as.data.frame(events)$value)))
    symbols <- unique(c(states, names(gobs),
                        getSymbols(as.character(c(fdyn, gobs))), extraSyms))
    # states with a right-hand side that evaluates to 0 are constant in time and
    # substituted by their per-condition value (boolean switches, held inputs)
    isZeroRHS <- vapply(as.character(fdyn), function(r)
      isTRUE(suppressWarnings(tryCatch(eval(parse(text = r)) == 0,
                                       error = function(e) FALSE))), logical(1))
    constStates <- names(fdyn)[isZeroRHS]
    # states forced to zero at the resting state are held at zero in the f = 0
    # solve but stay dynamic states in the observability tape
    equilZeroStates <- if (equilibrate && !is.null(feqnlist))
      intersect(.equil_zero_states(feqnlist, forcings), states) else character(0)
    res <- .sym_resolve_conditions(conditions, events, initial, symbols, states,
                                   constStates, forcings, t0, equilibrate)
    multi <- sd$compileObservabilityTapeMulti(
      model = toLines(fdyn), observation = toLines(gobs),
      conditionSubs = res$subs, conditionIC0 = res$ic0,
      fixed = if (length(fixed)) fixed else NULL,
      parameters = if (length(parameters)) parameters else NULL,
      equilibrate = equilibrate, segEquilibrate = as.list(res$segEquil),
      forcings = if (length(forcings)) forcings else NULL,
      conditionEvents = res$segEvents, conditionT0Events = res$events0)
    if (!isTRUE(multi$ok))
      stop("method = \"observability\" requires rational right-hand sides, ",
           "observables and initial conditions (built from +, -, *, / and ",
           "integer powers).\n  ",
           paste(unlist(multi$nonrational), collapse = "\n  "),
           "\nA logarithmic observable log10(h) + offset equals the rational ",
           "observable scale * h; supply it in that form, or use ",
           "method = \"liesym\".", call. = FALSE)
    spy <- tryCatch(reticulate::import("sympy", convert = TRUE),
                    error = function(err) NULL)
    return(.observability_analytic_multi(multi, spy = spy,
                                         closedForm = closedForm, sd = sd,
                                         cores = cores,
                                         equilZeroStates = equilZeroStates,
                                         t0events = res$events0,
                                         nConditions = res$nConditions,
                                         chainOf = res$chainOf,
                                         nGaps = res$nGaps,
                                         control = control))
  }

  # liesym / scaling: symbolic engines on the (trafo-substituted) f and g;
  # forcings are the externally driven (input) states. The scaling control carries
  # only backend and point, so the liesym-only fields fall back to their defaults.
  ctl <- if (method == "liesym") liesym else scaling
  fld <- function(nm, default) if (is.null(ctl[[nm]])) default else ctl[[nm]]
  # the Python engine writes its own progress and report to stdout; capture it so
  # the R print/summary methods are the single display path.
  reticulate::py_capture_output(
    res <- sd$symmetryDetectiondMod(
      model       = toLines(fdyn),
      observation = toLines(gobs),
      ansatz      = fld("ansatz", "uni"),
      pMax        = as.integer(fld("pMax", 2L)),
      inputs      = if (length(forcings)) forcings else NULL,
      fixed       = fixed,
      parallel    = as.integer(cores),
      allTrafos   = fld("allTrafos", FALSE),
      lieOrder    = as.integer(fld("lieOrder", 0L)),
      exact       = fld("exact", TRUE),
      verify      = fld("verify", TRUE),
      backend     = ctl$backend,
      parameters  = if (length(parameters)) parameters else NULL,
      method      = method,
      point       = ctl$point,
      symbolic    = closedForm
    ))
  structure(res, class = "symmetryDetection")
}


# `expr` identically zero once the `dead` symbols are set to zero, tested at
# random positive values for the surviving symbols
.expr_is_zero <- function(expr, dead, ntry = 4, tol = 1e-9) {
  free <- setdiff(getSymbols(expr), dead)
  for (i in seq_len(ntry)) {
    env <- as.list(setNames(rep(0, length(dead)), dead))
    if (length(free))
      env <- c(env, as.list(setNames(runif(length(free), 0.5, 1.5), free)))
    val <- tryCatch(eval(parse(text = expr), env), error = function(e) NA_real_)
    if (is.na(val) || abs(val) > tol) return(FALSE)
  }
  TRUE
}

# A sink cluster among the live reactions: a set of species admitting a nonneg
# combination whose net production over every live reaction is <= 0 and strictly
# negative overall, so mass can only leave and the cluster vanishes at steady
# state even when each member has a live producer. Structural (stoichiometry
# only), solved as an LP per candidate species (mirrors steadyStates).
.equil_sink_cluster <- function(M, eps = 1e-8, Mbig = 1e4) {
  nF <- nrow(M); nS <- ncol(M)
  if (nF == 0L || nS == 0L) return(integer(0))
  c_obj <- colSums(M); id <- diag(nS)
  for (i in seq_len(nS)) {
    lb <- rep(0, nS); ub <- rep(Mbig, nS); lb[i] <- 1; ub[i] <- 1
    res <- tryCatch(
      lpSolve::lp("min", c_obj, rbind(M, id, id),
                  c(rep("<=", nF), rep(">=", nS), rep("<=", nS)),
                  c(rep(0, nF), lb, ub)),
      error = function(e) NULL)
    if (!is.null(res) && res$status == 0 && res$objval < -eps)
      return(which(res$solution > eps))
  }
  integer(0)
}

# states forced to zero at a steady state once the forcings are held at zero: a
# species is dead if every producing reaction has a rate that vanishes under the
# already-dead symbols, iterated to a fixpoint over the reaction graph. When this
# propagation stalls, an LP sink-cluster step catches collectively-draining
# clusters that no single-species rule sees.
.equil_zero_states <- function(eqnlist, forcings) {
  S <- eqnlist$smatrix
  if (is.null(S) || !length(S)) return(character(0))
  S[is.na(S)] <- 0
  species <- colnames(S)
  rates <- as.character(eqnlist$rates)
  seed <- intersect(forcings, c(species, getSymbols(paste(rates, collapse = "+"))))
  dead <- seed
  repeat {
    live <- !vapply(rates, .expr_is_zero, logical(1), dead = dead)
    newdead <- character(0)
    for (sp in setdiff(species, dead)) {
      prod <- which(S[, sp] > 0)
      if (length(prod) && !any(live[prod])) newdead <- c(newdead, sp)
    }
    if (!length(setdiff(newdead, dead))) {
      alive <- setdiff(species, dead)
      cl <- .equil_sink_cluster(S[live, alive, drop = FALSE])
      if (length(cl)) newdead <- alive[cl] else break
    }
    if (!length(setdiff(newdead, dead))) break
    dead <- union(dead, newdead)
  }
  setdiff(dead, forcings)
}


# Analytic (closed-form) observability path for symmetryDetection().
#
# Builds the observability-identifiability matrix exactly over GF(p) through the
# C++ kernel (src/symmetry_kernel.cpp), certifies the rank and the
# non-identifiable directions, peels off the scaling directions (already exact
# from the integer kernel), and reconstructs each remaining direction as an
# exact rational function of the parameters by multi-point interpolation:
# sampling the modular nullspace at many generic points, fitting a rational
# function per nullspace entry over each prime, and lifting the coefficients to
# the rationals by Chinese remaindering and rational reconstruction. The whole
# path is exact; no floating-point arithmetic enters.

# four primes < 2^31 (their product < 2^124, within unsigned __int128)
.symPrimes <- c(2147483647, 2147483629, 2147483587, 2147483579)

# a fifth prime, disjoint from the reconstruction primes, used only to certify a
# reconstructed direction against the nullspace at a fresh evaluation
.symVerifyPrime <- 2147483563

#' Tuning controls for `symmetryDetection(method = "observability")`
#'
#' Bundles the saturation and closed-form-reconstruction parameters of the
#' observability engine into one object, passed as `symmetryDetection(...,
#' control = reconstControl(...))`. Raise the caps for models whose
#' non-identifiability directions are wide or high-degree rationals (at the cost of
#' more finite-field samples).
#'
#' @param relevanceCap Max variables a single nullspace entry may couple before it
#'   is fit by sparse interpolation instead of the dense fit.
#' @param relevanceCapDir Max variables a whole direction may couple before it is
#'   returned as support-only.
#' @param relevanceCapSparse Max variables a single entry may couple for the
#'   sparse (Ben-Or-Tiwari) path; beyond this the entry is support-only.
#' @param degreeCap Total-degree bound of the dense rational fit per entry.
#' @param sampleSlack Extra evaluation points beyond the minimum the fit needs.
#' @param probeRetries Fresh-randomness retries for the per-leaf relevance probe
#'   when a perturbation shifts the pivot set.
#' @param laurentDegNum,laurentDegDen Numerator and single-monomial denominator
#'   degree bounds of the sparse Laurent path.
#' @param laurentCandCap Cap on the candidate-monomial product enumerated by the
#'   Laurent / general path (guards memory).
#' @param termCap Max term count (Ben-Or-Tiwari order) of a sparse entry.
#' @param generalDegNum,generalDegDen Numerator and (multi-term) denominator
#'   degree bounds of the general sparse-rational (Cauchy + Ben-Or-Tiwari) path.
#' @param gapOrderCap Cap on the gap power-series order raised when propagating the
#'   state exactly across post-`t0` event boundaries.
#' @return A `reconstControl` list.
#' @seealso [symmetryDetection()]
#' @export
reconstControl <- function(relevanceCap       = 6L,
                           relevanceCapDir    = 24L,
                           relevanceCapSparse = 30L,
                           degreeCap          = 4L,
                           sampleSlack        = 5L,
                           probeRetries       = 8L,
                           laurentDegNum      = 4L,
                           laurentDegDen      = 2L,
                           laurentCandCap     = 200000L,
                           termCap            = 60L,
                           generalDegNum      = 4L,
                           generalDegDen      = 3L,
                           gapOrderCap        = 8L) {
  stopifnot(relevanceCap >= 0L, relevanceCapDir >= 1L,
            relevanceCapSparse >= relevanceCap, degreeCap >= 0L,
            sampleSlack >= 0L, probeRetries >= 1L, termCap >= 1L,
            laurentDegNum >= 1L, laurentDegDen >= 0L, laurentCandCap >= 1L,
            generalDegNum >= 1L, generalDegDen >= 1L, gapOrderCap >= 0L)
  structure(list(relevanceCap       = as.integer(relevanceCap),
                 relevanceCapDir    = as.integer(relevanceCapDir),
                 relevanceCapSparse = as.integer(relevanceCapSparse),
                 degreeCap          = as.integer(degreeCap),
                 sampleSlack        = as.integer(sampleSlack),
                 probeRetries       = as.integer(probeRetries),
                 laurentDegNum      = as.integer(laurentDegNum),
                 laurentDegDen      = as.integer(laurentDegDen),
                 laurentCandCap     = as.integer(laurentCandCap),
                 termCap            = as.integer(termCap),
                 generalDegNum      = as.integer(generalDegNum),
                 generalDegDen      = as.integer(generalDegDen),
                 gapOrderCap        = as.integer(gapOrderCap)),
            class = c("reconstControl", "list"))
}


#' Tuning control for `symmetryDetection(method = "liesym")`
#'
#' Bundles the polynomial Lie-symmetry engine's parameters into one object, passed
#' as `symmetryDetection(..., liesym = liesymControl(...))`.
#'
#' @param ansatz Type of infinitesimal ansatz: `"uni"`, `"par"` or `"multi"`.
#' @param pMax Maximal degree of the infinitesimal ansatz (integer `>= 1`).
#' @param lieOrder Integer `N >= 0`. Also require the generator to annihilate the
#'   `k`-th Lie derivative `L^k g` for `k = 1..N`.
#' @param exact Logical. Use exact modular linear algebra (`TRUE`) or the
#'   floating-point path (`FALSE`).
#' @param verify Logical. Symbolically verify each generator.
#' @param allTrafos Logical. Keep transformations that share a common parameter
#'   factor instead of dropping them.
#' @param backend `"symengine"` (faster, falls back to sympy) or `"sympy"`.
#' @param point Optional named numeric vector: the evaluation point for the
#'   symbolic search (default: a deterministic generic point).
#' @return A `liesymControl` list.
#' @seealso [symmetryDetection()], [scalingControl()]
#' @export
liesymControl <- function(ansatz    = c("uni", "par", "multi"),
                          pMax      = 2L,
                          lieOrder  = 0L,
                          exact     = TRUE,
                          verify    = TRUE,
                          allTrafos = FALSE,
                          backend   = c("symengine", "sympy"),
                          point     = NULL) {
  ansatz  <- match.arg(ansatz)
  backend <- match.arg(backend)
  stopifnot(pMax >= 1L, lieOrder >= 0L,
            is.logical(exact), is.logical(verify), is.logical(allTrafos),
            is.null(point) || is.numeric(point))
  structure(list(ansatz = ansatz, pMax = as.integer(pMax),
                 lieOrder = as.integer(lieOrder), exact = isTRUE(exact),
                 verify = isTRUE(verify), allTrafos = isTRUE(allTrafos),
                 backend = backend, point = point),
            class = c("liesymControl", "list"))
}


#' Tuning control for `symmetryDetection(method = "scaling")`
#'
#' Bundles the scaling (toric) symmetry engine's parameters into one object, passed
#' as `symmetryDetection(..., scaling = scalingControl(...))`.
#'
#' @param backend `"symengine"` (faster, falls back to sympy) or `"sympy"`.
#' @param point Optional named numeric vector: the evaluation point for the
#'   symbolic search (default: a deterministic generic point).
#' @return A `scalingControl` list.
#' @seealso [symmetryDetection()], [liesymControl()]
#' @export
scalingControl <- function(backend = c("symengine", "sympy"), point = NULL) {
  backend <- match.arg(backend)
  stopifnot(is.null(point) || is.numeric(point))
  structure(list(backend = backend, point = point),
            class = c("scalingControl", "list"))
}


.sym_sieve <- function(n) {
  limit <- 200L
  repeat {
    is_p <- rep(TRUE, limit)
    is_p[1] <- FALSE
    for (i in 2:floor(sqrt(limit))) if (is_p[i]) is_p[seq(i * i, limit, by = i)] <- FALSE
    pr <- which(is_p)
    if (length(pr) >= n) return(pr[seq_len(n)])
    limit <- limit * 2L
  }
}


# A stream of distinct primes used as generic evaluation coordinates; grows on
# demand so the interpolation never runs out of sample points.
.sym_pool <- function() {
  cache <- .sym_sieve(1000L)
  function(k) {
    if (max(k) > length(cache)) cache <<- .sym_sieve(max(2L * length(cache), max(k)))
    cache[k]
  }
}


.sym_null_residues <- function(res, freeCol, p) {
  nz <- res$dim
  piv <- res$pivots
  v <- integer(nz)
  v[freeCol + 1L] <- 1L
  for (ri in seq_along(piv))
    v[piv[ri] + 1L] <- as.integer((p - res$R[ri, freeCol + 1L]) %% p)
  v
}


# Modular arithmetic over GF(p) for the residual-direction gauge, p < 2^31. A
# double holds integers exactly only to 2^53, so a*b (up to 2^62) is split on a
# 15-bit boundary to keep every intermediate product below 2^47. Vectorised in a;
# b is a scalar (the elimination multiplier).
.sym_mulmod <- function(a, b, p) {
  a <- a %% p; b <- b %% p
  hi <- (a * (b %/% 32768)) %% p
  ((hi * 32768) %% p + a * (b %% 32768)) %% p
}
# inverse by Fermat, computed with the split multiply so the squarings stay exact
.sym_invmod <- function(a, p) {
  r <- 1; b <- a %% p; e <- p - 2
  while (e > 0) {
    if (e %% 2 == 1) r <- .sym_mulmod(r, b, p)
    e <- e %/% 2
    if (e > 0) b <- .sym_mulmod(b, b, p)
  }
  r
}

# Reduced row echelon form over GF(p); returns the reduced rows, the 0-based pivot
# columns and the rank. Sizes here are tiny (the nullspace dimension by nz).
.sym_rref_modp <- function(M, p) {
  M <- matrix(as.numeric(M) %% p, nrow(M), ncol(M))
  nr <- nrow(M); nc <- ncol(M); piv <- integer(0); row <- 1L
  for (col in seq_len(nc)) {
    if (row > nr) break
    nz <- which(M[row:nr, col] %% p != 0)
    if (!length(nz)) next
    pr <- row + nz[1] - 1L
    if (pr != row) { tmp <- M[row, ]; M[row, ] <- M[pr, ]; M[pr, ] <- tmp }
    M[row, ] <- .sym_mulmod(M[row, ], .sym_invmod(M[row, col], p), p)
    for (r2 in seq_len(nr)) if (r2 != row) {
      f2 <- M[r2, col] %% p
      if (f2 != 0) M[r2, ] <- (M[r2, ] - .sym_mulmod(M[row, ], f2, p)) %% p
    }
    piv <- c(piv, col - 1L); row <- row + 1L
  }
  list(R = M, piv = piv, rank = length(piv))
}

# k-by-k matrix inverse over GF(p) by Gauss-Jordan; NULL if singular.
.sym_matinv_modp <- function(B, p) {
  k <- nrow(B)
  A <- cbind(matrix(as.numeric(B) %% p, k, k), diag(k))
  for (col in seq_len(k)) {
    nz <- which(A[col:k, col] %% p != 0)
    if (!length(nz)) return(NULL)
    pr <- col + nz[1] - 1L
    if (pr != col) { tmp <- A[col, ]; A[col, ] <- A[pr, ]; A[pr, ] <- tmp }
    A[col, ] <- .sym_mulmod(A[col, ], .sym_invmod(A[col, col], p), p)
    for (r in seq_len(k)) if (r != col) {
      f <- A[r, col] %% p
      if (f != 0) A[r, ] <- (A[r, ] - .sym_mulmod(A[col, ], f, p)) %% p
    }
  }
  A[, (k + 1L):(2L * k), drop = FALSE]
}

# Decouple the residual directions the free-column gauge entangles. Each residual
# residue is first reduced modulo the (exact integer) scaling span, then a set of
# distinguished columns (the pivots of the reduced residuals at the base point) is
# fixed. Direction i is the unique span representative that is 1 on its own
# distinguished column and 0 on the others' (a k-by-k solve), which for genuinely
# independent directions is exactly the physical sparse direction. The anchor entry
# stays an exact 1, so the rational reconstruction is clean. Returns one (anchor,
# residue function) pair per direction, or the raw free-column gauge if the
# directions cannot be separated. The residue function returns the decoupled
# nullspace vector at a kernel result (NULL when degenerate at that sample).
.sym_canon_gauge <- function(residualFree, scalRows, P, nz, sc) {
  k <- length(residualFree)
  raw <- list(anchors = residualFree, residueFns = vector("list", k))
  if (k <= 1L) return(raw)

  reduceScal <- function(R, p) {
    if (nrow(scalRows) == 0L) return(R %% p)
    Sr <- .sym_rref_modp(scalRows, p)
    for (j in seq_along(Sr$piv)) {
      fac <- R[, Sr$piv[j] + 1L] %% p
      for (i in seq_len(nrow(R)))
        if (fac[i] != 0) R[i, ] <- (R[i, ] - .sym_mulmod(Sr$R[j, ], fac[i], p)) %% p
    }
    R
  }
  resid <- function(rp, p) reduceScal(
    t(vapply(residualFree, function(fc) .sym_null_residues(rp, fc, p),
             integer(nz))), p)

  R0 <- resid(sc$ref, P)
  rr <- .sym_rref_modp(R0, P)
  if (rr$rank < k) return(raw)
  dp <- rr$piv

  makeFn <- function(i) function(rp, p, zvals = NULL) {
    R <- resid(rp, p)
    Binv <- .sym_matinv_modp(R[, dp + 1L, drop = FALSE], p)
    if (is.null(Binv)) return(NULL)
    w <- numeric(nz)
    for (l in seq_len(k)) w <- (w + .sym_mulmod(R[l, ], Binv[i, l], p)) %% p
    as.integer(w)
  }
  list(anchors = dp, residueFns = lapply(seq_len(k), makeFn))
}


# Decouple the residual directions in LOGARITHMIC coordinates. A weighted scaling
# xi_i = c_i * z_i (e.g. a Hill-exponent feedback, c_i = -nhill) has a constant or
# low-degree log-residue eta_i = xi_i / z_i, so normalising "1 on a pivot column"
# introduces no rational denominator and the reconstruction stays sparse. The
# free-column / canonical gauges instead divide by a rational function of the
# leaves and blow the direction up across the whole feedback loop. Mirrors
# .sym_canon_gauge but on eta-rows; the residue functions divide each direction's
# free-column residue componentwise by the sample point's z-values (mod p). The
# reconstructed eta-direction is turned back into xi by .sym_logcoord_backsub.
.sym_logcoord_gauge <- function(residualFree, scalRows, P, nz, sc, zvals0) {
  k <- length(residualFree)
  raw <- list(anchors = residualFree, residueFns = vector("list", k))
  if (k == 0L) return(raw)

  # integer scaling weights in log coordinates: a scaling tangent is weight * z,
  # so its log-residue is the (point-independent) integer weight, recovered from
  # the base point by a centred lift (toric weights are small). Reducing the
  # eta-rows modulo these strips any scaling admixture the free-column gauge folded
  # in, which is what couples a feedback direction to its whole upstream chain.
  W <- if (nrow(scalRows) == 0L) matrix(0L, 0L, nz) else {
    Wm <- matrix(0L, nrow(scalRows), nz)
    for (c in seq_len(nz)) {
      zc <- as.numeric(zvals0[c]) %% P
      if (zc == 0) next
      wc <- .sym_mulmod(scalRows[, c], .sym_invmod(zc, P), P)
      Wm[, c] <- as.integer(ifelse(wc > P / 2, wc - P, wc))
    }
    Wm
  }
  reduceScal <- function(R, p) {
    if (nrow(W) == 0L) return(R %% p)
    Sr <- .sym_rref_modp(W, p)
    for (j in seq_along(Sr$piv)) {
      fac <- R[, Sr$piv[j] + 1L] %% p
      for (i in seq_len(nrow(R)))
        if (fac[i] != 0) R[i, ] <- (R[i, ] - .sym_mulmod(Sr$R[j, ], fac[i], p)) %% p
    }
    R
  }

  etaRows <- function(rp, p, zvals) {
    R <- t(vapply(residualFree, function(fc) .sym_null_residues(rp, fc, p),
                  integer(nz)))
    for (c in seq_len(nz)) {
      zc <- as.numeric(zvals[c]) %% p
      if (zc == 0) { if (any(R[, c] %% p != 0)) return(NULL); next }
      R[, c] <- .sym_mulmod(R[, c], .sym_invmod(zc, p), p)
    }
    reduceScal(R, p)
  }

  E0 <- etaRows(sc$ref, P, zvals0)
  if (is.null(E0)) return(raw)
  rr <- .sym_rref_modp(E0, P)
  if (rr$rank < k) return(raw)
  dp <- rr$piv

  makeFn <- function(i) function(rp, p, zvals = NULL) {
    E <- etaRows(rp, p, zvals)
    if (is.null(E)) return(NULL)
    Binv <- .sym_matinv_modp(E[, dp + 1L, drop = FALSE], p)
    if (is.null(Binv)) return(NULL)
    w <- numeric(nz)
    for (l in seq_len(k)) w <- (w + .sym_mulmod(E[l, ], Binv[i, l], p)) %% p
    as.integer(w)
  }
  list(anchors = dp, residueFns = lapply(seq_len(k), makeFn))
}


# Turn a direction reconstructed in log coordinates back into original
# coordinates: each entry is the log-residue eta_c of column znames[c], so the
# tangent is xi_c = eta_c * znames[c] (the anchor "1" becomes the column symbol
# itself). Every znames column is a multiplicative coordinate, so the symbol is
# just multiplied in and the product cancelled.
.sym_logcoord_backsub <- function(vector, spy) {
  out <- list()
  for (nm in names(vector)) {
    expr <- paste0("(", vector[[nm]], ")*", nm)
    out[[nm]] <- if (is.null(spy)) expr else .sym_simplify(expr, spy)
  }
  out
}


# exponent vectors of `nvar` variables with total degree <= `degree`, ordered by
# total degree. Enumerates the C(nvar+degree, nvar) monomials directly (never the
# full (degree+1)^nvar grid), so it stays feasible for many variables.
.sym_mono_table <- function(nvar, degree) {
  if (nvar == 0L) return(matrix(0L, 1L, 0L))
  gen <- function(n, d) {
    if (n == 1L) return(lapply(0:d, function(e) e))
    unlist(lapply(0:d, function(e) lapply(gen(n - 1L, d - e), function(s) c(e, s))),
           recursive = FALSE)
  }
  m <- do.call(rbind, gen(nvar, degree))
  m <- m[order(rowSums(m)), , drop = FALSE]
  storage.mode(m) <- "integer"
  dimnames(m) <- NULL
  m
}


.sym_mono_string <- function(expo, vars) {
  terms <- character(0)
  for (k in seq_along(vars)) if (expo[k] != 0)
    terms <- c(terms, if (expo[k] == 1) vars[k] else paste0(vars[k], "^", expo[k]))
  if (!length(terms)) "1" else paste(terms, collapse = "*")
}


.sym_poly_string <- function(numc, denc, mons, vars) {
  parts <- character(0)
  for (j in seq_len(nrow(mons))) {
    n <- numc[j]; d <- denc[j]
    if (n == "0") next
    coef <- if (d == "1") n else paste0("(", n, "/", d, ")")
    ms <- .sym_mono_string(mons[j, ], vars)
    parts <- c(parts, if (ms == "1") coef else paste0(coef, "*", ms))
  }
  if (!length(parts)) "0" else paste(parts, collapse = " + ")
}


# Reconstruct one nullspace entry as a rational function of the relevant
# variables. The per-prime fit returns the kernel vector (num and den
# coefficients) in a gauge fixed by its free column; reconstruction needs the
# same free column at every prime. Returns coefficient strings, or NULL if no
# closed form of bounded degree fits.
.sym_reconstruct_entry <- function(sampleU, mons, residues, primes) {
  nMon <- nrow(mons)
  coefRes <- matrix(0L, 2L * nMon, length(primes))
  freeCol <- NULL
  for (j in seq_along(primes)) {
    fit <- symFitRational(sampleU, mons, as.integer(residues[, j]), primes[j])
    if (!identical(fit$status, "ok")) return(NULL)
    if (is.null(freeCol)) freeCol <- fit$freeCol
    else if (fit$freeCol != freeCol) return(NULL)
    coefRes[, j] <- fit$coeffs
  }
  rec <- symRatRecon(coefRes, as.integer(primes))
  if (any(rec$den == "0")) return(NULL)
  num <- seq_len(nMon); den <- nMon + seq_len(nMon)
  list(numCoefN = rec$num[num], numCoefD = rec$den[num],
       denCoefN = rec$num[den], denCoefD = rec$den[den])
}


.sym_simplify <- function(expr, spy) {
  out <- tryCatch(as.character(spy$cancel(spy$sympify(expr))),
                  error = function(e) expr)
  if (length(out) != 1L || is.na(out)) expr else out
}


# Replace the recast coordinates in each entry of a reconstructed direction by
# their meaning: E -> base^exp, L -> log(base). The substitution and cancellation
# run in Python (exact symbolic arithmetic).
.sym_recast_backsub <- function(vector, recast, sd) {
  eN <- as.list(vapply(recast, function(r) as.character(r$E), ""))
  lN <- as.list(vapply(recast, function(r) as.character(r$L), ""))
  bN <- as.list(vapply(recast, function(r) as.character(r$base), ""))
  xN <- as.list(vapply(recast, function(r) as.character(r$exp), ""))
  lapply(vector, function(x) {
    out <- tryCatch(sd$recastBacksub(as.character(x), eN, lN, bN, xN),
                    error = function(e) as.character(x))
    if (length(out) != 1L || is.na(out)) as.character(x) else out
  })
}


# Exact value of a rational expression at integer coordinates `env`, reduced
# modulo q (computed in Python to keep the big-integer arithmetic exact). Returns
# NA when the denominator vanishes mod q (verification then inconclusive).
.sym_eval_modq <- function(expr, env, q, sd) {
  v <- tryCatch(sd$evalRationalMod(as.character(expr), as.list(names(env)),
                                   as.list(as.numeric(unlist(env))),
                                   as.integer(q)),
                error = function(err) NULL)
  if (is.null(v)) NA_integer_ else as.integer(v)
}


# Certify a reconstructed direction (free column `f`) against the nullspace at a
# fresh prime: the closed-form entries, evaluated at the base point modulo
# `.symVerifyPrime`, must equal the kernel's null vector there. Returns FALSE only
# on a definite mismatch; an unavailable evaluation leaves the verdict inconclusive
# (TRUE) so verification never rejects a correct direction it cannot re-check.
.sym_verify_direction <- function(entry, f, znames, leafNames, point0, NtUsed,
                                  kcall, sd, residueFn = NULL) {
  if (is.null(sd) || is.null(entry$vector)) return(TRUE)
  q <- .symVerifyPrime
  rq <- kcall(point0, q, as.integer(NtUsed))
  if (!isTRUE(rq$ok)) return(TRUE)
  aq <- if (is.null(residueFn)) .sym_null_residues(rq, f, q) else residueFn(rq, q, NULL)
  if (is.null(aq)) return(TRUE)
  env <- setNames(as.list(as.numeric(point0[seq_along(leafNames)])), leafNames)
  for (nm in names(entry$vector)) {
    col <- match(nm, znames)
    if (is.na(col)) next
    got <- .sym_eval_modq(entry$vector[[nm]], env, q, sd)
    if (is.na(got)) return(TRUE)
    if (as.integer(got) %% q != as.integer(aq[col]) %% q) return(FALSE)
  }
  TRUE
}


# Gauge-independent certificate that a reconstructed direction genuinely lies in
# the nullspace: at a fresh evaluation point, build the closed-form vector and
# check it is reproduced by the kernel's free-column residues (v in nullspace iff
# v = sum over free columns of v[free] times that column's null vector). This
# rejects a canonical representative that self-verifies against its own residue
# function but has baked in base-point values (it then fails away from the base
# point). Inconclusive evaluations leave the verdict TRUE so a correct direction is
# never rejected.
.sym_verify_in_nullspace <- function(entry, f, znames, leafNames, point0, NtUsed,
                                     kcall, pool, poolNext, nz, sd) {
  if (is.null(sd) || is.null(entry$vector)) return(TRUE)
  q <- .symVerifyPrime
  pt <- as.numeric(point0)
  pt[] <- pool(poolNext + seq_along(pt) - 1L)
  rq <- kcall(pt, q, as.integer(NtUsed))
  if (!isTRUE(rq$ok)) return(TRUE)
  env <- setNames(as.list(pt[seq_along(leafNames)]), leafNames)
  v <- integer(nz); v[f + 1L] <- 1L
  for (nm in names(entry$vector)) {
    col <- match(nm, znames)
    if (is.na(col)) next
    got <- .sym_eval_modq(entry$vector[[nm]], env, q, sd)
    if (is.na(got)) return(TRUE)
    v[col] <- as.integer(got %% q)
  }
  freeCols <- setdiff(seq_len(nz) - 1L, as.integer(rq$pivots))
  recon <- numeric(nz)
  for (fcol in freeCols) {
    vf <- v[fcol + 1L] %% q
    if (vf != 0) recon <- (recon + .sym_mulmod(.sym_null_residues(rq, fcol, q), vf, q)) %% q
  }
  all((recon - v) %% q == 0)
}


# Find a generic base point where the rank is maximal over the primes and
# saturate the Lie order (and, with event gaps, the gap power-series order Mtot)
# until the rank stops growing. kcall(point, p, Nt, Mtot) must return the kernel
# list (ok, R, pivots, rank, dim). `maxM` caps the gap order (0 for the no-gap
# path). Returns NULL if no usable point is found, else the reference reduction,
# the certified rank, and the Lie / gap orders used.
.sym_saturate_certify <- function(kcall, nLeaves, nz, maxM = 0L) {
  P <- .symPrimes[1]
  pool <- .sym_pool()
  point0 <- pool(seq_len(nLeaves))
  poolNext <- nLeaves + 1L

  saturateNt <- function(point, Mtot) {
    prev <- -1L; Nt <- 1L; res <- NULL
    repeat {
      r <- kcall(point, P, Nt, Mtot)
      if (!isTRUE(r$ok)) return(NULL)
      res <- r
      if (r$rank >= nz || (r$rank == prev && Nt >= 2L) || Nt > nz + 1L) break
      prev <- r$rank; Nt <- Nt + 1L
    }
    list(res = res, Nt = Nt)
  }

  # several generic points may be tried before one admits a steady-state point
  # over GF(p) for every condition (saturate the Lie order at gap order 0 first)
  sat <- NULL
  for (attempt in 1:50) {
    sat <- saturateNt(point0, 0L)
    if (!is.null(sat)) break
    point0 <- pool(poolNext + seq_len(nLeaves) - 1L); poolNext <- poolNext + nLeaves
  }
  if (is.null(sat)) return(NULL)
  NtUsed <- sat$Nt; MtotUsed <- 0L; saturatedM <- TRUE

  # raise the gap order until the rank stops growing (exact generic-timing rank);
  # if the cap is hit while still growing, the truncated rank is conservative
  if (maxM > 0L) repeat {
    if (sat$res$rank >= nz) break
    if (MtotUsed >= maxM) { saturatedM <- FALSE; break }
    satM <- saturateNt(point0, MtotUsed + 1L)
    if (is.null(satM)) break
    MtotUsed <- MtotUsed + 1L
    if (satM$res$rank <= sat$res$rank) break
    sat <- satM; NtUsed <- sat$Nt
  }

  rankMax <- sat$res$rank
  for (pj in .symPrimes[-1]) {
    rj <- kcall(point0, pj, NtUsed, MtotUsed)
    if (isTRUE(rj$ok)) rankMax <- max(rankMax, rj$rank)
  }
  while (sat$res$rank < rankMax) {
    point0 <- pool(poolNext + seq_len(nLeaves) - 1L); poolNext <- poolNext + nLeaves
    sat <- saturateNt(point0, MtotUsed)
    if (is.null(sat)) return(NULL)
    NtUsed <- sat$Nt
  }
  list(ref = sat$res, NtUsed = NtUsed, MtotUsed = MtotUsed, saturatedM = saturatedM,
       point0 = point0, pool = pool, poolNext = poolNext,
       rank = sat$res$rank, pivots = sat$res$pivots)
}


# Nullspace basis of the reference reduction: one column per free coordinate.
.sym_nullspace_basis <- function(ref, freeCols, P) {
  if (!length(freeCols)) return(matrix(0L, ref$dim, 0L))
  vapply(freeCols, function(fc) .sym_null_residues(ref, fc, P), integer(ref$dim))
}


# Peel the scaling directions common to all conditions. Each generator is
# projected onto the unknown coordinates z (state weights outside z, e.g. under a
# steady-state constraint, drop out), turned into its tangent w_c * z_c at the
# base point, verified to lie in the nullspace, and kept only if independent of
# the scalings already taken. Returns the scaling entries and their tangent span.
.sym_peel_scalings <- function(scalRes, znames, nz, zval, P, N) {
  scaling <- list()
  Bmat <- matrix(0L, nz, 0L)
  inSpan <- function(M, x) ncol(M) > 0L &&
    !is.null(symSolveMod(matrix(as.integer(M), nz), as.integer(x), P))
  for (d in scalRes$nonIdentifiable) {
    v <- integer(nz)
    cols <- match(names(d$vector), znames)
    keep <- !is.na(cols)
    if (!any(keep)) next
    v[cols[keep]] <- as.integer(unlist(d$vector))[keep]
    if (all(v == 0)) next
    tangent <- as.integer((as.numeric(v) * zval) %% P)
    if (all(tangent == 0) || !inSpan(N, tangent) || inSpan(Bmat, tangent)) next
    Bmat <- cbind(Bmat, tangent)
    scaling[[length(scaling) + 1L]] <- list(
      support = sort(znames[v != 0]),
      vector = setNames(as.list(as.character(v[v != 0])), znames[v != 0]),
      type = "scaling", closedForm = TRUE)
  }
  list(scaling = scaling, Bmat = Bmat)
}


# Multi-condition analytic observability. `multi` is the list returned by
# compileObservabilityTapeMulti: per-condition tapes over a shared coordinate
# space. The observability rows of all conditions are stacked and reduced once,
# so the verdict and every reconstructed direction reflect the intersection
# nullspace. Known inputs are baked into each condition and are not leaves, so a
# direction can only involve genuine parameters and free initial values.
.observability_analytic_multi <- function(multi, spy = NULL,
                                          closedForm = TRUE, sd = NULL,
                                          cores = 1, equilZeroStates = character(0),
                                          t0events = list(), nConditions = NULL,
                                          chainOf = NULL, nGaps = 0L,
                                          control = reconstControl()) {
  ctrl <- control
  cores <- as.integer(max(1L, cores))
  nLeaves <- as.integer(multi$nLeaves)
  nStates <- as.integer(multi$nStates)
  zSlots <- as.integer(multi$zSlots)
  znames <- as.character(multi$znames)
  nz <- length(znames)
  leafNames <- as.character(multi$leafNames)
  ssConstraint <- isTRUE(multi$equilibrate)
  nSegmentTapes <- length(multi$tapes)
  # with post-t0 event gaps the state is propagated exactly across each gap by the
  # kernel (generic gap length as a formal power series); segments group into one
  # chain per condition. maxM caps the gap power-series order raised at saturation.
  hasGaps <- isTRUE(nGaps > 0L) && !is.null(chainOf)
  chainGroups <- if (hasGaps) split(seq_along(multi$tapes), chainOf) else NULL
  firstOfChain <- if (hasGaps) !duplicated(chainOf) else rep(TRUE, length(multi$tapes))
  maxM <- if (hasGaps) ctrl$gapOrderCap else 0L
  # the recast coordinates E = base^exp, L = log base are appended to the leaf
  # space as extra sampled coordinates; reconstruction fits the rational
  # dependence on them and back-substitutes, recovering closed forms that involve
  # a free exponent. With no recast, the augmented space equals the leaf space.
  recast <- list()
  nAug <- nLeaves
  leafNamesAug <- leafNames

  tapeFields <- function(t) {
    out <- list(
      op = as.integer(t$op), a = as.integer(t$a), b = as.integer(t$b),
      cnum = as.character(t$cnum), cden = as.character(t$cden),
      stateSlots = as.integer(t$stateSlots), fOut = as.integer(t$fOut),
      gOut = as.integer(t$gOut), icLeaf = as.integer(t$icLeaf),
      icNum = as.character(t$icNum), icDen = as.character(t$icDen))
    # a segment seeded from an IC tape (free/carry initial values, doses, resets)
    # carries its IC instructions; a steady-state-seeded segment carries none and
    # gets an icSeed instead. This is decided per segment, not globally, so the
    # two seed kinds can coexist (equilibrate S0 + carried later segments).
    if (!is.null(t$icOp))
      out <- c(out, list(
        icOp = as.integer(t$icOp), icA = as.integer(t$icA),
        icB = as.integer(t$icB), icCnum = as.character(t$icCnum),
        icCden = as.character(t$icCden), icOut = as.integer(t$icOut)))
    # a later segment's left-boundary state-dose event map (replace/add/multiply by
    # a parametric value), applied by the kernel to the propagated state
    if (!is.null(t$evVarIdx))
      out <- c(out, list(
        evOp = as.integer(t$evOp), evA = as.integer(t$evA),
        evB = as.integer(t$evB), evCnum = as.character(t$evCnum),
        evCden = as.character(t$evCden), evOut = as.integer(t$evOut),
        evVarIdx = as.integer(t$evVarIdx), evMethod = as.integer(t$evMethod)))
    out
  }
  tapes <- lapply(multi$tapes, tapeFields)

  if (ssConstraint) {
    # constraint mode: states are not free coordinates. Per (point, prime) solve
    # f = 0 on the interior component and seed each tape with the modular steady
    # state plus its IFT parameter-duals (icSeed), then stack as usual.
    stateNames <- as.character(multi$stateNames)
    paramNames <- as.character(multi$paramNames)
    forcings <- if (is.null(multi$forcings)) character(0) else as.character(multi$forcings)
    # forcings and forced-zero states are held at zero in the f = 0 solve
    solveHeld <- union(forcings, as.character(equilZeroStates))
    models <- lapply(multi$tapes, function(t) as.character(t$constraintModel))
    # only the equilibrate-seeded segments are solved for a steady state; with
    # gaps only the first segment of each chain anchors (the rest are propagated)
    isEquilTape <- vapply(multi$tapes,
                          function(t) length(t$constraintModel) > 0L, logical(1)) &
                   firstOfChain
    w <- nz + 1L
    # power/Hill recast atoms. Each base contributes a generic coordinate held at a
    # fixed distinct prime (reduced per prime) plus L = log base. A normal base is
    # solved with E = base^exp generic; an inverted base is generic with E solved.
    # realStateNames lists the names actually solved for the steady state.
    recast <- multi$powerRecast
    if (is.null(recast)) recast <- list()
    realStateNames <- if (length(recast)) as.character(multi$realStateNames) else stateNames
    genName <- function(r) if (isTRUE(r$inverted)) as.character(r$base) else as.character(r$E)
    genNames <- vapply(recast, genName, "")
    lNames <- unique(vapply(recast, function(r) as.character(r$L), ""))
    solveParamNames <- c(paramNames, genNames)
    if (length(recast)) {
      leafNamesAug <- c(leafNames, genNames, lNames)
      nAug <- length(leafNamesAug)
    }
    # per-condition seed plan: the steady-state value/Jacobian/event tapes, so the
    # C++ kernel produces the icSeed in integer arithmetic. A non-generically-linear
    # condition has no plan and is seeded by the symbolic solve instead.
    useCppSeed <- isTRUE(getOption("dMod.symSeedCpp", TRUE))
    seedPlans <- lapply(seq_along(tapes), function(ci) {
      if (!isEquilTape[ci] || !useCppSeed) return(NULL)
      evC <- if (ci <= length(t0events) && length(t0events[[ci]]))
        t0events[[ci]] else NULL
      plan <- tryCatch(sd$compileSeedPlan(
        model = models[[ci]], stateNames = realStateNames,
        paramNames = solveParamNames,
        forcings = if (length(solveHeld)) solveHeld else NULL,
        recast = if (length(recast)) recast else NULL, t0events = evC),
        error = function(e) NULL)
      if (is.list(plan) && isTRUE(plan$genericLinear)) plan else NULL
    })
    ssCache <- new.env(parent = emptyenv())
    ssWhy <- NULL              # last steady-state failure reason, for diagnostics
    seedTapes <- function(point, p) {
      pmod <- as.numeric(point) %% p
      paramVals <- as.list(setNames(pmod[seq_along(leafNames)], leafNames))
      lValsList <- NULL
      if (length(recast)) {
        for (i in seq_along(genNames)) paramVals[[genNames[i]]] <- pmod[nLeaves + i]
        lValsList <- as.list(setNames(
          pmod[nLeaves + length(genNames) + seq_along(lNames)], lNames))
      }
      # seed each condition's first segment from its steady state: the compiled
      # C++ kernel for a generically-linear resting state, the symbolic solve when
      # a coupled residual remains
      solveOne <- function(ci) {
        if (!isEquilTape[ci]) return(NULL)
        plan <- seedPlans[[ci]]
        if (!is.null(plan)) {
          pvec <- as.integer(vapply(plan$paramNames, function(nm) {
            v <- paramVals[[nm]]
            if (is.null(v)) 0L else as.integer(as.numeric(v) %% p)
          }, integer(1)))
          lvec <- if (length(plan$recast))
            as.integer(vapply(plan$recast, function(rc) {
              v <- lValsList[[rc$L]]
              if (is.null(v)) 0L else as.integer(as.numeric(v) %% p)
            }, integer(1))) else integer(0)
          return(symSteadyStateSeed(plan, pvec, lvec, p))
        }
        evC <- if (ci <= length(t0events) && length(t0events[[ci]]))
          t0events[[ci]] else NULL
        sd$solveSteadyStateModular(
          model = models[[ci]], stateNames = realStateNames,
          paramNames = solveParamNames, paramVals = paramVals, prime = p,
          forcings = if (length(solveHeld)) solveHeld else NULL,
          t0events = evC,
          recast = if (length(recast)) recast else NULL, lVals = lValsList)
      }
      sols <- lapply(seq_along(tapes), solveOne)
      out <- vector("list", length(tapes))
      for (ci in seq_along(tapes)) {
        # a carried/IC-tape segment is already fully specified; pass it through
        if (!isEquilTape[ci]) { out[[ci]] <- tapes[[ci]]; next }
        sol <- sols[[ci]]
        if (inherits(sol, "try-error") || !is.list(sol) || !isTRUE(sol$ok)) {
          ssWhy <<- if (!is.list(sol)) "parallel solve failed"
                    else if (is.null(sol$why)) "unknown" else as.character(sol$why)
          return(NULL)
        }
        icSeed <- matrix(0L, nStates, w)
        # solved rows are keyed by name: an inverted base solves its E, whose row
        # in stateNames is not the base's own position.
        xs <- as.numeric(sol$xstar)
        for (k in seq_along(realStateNames)) {
          row <- match(realStateNames[k], stateNames)
          icSeed[row, 1L] <- as.integer(xs[k])
          for (c in seq_len(nz))
            icSeed[row, c + 1L] <- as.integer(as.numeric(sol$dx[[znames[c]]])[k])
        }
        # each recast entry seeds its generic partner (gen) and L with their duals
        for (rc in sol$recast) {
          gr <- match(rc$gen, stateNames); lr <- match(rc$L, stateNames)
          icSeed[gr, 1L] <- as.integer(rc$gen0)
          icSeed[lr, 1L] <- as.integer(rc$L0)
          for (c in seq_len(nz)) {
            gv <- rc$genDual[[znames[c]]]; lv <- rc$lDual[[znames[c]]]
            if (!is.null(gv)) icSeed[gr, c + 1L] <- as.integer(gv)
            if (!is.null(lv)) icSeed[lr, c + 1L] <- as.integer(lv)
          }
        }
        tp <- tapes[[ci]]
        tp$icSeed <- icSeed
        out[[ci]] <- tp
      }
      out
    }
    # cache both successful seeds and failed (point, prime) solves, so a point
    # that has no steady state over GF(p) is not solved again on a later visit
    kcall4 <- function(point, p, Nt, Mtot = 0L) {
      key <- paste(c(format(as.numeric(point), scientific = FALSE), p),
                   collapse = ",")
      seeded <- ssCache[[key]]
      if (is.null(seeded)) {
        seeded <- seedTapes(point, p)
        ssCache[[key]] <- if (is.null(seeded)) list(ok = FALSE) else seeded
      }
      if (is.list(seeded) && identical(seeded$ok, FALSE)) return(list(ok = FALSE))
      pt <- as.integer(point[seq_len(nLeaves)])
      if (hasGaps)
        symObsNullChain(lapply(chainGroups, function(idx) seeded[idx]), nLeaves,
                        nStates, zSlots, pt, p, as.integer(Nt), as.integer(Mtot), cores)
      else
        symObsNullMulti(seeded, nLeaves, nStates, zSlots, pt, p, as.integer(Nt), cores)
    }
  } else {
    kcall4 <- function(point, p, Nt, Mtot = 0L) {
      pt <- as.integer(point)
      if (hasGaps)
        symObsNullChain(lapply(chainGroups, function(idx) tapes[idx]), nLeaves,
                        nStates, zSlots, pt, p, as.integer(Nt), as.integer(Mtot), cores)
      else
        symObsNullMulti(tapes, nLeaves, nStates, zSlots, pt, p, as.integer(Nt), cores)
    }
  }

  sc <- .sym_saturate_certify(kcall4, nAug, nz, maxM)
  if (is.null(sc)) {
    if (ssConstraint && !is.null(ssWhy))
      warning("symmetryDetection(): no steady-state point over the finite field ",
              "after the generic-point retries (", ssWhy, "). The equilibrate ",
              "constraint could not be evaluated.", call. = FALSE)
    return(NULL)
  }
  # the reconstruction samples at the saturated gap order, baked into a 3-arg kcall
  MtotUsed <- if (is.null(sc$MtotUsed)) 0L else as.integer(sc$MtotUsed)
  kcall <- function(point, p, Nt) kcall4(point, p, Nt, MtotUsed)
  # batched solve over many (point, prime) pairs for the reconstruction's sample
  # bank. The single-segment, non-constraint path shares one tape set across all
  # points, so the whole batch runs in one OpenMP-over-points kernel call; the
  # constraint and gap paths need per-point seeding, so they fall back to looping
  # the per-call kcall (which still threads cores over conditions/segments).
  canBatch <- !ssConstraint && !hasGaps
  kbatch <- function(pointList, primeVec, Nt) {
    if (canBatch && length(pointList) > 0L) {
      M <- do.call(rbind, lapply(pointList, function(pp)
        as.integer(pp[seq_len(nLeaves)])))
      symObsNullBatch(tapes, nLeaves, nStates, zSlots, M, as.numeric(primeVec),
                      as.integer(Nt), cores)
    } else {
      Map(function(pp, pr) kcall(pp, pr, Nt), pointList, primeVec)
    }
  }
  if (hasGaps && isFALSE(sc$saturatedM))
    warning("symmetryDetection(): the gap power-series order hit the cap (",
            maxM, ") before the rank stabilised; the reported rank is a sound ",
            "(conservative) lower bound. Raise reconstControl(gapOrderCap=) for ",
            "the exact generic-timing verdict.", call. = FALSE)

  result <- list(method = "observability", engine = "analytic",
                 conditions = if (is.null(nConditions)) nSegmentTapes
                              else as.integer(nConditions),
                 segments = nSegmentTapes, gapOrderUsed = MtotUsed,
                 identifiable = (sc$rank == nz), rank = as.integer(sc$rank),
                 dim = as.integer(nz), lieOrderUsed = as.integer(sc$NtUsed),
                 nonIdentifiable = list())
  class(result) <- "symmetryDetection"
  if (sc$rank == nz)
    return(result)

  P <- .symPrimes[1]
  freeCols <- setdiff(0:(nz - 1L), sc$pivots)
  N <- .sym_nullspace_basis(sc$ref, freeCols, P)

  # scalings common to every condition are exact (integer kernel) and always
  # reported in closed form; their span is excluded before any reconstruction
  scaling <- list(); Bmat <- matrix(0L, nz, 0L)
  if (!is.null(sd)) {
    # each segment's regime dynamics/observation feed the scaling-candidate search;
    # candidates are validated against the nullspace N (which reflects propagation)
    modelLines <- lapply(multi$tapes, function(t) as.character(t$modelLines))
    obsLines   <- lapply(multi$tapes, function(t) as.character(t$obsLines))
    inputs <- if (!is.null(multi$forcings)) as.character(multi$forcings) else NULL
    fixed  <- setdiff(leafNames, znames)
    scalRes <- tryCatch(sd$scalingSymmetriesMulti(
      perCondModel = modelLines, perCondObs = obsLines,
      inputs = if (length(inputs)) inputs else NULL,
      fixed = if (length(fixed)) fixed else NULL), error = function(e) NULL)
    if (!is.null(scalRes)) {
      peel <- .sym_peel_scalings(scalRes, znames, nz, sc$point0[zSlots + 1L], P, N)
      scaling <- peel$scaling; Bmat <- peel$Bmat
    }
  }

  # free directions not spanned by the scalings are the residual non-scaling ones
  scalCols <- ncol(Bmat)
  residualFree <- integer(0)
  for (fc in freeCols) {
    bf <- .sym_null_residues(sc$ref, fc, P)
    inSpan <- ncol(Bmat) > 0L &&
      !is.null(symSolveMod(matrix(as.integer(Bmat), nz), as.integer(bf), P))
    if (!inSpan) { residualFree <- c(residualFree, fc); Bmat <- cbind(Bmat, bf) }
  }
  scalRows <- if (scalCols > 0L) t(Bmat[, seq_len(scalCols), drop = FALSE])
              else matrix(0L, 0L, nz)

  if (isTRUE(closedForm)) {
    interp <- list()
    # shared relevance probe: a single solve at a one-leaf perturbation yields the
    # residues of every direction, so the per-leaf scan is done once here rather
    # than once per direction (the dominant cost at scale).
    # a perturbation that shifts the pivot set yields an uncomparable probe, which
    # would mark the leaf relevant to every entry and inflate the fit past its caps;
    # retry with fresh random values until the pivots match before accepting the
    # pessimistic last candidate, so genuinely irrelevant leaves stay out of the fit
    probeNext <- sc$poolNext
    relProbe <- lapply(seq_len(nAug), function(li) {
      rp <- NULL; pert <- sc$point0
      for (att in seq_len(ctrl$probeRetries)) {
        pert <- sc$point0
        pert[li] <- sc$pool(probeNext); probeNext <<- probeNext + 1L
        cand <- kcall(pert, P, sc$NtUsed)
        rp <- cand
        if (!is.null(cand) && isTRUE(cand$ok) &&
            identical(as.integer(cand$pivots), as.integer(sc$pivots))) break
      }
      list(rp = rp, zvals = if (length(zSlots)) pert[zSlots + 1L] else NULL)
    })
    poolNext <- probeNext

    # reconstruct one direction in a given gauge (free-column when residueFn is
    # NULL, decoupled-canonical otherwise): interpolate, certify against the
    # nullspace at a fresh prime, and on failure downgrade to support-only. Also
    # certifies that the reconstructed vector actually lies in the nullspace at a
    # fresh point (gauge-independent), which rejects a self-consistent but
    # base-point-contaminated canonical representative.
    zvals0 <- if (length(zSlots)) sc$point0[zSlots + 1L] else NULL
    reconstructOne <- function(fc, rfn, logCoords = FALSE, sharedBank = NULL) {
      dir <- .sym_interpolate_direction(fc, sc$ref, sc$pivots, znames, nz, zSlots,
                                        leafNamesAug, nAug, sc$point0, sc$pool,
                                        poolNext, sc$NtUsed, kcall, spy, relProbe,
                                        rfn, ctrl, kbatch, sharedBank)
      poolNext <<- dir$poolNext
      e <- dir$entry
      # log-gauge entries are eta = xi / z; turn them back into the tangent xi
      # before verifying, then certify in original coordinates (gauge-independent)
      if (logCoords && isTRUE(e$closedForm))
        e$vector <- .sym_logcoord_backsub(e$vector, spy)
      ok <- isTRUE(e$closedForm) && (
        if (logCoords)
          .sym_verify_in_nullspace(e, fc, znames, leafNamesAug, sc$point0,
                                   sc$NtUsed, kcall, sc$pool, poolNext, nz, sd)
        else
          .sym_verify_direction(e, fc, znames, leafNamesAug, sc$point0,
                                sc$NtUsed, kcall, sd, rfn) &&
          (is.null(rfn) ||
           .sym_verify_in_nullspace(e, fc, znames, leafNamesAug, sc$point0,
                                    sc$NtUsed, kcall, sc$pool, poolNext, nz, sd)))
      if (isTRUE(e$closedForm) && !ok) {
        v <- if (!is.null(rfn)) rfn(sc$ref, P, zvals0) else .sym_null_residues(sc$ref, fc, P)
        if (is.null(v)) v <- .sym_null_residues(sc$ref, fc, P)
        e <- list(support = sort(znames[v != 0]), type = "general",
                  closedForm = FALSE,
                  reason = paste("a closed form was reconstructed but failed",
                                 "verification at a fresh prime"))
      }
      if (isTRUE(e$closedForm) && length(recast) && !is.null(sd))
        e$vector <- .sym_recast_backsub(e$vector, recast, sd)
      e
    }

    # raw free-column gauge first; this is exact for every model whose directions
    # the gauge does not entangle. One shared dense sample bank over the union of
    # all directions' relevant leaves serves every direction: the expensive kernel
    # is evaluated once per point and reused across directions.
    metas <- lapply(residualFree, function(fc)
      .sym_direction_relevance(fc, sc$ref, sc$pivots, zSlots, nAug, sc$point0,
                               relProbe, NULL))
    needs <- vapply(metas, function(m) .sym_dense_need(m$relByEntry, ctrl)$maxNeed,
                    integer(1))
    denseDir <- which(needs > 0L)
    bank <- NULL
    if (length(denseDir)) {
      unionLeaves <- sort(unique(unlist(lapply(metas[denseDir], function(m) m$relevant))))
      bk <- .sym_build_shared_bank(unionLeaves, max(needs[denseDir]), sc$point0,
                                   sc$pool, poolNext, kbatch, sc$pivots, sc$NtUsed,
                                   length(.symPrimes))
      poolNext <- bk$poolNext
      if (isTRUE(bk$ok)) bank <- bk
    }
    interp <- lapply(residualFree, function(fc) reconstructOne(fc, NULL, sharedBank = bank))

    # rescue: if a direction failed to close, the free-column gauge may be
    # entangling a few genuinely independent directions (e.g. two Hill feedbacks).
    # Retry the whole residual set in the decoupled-canonical gauge and adopt it
    # only if every direction then closes; otherwise keep the raw result.
    if (any(!vapply(interp, function(e) isTRUE(e$closedForm), logical(1)))) {
      gauge <- .sym_canon_gauge(residualFree, scalRows, P, nz, sc)
      if (length(gauge$anchors) == length(residualFree) &&
          !all(vapply(gauge$residueFns, is.null, logical(1)))) {
        canon <- lapply(seq_along(gauge$anchors), function(gi)
          reconstructOne(gauge$anchors[gi], gauge$residueFns[[gi]]))
        if (all(vapply(canon, function(e) isTRUE(e$closedForm), logical(1))))
          interp <- canon
      }
    }

    # second rescue: weighted scalings whose weight is a parameter (e.g. a Hill
    # exponent) slip through the integer scaling peel and the free-column gauge
    # blows them up. They are sparse in log coordinates, so decouple there and
    # back-substitute; adopt only if every residual direction then closes.
    if (any(!vapply(interp, function(e) isTRUE(e$closedForm), logical(1)))) {
      lg <- .sym_logcoord_gauge(residualFree, scalRows, P, nz, sc, zvals0)
      if (length(lg$anchors) == length(residualFree) &&
          !all(vapply(lg$residueFns, is.null, logical(1)))) {
        logrec <- lapply(seq_along(lg$anchors), function(gi)
          reconstructOne(lg$anchors[gi], lg$residueFns[[gi]], logCoords = TRUE))
        if (all(vapply(logrec, function(e) isTRUE(e$closedForm), logical(1))))
          interp <- logrec
      }
    }
    result$nonIdentifiable <- c(scaling, interp)
  } else {
    support <- lapply(residualFree, function(fc) {
      v <- .sym_null_residues(sc$ref, fc, P)
      list(support = sort(znames[v != 0]), type = "general", closedForm = FALSE)
    })
    result$nonIdentifiable <- c(scaling, support)
  }
  result
}


# Laurent candidate exponents: every numerator monomial (degree <= dNum) minus
# every single denominator monomial (degree <= dDen). The recovered entry is the
# Laurent polynomial supported on a subset of these, i.e. a rational whose
# denominator is a single monomial.
.sym_laurent_candidates <- function(nvar, dNum, dDen) {
  num <- .sym_mono_table(nvar, dNum)
  den <- .sym_mono_table(nvar, dDen)
  cand <- do.call(rbind, lapply(seq_len(nrow(den)),
                                function(d) sweep(num, 2L, den[d, ], "-")))
  cand <- unique(cand)
  storage.mode(cand) <- "integer"
  cand
}


# Lift one Ben-Or-Tiwari polynomial recovered per prime to the rationals: require
# the supports (exponent sets) to agree across primes, then Chinese-remainder the
# coefficients. Returns list(exps, num, den) of coefficient strings, or NULL.
.sym_bot_reconcile <- function(perPrime) {
  np <- length(perPrime)
  key <- function(m) apply(m, 1L, paste, collapse = ",")
  k1 <- key(perPrime[[1]]$exps)
  t <- length(k1)
  coefMat <- matrix(0L, t, np)
  coefMat[, 1L] <- perPrime[[1]]$coeffs
  for (pj in 2:np) {
    kp <- key(perPrime[[pj]]$exps)
    if (length(kp) != t || !setequal(kp, k1)) return(NULL)
    coefMat[, pj] <- perPrime[[pj]]$coeffs[match(k1, kp)]
  }
  rec <- symRatRecon(coefMat, as.integer(.symPrimes))
  if (any(rec$den == "0")) return(NULL)
  list(exps = perPrime[[1]]$exps, num = rec$num, den = rec$den)
}


# Assemble the Laurent terms into one entry string, clearing the common
# single-monomial denominator by shifting every exponent up by its largest
# negative part.
.sym_laurent_assemble <- function(perPrime, reli, leafNames, nvar) {
  rc <- .sym_bot_reconcile(perPrime)
  if (is.null(rc)) return(NULL)
  expRows <- rc$exps
  t <- nrow(expRows)
  mu <- apply(expRows, 2L, function(col) max(0L, -min(col)))
  vars <- leafNames[reli]
  numStr <- .sym_poly_string(rc$num, rc$den,
                             expRows + matrix(mu, t, nvar, byrow = TRUE), vars)
  denStr <- .sym_mono_string(mu, vars)
  if (denStr == "1") numStr else paste0("(", numStr, ")/(", denStr, ")")
}


# Reconstruct one wide entry by sparse Laurent interpolation: sample it on a
# geometric schedule of distinct prime bases, grow the sample length until the
# Ben-Or-Tiwari term count stabilises, then identify the term monomials with the
# smallest numerator/denominator degrees that fit (so a low-degree entry needs
# only a small candidate set). Coefficients are lifted across primes. Returns the
# entry string, or NULL when it is not a bounded-degree Laurent polynomial.
.sym_sparse_entry <- function(reli, supportCol, f, point0, leafNames, NtUsed,
                              kcall, pivots, residueFn = NULL,
                              ctrl = reconstControl(), zSlots = NULL) {
  nvar <- length(reli)
  bases <- .sym_sieve(nvar)
  np <- length(.symPrimes)
  maxLen <- 2L * ctrl$termCap + 2L

  # grow geometric samples until the term count (BM order) stabilises
  seqs <- lapply(seq_len(np), function(.) integer(0))
  cur <- lapply(seq_len(np), function(.) rep(1, nvar))
  have <- 0L
  repeat {
    target <- min(if (have == 0L) 16L else 2L * have, maxLen)
    for (pj in seq_len(np)) {
      p <- .symPrimes[pj]
      for (k in (have + 1L):target) {
        pt <- point0; pt[reli] <- cur[[pj]]
        rp <- kcall(pt, p, NtUsed)
        if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots)))
          return(NULL)
        zv <- if (is.null(zSlots)) NULL else pt[zSlots + 1L]
        nvp <- if (is.null(residueFn)) .sym_null_residues(rp, f, p) else residueFn(rp, p, zv)
        if (is.null(nvp)) return(NULL)
        seqs[[pj]][k] <- nvp[supportCol + 1L]
        cur[[pj]] <- (cur[[pj]] * bases) %% p
      }
    }
    have <- target
    if (2L * symBMorder(seqs[[1]], .symPrimes[1]) < have) break
    # the term count never stabilised within budget: the entry is not a
    # bounded-term Laurent polynomial (e.g. a genuine multi-term denominator)
    if (have >= maxLen) return(NULL)
  }

  # identify the terms at the smallest degrees that fit, smallest candidate set
  # first; the denominator is a single monomial (dDen, including the polynomial
  # case dDen = 0)
  for (dDen in 0:ctrl$laurentDegDen) for (dNum in seq_len(ctrl$laurentDegNum)) {
    if (choose(nvar + dNum, nvar) * choose(nvar + dDen, nvar) > ctrl$laurentCandCap)
      next
    candM <- matrix(as.integer(.sym_laurent_candidates(nvar, dNum, dDen)), ncol = nvar)
    perPrime <- vector("list", np)
    okAll <- TRUE
    for (pj in seq_len(np)) {
      monoRes <- symMonoResidues(candM, as.integer(bases), .symPrimes[pj])
      res <- symSparsePoly(seqs[[pj]], candM, as.integer(monoRes), .symPrimes[pj])
      if (!identical(res$status, "ok") || res$nterms == 0L) { okAll <- FALSE; break }
      perPrime[[pj]] <- res
    }
    if (!okAll) next
    out <- .sym_laurent_assemble(perPrime, reli, leafNames, nvar)
    if (!is.null(out)) return(out)
  }
  NULL
}


# Reconstruct one wide entry with a general (multi-term) denominator. From a
# generic shift s, the entry along a ray s + t*(u - s) is a univariate rational
# in t; a Cauchy fit normalised to B(0) = 1 evaluates N(u)/D(s) and D(u)/D(s) at
# any u (the gauge D(s) is a constant that cancels in N/D). Both are sparse
# polynomials, recovered by Ben-Or-Tiwari on the geometric schedule and lifted
# across primes. Returns the entry string, or NULL when no bounded-degree
# rational form fits. The common D(s) factor is cancelled by the caller's
# simplification.
.sym_general_rational_entry <- function(reli, supportCol, f, point0, leafNames,
                                        NtUsed, kcall, pivots, residueFn = NULL,
                                        ctrl = reconstControl(), zSlots = NULL) {
  nvar <- length(reli)
  np <- length(.symPrimes)
  bases <- .sym_sieve(nvar)
  s <- as.numeric(point0[reli])

  sampleR <- function(relvals, p) {
    pt <- point0; pt[reli] <- relvals %% p
    rp <- kcall(pt, p, NtUsed)
    if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots)))
      return(NULL)
    zv <- if (is.null(zSlots)) NULL else pt[zSlots + 1L]
    nvp <- if (is.null(residueFn)) .sym_null_residues(rp, f, p) else residueFn(rp, p, zv)
    if (is.null(nvp)) return(NULL)
    nvp[supportCol + 1L]
  }
  # value of (N(u)/D(s), D(u)/D(s)) at u over prime p, via the ray Cauchy fit
  evalND <- function(uvec, p, dN, dD) {
    v <- (uvec - s) %% p
    tn <- 0:(dN + dD + 4L)
    rv <- lapply(tn, function(tj) sampleR((s + tj * v) %% p, p))
    if (any(vapply(rv, is.null, logical(1)))) return(NULL)
    res <- symCauchyEval(as.integer(tn), as.integer(unlist(rv)), dN, dD, p)
    if (!identical(res$status, "ok")) return(NULL)
    c(res$N, res$D)
  }

  # numerator/denominator degrees, probed once at a generic point
  P1 <- .symPrimes[1]
  u0 <- vapply(seq_len(nvar), function(i) (bases[i] * bases[i] + 3) %% P1, numeric(1))
  deg <- NULL
  for (tot in 2:(ctrl$generalDegNum + ctrl$generalDegDen)) {
    for (a in seq.int(max(1L, tot - ctrl$generalDegDen), min(tot - 1L, ctrl$generalDegNum))) {
      bD <- tot - a
      if (bD < 1L || bD > ctrl$generalDegDen) next
      if (!is.null(evalND(u0, P1, a, bD))) { deg <- c(a, bD); break }
    }
    if (!is.null(deg)) break
  }
  if (is.null(deg)) return(NULL)
  dN <- deg[1]; dD <- deg[2]
  if (choose(nvar + dN, nvar) + choose(nvar + dD, nvar) > ctrl$laurentCandCap)
    return(NULL)
  candN <- matrix(as.integer(.sym_mono_table(nvar, dN)), ncol = nvar)
  candD <- matrix(as.integer(.sym_mono_table(nvar, dD)), ncol = nvar)

  # geometric Ben-Or-Tiwari sampling of both N/D(s) and D/D(s), grown until both
  # term counts stabilise
  maxLen <- 2L * ctrl$termCap + 2L
  Nseq <- lapply(seq_len(np), function(.) integer(0)); Dseq <- Nseq
  cur <- lapply(seq_len(np), function(.) rep(1, nvar))
  have <- 0L
  repeat {
    target <- min(if (have == 0L) 16L else 2L * have, maxLen)
    for (pj in seq_len(np)) {
      p <- .symPrimes[pj]
      for (k in (have + 1L):target) {
        e <- evalND(cur[[pj]], p, dN, dD)
        if (is.null(e)) return(NULL)
        Nseq[[pj]][k] <- e[1]; Dseq[[pj]][k] <- e[2]
        cur[[pj]] <- (cur[[pj]] * bases) %% p
      }
    }
    have <- target
    if (2L * symBMorder(Nseq[[1]], P1) < have &&
        2L * symBMorder(Dseq[[1]], P1) < have) break
    if (have >= maxLen) return(NULL)
  }

  fitPoly <- function(seqList, candM) {
    perPrime <- vector("list", np)
    for (pj in seq_len(np)) {
      mr <- symMonoResidues(candM, as.integer(bases), .symPrimes[pj])
      res <- symSparsePoly(seqList[[pj]], candM, as.integer(mr), .symPrimes[pj])
      if (!identical(res$status, "ok") || res$nterms == 0L) return(NULL)
      perPrime[[pj]] <- res
    }
    .sym_bot_reconcile(perPrime)
  }
  Nrc <- fitPoly(Nseq, candN); if (is.null(Nrc)) return(NULL)
  Drc <- fitPoly(Dseq, candD); if (is.null(Drc)) return(NULL)
  vars <- leafNames[reli]
  numS <- .sym_poly_string(Nrc$num, Nrc$den, Nrc$exps, vars)
  denS <- .sym_poly_string(Drc$num, Drc$den, Drc$exps, vars)
  if (denS == "0") return(NULL)
  if (denS == "1") numS else paste0("(", numS, ")/(", denS, ")")
}


# Per-direction relevance from the shared probe: the direction's nullspace vector
# at the base point (raw free-column residue, or canonicalised row when residueFn
# is set), the support coordinates to fit, and which leaves move the direction and
# each entry. A leaf is relevant to an entry if perturbing it (a probe result,
# shared across directions) moves that entry; each entry is then fit over only its
# own variables, so a wide direction with narrow entries stays cheap.
.sym_direction_relevance <- function(f, ref, pivots, zSlots, nLeaves, point0,
                                     relProbe, residueFn) {
  P <- .symPrimes[1]
  zvals0 <- if (is.null(zSlots)) NULL else point0[zSlots + 1L]
  nv <- function(rp, p, zvals = zvals0) if (is.null(residueFn)) .sym_null_residues(rp, f, p)
                        else residueFn(rp, p, zvals)
  base_nv <- nv(ref, P)
  if (is.null(base_nv)) base_nv <- .sym_null_residues(ref, f, P)
  supportCols <- setdiff(which(base_nv != 0) - 1L, f)
  relevant <- integer(0)
  relByEntry <- replicate(length(supportCols), integer(0), simplify = FALSE)
  for (li in seq_len(nLeaves)) {
    rp <- relProbe[[li]]$rp; zvp <- relProbe[[li]]$zvals
    nvp <- if (!is.null(rp) && isTRUE(rp$ok) &&
               identical(as.integer(rp$pivots), as.integer(pivots))) nv(rp, P, zvp)
           else NULL
    if (is.null(nvp)) {
      relevant <- c(relevant, li)
      relByEntry <- lapply(relByEntry, function(s) c(s, li))
      next
    }
    if (any(nvp != base_nv)) relevant <- c(relevant, li)
    for (i in which(nvp[supportCols + 1L] != base_nv[supportCols + 1L]))
      relByEntry[[i]] <- c(relByEntry[[i]], li)
  }
  list(base_nv = base_nv, supportCols = supportCols, relevant = relevant,
       relByEntry = relByEntry)
}


# Points (over the union of relevant leaves) needed for the dense fit of the widest
# entry, and whether any entry is dense at all.
.sym_dense_need <- function(relByEntry, ctrl) {
  isDense <- function(r) length(r) >= 1L && length(r) <= ctrl$relevanceCap
  denseRel <- max(0L, vapply(relByEntry,
                             function(r) if (isDense(r)) length(r) else 0L, integer(1)))
  anyDense <- any(vapply(relByEntry, isDense, logical(1)))
  maxNeed <- if (anyDense)
    2L * choose(denseRel + ctrl$degreeCap, denseRel) - 1L + ctrl$sampleSlack else 0L
  list(maxNeed = as.integer(maxNeed), anyDense = anyDense)
}


# One dense sample bank shared across all free-column directions: points over the
# union of their relevant leaves, sized for the widest entry, each evaluated once
# per prime (a single kernel evaluation yields every direction's residues). Only
# pivot-consistent points are kept; the full kernel result is stored per point and
# prime so each direction extracts its own residues later. Returns NULL when the
# bank could not be filled (the caller then samples per direction).
.sym_build_shared_bank <- function(union, needPts, point0, pool, poolNext, kbatch,
                                   pivots, NtUsed, nP) {
  nu <- length(union)
  U <- matrix(0L, 0L, nu); points <- list(); rps <- list()
  tries <- 0L
  while (length(points) < needPts && tries < 30L * needPts) {
    chunk <- max(1L, min(needPts - length(points) + 5L, 30L * needPts - tries))
    tries <- tries + chunk
    uv <- lapply(seq_len(chunk), function(ci) {
      u <- pool(poolNext + seq_len(nu) - 1L); poolNext <<- poolNext + nu; u
    })
    cand <- lapply(uv, function(u) { pt <- point0; pt[union] <- u; pt })
    res <- kbatch(rep(cand, each = nP), rep(.symPrimes, times = chunk), NtUsed)
    for (ci in seq_len(chunk)) {
      if (length(points) >= needPts) break
      good <- TRUE; rpc <- vector("list", nP)
      for (j in seq_len(nP)) {
        rp <- res[[(ci - 1L) * nP + j]]
        if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots))) {
          good <- FALSE; break
        }
        rpc[[j]] <- rp
      }
      if (!good) next
      U <- rbind(U, uv[[ci]])
      points[[length(points) + 1L]] <- cand[[ci]]
      rps[[length(rps) + 1L]] <- rpc
    }
  }
  if (length(points) < needPts) return(list(ok = FALSE, poolNext = poolNext))
  list(ok = TRUE, U = U, points = points, rps = rps, union = union,
       poolNext = poolNext)
}


.sym_interpolate_direction <- function(f, ref, pivots, znames, nz, zSlots, leafNames,
                                       nLeaves, point0, pool, poolNext, NtUsed,
                                       kcall, spy, relProbe = NULL,
                                       residueFn = NULL, ctrl = reconstControl(),
                                       kbatch = NULL, sharedBank = NULL) {
  if (is.null(kbatch))
    kbatch <- function(pl, pv, Nt) Map(function(pp, pr) kcall(pp, pr, Nt), pl, pv)
  P <- .symPrimes[1]
  zvalsOf <- function(pt) if (is.null(zSlots)) NULL else pt[zSlots + 1L]
  zvals0 <- zvalsOf(point0)
  # `nv` is the direction's nullspace vector at a kernel result: the raw
  # free-column residue, or the canonicalised residual row when residueFn is set.
  # log-gauge residue functions also need the sample point's z-values. `f` is the
  # gauge column whose entry is fixed to 1 (the free column, or the canonical
  # pivot). supportCols are the remaining nonzero coordinates to fit.
  nv <- function(rp, p, zvals = zvals0) if (is.null(residueFn)) .sym_null_residues(rp, f, p)
                        else residueFn(rp, p, zvals)
  rel <- .sym_direction_relevance(f, ref, pivots, zSlots, nLeaves, point0,
                                  relProbe, residueFn)
  base_nv <- rel$base_nv; supportCols <- rel$supportCols
  relevant <- rel$relevant; relByEntry <- rel$relByEntry
  fallback <- function(reason = NULL) {
    list(poolNext = poolNext,
         entry = list(support = sort(znames[base_nv != 0]), type = "general",
                      closedForm = FALSE, reason = reason))
  }

  # a constant entry (no relevant leaf) is reconstructed from the base point
  constEntry <- function(col) {
    res <- vapply(.symPrimes, function(pj) {
      rp <- kcall(point0, pj, NtUsed)
      if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots)))
        return(NA_real_)
      nvp <- nv(rp, pj)
      if (is.null(nvp)) return(NA_real_)
      nvp[col + 1L]
    }, numeric(1))
    if (anyNA(res)) return(NULL)
    rec <- symRatRecon(matrix(as.integer(res), 1L), as.integer(.symPrimes))
    if (rec$den[1] == "0") return(NULL)
    if (rec$den[1] == "1") rec$num[1] else paste0(rec$num[1], "/", rec$den[1])
  }

  if (!length(relevant)) {
    entries <- list()
    for (i in seq_along(supportCols)) {
      e <- constEntry(supportCols[i])
      if (is.null(e)) return(fallback("a constant entry could not be lifted from its residues"))
      entries[[znames[supportCols[i] + 1L]]] <- e
    }
    entries[[znames[f + 1L]]] <- "1"
    return(list(poolNext = poolNext,
                entry = list(support = sort(names(entries)), vector = entries,
                             type = "general", closedForm = TRUE)))
  }

  maxRel <- max(0L, vapply(relByEntry, length, integer(1)))
  if (length(relevant) > ctrl$relevanceCapDir || maxRel > ctrl$relevanceCapSparse)
    return(fallback(sprintf(
      paste("direction couples %d variables (an entry up to %d), above",
            "relevanceCapDir=%d / relevanceCapSparse=%d; raise the relevant cap"),
      length(relevant), maxRel, ctrl$relevanceCapDir, ctrl$relevanceCapSparse)))
  nrel <- length(relevant)
  nP <- length(.symPrimes)
  dn <- .sym_dense_need(relByEntry, ctrl)
  sampleU <- matrix(0L, 0L, nrel)
  rstore <- lapply(seq_along(supportCols), function(i) matrix(0L, 0L, nP))
  if (dn$anyDense) {
    maxNeed <- dn$maxNeed
    if (!is.null(sharedBank)) {
      # draw the dense fit from the shared bank: project its union points to this
      # direction's relevant leaves and extract this direction's residues from the
      # stored kernel results, so one kernel evaluation per point serves every
      # direction
      nAcc <- nrow(sharedBank$U)
      if (nAcc < maxNeed)
        return(fallback("the shared bank holds fewer points than this direction needs"))
      sampleU <- sharedBank$U[, match(relevant, sharedBank$union), drop = FALSE]
      rstore <- lapply(seq_along(supportCols), function(i) matrix(0L, nAcc, nP))
      for (a in seq_len(nAcc)) {
        zvc <- zvalsOf(sharedBank$points[[a]])
        for (j in seq_len(nP)) {
          nvp <- nv(sharedBank$rps[[a]][[j]], .symPrimes[j], zvc)
          if (is.null(nvp)) return(fallback("shared-bank residue extraction failed"))
          for (i in seq_along(supportCols)) rstore[[i]][a, j] <- nvp[supportCols[i] + 1L]
        }
      }
    } else {
      # accumulate accepted samples in chunks: draw a batch of candidate points,
      # solve every (candidate, prime) pair at once, then keep the candidates whose
      # pivots match at all primes. A pivot-shifting candidate is discarded and the
      # loop refills from the next chunk.
      tries <- 0L
      while (nrow(sampleU) < maxNeed && tries < 30L * maxNeed) {
        chunk <- max(1L, min(maxNeed - nrow(sampleU) + ctrl$sampleSlack,
                             30L * maxNeed - tries))
        tries <- tries + chunk
        uv <- lapply(seq_len(chunk), function(ci) {
          u <- pool(poolNext + seq_len(nrel) - 1L); poolNext <<- poolNext + nrel
          u
        })
        cand <- lapply(uv, function(u) { pt <- point0; pt[relevant] <- u; pt })
        res <- kbatch(rep(cand, each = nP), rep(.symPrimes, times = chunk), NtUsed)
        for (ci in seq_len(chunk)) {
          if (nrow(sampleU) >= maxNeed) break
          good <- TRUE
          rrow <- matrix(0L, length(supportCols), nP)
          zvc <- zvalsOf(cand[[ci]])
          for (j in seq_len(nP)) {
            rp <- res[[(ci - 1L) * nP + j]]
            if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots))) {
              good <- FALSE; break
            }
            nvp <- nv(rp, .symPrimes[j], zvc)
            if (is.null(nvp)) { good <- FALSE; break }
            rrow[, j] <- nvp[supportCols + 1L]
          }
          if (!good) next
          sampleU <- rbind(sampleU, uv[[ci]])
          for (i in seq_along(supportCols)) rstore[[i]] <- rbind(rstore[[i]], rrow[i, ])
        }
      }
      if (nrow(sampleU) < maxNeed)
        return(fallback(paste("dense sampling could not collect enough",
                              "pivot-consistent points; raise sampleSlack or probeRetries")))
    }
  }

  # reconstruct each entry over its own relevant variables: a constant from the
  # base point, a narrow entry by the dense fit, a wide one by sparse Laurent
  entries <- list()
  for (i in seq_along(supportCols)) {
    reli <- relByEntry[[i]]
    if (!length(reli)) {
      e <- constEntry(supportCols[i])
      if (is.null(e)) return(fallback("a constant entry could not be lifted from its residues"))
      entries[[znames[supportCols[i] + 1L]]] <- e
      next
    }
    if (length(reli) > ctrl$relevanceCap) {
      # single-monomial denominator (Laurent) first, then a general denominator
      e <- .sym_sparse_entry(reli, supportCols[i], f, point0, leafNames, NtUsed,
                             kcall, pivots, residueFn, ctrl, zSlots)
      if (is.null(e))
        e <- .sym_general_rational_entry(reli, supportCols[i], f, point0,
                                         leafNames, NtUsed, kcall, pivots,
                                         residueFn, ctrl, zSlots)
      if (is.null(e)) return(fallback(sprintf(
        paste("an entry couples %d variables and the sparse fit hit its caps",
              "(laurentDegNum=%d, generalDegNum=%d, generalDegDen=%d, termCap=%d);",
              "raise these, or raise relevanceCap=%d to use the dense fit"),
        length(reli), ctrl$laurentDegNum, ctrl$generalDegNum, ctrl$generalDegDen,
        ctrl$termCap, ctrl$relevanceCap)))
      if (!is.null(spy)) e <- .sym_simplify(e, spy)
      entries[[znames[supportCols[i] + 1L]]] <- e
      next
    }
    cols_i <- match(reli, relevant)
    vars_i <- leafNames[reli]
    fitted <- NULL
    for (degree in 0:ctrl$degreeCap) {
      mons <- .sym_mono_table(length(reli), degree)
      need <- 2L * nrow(mons) - 1L + ctrl$sampleSlack
      rec <- .sym_reconstruct_entry(
        matrix(as.integer(sampleU[seq_len(need), cols_i, drop = FALSE]), need),
        matrix(as.integer(mons), nrow(mons)),
        rstore[[i]][seq_len(need), , drop = FALSE], .symPrimes)
      if (!is.null(rec)) { fitted <- list(rec = rec, mons = mons); break }
    }
    if (is.null(fitted))
      return(fallback(sprintf("an entry exceeded degreeCap=%d; raise it",
                              ctrl$degreeCap)))
    numStr <- .sym_poly_string(fitted$rec$numCoefN, fitted$rec$numCoefD,
                               fitted$mons, vars_i)
    denStr <- .sym_poly_string(fitted$rec$denCoefN, fitted$rec$denCoefD,
                               fitted$mons, vars_i)
    expr <- if (denStr == "1") numStr else paste0("(", numStr, ")/(", denStr, ")")
    if (!is.null(spy)) expr <- .sym_simplify(expr, spy)
    entries[[znames[supportCols[i] + 1L]]] <- expr
  }

  entries[[znames[f + 1L]]] <- "1"
  list(poolNext = poolNext,
       entry = list(support = sort(names(entries)), vector = entries,
                    type = "general", closedForm = TRUE))
}


# Resolve the condition grid, event list and initial values into the per-segment
# inputs of compileObservabilityTapeMulti, ordered into one chain per condition.
# The post-t0 timeline of each condition is split into analytic segments at the
# distinct post-t0 event times; the state is then propagated exactly across each
# inter-event gap by the kernel (generic gap length as a formal power-series
# variable), so later segments are seeded from the propagated state, not from free
# coordinates. The expansion start point t0 (default: earliest event time, else 0)
# buckets events: events at or before t0 set the pre-stimulus regime and compose
# the first segment's start-point initial condition (or, under equilibrate, the
# resting steady state with the t0 dose on top); events after t0 open further
# segments, each carrying the regime in force at its left endpoint.
#
# Returns flat per-segment lists aligned with the tapes, plus chainOf/posInChain
# grouping them into per-condition chains for the propagation kernel. Only the
# first segment of a chain carries a real seed; later segments get a dummy ic0 (so
# no spurious free coordinates) and are seeded by propagation.
.sym_resolve_conditions <- function(conditions, events, initial, symbols, states,
                                    constStates = character(0),
                                    forcings = character(0), t0 = NULL,
                                    equilibrate = FALSE) {
  if (is.null(conditions)) conditions <- data.frame(row.names = "c1")
  else conditions <- as.data.frame(conditions, stringsAsFactors = FALSE)
  conds <- rownames(conditions)
  if (is.null(conds) || !length(conds))
    conds <- as.character(seq_len(max(1L, nrow(conditions))))
  cols <- colnames(conditions)
  subCols <- intersect(cols, symbols)
  initial <- if (is.null(initial)) NULL else as.eqnvec(initial)
  dynStates <- setdiff(states, constStates)

  ev <- if (is.null(events)) NULL else as.data.frame(events, stringsAsFactors = FALSE)
  tnum <- if (is.null(ev) || !nrow(ev)) numeric(0)
          else suppressWarnings(as.numeric(as.character(ev$time)))
  if (is.null(t0)) t0 <- if (length(tnum) && any(!is.na(tnum)))
    min(tnum, na.rm = TRUE) else 0
  t0 <- as.numeric(t0)
  evVar    <- if (is.null(ev) || !nrow(ev)) character(0) else as.character(ev$var)
  evMethod <- if (is.null(ev) || !nrow(ev)) character(0) else as.character(ev$method)
  isSwitch <- evVar %in% constStates

  cell <- function(ci, col) as.character(conditions[ci, col])
  resolve <- function(val, ci) {
    val <- as.character(val)
    if (val %in% cols) cell(ci, val) else val
  }

  # post-t0 segment boundaries: the distinct event times strictly after t0
  postT <- if (length(tnum)) sort(unique(tnum[!is.na(tnum) & tnum > t0])) else numeric(0)
  segTimes <- c(t0, postT)               # left endpoints of S0, S1, ...
  nSeg <- length(segTimes)
  # dummy seed for a propagation-seeded later segment: not a free unknown
  dummyIc0 <- setNames(as.list(rep("0", length(dynStates))), dynStates)

  subsList <- list(); ic0List <- list(); ev0List <- list(); segEvList <- list()
  equilList <- logical(0); chainOf <- integer(0); posInChain <- integer(0)
  for (k in seq_along(conds)) {
    subsBase <- list()
    for (col in subCols) subsBase[[col]] <- cell(k, col)
    for (j in seq_len(nSeg)) {
      tau <- segTimes[j]
      # regime in force during this segment: the latest switch value <= tau for
      # each constant state. A forcing defaults to 0 before its first event fires,
      # so it stays a baked constant (never a free coordinate) in early segments.
      subs <- subsBase
      for (cs in constStates) {
        si <- which(isSwitch & evVar == cs & evMethod == "replace" &
                    !is.na(tnum) & tnum <= tau)
        if (length(si)) subs[[cs]] <- resolve(ev$value[si[which.max(tnum[si])]], k)
        else if (cs %in% forcings) subs[[cs]] <- "0"
      }
      if (j == 1L) {
        # S0: t0 events compose the start-point initial condition (and feed the
        # equilibrate solver as the dose applied on top of the resting state)
        leIdx <- which(!isSwitch & !is.na(tnum) & tnum <= t0)
        ev0 <- lapply(leIdx, function(i) list(
          var = evVar[i], method = evMethod[i], value = resolve(ev$value[i], k)))
        ic0 <- list()
        for (X in dynStates) {
          isForcing <- X %in% forcings
          hasInit <- !is.null(initial) && X %in% names(initial)
          ops <- leIdx[evVar[leIdx] == X]
          if (!isForcing && !hasInit && !length(ops)) next
          cur <- if (isForcing) "0" else if (hasInit) as.character(initial[[X]]) else X
          for (i in ops) cur <- .sym_compose_event(cur, evMethod[i], resolve(ev$value[i], k))
          ic0[[X]] <- cur
        }
      } else {
        # later segment: propagation-seeded, dummy ic0 (no free coordinates). State
        # doses at this boundary are applied to the propagated state by the kernel.
        ev0 <- list(); ic0 <- dummyIc0
        evIdx <- which(!isSwitch & !is.na(tnum) & tnum == tau)
        segEv <- lapply(evIdx, function(i) list(
          var = evVar[i], method = evMethod[i], value = resolve(ev$value[i], k)))
      }
      subsList[[length(subsList) + 1L]] <- subs
      ic0List[[length(ic0List) + 1L]] <- ic0
      ev0List[[length(ev0List) + 1L]] <- ev0
      segEvList[[length(segEvList) + 1L]] <- if (j == 1L) list() else segEv
      equilList <- c(equilList, isTRUE(equilibrate))
      chainOf <- c(chainOf, k); posInChain <- c(posInChain, j)
    }
  }
  list(subs = subsList, ic0 = ic0List, segEquil = equilList, events0 = ev0List,
       segEvents = segEvList, chainOf = chainOf, posInChain = posInChain,
       conditions = conds, nConditions = length(conds), nGaps = nSeg - 1L)
}

# compose one event onto a current initial-value expression string
.sym_compose_event <- function(cur, method, value) {
  switch(as.character(method),
    replace  = value,
    add      = paste0("(", cur, ") + (", value, ")"),
    multiply = paste0("(", cur, ") * (", value, ")"),
    stop("unsupported event method '", method, "'"))
}


# which engine produced a result, from the fields it carries
.sym_engine <- function(x) {
  if (!is.null(x$rank)) "observability"
  else if (!is.null(x$count)) "scaling"
  else "liesym"
}

# one formatted line for a non-identifiability / symmetry direction
.sym_direction_line <- function(d) {
  if (!is.null(d$vector)) {
    txt <- paste(vapply(names(d$vector),
                        function(k) paste0(k, " : ", d$vector[[k]]),
                        character(1)), collapse = ",  ")
    paste0("[", d$type, "] ", txt)
  } else {
    paste0("[", d$type, ", closed form not found] involves: ",
           paste(d$support, collapse = ", "))
  }
}

# print one line indented and wrapped to the console width
.sym_cat_wrapped <- function(line, indent = 2L)
  cat(strwrap(line, indent = indent, exdent = indent + 4L), "", sep = "\n")

# print a direction line and, when it could not be closed, the reason below it
.sym_print_direction <- function(d) {
  out <- strwrap(.sym_direction_line(d), indent = 2L, exdent = 6L)
  if (!isTRUE(d$closedForm) && !is.null(d$reason))
    out <- c(out, strwrap(paste0("reason: ", d$reason), indent = 6L, exdent = 10L))
  cat(out, "", sep = "\n")
}

#' @export
print.symmetryDetection <- function(x, ...) {
  engine <- .sym_engine(x)
  if (engine == "observability") {
    head <- if (isTRUE(x$identifiable))
      sprintf("structurally identifiable (rank %d/%d)", x$rank, x$dim)
    else
      sprintf("rank %d/%d, %d non-identifiable direction(s)",
              x$rank, x$dim, length(x$nonIdentifiable))
    n <- length(x$nonIdentifiable)
  } else if (engine == "scaling") {
    head <- sprintf("%d scaling symmetry/ies", x$count)
    n <- x$count
  } else {
    head <- sprintf("%d transformation(s)", length(x))
    n <- length(x)
  }
  cat(sprintf("<symmetryDetection: %s> %s\n", engine, head))
  nopen <- if (engine == "observability")
    sum(vapply(x$nonIdentifiable,
               function(d) !isTRUE(d$closedForm) && !is.null(d$reason), logical(1)))
  else 0L
  if (nopen > 0L)
    cat(sprintf("  %d direction(s) without a closed form; summary() shows why\n", nopen))
  else if (n > 0L)
    cat("  (call summary() for the directions)\n")
  invisible(x)
}

#' @export
summary.symmetryDetection <- function(object, ...) {
  engine <- .sym_engine(object)
  cat(strrep("-", 60), "\n", sep = "")
  if (engine == "observability") {
    cat(sprintf("Observability (exact analytic): rank %d / %d  (Lie order %d)\n",
                object$rank, object$dim, object$lieOrderUsed))
    if (!is.null(object$conditions) && object$conditions > 1L)
      cat(sprintf("  %d conditions, %d segments, gap order %d\n",
                  object$conditions, object$segments, object$gapOrderUsed))
    if (isTRUE(object$identifiable)) {
      cat("Model is structurally locally identifiable (full rank).\n")
      return(invisible(object))
    }
    cat(sprintf("%d structural non-identifiability direction(s):\n",
                length(object$nonIdentifiable)))
    for (d in object$nonIdentifiable) .sym_print_direction(d)
  } else if (engine == "scaling") {
    cat(sprintf("%d scaling symmetry/ies (exact integer kernel):\n", object$count))
    for (d in object$nonIdentifiable) .sym_print_direction(d)
  } else {
    cat(sprintf("%d Lie-point transformation(s):\n", length(object)))
    for (tr in object) {
      vmap <- tr$infinitesimals
      txt <- paste(vapply(names(vmap),
                          function(k) paste0(k, " : ", vmap[[k]]),
                          character(1)), collapse = ",  ")
      flag <- if (isTRUE(tr$verified)) "" else " (unverified)"
      .sym_cat_wrapped(paste0("[", sub("^Type: ", "", tr$type), flag, "] ", txt))
    }
  }
  invisible(object)
}
