
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
#' @param events Optional [eventlist]. For `"observability"` with `conditions`,
#'   a `replace` event at a known time `<= 0` defines the initial value of its
#'   target state per condition; the event value is looked up in `conditions`
#'   when it names a grid column. This is how the boolean switches and the
#'   stimulus level of a condition grid reach the analysis.
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
#' @param t0 Optional numeric start time of the local Taylor expansion (default:
#'   the earliest event time, else 0). Events are bucketed relative to `t0`: a
#'   dynamic state's event at `t0` is the dose that perturbs its initial value;
#'   constant-state (switch) events are substitutions at any time.
#' @param equilibrate Logical. For `"observability"`, start the states at a steady
#'   state of `f` with forcings held at 0, solved exactly over a finite field
#'   without a symbolic solution (parameters are never solved for); t0 events such
#'   as a dose are applied on top. A free power or Hill exponent is supported: the
#'   term `base^exp` is recast to a rational coordinate `E = base^exp` with a
#'   companion `L = log(base)`, so the verdict (which parameters, including the
#'   exponent, are identifiable) is exact. With `closedForm`, `E` and `L` are
#'   sampled as coordinates and back-substituted, so a non-identifiable direction
#'   is returned in closed form even when it involves `base^exp` or `log(base)`.
#'   State initial conditions given in `trafo` are ignored when equilibrating.
#' @param ansatz Type of infinitesimal ansatz: `"uni"`, `"par"` or `"multi"`
#'   (`"liesym"`).
#' @param pMax Maximal degree of the infinitesimal ansatz (`"liesym"`).
#' @param lieOrder Integer `N >= 0`. Also require the symmetry generator to
#'   annihilate the `k`-th Lie derivative `L^k g` for `k = 1..N` (`"liesym"`;
#'   `"observability"` auto-saturates instead).
#' @param reduceCQ Logical. If the model is an [eqnlist] with conserved
#'   quantities, eliminate one pivot species per conserved quantity first.
#' @param backend `"symengine"` (faster, falls back to sympy) or `"sympy"`.
#' @param exact Logical. Use exact modular linear algebra (`TRUE`) or the legacy
#'   floating-point path (`FALSE`) for `"liesym"`.
#' @param verify Logical. Symbolically verify each `"liesym"` generator.
#' @param point Optional named numeric vector: the evaluation point for the
#'   symbolic `"liesym"` and `"scaling"` engines (default: a deterministic
#'   generic point). It does not apply to `"observability"`, which generates and
#'   certifies its own generic points over a finite field internally.
#' @param closedForm Logical. For `"observability"`, reconstruct the non-scaling
#'   directions as exact rational functions of the parameters, not just their
#'   support. Scaling directions are peeled exactly from the integer kernel and
#'   reported in closed form regardless of this flag. A narrow entry is fit
#'   densely; an entry coupling many variables is recovered by sparse
#'   (Ben-Or-Tiwari) interpolation. Every reconstructed direction is certified
#'   against the nullspace at a fresh prime; one that cannot be reconstructed or
#'   certified is returned with its support and `closedForm = FALSE` instead.
#' @param cores Number of threads. For `"observability"` the per-condition Taylor
#'   build is parallelised over conditions; for `"liesym"` and `"scaling"` it is
#'   the worker count of the symbolic search.
#' @param allTrafos Do not drop transformations with a common parameter factor
#'   (`"liesym"`).
#'
#' @return For `"observability"`: a list with `identifiable`, `rank`, `dim`,
#'   `lieOrderUsed` and `nonIdentifiable`. Each direction carries a `type`:
#'   `"scaling"` directions give the integer toric weights `c_i` (the symmetry is
#'   `z_i -> lambda^{c_i} z_i`); `"general"` directions give the reconstructed
#'   tangent components as rational expressions (or only the `support` when
#'   `closedForm = FALSE`). The multi-condition path adds `conditions`, the number
#'   of distinct conditions stacked. For `"liesym"`: a list of generators, each with
#'   `infinitesimals`, `transformation`, `type` and `verified`. Both also print
#'   a summary.
#'
#' @note This function is under active development. Its interface, defaults and
#'   output structure may still change, and some paths (notably the closed-form
#'   reconstruction of non-scaling directions) are being extended.
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
#' symmetryDetection(eq, eqnvec(Aobs = "alpha * A"))                   # observability
#' symmetryDetection(eq, eqnvec(Aobs = "alpha * A"), method = "liesym")
#' symmetryDetection(eq, eqnvec(Aobs = "alpha * A"), method = "scaling")
#'
#' # A steady-state initial condition as an expression: a state-named trafo entry
#' # is the initial condition, seeded with its parameter sensitivities.
#' symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
#'                   trafo = eqnvec(x = "b/a"), closedForm = TRUE)
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
#' symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
#'                   method = "observability", events = dose,
#'                   conditions = data.frame(dose = 2, row.names = "stim"))
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
                              ansatz = "uni", pMax = 2, lieOrder = 0L,
                              reduceCQ = TRUE,
                              backend = c("symengine", "sympy"),
                              exact = TRUE, verify = TRUE, point = NULL,
                              closedForm = FALSE, cores = 1, allTrafos = FALSE) {

  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required for symmetryDetection().")
  method <- match.arg(method)
  backend <- match.arg(backend)
  equilibrate <- isTRUE(equilibrate)

  # warn about arguments that do not apply to the chosen engine instead of
  # silently ignoring them
  supplied <- setdiff(names(match.call())[-1], "")
  applies <- switch(method,
    observability = c("events", "conditions", "t0", "equilibrate"),
    liesym  = c("ansatz", "pMax", "lieOrder", "exact", "verify", "allTrafos",
                "point", "backend"),
    scaling = c("point", "backend"))
  methodSpecific <- c("events", "conditions", "t0",
                      "equilibrate", "ansatz", "pMax", "lieOrder",
                      "exact", "verify", "allTrafos", "point", "backend")
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
                                   constStates, forcings, t0)
    multi <- sd$compileObservabilityTapeMulti(
      model = toLines(fdyn), observation = toLines(gobs),
      conditionSubs = res$subs, conditionIC0 = res$ic0,
      fixed = if (length(fixed)) fixed else NULL,
      parameters = if (length(parameters)) parameters else NULL,
      equilibrate = equilibrate,
      forcings = if (length(forcings)) forcings else NULL)
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
    return(.observability_analytic_multi(multi, spy = spy, verbose = TRUE,
                                         closedForm = closedForm, sd = sd,
                                         cores = cores,
                                         equilZeroStates = equilZeroStates,
                                         t0events = res$events0))
  }

  # liesym / scaling: symbolic engines on the (trafo-substituted) f and g;
  # forcings are the externally driven (input) states
  sd$symmetryDetectiondMod(
    model       = toLines(fdyn),
    observation = toLines(gobs),
    ansatz      = ansatz,
    pMax        = as.integer(pMax),
    inputs      = if (length(forcings)) forcings else NULL,
    fixed       = fixed,
    parallel    = as.integer(cores),
    allTrafos   = allTrafos,
    lieOrder    = as.integer(lieOrder),
    exact       = exact,
    verify      = verify,
    backend     = backend,
    parameters  = if (length(parameters)) parameters else NULL,
    method      = method,
    point       = point,
    symbolic    = closedForm
  )
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

# .symRelevanceCap bounds the variables a single nullspace entry may couple (the
# per-entry interpolation cost); .symRelevanceCapDir bounds the variables a whole
# direction may couple. A direction can be wide while each of its entries stays
# narrow, so the two caps differ and each entry is fit over its own variables.
.symRelevanceCap <- 6L
.symRelevanceCapDir <- 24L
.symDegreeCap <- 4L
.symSampleSlack <- 5L

# an entry coupling more than .symRelevanceCap variables is reconstructed by sparse
# Laurent interpolation (Ben-Or-Tiwari) instead of the dense fit, up to
# .symRelevanceCapSparse variables; .symLaurentDegNum / .symLaurentDegDen bound the
# numerator and (single-monomial) denominator degrees of such an entry.
.symRelevanceCapSparse <- 30L
.symLaurentDegNum <- 4L
.symLaurentDegDen <- 2L
.symTermCap <- 60L
.symLaurentCandCap <- 200000L
# general (multi-term denominator) sparse-rational path: numerator/denominator
# degree bounds for the shift + Cauchy + Ben-Or-Tiwari reconstruction
.symGeneralDegNum <- 4L
.symGeneralDegDen <- 3L


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
                                  kcall, sd) {
  if (is.null(sd) || is.null(entry$vector)) return(TRUE)
  q <- .symVerifyPrime
  rq <- kcall(point0, q, as.integer(NtUsed))
  if (!isTRUE(rq$ok)) return(TRUE)
  aq <- .sym_null_residues(rq, f, q)
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


# Find a generic base point where the rank is maximal over the primes and
# saturate the Lie order until the rank stops growing. kcall(point, p, Nt) must
# return the kernel list (ok, R, pivots, rank, dim). Returns NULL if no usable
# point is found, else the reference reduction and the certified rank.
.sym_saturate_certify <- function(kcall, nLeaves, nz) {
  P <- .symPrimes[1]
  pool <- .sym_pool()
  point0 <- pool(seq_len(nLeaves))
  poolNext <- nLeaves + 1L

  saturate <- function(point) {
    prev <- -1L; Nt <- 1L; res <- NULL
    repeat {
      r <- kcall(point, P, Nt)
      if (!isTRUE(r$ok)) return(NULL)
      res <- r
      if (r$rank >= nz || (r$rank == prev && Nt >= 2L) || Nt > nz + 1L) break
      prev <- r$rank; Nt <- Nt + 1L
    }
    list(res = res, Nt = Nt)
  }

  # several generic points may be tried before one admits a steady-state point
  # over GF(p) for every condition
  sat <- NULL
  for (attempt in 1:50) {
    sat <- saturate(point0)
    if (!is.null(sat)) break
    point0 <- pool(poolNext + seq_len(nLeaves) - 1L); poolNext <- poolNext + nLeaves
  }
  if (is.null(sat)) return(NULL)
  NtUsed <- sat$Nt

  rankMax <- sat$res$rank
  for (pj in .symPrimes[-1]) {
    rj <- kcall(point0, pj, NtUsed)
    if (isTRUE(rj$ok)) rankMax <- max(rankMax, rj$rank)
  }
  while (sat$res$rank < rankMax) {
    point0 <- pool(poolNext + seq_len(nLeaves) - 1L); poolNext <- poolNext + nLeaves
    sat <- saturate(point0)
    if (is.null(sat)) return(NULL)
    NtUsed <- sat$Nt
  }
  list(ref = sat$res, NtUsed = NtUsed, point0 = point0, pool = pool,
       poolNext = poolNext, rank = sat$res$rank, pivots = sat$res$pivots)
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
.observability_analytic_multi <- function(multi, spy = NULL, verbose = TRUE,
                                          closedForm = TRUE, sd = NULL,
                                          cores = 1, equilZeroStates = character(0),
                                          t0events = list()) {
  cores <- as.integer(max(1L, cores))
  nLeaves <- as.integer(multi$nLeaves)
  nStates <- as.integer(multi$nStates)
  zSlots <- as.integer(multi$zSlots)
  znames <- as.character(multi$znames)
  nz <- length(znames)
  leafNames <- as.character(multi$leafNames)
  ssConstraint <- isTRUE(multi$equilibrate)
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
    if (!ssConstraint)
      out <- c(out, list(
        icOp = as.integer(t$icOp), icA = as.integer(t$icA),
        icB = as.integer(t$icB), icCnum = as.character(t$icCnum),
        icCden = as.character(t$icCden), icOut = as.integer(t$icOut)))
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
    w <- nz + 1L
    # power/Hill recast atoms (E = base^exp, L = log base) are extra states held
    # generic in the f = 0 solve; their generic values are fixed distinct primes,
    # reduced per prime. realStateNames are the genuine states solved for.
    recast <- multi$powerRecast
    if (is.null(recast)) recast <- list()
    realStateNames <- if (length(recast)) as.character(multi$realStateNames) else stateNames
    eNames <- vapply(recast, function(r) as.character(r$E), "")
    lNames <- unique(vapply(recast, function(r) as.character(r$L), ""))
    solveParamNames <- c(paramNames, eNames)
    if (length(recast)) {
      leafNamesAug <- c(leafNames, eNames, lNames)
      nAug <- length(leafNamesAug)
    }
    ssCache <- new.env(parent = emptyenv())
    ssWhy <- NULL              # last steady-state failure reason, for diagnostics
    seedTapes <- function(point, p) {
      pmod <- as.numeric(point) %% p
      paramVals <- as.list(setNames(pmod[seq_along(leafNames)], leafNames))
      lValsList <- NULL
      if (length(recast)) {
        for (i in seq_along(eNames)) paramVals[[eNames[i]]] <- pmod[nLeaves + i]
        lValsList <- as.list(setNames(
          pmod[nLeaves + length(eNames) + seq_along(lNames)], lNames))
      }
      # the conditions are independent steady-state solves; run them across cores
      # (forked workers, each its own Python) since this is the dominant cost
      solveOne <- function(ci) {
        evC <- if (ci <= length(t0events) && length(t0events[[ci]]))
          t0events[[ci]] else NULL
        sd$solveSteadyStateModular(
          model = models[[ci]], stateNames = realStateNames,
          paramNames = solveParamNames, paramVals = paramVals, prime = p,
          forcings = if (length(solveHeld)) solveHeld else NULL,
          t0events = evC,
          recast = if (length(recast)) recast else NULL, lVals = lValsList)
      }
      # conditions are solved serially: the steady-state solve is Python/sympy
      # (GIL-bound) and the process is large, so forking per call costs more than
      # it saves; running in-process instead lets the solve caches persist.
      sols <- lapply(seq_along(tapes), solveOne)
      out <- vector("list", length(tapes))
      for (ci in seq_along(tapes)) {
        sol <- sols[[ci]]
        if (inherits(sol, "try-error") || !is.list(sol) || !isTRUE(sol$ok)) {
          ssWhy <<- if (!is.list(sol)) "parallel solve failed"
                    else if (is.null(sol$why)) "unknown" else as.character(sol$why)
          return(NULL)
        }
        icSeed <- matrix(0L, nStates, w)
        nReal <- length(realStateNames)
        icSeed[seq_len(nReal), 1L] <- as.integer(as.numeric(sol$xstar))
        for (c in seq_len(nz))
          icSeed[seq_len(nReal), c + 1L] <- as.integer(as.numeric(sol$dx[[znames[c]]]))
        # E, L atoms (recast): generic value with their parameter-duals
        for (rc in sol$recast) {
          er <- match(rc$E, stateNames); lr <- match(rc$L, stateNames)
          icSeed[er, 1L] <- as.integer(rc$E0)
          icSeed[lr, 1L] <- as.integer(rc$L0)
          for (c in seq_len(nz)) {
            ev <- rc$eDual[[znames[c]]]; lv <- rc$lDual[[znames[c]]]
            if (!is.null(ev)) icSeed[er, c + 1L] <- as.integer(ev)
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
    kcall <- function(point, p, Nt) {
      key <- paste(c(format(as.numeric(point), scientific = FALSE), p),
                   collapse = ",")
      seeded <- ssCache[[key]]
      if (is.null(seeded)) {
        seeded <- seedTapes(point, p)
        ssCache[[key]] <- if (is.null(seeded)) list(ok = FALSE) else seeded
      }
      if (is.list(seeded) && identical(seeded$ok, FALSE)) return(list(ok = FALSE))
      symObsNullMulti(seeded, nLeaves, nStates, zSlots,
                      as.integer(point[seq_len(nLeaves)]), p, as.integer(Nt), cores)
    }
  } else {
    kcall <- function(point, p, Nt)
      symObsNullMulti(tapes, nLeaves, nStates, zSlots, as.integer(point), p,
                      as.integer(Nt), cores)
  }

  sc <- .sym_saturate_certify(kcall, nAug, nz)
  if (is.null(sc)) {
    if (ssConstraint && !is.null(ssWhy))
      warning("symmetryDetection(): no steady-state point over the finite field ",
              "after the generic-point retries (", ssWhy, "). The equilibrate ",
              "constraint could not be evaluated.", call. = FALSE)
    return(NULL)
  }

  result <- list(method = "observability", engine = "analytic",
                 conditions = length(tapes),
                 identifiable = (sc$rank == nz), rank = as.integer(sc$rank),
                 dim = as.integer(nz), lieOrderUsed = as.integer(sc$NtUsed),
                 nonIdentifiable = list())
  if (sc$rank == nz) {
    if (verbose) .sym_print_analytic(result, sc$NtUsed)
    return(result)
  }

  P <- .symPrimes[1]
  freeCols <- setdiff(0:(nz - 1L), sc$pivots)
  N <- .sym_nullspace_basis(sc$ref, freeCols, P)

  # scalings common to every condition are exact (integer kernel) and always
  # reported in closed form; their span is excluded before any reconstruction
  scaling <- list(); Bmat <- matrix(0L, nz, 0L)
  if (!is.null(sd)) {
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
  residualFree <- integer(0)
  for (fc in freeCols) {
    bf <- .sym_null_residues(sc$ref, fc, P)
    inSpan <- ncol(Bmat) > 0L &&
      !is.null(symSolveMod(matrix(as.integer(Bmat), nz), as.integer(bf), P))
    if (!inSpan) { residualFree <- c(residualFree, fc); Bmat <- cbind(Bmat, bf) }
  }

  if (isTRUE(closedForm)) {
    interp <- list()
    # shared relevance probe: a single solve at a one-leaf perturbation yields the
    # residues of every direction, so the per-leaf scan is done once here rather
    # than once per direction (the dominant cost at scale).
    probeBase <- sc$poolNext
    relProbe <- lapply(seq_len(nAug), function(li) {
      pert <- sc$point0
      pert[li] <- sc$pool(probeBase + li - 1L)
      kcall(pert, P, sc$NtUsed)
    })
    poolNext <- probeBase + nAug
    for (fc in residualFree) {
      dir <- .sym_interpolate_direction(fc, sc$ref, sc$pivots, znames, nz,
                                        leafNamesAug, nAug, sc$point0, sc$pool,
                                        poolNext, sc$NtUsed, kcall, spy, relProbe)
      poolNext <- dir$poolNext
      e <- dir$entry
      # certify the closed form against the nullspace at a fresh prime; a failed
      # reconstruction is downgraded to support-only, never reported as wrong
      if (isTRUE(e$closedForm) &&
          !.sym_verify_direction(e, fc, znames, leafNamesAug, sc$point0,
                                 sc$NtUsed, kcall, sd)) {
        v <- .sym_null_residues(sc$ref, fc, P)
        e <- list(support = sort(znames[v != 0]), type = "general",
                  closedForm = FALSE)
      }
      # back-substitute the recast coordinates: E -> base^exp, L -> log(base), so
      # a direction fit over E, L is returned as a closed form in the parameters
      if (isTRUE(e$closedForm) && length(recast) && !is.null(sd))
        e$vector <- .sym_recast_backsub(e$vector, recast, sd)
      interp[[length(interp) + 1L]] <- e
    }
    result$nonIdentifiable <- c(scaling, interp)
  } else {
    support <- lapply(residualFree, function(fc) {
      v <- .sym_null_residues(sc$ref, fc, P)
      list(support = sort(znames[v != 0]), type = "general", closedForm = FALSE)
    })
    result$nonIdentifiable <- c(scaling, support)
  }
  if (verbose) .sym_print_analytic(result, sc$NtUsed)
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
                              kcall, pivots) {
  nvar <- length(reli)
  bases <- .sym_sieve(nvar)
  np <- length(.symPrimes)
  maxLen <- 2L * .symTermCap + 2L

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
        seqs[[pj]][k] <- .sym_null_residues(rp, f, p)[supportCol + 1L]
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
  for (dDen in 0:.symLaurentDegDen) for (dNum in seq_len(.symLaurentDegNum)) {
    if (choose(nvar + dNum, nvar) * choose(nvar + dDen, nvar) > .symLaurentCandCap)
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
                                        NtUsed, kcall, pivots) {
  nvar <- length(reli)
  np <- length(.symPrimes)
  bases <- .sym_sieve(nvar)
  s <- as.numeric(point0[reli])

  sampleR <- function(relvals, p) {
    pt <- point0; pt[reli] <- relvals %% p
    rp <- kcall(pt, p, NtUsed)
    if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots)))
      return(NULL)
    .sym_null_residues(rp, f, p)[supportCol + 1L]
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
  for (tot in 2:(.symGeneralDegNum + .symGeneralDegDen)) {
    for (a in seq.int(max(1L, tot - .symGeneralDegDen), min(tot - 1L, .symGeneralDegNum))) {
      bD <- tot - a
      if (bD < 1L || bD > .symGeneralDegDen) next
      if (!is.null(evalND(u0, P1, a, bD))) { deg <- c(a, bD); break }
    }
    if (!is.null(deg)) break
  }
  if (is.null(deg)) return(NULL)
  dN <- deg[1]; dD <- deg[2]
  if (choose(nvar + dN, nvar) + choose(nvar + dD, nvar) > .symLaurentCandCap)
    return(NULL)
  candN <- matrix(as.integer(.sym_mono_table(nvar, dN)), ncol = nvar)
  candD <- matrix(as.integer(.sym_mono_table(nvar, dD)), ncol = nvar)

  # geometric Ben-Or-Tiwari sampling of both N/D(s) and D/D(s), grown until both
  # term counts stabilise
  maxLen <- 2L * .symTermCap + 2L
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


.sym_interpolate_direction <- function(f, ref, pivots, znames, nz, leafNames,
                                       nLeaves, point0, pool, poolNext, NtUsed,
                                       kcall, spy, relProbe = NULL) {
  P <- .symPrimes[1]
  base_nv <- .sym_null_residues(ref, f, P)
  supportCols <- pivots[base_nv[pivots + 1L] != 0]
  fallback <- function() {
    supp <- sort(znames[base_nv != 0])
    list(poolNext = poolNext,
         entry = list(support = supp, type = "general", closedForm = FALSE))
  }

  # per-entry relevance: a leaf is relevant to an entry if perturbing it moves
  # that entry. The direction's relevant set is the union; each entry is fit over
  # only its own variables, so a wide direction with narrow entries stays cheap.
  # Known inputs are leaves too.
  # `relProbe`, when supplied, holds the kcall result of perturbing each leaf once,
  # shared across all directions (a single kcall yields every direction's
  # residues), so the relevance scan is not repeated per direction.
  relevant <- integer(0)
  relByEntry <- replicate(length(supportCols), integer(0), simplify = FALSE)
  for (li in seq_len(nLeaves)) {
    if (!is.null(relProbe)) {
      rp <- relProbe[[li]]
    } else {
      pert <- point0
      pert[li] <- pool(poolNext); poolNext <- poolNext + 1L
      rp <- kcall(pert, P, NtUsed)
    }
    if (is.null(rp) || !isTRUE(rp$ok) ||
        !identical(as.integer(rp$pivots), as.integer(pivots))) {
      relevant <- c(relevant, li)
      relByEntry <- lapply(relByEntry, function(s) c(s, li))
      next
    }
    nvp <- .sym_null_residues(rp, f, P)
    if (any(nvp != base_nv)) relevant <- c(relevant, li)
    for (i in which(nvp[supportCols + 1L] != base_nv[supportCols + 1L]))
      relByEntry[[i]] <- c(relByEntry[[i]], li)
  }

  # a constant entry (no relevant leaf) is reconstructed from the base point
  constEntry <- function(col) {
    res <- vapply(.symPrimes, function(pj) {
      rp <- kcall(point0, pj, NtUsed)
      if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots)))
        return(NA_real_)
      .sym_null_residues(rp, f, pj)[col + 1L]
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
      if (is.null(e)) return(fallback())
      entries[[znames[supportCols[i] + 1L]]] <- e
    }
    entries[[znames[f + 1L]]] <- "1"
    return(list(poolNext = poolNext,
                entry = list(support = sort(names(entries)), vector = entries,
                             type = "general", closedForm = TRUE)))
  }

  maxRel <- max(0L, vapply(relByEntry, length, integer(1)))
  if (length(relevant) > .symRelevanceCapDir || maxRel > .symRelevanceCapSparse)
    return(fallback())
  nrel <- length(relevant)

  # dense sample bank over the union of relevant leaves, sized for the widest
  # entry the dense fit handles; entries wider than .symRelevanceCap go sparse
  isDense <- function(r) length(r) >= 1L && length(r) <= .symRelevanceCap
  denseRel <- max(0L, vapply(relByEntry,
                             function(r) if (isDense(r)) length(r) else 0L,
                             integer(1)))
  sampleU <- matrix(0L, 0L, nrel)
  rstore <- lapply(seq_along(supportCols),
                   function(i) matrix(0L, 0L, length(.symPrimes)))
  if (any(vapply(relByEntry, isDense, logical(1)))) {
    maxMon <- choose(denseRel + .symDegreeCap, denseRel)
    maxNeed <- 2L * maxMon - 1L + .symSampleSlack
    tries <- 0L
    while (nrow(sampleU) < maxNeed && tries < 30L * maxNeed) {
      tries <- tries + 1L
      pt <- point0
      uvals <- pool(poolNext + seq_len(nrel) - 1L); poolNext <- poolNext + nrel
      pt[relevant] <- uvals
      good <- TRUE
      rrow <- matrix(0L, length(supportCols), length(.symPrimes))
      for (j in seq_along(.symPrimes)) {
        rp <- kcall(pt, .symPrimes[j], NtUsed)
        if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots))) {
          good <- FALSE; break
        }
        rrow[, j] <- .sym_null_residues(rp, f, .symPrimes[j])[supportCols + 1L]
      }
      if (!good) next
      sampleU <- rbind(sampleU, uvals)
      for (i in seq_along(supportCols)) rstore[[i]] <- rbind(rstore[[i]], rrow[i, ])
    }
    if (nrow(sampleU) < maxNeed) return(fallback())
  }

  # reconstruct each entry over its own relevant variables: a constant from the
  # base point, a narrow entry by the dense fit, a wide one by sparse Laurent
  entries <- list()
  for (i in seq_along(supportCols)) {
    reli <- relByEntry[[i]]
    if (!length(reli)) {
      e <- constEntry(supportCols[i])
      if (is.null(e)) return(fallback())
      entries[[znames[supportCols[i] + 1L]]] <- e
      next
    }
    if (length(reli) > .symRelevanceCap) {
      # single-monomial denominator (Laurent) first, then a general denominator
      e <- .sym_sparse_entry(reli, supportCols[i], f, point0, leafNames, NtUsed,
                             kcall, pivots)
      if (is.null(e))
        e <- .sym_general_rational_entry(reli, supportCols[i], f, point0,
                                         leafNames, NtUsed, kcall, pivots)
      if (is.null(e)) return(fallback())
      if (!is.null(spy)) e <- .sym_simplify(e, spy)
      entries[[znames[supportCols[i] + 1L]]] <- e
      next
    }
    cols_i <- match(reli, relevant)
    vars_i <- leafNames[reli]
    fitted <- NULL
    for (degree in 0:.symDegreeCap) {
      mons <- .sym_mono_table(length(reli), degree)
      need <- 2L * nrow(mons) - 1L + .symSampleSlack
      rec <- .sym_reconstruct_entry(
        matrix(as.integer(sampleU[seq_len(need), cols_i, drop = FALSE]), need),
        matrix(as.integer(mons), nrow(mons)),
        rstore[[i]][seq_len(need), , drop = FALSE], .symPrimes)
      if (!is.null(rec)) { fitted <- list(rec = rec, mons = mons); break }
    }
    if (is.null(fitted)) return(fallback())
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


# Resolve the condition grid, event list and initial values into the
# per-condition inputs of compileObservabilityTapeMulti. The expansion start
# point t0 (default: earliest event time, else 0) buckets events: events before
# t0 set the pre-equilibration regime (their target is substituted as a
# constant); events at t0 compose the start-point initial condition by method
# (replace / add / multiply); events after t0 do not affect the local t0 jet.
# A grid column whose name is a model symbol is a per-condition substitution; an
# event value naming a grid column is looked up there. Each state's t0+ initial
# condition is `initial[X]` (if given) composed with its t0 events, else the bare
# state symbol (free). Identical resolved conditions are de-duplicated.
.sym_resolve_conditions <- function(conditions, events, initial, symbols, states,
                                    constStates = character(0),
                                    forcings = character(0), t0 = NULL) {
  if (is.null(conditions)) conditions <- data.frame(row.names = "c1")
  else conditions <- as.data.frame(conditions, stringsAsFactors = FALSE)
  conds <- rownames(conditions)
  if (is.null(conds) || !length(conds))
    conds <- as.character(seq_len(max(1L, nrow(conditions))))
  cols <- colnames(conditions)
  subCols <- intersect(cols, symbols)
  initial <- if (is.null(initial)) NULL else as.eqnvec(initial)

  ev <- if (is.null(events)) NULL else as.data.frame(events, stringsAsFactors = FALSE)
  tnum <- if (is.null(ev) || !nrow(ev)) numeric(0)
          else suppressWarnings(as.numeric(as.character(ev$time)))
  if (is.null(t0)) t0 <- if (length(tnum) && any(!is.na(tnum)))
    min(tnum, na.rm = TRUE) else 0
  t0 <- as.numeric(t0)
  # a RHS = 0 state outside `forcings` is substituted by its value; every other
  # state composes its t0+ initial condition from the events up to t0. A forcing
  # defaults to 0; a non-forcing to its `initial` expression, else a free unknown.
  evVar <- if (is.null(ev) || !nrow(ev)) character(0) else as.character(ev$var)
  subIdx <- if (length(evVar)) which(evVar %in% constStates &
                                     as.character(ev$method) == "replace") else integer(0)
  leIdx  <- if (length(tnum)) which(!is.na(tnum) & tnum <= t0 &
                                    !(evVar %in% constStates)) else integer(0)

  cell <- function(ci, col) as.character(conditions[ci, col])
  resolve <- function(val, ci) {
    val <- as.character(val)
    if (val %in% cols) cell(ci, val) else val
  }

  subsList <- vector("list", length(conds))
  ic0List <- vector("list", length(conds))
  ev0List <- vector("list", length(conds))
  for (k in seq_along(conds)) {
    subs <- list()
    for (col in subCols) subs[[col]] <- cell(k, col)
    for (j in subIdx) subs[[evVar[j]]] <- resolve(ev$value[j], k)
    subsList[[k]] <- subs
    # the t0 events (a dose), passed to the equilibrate solver to apply on top of
    # the resting steady state
    ev0List[[k]] <- lapply(leIdx, function(j) list(
      var = evVar[j], method = as.character(ev$method[j]),
      value = resolve(ev$value[j], k)))

    ic0 <- list()
    for (X in setdiff(states, constStates)) {
      isForcing <- X %in% forcings
      hasInit <- !is.null(initial) && X %in% names(initial)
      ops <- leIdx[evVar[leIdx] == X]
      if (!isForcing && !hasInit && !length(ops)) next
      cur <- if (isForcing) "0" else if (hasInit) as.character(initial[[X]]) else X
      for (j in ops) {
        v <- resolve(ev$value[j], k)
        cur <- switch(as.character(ev$method[j]),
          replace  = v,
          add      = paste0("(", cur, ") + (", v, ")"),
          multiply = paste0("(", cur, ") * (", v, ")"),
          stop("unsupported event method '", ev$method[j], "'"))
      }
      ic0[[X]] <- cur
    }
    ic0List[[k]] <- ic0
  }

  # deparse() of a large subs/ic0 list returns several lines; collapse each to a
  # single string so the dedup key stays length one (matters for big models)
  keys <- vapply(seq_along(conds), function(k)
    paste(paste(deparse(subsList[[k]]), collapse = ""),
          paste(deparse(ic0List[[k]]), collapse = ""), sep = "|"), character(1))
  keep <- !duplicated(keys)
  list(subs = subsList[keep], ic0 = ic0List[keep], conditions = conds[keep],
       events0 = ev0List[keep])
}


.sym_print_analytic <- function(result, NtUsed) {
  cat(strrep("-", 60), "\n", sep = "")
  cat(sprintf("Observability (exact analytic): rank %d / %d  (Lie order %d)\n",
              result$rank, result$dim, NtUsed))
  if (result$identifiable) {
    cat("Model is structurally locally identifiable (full rank).\n")
    return(invisible())
  }
  cat(sprintf("%d structural non-identifiability direction(s):\n",
              length(result$nonIdentifiable)))
  for (d in result$nonIdentifiable) {
    if (isTRUE(d$closedForm)) {
      txt <- paste(vapply(names(d$vector),
                          function(k) paste0(k, " : ", d$vector[[k]]),
                          character(1)), collapse = ",  ")
      cat("  [", d$type, "] ", txt, "\n", sep = "")
    } else {
      cat("  [", d$type, ", closed form not found] involves: ",
          paste(d$support, collapse = ", "), "\n", sep = "")
    }
  }
}
