
#' Search for structural non-identifiabilities of a model
#'
#' @description Detects structural non-identifiabilities of a reaction network and
#'   its observation map, via the Python module `symmetryDetectionVersion2`
#'   (`reticulate`). The model `f`, observables `g` and optional `trafo` are given as
#'   equations (anything [as.eqnvec] accepts). The engine is chosen with `method`:
#'
#'   * `"observability"` (default): rank of the observability-identifiability matrix
#'     over a finite field; the nullspace gives the non-identifiable directions.
#'     Exact and scalable; needs rational right-hand sides and observables. A free
#'     power/Hill exponent (`base^exp`) is supported, recast internally.
#'   * `"polynomial"`: the polynomial Lie-symmetry ansatz of Merkt et al. (2015),
#'     returning the explicit generator and finite transformation (grows with `pMax`).
#'   * `"scaling"`: the scaling symmetries only, from the integer kernel of the
#'     monomial-exponent conditions (`equilibrate` does not apply).
#'   * `"translation"`: the exact additive (constant-shift) symmetries, the additive
#'     analogue of `"scaling"`.
#'
#' @param f The model right-hand sides: an [eqnlist], an [eqnvec], or a named
#'   character vector keyed by state name (anything [as.eqnvec] accepts). With an
#'   [eqnlist], conserved quantities can be reduced first (`reduceCQ`).
#' @param g The observation functions, as an [eqnvec] or named character vector.
#' @param trafo Optional parameter transformation: one [eqnvec] (applied to every
#'   condition) or, for `"observability"`, a *list* of [eqnvec]s (one per condition).
#'   A parameter-named entry is substituted into `f` and `g`; a state-named entry is
#'   that state's initial condition. A solved steady state from [steadyStates] can be
#'   passed whole, and a per-condition list is the explicit route for per-condition
#'   steady states.
#' @param method One of `"observability"` (default), `"polynomial"`, `"scaling"`
#'   or `"translation"`; see Description.
#' @param parameters Character vector of extra symbols to treat as parameters.
#' @param forcings Character vector of externally driven (input) state names. For
#'   `"observability"` a forcing is an integrated state (default initial value 0)
#'   excluded from the `f = 0` steady-state constraint; for `"polynomial"`/`"scaling"`
#'   it is an input whose infinitesimal ansatz is fixed.
#' @param events Optional [eventlist]. For `"observability"` an event value naming a
#'   `conditions` grid column is read from that grid. The Taylor expansion starts at
#'   the earliest event time; later events split the timeline into segments that are
#'   stacked, the state propagated exactly across each gap (Details).
#' @param conditions Optional data frame of experimental conditions (one row per
#'   condition; columns named by model symbols or event-value placeholders). For
#'   `"observability"` each condition is compiled with its own substitutions (a
#'   numeric cell bakes a symbol to a constant, a symbol cell renames it) and their
#'   observability matrices are stacked: a direction is non-identifiable only when
#'   unobservable in every condition (the intersection nullspace).
#' @param fixed Character vector of symbols that are known and therefore not
#'   unknowns. For `"observability"` they are excluded from the coordinates `z`:
#'   a fixed parameter is a known constant, and a fixed state keeps its dynamics
#'   but carries no unknown initial value. For `"polynomial"` their infinitesimal
#'   ansatz is set to zero.
#' @param equilibrate Logical. For `"observability"`, start the states at a steady
#'   state of `f` (forcings held at 0), solved over a finite field; the earliest
#'   events are applied on top. State initial conditions in `trafo` are then ignored.
#'   A free power/Hill exponent is supported with or without `equilibrate`.
#' @param reduceCQ Logical. Controls how a conserved moiety (which makes the
#'   `equilibrate` steady-state system `f = 0` rank-deficient) is handled for an
#'   [eqnlist]. `FALSE` (default) keeps every species a coordinate and, under
#'   `equilibrate`, holds one pivot species' resting value as its initial-value
#'   parameter (the "held-variable" parameterisation): the steady state becomes a
#'   function of that value with the other moiety species solved from `f = 0`.
#'   Following the dMod / deSolve convention that a parameter named like a state is
#'   that state's initial value, the moiety freedom is reported under the pivot
#'   species' own name. `TRUE` instead eliminates one pivot species per conserved
#'   quantity and reports the moiety freedom on a `total` parameter (the pivot
#'   species is then no longer a coordinate); it trades a simpler (often all-scaling)
#'   basis for fewer coordinates, which helps large models. Both give the same
#'   identifiability; they differ only in the moiety coordinate. If a state-named
#'   `trafo` initial condition names a species that participates in a conserved
#'   quantity, `reduceCQ` is forced to `FALSE` (with a warning): `TRUE` would
#'   otherwise eliminate a moiety species and silently discard the supplied
#'   steady-state relation, so the model would no longer start at the given steady
#'   state.
#' @param freeInitial Character vector of state names. Under the held-variable
#'   parameterisation (`equilibrate = TRUE`, `reduceCQ = FALSE`), `f = 0` pins every
#'   moiety species but one, and that one (the pivot) keeps its free resting value as
#'   its initial-value parameter. `freeInitial` chooses that pivot per conserved
#'   quantity (one species each) instead of the automatic choice, so the moiety
#'   freedom is reported on the species you name. A name that is not a valid moiety
#'   pivot is ignored, and the argument is ignored (with a warning) outside the
#'   held-variable case, where every species already carries a free initial value.
#' @param reconstruct Logical. For `"observability"`, reconstruct the non-scaling
#'   directions as exact rational functions (a free exponent via its `base^exp` /
#'   `log(base)` recast), not just their support; scalings are always exact. A
#'   direction that cannot be reconstructed or certified is returned support-only
#'   (its `explicit` field is then `FALSE`).
#' @param verify Logical (default `TRUE`). For `"observability"`, run a fast
#'   Schwartz-Zippel saturation guard that extends the Lie order past the plateau
#'   heuristic and warns if the rank is still growing (a premature stop would
#'   over-report non-identifiability); the verdict is attached as `$verification`.
#'   The margin of extra orders is set by `DMOD_SYM_VERIFY_MARGIN` (default 6).
#' @param cores Number of threads for `"observability"`. The reconstruction's
#'   finite-field steady-state solves are parallelised over the sample points (fork
#'   on Linux, PSOCK pool on Windows), and the observability kernel threads over
#'   conditions/segments; `cores` is split across these two (nested) axes so they
#'   never oversubscribe. The `"polynomial"` and `"scaling"` engines are exact and
#'   serial and ignore it.
#' @param control A [reconstControl()] list tuning the `"observability"` engine's
#'   saturation and closed-form reconstruction (relevance caps, fit degrees, term
#'   and gap-order caps). Raise the caps to recover wide or high-degree directions.
#' @param polynomial A [polynomialControl()] list tuning the `"polynomial"` engine: the
#'   infinitesimal ansatz and degree, the extra Lie-derivative order, the symbolic
#'   backend and verification.
#' @param scaling A [scalingControl()] list tuning the `"scaling"` engine: the
#'   symbolic backend.
#' @param symEngine For `"observability"`, the engine computing the matrix:
#'   `"modular"` (default; finite fields GF(p) + CRT, fast and scalable, drives the
#'   closed-form reconstruction and `verify`) or `"symbolic"` (exact sympy
#'   rational-function field, an independent cross-check for SMALL models; no recast,
#'   no `equilibrate`, no event gaps, no closed-form reconstruction).
#' @param verbose Logical (default `TRUE`, overridable with
#'   `options(dMod.sym.verbose =)`). Print the result report on return (as
#'   `print()` renders it); set `FALSE` to compute silently.
#'
#' @return An object of class `symmetryDetection`, the same shape for every
#'   `method`, holding the *verdict* at the top level and the *how* under `$info`:
#'   \describe{
#'     \item{`method`}{the engine that ran (`"observability"`, `"translation"`,
#'       `"scaling"` or `"polynomial"`).}
#'     \item{`identifiable`}{`TRUE`/`FALSE` for `"observability"`/`"translation"`
#'       (a full-rank verdict); `NA` for `"scaling"`/`"polynomial"`, which search
#'       non-exhaustively so "nothing found" is no proof of identifiability.}
#'     \item{`rank`, `dim`}{the observability-matrix rank and coordinate-space
#'       dimension (both `NA` for the scaling/polynomial engines).}
#'     \item{`symmetries`}{a list of the found generators/directions (empty when
#'       identifiable), each a differential generator `sum_i xi_i d/dz_i` in a
#'       canonical gauge. Fields: `generator` (components `xi_i`, keyed by
#'       coordinate); `weights` (a scaling's integer weights, else `NULL`); `type`
#'       (`"scaling"`, `"translation"`, `"affine"`, `"polynomial"` or `"general"`);
#'       `degree`; `support` (the coordinates involved, the only field set when no
#'       closed form was reached); `explicit`; `reason`; `certified`;
#'       `transformation` (the finite map, polynomial engine); `verified`.}
#'     \item{`info`}{`engine`, `lieOrderUsed`, `gapOrderUsed`, `conditions`,
#'       `segments`, the `settings` used, `elapsed` seconds and the `verification`
#'       guard.}
#'     \item{`call`}{the matched call.}
#'   }
#'   `print()` shows the verdict and generators; `summary()` adds the computation
#'   block, and `summary(x, verbose = TRUE)` the settings and guard detail.
#'
#' @details With events after the earliest one, `"observability"` splits the timeline
#'   at the event times into segments and propagates the state exactly across each
#'   inter-event gap (for generic timing), so a parameter entering only through a
#'   transient between events is identified -- place such an event at the earliest
#'   time. A `replace`/`add`/`multiply` event acts on the propagated state; a dose on
#'   a species eliminated by `reduceCQ` is not seen (keep it with `reduceCQ = FALSE`).
#'   The method internals are covered in `vignette("symmetryDetection")`.
#'
#' @note The interface, defaults and output structure may still change between
#'   releases.
#'
#' @references \[1\]
#' <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.012920>
#'
#' @examples
#' \dontrun{
#' # The canonical scaling symmetry: a reversible reaction observed only through
#' # alpha * A leaves the absolute scale free.
#' eq <- eqnlist() |>
#'   addReaction("A", "B", "k1 * A") |>
#'   addReaction("B", "A", "k2 * B")
#'
#' # Assigning still shows the result: symmetryDetection() prints its report on
#' # return. print(out) is terse (verdict + generators); summary(out) adds the
#' # computation block (engine, Lie order, saturation guard, timing).
#' out <- symmetryDetection(eq, eqnvec(Aobs = "alpha * A"))   # observability
#' print(out)
#' summary(out)
#' out <- symmetryDetection(eq, eqnvec(Aobs = "alpha * A"), method = "polynomial")
#' out <- symmetryDetection(eq, eqnvec(Aobs = "alpha * A"), method = "scaling")
#'
#' # A steady-state initial condition as an expression: a state-named trafo entry
#' # is the initial condition, seeded with its parameter sensitivities.
#' out <- symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
#'                          trafo = eqnvec(x = "b/a"), reconstruct = TRUE)
#'
#' # Several experimental conditions: a switch held at one generic value makes
#' # k1 and k2 look unidentifiable, but two values (set per condition by an
#' # event, read from the grid) identify both.
#' fu <- eqnvec(A = "-(k1 + u * k2) * A", u = "0")
#' events <- addEvent(eventlist(), var = "u", time = -1, value = "var_u",
#'                    method = "replace")
#' grid <- data.frame(var_u = c(0, 1), row.names = c("ctrl", "stim"))
#' out <- symmetryDetection(fu, eqnvec(y = "A"), method = "observability",
#'                          events = events, conditions = grid)
#' out$identifiable   # the object's fields are there for programmatic use
#'
#' # Pre-equilibrated model with a dose event: the steady state b/a is the
#' # relaxation attractor; a known dose makes (a, b, s) identifiable.
#' dose <- addEvent(eventlist(), var = "x", time = 0, value = "dose",
#'                  method = "replace")
#' out <- symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"),
#'                          method = "observability", events = dose,
#'                          conditions = data.frame(dose = 2, row.names = "stim"))
#' summary(out)
#'
#' # Full worked script (scaling, closed-form directions, an enzyme assay, a
#' # transcription-translation rate curve, conditions and steady states):
#' file.edit(system.file("examples", "symmetryDetection.R", package = "dMod"))
#' }
#' @export
symmetryDetection <- function(f = NULL, g = NULL, trafo = NULL,
                              method = c("observability", "polynomial", "scaling",
                                         "translation"),
                              parameters = NULL, fixed = NULL, forcings = NULL,
                              events = NULL, conditions = NULL,
                              equilibrate = FALSE,
                              reduceCQ = FALSE, freeInitial = NULL,
                              reconstruct = FALSE, verify = TRUE, cores = 1,
                              control = reconstControl(),
                              polynomial = polynomialControl(),
                              scaling = scalingControl(),
                              symEngine = c("modular", "symbolic"),
                              verbose = getOption("dMod.sym.verbose", TRUE)) {

  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required for symmetryDetection().")
  symEngine <- match.arg(symEngine)
  method <- match.arg(method)
  equilibrate <- isTRUE(equilibrate)
  # captured for the summary header (reproducibility) and to time the whole run
  .sym_call <- match.call()
  .sym_t0 <- Sys.time()
  # the engine settings surfaced in summary()'s computation report; summary()
  # shows only the entries relevant to the chosen method
  .sym_settings <- list(reduceCQ = isTRUE(reduceCQ), equilibrate = isTRUE(equilibrate),
                        reconstruct = isTRUE(reconstruct), verify = isTRUE(verify),
                        symEngine = symEngine, degreeCap = control$degreeCap,
                        certifyPoly = isTRUE(control$certifyPoly),
                        ansatz = polynomial$ansatz, pMax = polynomial$pMax,
                        polyBackend = if (is.null(polynomial$backend)) "symengine"
                                      else polynomial$backend)
  # every engine return funnels through here: normalise to the public object,
  # print the report unless verbose = FALSE (so `out <- symmetryDetection(...)`
  # still shows it), and return it invisibly to avoid a double print at top level.
  deliver <- function(raw, method) {
    res <- .sym_finalize(raw, method, .sym_settings, .sym_call,
                         elapsed = as.numeric(Sys.time() - .sym_t0, units = "secs"))
    if (isTRUE(verbose)) print(res)
    invisible(res)
  }

  # warn about arguments that do not apply to the chosen engine instead of
  # silently ignoring them
  supplied <- setdiff(names(match.call())[-1], "")
  applies <- switch(method,
    observability = c("events", "conditions", "equilibrate", "control"),
    translation   = c("events", "conditions", "equilibrate", "control"),
    polynomial = "polynomial",
    scaling = c("scaling", "events", "conditions"))
  methodSpecific <- c("events", "conditions", "equilibrate",
                      "control", "polynomial", "scaling")
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

  # freeInitial names the moiety species that keep a free resting initial value (the
  # held pivot) under the held-variable parameterisation. It only applies to
  # observability/translation with equilibrate = TRUE and reduceCQ = FALSE; elsewhere
  # every species already carries a free initial (reduceCQ = FALSE, no equilibrate) or
  # the moiety is reduced to a `total` (reduceCQ = TRUE), so it is warned and ignored.
  freeInitial <- freeInitial %||% character(0)
  if (length(freeInitial)) {
    if (length(bad <- setdiff(freeInitial, states)))
      warning("symmetryDetection(): freeInitial names non-state(s) ",
              paste(bad, collapse = ", "), "; ignored.", call. = FALSE)
    freeInitial <- intersect(freeInitial, states)
    if (length(freeInitial) && !(method %in% c("observability", "translation") &&
                                 equilibrate && !isTRUE(reduceCQ))) {
      warning("symmetryDetection(): freeInitial applies only to ",
              "method = \"observability\"/\"translation\" with equilibrate = TRUE and ",
              "reduceCQ = FALSE (the held-variable moiety parameterisation); ignored.",
              call. = FALSE)
      freeInitial <- character(0)
    }
  }

  # A `trafo` is either one eqnvec (applied to every condition) or a LIST of eqnvecs
  # (one per condition). In both cases a parameter-named entry is a substitution and
  # a state-named entry is an initial condition. A single trafo is substituted into
  # f and g up front; a per-condition list defers the substitution to each condition
  # (so a parameter may be baked to a different value or expression per condition,
  # exactly like a condition grid, and the two compose).
  initial <- NULL
  condSubs <- NULL          # per-condition parameter substitutions (trafo list)
  condInitial <- NULL       # per-condition initial conditions (trafo list)
  trafoSyms <- character(0) # symbols the trafo introduces (substitution targets/params)
  trafoList <- !is.null(trafo) && is.list(trafo) && !inherits(trafo, "eqnvec")
  if (trafoList) {
    trafos <- lapply(trafo, as.eqnvec)
    condSubs <- lapply(trafos, function(tr) {
      se <- setdiff(names(tr), states)
      as.list(setNames(as.character(tr[se]), se))
    })
    condInitial <- lapply(trafos, function(tr) {
      ic <- intersect(names(tr), states)
      if (length(ic)) as.eqnvec(tr[ic]) else NULL
    })
    trafoSyms <- unique(unlist(lapply(trafos,
                        function(tr) getSymbols(as.character(tr)))))
    if (method == "polynomial")
      warning("symmetryDetection(): a per-condition `trafo` list is not yet ",
              "supported for method = \"polynomial\"; pass a single trafo.", call. = FALSE)
  } else if (!is.null(trafo)) {
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
  # (single or per-condition) are dropped with a warning
  icNames <- unique(c(names(initial),
                      unlist(lapply(condInitial, names))))
  if (equilibrate && length(icNames)) {
    warning("symmetryDetection(): equilibrate solves the steady state from f = 0; ",
            "the initial condition(s) for ", paste(icNames, collapse = ", "),
            " in `trafo` are ignored.", call. = FALSE)
    initial <- NULL
    condInitial <- if (!is.null(condInitial))
      lapply(condInitial, function(x) NULL) else NULL
  }

  # A `trafo` initial condition for a species that participates in a conserved
  # quantity is incompatible with reduceCQ = TRUE: the CQ reduction expresses one
  # moiety species through its `total` and silently discards the supplied
  # steady-state relation for the eliminated species, so the model no longer starts
  # at the given steady state (it relaxes, over-reporting identifiability). Force
  # reduceCQ = FALSE with a warning so the explicit steady state is used as given.
  # (equilibrate has its own moiety handling and has already dropped trafo ICs above,
  # so it is exempt.)
  if (isTRUE(reduceCQ) && !equilibrate && !is.null(feqnlist) && length(icNames)) {
    moietyStates <- getSymbols(as.character(getTotals(feqnlist)))
    clash <- intersect(icNames, moietyStates)
    if (length(clash)) {
      warning("symmetryDetection(): a `trafo` initial condition was supplied for ",
              paste(clash, collapse = ", "), ", which participate(s) in a conserved ",
              "quantity; reduceCQ = TRUE would eliminate a moiety species and discard ",
              "that steady-state relation. Forcing reduceCQ = FALSE so the supplied ",
              "steady state is used as given.", call. = FALSE)
      reduceCQ <- FALSE
      .sym_settings$reduceCQ <- FALSE
    }
  }

  # conserved-quantity reduction (eqnlist input only)
  if (!is.null(feqnlist) && isTRUE(reduceCQ)) {
    totals <- getTotals(feqnlist)
    if (length(totals)) {
      # a moiety species under a free exponent must survive the reduction as a
      # bare symbol; keep it and eliminate another species of its total instead
      avoidCQ <- .sym_free_exponent_bases(as.character(c(fdyn, gobs)), names(fdyn))
      cq <- .detect_and_substitute_cq(totals, TRUE, fdyn, names(fdyn),
                                      parameters, expressInTotals = TRUE,
                                      avoid = avoidCQ)
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
        if (!is.null(condInitial))
          condInitial <- lapply(condInitial,
                                function(x) if (is.null(x)) NULL else recon(x))
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

  if (method %in% c("observability", "translation")) {
    # method = "translation" is the observability engine restricted to the exact
    # constant (additive) symmetry lattice: it peels the translation directions and
    # returns only those (the additive analogue of method = "scaling"). It always
    # uses the modular kernel (the peel needs it).
    translationsOnly <- (method == "translation")
    # shared condition/event resolution -- both the symbolic and modular engines use
    # the same per-condition substitutions and initial conditions, so it is computed
    # once here (the modular path below reuses `res`, `spy`, `constStates`, ...).
    # grid-substitution targets include symbols that appear only in initial
    # values, event values or a per-condition trafo, so they can be fixed too
    extraSyms <- c(if (!is.null(initial)) getSymbols(as.character(as.eqnvec(initial))),
                   if (!is.null(events)) getSymbols(as.character(as.data.frame(events)$value)),
                   if (length(condInitial)) unlist(lapply(condInitial,
                     function(x) if (is.null(x)) NULL else getSymbols(as.character(x)))),
                   trafoSyms)
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
    # held-variable conserved-moiety parameterisation (equilibrate without reduceCQ):
    # f = 0 is rank-deficient by one equation per moiety, so one pivot species per
    # moiety keeps its resting value free and the rest are solved. The pivots stay
    # dynamic states (all species remain coordinates); only the point solve freezes
    # them. Reuse the reduceCQ pivot choice, but never a dead (zero) or forced state,
    # or a species carrying a free exponent.
    heldStateParams <- character(0)   # named: pivot state -> initial-value parameter
    if (equilibrate && !isTRUE(reduceCQ) && !is.null(feqnlist)) {
      totalsFV <- getTotals(feqnlist)
      if (length(totalsFV)) {
        avoidFV <- unique(c(equilZeroStates, forcings,
          .sym_free_exponent_bases(as.character(c(fdyn, gobs)), names(fdyn))))
        decFV <- .cq_pivot_decomposition(totalsFV, states, parameters, avoid = avoidFV,
                                         prefer = freeInitial)
        piv <- intersect(decFV$pivots[!is.na(decFV$pivots)], states)
        piv <- setdiff(piv, c(equilZeroStates, forcings))
        if (length(piv)) {
          # dMod convention (as in deSolve): a parameter named like a state IS that
          # state's initial value. Hold each pivot's resting value under the pivot's
          # own name, so the moiety freedom is reported as that initial value, not a
          # `total`. (The tape adds these as parameter coordinates; the solver gets the
          # pivot value through the held-state channel, not as a parameter.)
          heldStateParams <- setNames(piv, piv)
        }
      }
    }
    res <- .sym_resolve_conditions(conditions, events, initial, symbols, states,
                                   constStates, forcings, equilibrate = equilibrate,
                                   condSubs = condSubs, condInitial = condInitial)
    spy <- tryCatch(reticulate::import("sympy", convert = TRUE),
                    error = function(err) NULL)

    # pure-symbolic cross-check engine: builds the observability-identifiability matrix
    # d/dz[g, L_f g, ...] and reduces its rank/nullspace with sympy over the exact
    # rational-function field -- no finite fields, no power/Hill recast (base^exp stays a
    # symbolic atom), and the Lie order is carried to the exact saturation bound so no
    # verify guard is needed. Independent from the modular kernel, hence a strong
    # cross-check, but only for SMALL models. Handles multiple conditions and
    # single-time events by stacking each condition's observability rows over one
    # shared coordinate space (the intersection of the per-condition codistributions);
    # later events (gaps) and equilibrate stay a modular feature.
    if (symEngine == "symbolic" && !translationsOnly) {
      if (isTRUE(equilibrate))
        stop("symEngine = \"symbolic\" does not support equilibrate; supply the ",
             "steady state explicitly through `trafo` (e.g. from steadyStates()), ",
             "which is an exact substitution, or use symEngine = \"modular\".",
             call. = FALSE)
      if (res$nGaps > 0L)
        stop("symEngine = \"symbolic\" handles single-segment conditions only; ",
             "later events (gaps) need symEngine = \"modular\".", call. = FALSE)
      sr <- sd$observabilitySympyMulti(
        model = toLines(fdyn), observation = toLines(gobs),
        conditionSubs = res$subs, conditionIC0 = res$ic0,
        fixed = if (length(fixed)) fixed else NULL,
        parameters = if (length(parameters)) parameters else NULL,
        inputs = if (length(forcings)) forcings else NULL)
      if (!isTRUE(sr$ok))
        stop("symEngine = \"symbolic\": ",
             if (!is.null(sr$why)) sr$why else "could not build the symbolic system.",
             call. = FALSE)
      sr$method <- "observability"; sr$engine <- "symbolic"
      sr$lieOrderUsed <- as.integer(sr$lieOrder)
      sr$conditions <- as.integer(res$nConditions)
      sr$segments <- as.integer(res$nConditions)   # single-segment: one per condition
      sr$gapOrderUsed <- 0L
      sr$nonIdentifiable <- .sym_relabel_directions(sr$nonIdentifiable, sd,
        if (is.null(control$degreeCap)) 4L else control$degreeCap)
      if (isTRUE(control$certifyPoly))
        sr$nonIdentifiable <- .sym_certify_poly(sr$nonIdentifiable,
          toLines(fdyn), toLines(gobs), forcings, fixed, parameters, control, sd)
      return(deliver(sr, method))
    }
    # equilibrate always uses the implicit determining system: the states stay
    # coordinates and the steady state enters as the f=0 tangency constraint df.xi=0
    # (never forming or eliminating a symbolic x*), with each condition carrying its
    # own resting state, so non-scaling multi-condition directions and the recast
    # (free Hill/power exponent) case are both exact and closed-form. An explicit
    # steady state is instead supplied through `trafo` (from steadyStates()).
    useImplicit <- isTRUE(equilibrate)
    runObs <- function(ui) {
      multi <- sd$compileObservabilityTapeMulti(
        model = toLines(fdyn), observation = toLines(gobs),
        conditionSubs = res$subs, conditionIC0 = res$ic0,
        fixed = if (length(fixed)) fixed else NULL,
        parameters = if (length(parameters)) parameters else NULL,
        equilibrate = equilibrate, segEquilibrate = as.list(res$segEquil),
        forcings = if (length(forcings)) forcings else NULL,
        conditionEvents = res$segEvents, conditionT0Events = res$events0,
        jointSteadyState = isTRUE(ui),
        jointFixedStates = if (isTRUE(ui) && length(equilZeroStates))
          equilZeroStates else NULL,
        heldStateParams = if (isTRUE(ui) && length(heldStateParams))
          as.list(heldStateParams) else NULL)
      if (!isTRUE(multi$ok)) return(list(ok = FALSE, nonrational = multi$nonrational))
      list(ok = TRUE, result = .observability_analytic_multi(multi, spy = spy,
             closedForm = reconstruct, sd = sd, cores = cores,
             equilZeroStates = equilZeroStates, t0events = res$events0,
             nConditions = res$nConditions, chainOf = res$chainOf,
             nGaps = res$nGaps, implicitSteadyState = isTRUE(ui), control = control,
             verify = verify, translationsOnly = translationsOnly))
    }
    ro <- runObs(useImplicit)
    if (useImplicit && !isFALSE(ro$ok) && is.null(ro$result))
      stop("symmetryDetection(): the implicit steady-state path could not be ",
           "evaluated (e.g. a singular resting Jacobian from an unreduced conserved ",
           "moiety). Reduce conserved moieties (reduceCQ = TRUE) or supply an ",
           "explicit steady state through `trafo` (from steadyStates()).",
           call. = FALSE)
    if (isFALSE(ro$ok))
      stop("method = \"observability\" requires rational right-hand sides, ",
           "observables and initial conditions (built from +, -, *, / and ",
           "integer powers).\n  ",
           paste(unlist(ro$nonrational), collapse = "\n  "),
           "\nA logarithmic observable log10(h) + offset equals the rational ",
           "observable scale * h; supply it in that form, or use ",
           "method = \"polynomial\".", call. = FALSE)
    res <- ro$result
    if (is.list(res)) {
      res$nonIdentifiable <- .sym_relabel_directions(res$nonIdentifiable, sd,
        if (is.null(control$degreeCap)) 4L else control$degreeCap)
      if (!translationsOnly && isTRUE(control$certifyPoly))
        res$nonIdentifiable <- .sym_certify_poly(res$nonIdentifiable,
          toLines(fdyn), toLines(gobs), forcings, fixed, parameters, control, sd)
      if (translationsOnly) {
        # keep only the translation class (the peel returns just these; a fallback
        # full reconstruction may also carry other classes)
        res$nonIdentifiable <- Filter(function(d) isTRUE(d$type == "translation"),
                                      res$nonIdentifiable)
        res$method <- "translation"
      }
    }
    if (isTRUE(verify) && is.list(res) && is.list(res$verification) &&
        isFALSE(res$verification$ok))
      warning("symmetryDetection(verify = TRUE): the Schwartz-Zippel saturation guard ",
              "found the rank still growing past the reported Lie order -- the ",
              "directions may be over-reported; inspect $verification (",
              res$verification$reason, ").", call. = FALSE)
    return(deliver(res, if (!is.null(res$method)) res$method else method))
  }

  # scaling: exact integer-kernel engine. A known-value dose (replace/add) pins the
  # dosed state to an absolute value, so it cannot scale and its weight is forced to
  # 0 (fixed); a condition grid / per-condition trafo list intersects the
  # per-condition scaling lattices. Steady states are a no-op for scaling (a scaling
  # of f leaves f = 0 invariant), so `equilibrate` stays observability-only.
  if (method == "scaling") {
    fixedScal <- unique(c(fixed, .sym_event_pinned_states(events, conditions)))
    syms <- unique(c(states, names(gobs), getSymbols(as.character(c(fdyn, gobs)))))
    multiCond <- (!is.null(conditions) && nrow(as.data.frame(conditions)) > 1L) ||
                 length(condSubs) > 1L
    reticulate::py_capture_output(
      res <- if (multiCond) {
        pc <- .sym_percond_lines(fdyn, gobs, conditions, condSubs, syms)
        sd$scalingSymmetriesMulti(
          perCondModel = lapply(pc, `[[`, "f"),
          perCondObs   = lapply(pc, `[[`, "g"),
          inputs = if (length(forcings)) forcings else NULL,
          fixed  = if (length(fixedScal)) fixedScal else NULL)
      } else {
        sd$symmetryDetectiondMod(
          model = toLines(fdyn), observation = toLines(gobs),
          inputs = if (length(forcings)) forcings else NULL,
          fixed = fixedScal, backend = scaling$backend,
          parameters = if (length(parameters)) parameters else NULL, method = "scaling")
      })
    return(deliver(res, "scaling"))
  }

  # polynomial: the polynomial Lie-symmetry ansatz (Merkt et al. 2015) on the
  # (trafo-substituted) f and g; forcings are the externally driven (input) states.
  # (This engine was named "liesym" before; the Python method tag stays "liesym".)
  fld <- function(nm, default) if (is.null(polynomial[[nm]])) default else polynomial[[nm]]
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
      allTrafos   = fld("allTrafos", FALSE),
      lieOrder    = as.integer(fld("lieOrder", 0L)),
      exact       = fld("exact", TRUE),
      verify      = fld("verify", TRUE),
      backend     = if (is.null(polynomial$backend)) "symengine" else polynomial$backend,
      parameters  = if (length(parameters)) parameters else NULL,
      method      = "liesym"
    ))
  deliver(res, "polynomial")
}


# ---- Schwartz-Zippel saturation guard (verify = TRUE) ---------------------------
# A fast, single-point cross-check of the ONE thing the observability heuristic can get
# wrong silently: stopping the Lie order too early. The plateau rule declares the rank
# saturated once it stops growing for a few orders, but a rank can plateau then grow again
# (a Hill exponent only observable through a high-order derivative), so a premature stop
# over-reports non-identifiability. The guard re-uses the kernel that the main analysis
# already built and its base point -- whose modular steady-state solve is cached at the
# saturation prime -- and simply extends the Lie order past NtUsed. Only the Lie jet grows;
# the (expensive) f = 0 solve is not repeated, so this is a handful of matrix ranks, not a
# second analysis. If the rank climbs past the reported value the saturation was premature.
# By Schwartz-Zippel the modular rank at the point equals the generic rank almost surely.
# Directions themselves are already certified at a fresh prime during reconstruction and
# the scalings are exact, so the residual risk the guard closes is exactly the Lie stop.
.sym_sz_saturation_guard <- function(kcall, point0Solved, NtUsed, reportedRank,
                                     margin = as.integer(Sys.getenv("DMOD_SYM_VERIFY_MARGIN", "6"))) {
  P <- .symPrimes[1]                          # saturation prime: point0Solved's solve is cached
  r0 <- tryCatch(kcall(point0Solved, P, as.integer(NtUsed)), error = function(e) NULL)
  if (is.null(r0) || !isTRUE(r0$ok))
    return(list(ok = NA, method = "saturation guard",
                reason = "base point not re-evaluable"))
  base <- as.integer(r0$rank); maxR <- base; growAt <- NA_integer_
  for (k in seq_len(max(1L, margin))) {       # extend the Lie order; only the jet grows
    rk <- tryCatch(kcall(point0Solved, P, as.integer(NtUsed + k)), error = function(e) NULL)
    if (is.null(rk) || !isTRUE(rk$ok)) next
    if (as.integer(rk$rank) > maxR) maxR <- as.integer(rk$rank)
    if (as.integer(rk$rank) > base) { growAt <- as.integer(NtUsed + k); break }
  }
  ok <- is.na(growAt)
  list(ok = ok, method = "saturation guard", lieOrderUsed = as.integer(NtUsed),
       ordersChecked = as.integer(NtUsed + max(1L, margin)),
       kernelRank = base, kernelRankExtended = maxR, growAt = growAt,
       reason = if (ok)
         sprintf("rank %d stable through Lie order %d (%d orders beyond the reported saturation)",
                 base, NtUsed + max(1L, margin), max(1L, margin))
       else sprintf("rank grows %d -> %d at Lie order %d (the reported Lie order was premature)",
                    base, maxR, growAt))}


# States pinned to a known absolute value by an event (a numeric replace/add dose):
# such a state cannot scale, so for method = "scaling" its weight is fixed to 0. A
# multiply dose, or a dose to a parameter (symbolic) value, imposes no constraint. A
# grid-column value is resolved to its cells; a state pinned in ANY condition is
# pinned (the scaling common to all conditions must respect every condition).
.sym_event_pinned_states <- function(events, conditions) {
  if (is.null(events)) return(character(0))
  ev <- as.data.frame(events, stringsAsFactors = FALSE)
  if (!nrow(ev)) return(character(0))
  grid <- if (is.null(conditions)) NULL else as.data.frame(conditions, stringsAsFactors = FALSE)
  cols <- if (is.null(grid)) character(0) else colnames(grid)
  isKnown <- function(v) {
    v <- as.character(v)
    vals <- if (v %in% cols) as.character(unlist(grid[[v]])) else v
    length(vals) > 0L && all(!is.na(suppressWarnings(as.numeric(vals))))
  }
  pinned <- character(0)
  for (i in seq_len(nrow(ev)))
    if (as.character(ev$method[i]) %in% c("replace", "add") && isKnown(ev$value[i]))
      pinned <- c(pinned, as.character(ev$var[i]))
  unique(pinned)
}


# Per-condition (f, g) line lists for the scaling engine: each condition applies its
# own parameter substitutions (grid cells for symbol columns, plus a per-condition
# trafo) to f and g. The states are unchanged, so all conditions share them (as the
# multi-condition scaling kernel requires); only parameters are baked/renamed.
.sym_percond_lines <- function(fdyn, gobs, conditions, condSubs, symbols) {
  grid <- if (is.null(conditions)) NULL else as.data.frame(conditions, stringsAsFactors = FALSE)
  nGrid <- if (is.null(grid)) 0L else nrow(grid)
  K <- max(length(condSubs), nGrid, 1L)
  subCols <- if (is.null(grid)) character(0) else intersect(colnames(grid), symbols)
  lines <- function(e) { e <- as.eqnvec(e)
    as.character(vapply(seq_along(e),
                        function(i) paste(names(e)[i], "=", e[i]), character(1))) }
  lapply(seq_len(K), function(k) {
    subs <- list()
    for (col in subCols) subs[[col]] <- as.character(grid[k, col])
    if (length(condSubs) && !is.null(condSubs[[k]]))
      for (nm in names(condSubs[[k]])) subs[[nm]] <- condSubs[[k]][[nm]]
    sub <- function(e) { e <- as.eqnvec(e)
      if (!length(subs)) e else
        setNames(replaceSymbols(names(subs), unlist(subs), e), names(e)) }
    list(f = lines(sub(fdyn)), g = lines(sub(gobs)))
  })
}


# state species that appear as the base of a free (non-numeric) exponent, e.g.
# `C3` in `C3^nhill`. Such a species must stay a bare symbol through the
# conserved-quantity reduction (eliminating it by subtraction would put a sum
# under that exponent and break rationality), so it is kept out of the pivot set.
.sym_free_exponent_bases <- function(exprs, states) {
  if (!length(exprs) || !length(states)) return(character(0))
  hits <- regmatches(exprs, gregexpr(
    "[A-Za-z_][A-Za-z0-9_]*\\s*\\^\\s*(?![0-9])", exprs, perl = TRUE))
  bases <- sub("\\s*\\^.*$", "", unlist(hits, use.names = FALSE))
  intersect(unique(bases), states)
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

# TRUE once the closed-form reconstruction has run past its wall-clock budget
# (reconstControl(timeout=)). The deadline is a POSIXct stamped on the live control
# list at the start of the reconstruction, or NULL for an unbounded run.
.sym_expired <- function(ctrl) {
  d <- ctrl$deadline
  !is.null(d) && Sys.time() > d
}

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
#'   state exactly across later (post-stimulus) event boundaries.
#' @param minsupportCandCap Cap on the number of column subsets the minimal-support
#'   search enumerates when hunting narrow nullspace cocircuits. A genuinely wide
#'   direction has no small-support cocircuit, so an exhaustive scan would test
#'   `choose(support, s)` subsets fruitlessly; the search stops at this cap and lets
#'   any wide direction fall through to the free-column fit. Raise it only to peel
#'   cocircuits of unusually large support.
#' @param perprimeCap Max sample points per prime for the per-prime reconstruction
#'   used by the equilibrate/joint (coupled steady-state) path. There a random
#'   multi-parameter perturbation admits an interior modular steady state at only
#'   some primes, so a point that solves at every prime at once is vanishingly rare
#'   and the shared all-prime bank cannot fill; instead each prime collects its own
#'   solvable, pivot-consistent points and the per-prime entry fits are lifted by
#'   Chinese remaindering. It doubles as the width gate: a direction whose widest
#'   entry needs more than `perprimeCap` points already at degree 2 is too wide for
#'   the dense per-prime fit and is returned support-only immediately (no sampling),
#'   so an intrinsically wide/transcendental joint confound is reported by support
#'   quickly rather than grinding. The default keeps the coupled reconstruction
#'   bounded; raise it (with `timeout`) to attempt a genuinely wide joint direction.
#' @param perprimeMinPrimes Minimum number of primes that must fill for the
#'   per-prime path to lift a coefficient by CRT. A prime whose coupled solve never
#'   succeeds under the direction's perturbations is dropped; the rational
#'   reconstruction then runs over the remaining primes (their product still bounds
#'   the coefficient height). Below this many live primes the entry is support-only.
#' @param certifyPoly Logical. For each `"affine"`/`"polynomial"` direction, run the
#'   polynomial Lie-symmetry engine (`method = "polynomial"`) and set `$certified`
#'   to whether the direction is a nonzero constant combination of its generators --
#'   i.e. an exact polynomial Lie point symmetry, not merely an observability
#'   non-identifiability. `FALSE` by default; it is comparatively expensive.
#' @param certifyPolyDeg Ansatz degree for `certifyPoly`; `NULL` (default) uses each
#'   direction's own classified degree.
#' @param timeout Wall-clock budget in seconds for the closed-form reconstruction
#'   (`reconstruct = TRUE`). When it is exceeded, the directions still being
#'   reconstructed are returned support-only (`reconstruct = FALSE`) with a
#'   `reason`, so a hard general (non-scaling) direction aborts cleanly instead of
#'   running unbounded. `Inf` (the default) imposes no limit and reproduces the
#'   previous behaviour exactly.
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
                           gapOrderCap        = 8L,
                           minsupportCandCap  = 20000L,
                           perprimeCap        = 120L,
                           perprimeMinPrimes  = 3L,
                           certifyPoly        = FALSE,
                           certifyPolyDeg     = NULL,
                           timeout            = Inf) {
  stopifnot(relevanceCap >= 0L, relevanceCapDir >= 1L,
            relevanceCapSparse >= relevanceCap, degreeCap >= 0L,
            sampleSlack >= 0L, probeRetries >= 1L, termCap >= 1L,
            laurentDegNum >= 1L, laurentDegDen >= 0L, laurentCandCap >= 1L,
            generalDegNum >= 1L, generalDegDen >= 1L, gapOrderCap >= 0L,
            minsupportCandCap >= 1L, perprimeCap >= 1L, perprimeMinPrimes >= 2L,
            is.logical(certifyPoly), is.null(certifyPolyDeg) || certifyPolyDeg >= 1L,
            is.numeric(timeout), timeout > 0)
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
                 gapOrderCap        = as.integer(gapOrderCap),
                 minsupportCandCap  = as.integer(minsupportCandCap),
                 perprimeCap        = as.integer(perprimeCap),
                 perprimeMinPrimes  = as.integer(perprimeMinPrimes),
                 certifyPoly        = isTRUE(certifyPoly),
                 certifyPolyDeg     = if (is.null(certifyPolyDeg)) NULL
                                      else as.integer(certifyPolyDeg),
                 timeout            = timeout),
            class = c("reconstControl", "list"))
}


#' Tuning control for `symmetryDetection(method = "polynomial")`
#'
#' Bundles the polynomial Lie-symmetry engine's parameters into one object, passed
#' as `symmetryDetection(..., polynomial = polynomialControl(...))`.
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
#' @return A `polynomialControl` list.
#' @seealso [symmetryDetection()], [scalingControl()]
#' @export
polynomialControl <- function(ansatz    = c("uni", "par", "multi"),
                              pMax      = 2L,
                              lieOrder  = 0L,
                              exact     = TRUE,
                              verify    = TRUE,
                              allTrafos = FALSE,
                              backend   = c("symengine", "sympy")) {
  ansatz  <- match.arg(ansatz)
  backend <- match.arg(backend)
  stopifnot(pMax >= 1L, lieOrder >= 0L,
            is.logical(exact), is.logical(verify), is.logical(allTrafos))
  structure(list(ansatz = ansatz, pMax = as.integer(pMax),
                 lieOrder = as.integer(lieOrder), exact = isTRUE(exact),
                 verify = isTRUE(verify), allTrafos = isTRUE(allTrafos),
                 backend = backend),
            class = c("polynomialControl", "list"))
}


#' Tuning control for `symmetryDetection(method = "scaling")`
#'
#' Bundles the scaling (toric) symmetry engine's parameters into one object, passed
#' as `symmetryDetection(..., scaling = scalingControl(...))`.
#'
#' @param backend `"symengine"` (faster, falls back to sympy) or `"sympy"`.
#' @return A `scalingControl` list.
#' @seealso [symmetryDetection()], [polynomialControl()]
#' @export
scalingControl <- function(backend = c("symengine", "sympy")) {
  backend <- match.arg(backend)
  structure(list(backend = backend),
            class = c("scalingControl", "list"))
}


# ---- modular linear algebra over GF(p) -----------------------------------------------

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


# TRUE iff vector x lies in the column span of the nz-row matrix M over GF(p).
.sym_in_span <- function(M, x, nz, p)
  ncol(M) > 0L && !is.null(symSolveMod(matrix(as.integer(M), nz), as.integer(x), p))


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

# RREF of the rows of B, but with the columns in `physCols0` (0-based) taken as
# pivots first, so distinct residual directions receive distinct PHYSICAL pivots
# (the forward path anchors on them). Returns the reduced rows in the original
# column order and each row's 0-based pivot column.
.sym_phys_rref <- function(B, physCols0, nz, p) {
  if (nrow(B) == 0L) return(list(R = B, piv = integer(0)))
  ord <- c(sort(as.integer(physCols0)),
           sort(setdiff(0:(nz - 1L), as.integer(physCols0))))   # physical columns first
  rr <- .sym_rref_modp(B[, ord + 1L, drop = FALSE], p)
  Rback <- matrix(0, nrow(rr$R), nz)
  Rback[, ord + 1L] <- rr$R
  list(R = Rback, piv = ord[rr$piv + 1L])                        # 0-based original pivots
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
# Reduce rows R modulo the RREF of the loop-invariant scaling-lattice matrix S,
# memoised per prime (the RREF is recomputed once per prime, not per sample point).
.sym_reduce_mod_rows <- function(S) {
  cache <- new.env(parent = emptyenv())
  function(R, p) {
    if (nrow(S) == 0L) return(R %% p)
    key <- as.character(p)
    Sr <- cache[[key]]
    if (is.null(Sr)) { Sr <- .sym_rref_modp(S, p); cache[[key]] <- Sr }
    for (j in seq_along(Sr$piv)) {
      fac <- R[, Sr$piv[j] + 1L] %% p
      for (i in seq_len(nrow(R)))
        if (fac[i] != 0) R[i, ] <- (R[i, ] - .sym_mulmod(Sr$R[j, ], fac[i], p)) %% p
    }
    R
  }
}

# Shared tail of the canon/logcoord gauges: `rowsFn(rp, p, zvals)` maps a kernel
# result to the k reduced direction rows (or NULL). Pin the k pivot columns of the
# reference rows to the identity, so each direction's entry reads off directly.
.sym_gauge_from_rows <- function(residualFree, rowsFn, refRows, P, k, nz) {
  raw <- list(anchors = residualFree, residueFns = vector("list", k))
  if (is.null(refRows)) return(raw)
  rr <- .sym_rref_modp(refRows, P)
  if (rr$rank < k) return(raw)
  dp <- rr$piv
  makeFn <- function(i) function(rp, p, zvals = NULL) {
    R <- rowsFn(rp, p, zvals)
    if (is.null(R)) return(NULL)
    Binv <- .sym_matinv_modp(R[, dp + 1L, drop = FALSE], p)
    if (is.null(Binv)) return(NULL)
    w <- numeric(nz)
    for (l in seq_len(k)) w <- (w + .sym_mulmod(R[l, ], Binv[i, l], p)) %% p
    as.integer(w)
  }
  list(anchors = dp, residueFns = lapply(seq_len(k), makeFn))
}


.sym_canon_gauge <- function(residualFree, scalRows, P, nz, sc) {
  k <- length(residualFree)
  if (k <= 1L) return(list(anchors = residualFree, residueFns = vector("list", k)))
  reduceScal <- .sym_reduce_mod_rows(scalRows)
  resid <- function(rp, p, zvals = NULL) reduceScal(
    t(vapply(residualFree, function(fc) .sym_null_residues(rp, fc, p),
             integer(nz))), p)
  .sym_gauge_from_rows(residualFree, resid, resid(sc$ref, P), P, k, nz)
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
  if (k == 0L) return(list(anchors = residualFree, residueFns = vector("list", k)))

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
  reduceScal <- .sym_reduce_mod_rows(W)

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
  .sym_gauge_from_rows(residualFree, etaRows, etaRows(sc$ref, P, zvals0), P, k, nz)
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


# Decouple the residual directions by their MINIMAL SUPPORT. A weighted scaling
# whose weight is a free parameter (e.g. a Hill-exponent feedback, xi_kinh =
# -nhill * kinh) is a sparse circuit of the full nullspace, but .sym_canon_gauge
# reduces modulo the integer scaling lattice and so lifts it to a dense
# representative spread over the whole feedback loop. Here the scalings are kept
# in the gauge freedom instead of quotiented out: each residual direction is the
# unique nullspace vector supported on a minimal column set (a cocircuit of the
# nullspace matroid), which for a genuine sparse symmetry is exactly the physical
# direction. Found by enumerating small column subsets and testing whether the
# row space carries a vector vanishing off them; scaling cocircuits (already
# reported by the peel) are filtered out. Mirrors .sym_canon_gauge: returns one
# (anchor, residue function) pair per residual direction, or the raw gauge.
.sym_minsupport_gauge <- function(residualFree, scalRows, P, nz, sc, freeCols,
                                  supportCap = 6L, candCap = 20000L, maxSecs = 20) {
  k <- length(residualFree)
  raw <- list(anchors = residualFree, residueFns = vector("list", k))
  if (k == 0L) return(raw)
  deadline <- Sys.time() + maxSecs

  basisRows <- function(rp, p)
    t(vapply(freeCols, function(fc) .sym_null_residues(rp, fc, p), integer(nz)))
  rowComb <- function(c, B, p) {
    w <- numeric(nz)
    for (l in seq_along(c)) if (c[l] %% p != 0)
      w <- (w + .sym_mulmod(B[l, ], c[l], p)) %% p
    as.integer(w)
  }
  # left null vector of the column submatrix `keep` of B (rows = nullspace dim):
  # the combination of rows that vanishes on those columns, i.e. a row-space
  # vector supported off them. NULL when the kept columns have full row rank.
  leftNull <- function(B, keep, p) {
    nr <- nrow(B)
    M <- B[, keep, drop = FALSE]
    rr <- .sym_rref_modp(cbind(matrix(as.numeric(M) %% p, nr), diag(nr)), p)
    z <- which(rowSums(rr$R[, seq_len(ncol(M)), drop = FALSE] %% p != 0) == 0)
    if (!length(z)) return(NULL)
    rr$R[z[1], ncol(M) + seq_len(nr)] %% p
  }
  # the row-space vector supported within columns S (1-based), or NULL
  supported <- function(B, S, p) {
    c <- leftNull(B, setdiff(seq_len(nz), S), p)
    if (is.null(c)) return(NULL)
    rowComb(c, B, p)
  }
  isScaling <- function(v) .sym_in_span(t(scalRows), v, nz, P)

  B0 <- basisRows(sc$ref, P)
  # candidate columns: where any residual direction's free-column residue lives
  cols <- sort(unique(unlist(lapply(residualFree,
    function(fc) which(B0[match(fc, freeCols), ] %% P != 0)))))
  if (length(cols) < 2L) return(raw)

  # A genuinely wide direction has no small-support cocircuit, so exhaustively
  # scanning the union columns would test choose(|cols|, s) subsets (millions)
  # fruitlessly. Iterate the subsets in place (no combn materialisation) under a
  # global budget; a narrow cocircuit is found in the first few sizes, and once
  # the budget is spent the wide directions fall through to the free-column fit.
  nCols <- length(cols)
  budget <- as.integer(candCap)
  nextCombo <- function(idx, n, s) {
    i <- s
    while (i >= 1L && idx[i] == n - s + i) i <- i - 1L
    if (i < 1L) return(NULL)
    idx[i] <- idx[i] + 1L
    j <- i + 1L
    while (j <= s) { idx[j] <- idx[j - 1L] + 1L; j <- j + 1L }
    idx
  }
  found <- list(); sel <- matrix(0L, nz, 0L); iter <- 0L
  for (s in 2:min(nCols, supportCap)) {
    if (length(found) >= k || budget <= 0L) break
    idx <- seq_len(s)
    repeat {
      if (length(found) >= k || budget <= 0L) break
      # a wide direction has no small cocircuit, so the scan is fruitless; cap it by
      # wall clock too (the candidate budget alone can still be minutes on a wide
      # residual support). A genuine narrow cocircuit is found in the first sizes.
      iter <- iter + 1L
      if (iter %% 256L == 0L && Sys.time() > deadline) { budget <- 0L; break }
      S <- cols[idx]
      idx <- nextCombo(idx, nCols, s)
      budget <- budget - 1L
      spanned <- any(vapply(found, function(fc) all(fc$supp %in% S), logical(1)))
      if (!spanned) {
        v <- supported(B0, S, P)
        if (!is.null(v) && !isScaling(v)) {
          base <- cbind(if (nrow(scalRows)) t(scalRows) else matrix(0L, nz, 0L), sel)
          if (!.sym_in_span(base, v, nz, P)) {
            supp <- which(v %% P != 0)
            found[[length(found) + 1L]] <- list(supp = supp, anchor = supp[1] - 1L, S = S)
            sel <- cbind(sel, v)
          }
        }
      }
      if (is.null(idx)) break
    }
  }
  # partial gauge: return whatever minimal-support cocircuits were found (possibly
  # fewer than k). The caller reconstructs these cheaply and sends only the
  # directions no narrow cocircuit spans to the wide free-column fit.
  if (!length(found)) return(raw)

  fns <- lapply(found, function(fc) {
    S <- fc$S; anchor <- fc$anchor
    # log-coordinate residue eta_c = xi_c / z_c, normalised to 1 on the anchor. A
    # weighted scaling whose weight is a parameter (the Hill exponent) then has a
    # constant-leaf entry (eta = -nhill) instead of the rational xi/xi that couples
    # the feedback species' whole steady state, so the entry relevance collapses to
    # the weight's own leaf and the reconstruction needs no loop-coupled samples.
    fn <- function(rp, p, zvals = NULL) {
      v <- supported(basisRows(rp, p), S, p)
      if (is.null(v) || v[anchor + 1L] %% p == 0) return(NULL)
      v <- .sym_mulmod(v, .sym_invmod(v[anchor + 1L] %% p, p), p)
      if (is.null(zvals)) return(as.integer(v))
      za <- as.numeric(zvals[anchor + 1L]) %% p
      if (za == 0) return(NULL)
      eta <- integer(nz)
      for (c in S) {
        zc <- as.numeric(zvals[c]) %% p
        if (zc == 0) { if (v[c] %% p != 0) return(NULL) else next }
        eta[c] <- as.integer(.sym_mulmod(v[c], .sym_mulmod(za, .sym_invmod(zc, p), p), p))
      }
      eta
    }
    # support fixed to S, so a pivot-shifted probe (residue NULL) only marks that
    # sample unusable, not the leaf relevant: relevance is read optimistically and
    # the reconstructed form is certified at a fresh point
    attr(fn, "pinnedSupport") <- TRUE
    fn
  })
  list(anchors = vapply(found, function(fc) fc$anchor, integer(1)),
       residueFns = fns, vectors = sel)
}


# ---- symbolic reconstruction: monomials, per-prime lift, recast back-substitution ----

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


# FORWARD-sampling reconstruction of a joint residual direction whose entries depend on the
# resting state x*(theta) (the SMAD pool/exponent direction [9]). The backward per-prime path
# fails there because each prime's steady-state slice differs (no single base solves at every
# prime) so the coefficient CRT is inconsistent. Here f = 0 is solved by the LINEAR forward
# solve (choose the resting states + free params, solve a turnover-rate subset), which is valid
# at EVERY prime -- so the resting states become FREE, SHARED sample coordinates and the same
# points fill one bank across all primes (consistent CRT). The direction is reconstructed as a
# rational in (theta, resting states, log base, exponent); the resting states stay SYMBOLS in
# the report (substitute x* = g(theta) via the trafo downstream). Coefficients lift with the
# bignum CRT (no 4-prime cap). Returns a verified closed-form entry or a support-only fallback.
.sym_perprime_forward <- function(f, sc, kcall, kcallFwd, znames, zSlots, leafNames, nz,
                                  scaling, stateColNames, paramNames, recast, sd, spy, ctrl,
                                  physCols, models, realStateNames, solveParamNames, solveHeld,
                                  acIn = NULL, listAnchors = FALSE) {
  # three 31-bit primes give a ~2^93 CRT modulus -- ample for the (clean, small) coefficients of a
  # physical symmetry, while each extra prime is another full per-point kernel pass, so keep it low.
  rp8 <- unique(as.integer(c(.symPrimes, 2147483563, 2147483549, 2147483543, 2147483497)))
  nPrimeFwd <- max(2L, min(3L, length(rp8) - 1L))
  primes <- rp8[seq_len(nPrimeFwd)]; qv <- rp8[length(primes) + 1L]; P1 <- primes[1]
  dirSupport <- NULL   # this direction's support, set once its anchor is known
  fb <- function(reason) {
    supp <- if (!is.null(dirSupport)) dirSupport
            else { v <- .sym_null_residues(sc$ref, f, P1); sort(znames[v != 0]) }
    list(support = supp, type = "general", closedForm = FALSE, reason = reason) }
  if (is.null(sd) || is.null(kcallFwd) || .sym_expired(ctrl)) {
    if (listAnchors) return(integer(0)); return(fb("forward path unavailable")) }
  rdiag <- nzchar(Sys.getenv("DMOD_SYM_ROBUSTDIAG"))
  stripName <- function(nm) sub("\\|c[0-9]+$", "", nm)
  npt <- length(sc$point0); poolN <- sc$poolNext
  draw <- function(n) { u <- sc$pool(poolN + seq_len(n) - 1L); poolN <<- poolN + n; u }
  freeColsOf <- function(rp) setdiff(0:(nz - 1L), as.integer(rp$pivots))
  isStateCol <- logical(nz); isStateCol[which(znames %in% stateColNames)] <- TRUE
  physSet <- if (is.null(physCols)) (seq_len(nz) - 1L) else as.integer(physCols)
  nS <- length(scaling)
  Wn <- lapply(scaling, function(s) names(s$vector)); Wv <- lapply(scaling, function(s) as.character(unlist(s$vector)))
  tangentsAt <- function(zvals, p) {
    env <- as.list(setNames(as.numeric(zvals), znames)); M <- matrix(0L, nS, nz); drp <- logical(nS)
    for (j in seq_len(nS)) { cols <- match(Wn[[j]], znames)
      for (t in seq_along(cols)) { cc <- cols[t]; if (is.na(cc)) next
        w <- suppressWarnings(as.numeric(Wv[[j]][t]))
        if (is.na(w)) { wv <- .sym_eval_modq(Wv[[j]][t], env, p, sd)
          if (is.null(wv) || is.na(wv)) { drp[j] <- TRUE; break }; w <- as.numeric(wv) }
        zc <- if (isStateCol[cc]) 1 else (as.numeric(zvals[cc]) %% p)
        M[j, cc] <- as.integer((w %% p * zc) %% p) } }
    M[!drp, , drop = FALSE] }
  reduceRows <- function(B, S, p) { if (nrow(S) == 0L) return(B %% p)
    Sr <- .sym_rref_modp(S, p)
    for (j in seq_along(Sr$piv)) { fac <- B[, Sr$piv[j] + 1L] %% p
      for (i in seq_len(nrow(B))) if (fac[i] != 0) B[i, ] <- (B[i, ] - .sym_mulmod(Sr$R[j, ], fac[i], p)) %% p }
    B }
  # the residual direction normalised on `ac`, in a CONSISTENT gauge across sample
  # points: reduce the free-column null space by the scalings, physical-first RREF, and
  # return the direction with ac = 1 and every physical pivot = 0. This is well defined
  # whether `ac` is itself a pivot here or not (the forward kernel frees the states, so
  # its pivot set can differ from the backward one that chose the anchor), and gives a
  # distinct direction per anchor -- unlike "first row with ac != 0", which collapses
  # anchors that share a row.
  extract <- function(rp, zvals, p, ac) {
    B <- t(vapply(freeColsOf(rp), function(fc) .sym_null_residues(rp, fc, p), integer(nz)))
    Bred <- reduceRows(B, tangentsAt(zvals, p), p)
    Bred <- Bred[apply(Bred, 1L, function(x) any(x %% p != 0)), , drop = FALSE]
    if (nrow(Bred) == 0L) return(NULL)
    rr <- .sym_phys_rref(Bred, physSet, nz, p)
    pr <- match(ac, rr$piv)
    if (!is.na(pr)) return(as.integer(rr$R[pr, ] %% p))          # ac is a pivot: its row
    if ((ac + 1L) > nz) return(NULL)
    v <- integer(nz); v[ac + 1L] <- 1L                          # ac free: 1 on ac, 0 on pivots
    for (i in seq_along(rr$piv)) v[rr$piv[i] + 1L] <- as.integer((p - rr$R[i, ac + 1L]) %% p)
    v }
  kxF <- function(pt, p, ac, sr) { rp <- kcallFwd(pt, p, sc$NtUsed, sr)
    if (!isTRUE(rp$ok) || rp$rank != sc$rank) return(NULL); extract(rp, pt[zSlots + 1L], p, ac) }

  # 1. anchor set: RREF the residual free-column null vectors (mod scalings) at a fresh
  # backward-valid point, over the PHYSICAL columns, so each residual direction gets a
  # DISTINCT physical pivot. With `listAnchors` return the whole set; otherwise
  # reconstruct the one at `acIn` (or the first). This replaces "always the first
  # residual", which collapsed every call -- over a >1-dimensional residual space --
  # onto a single (often already-known) direction.
  ptbB <- NULL
  for (att in seq_len(200L)) { cand <- draw(npt)
    r <- kcall(cand, P1, sc$NtUsed); if (isTRUE(r$ok) && r$rank == sc$rank) { ptbB <- cand; break } }
  if (is.null(ptbB)) { if (listAnchors) return(integer(0)); return(fb("forward: no backward base for support")) }
  rp0 <- kcall(ptbB, P1, sc$NtUsed)
  B0 <- t(vapply(freeColsOf(rp0), function(fc) .sym_null_residues(rp0, fc, P1), integer(nz)))
  Bred0 <- reduceRows(B0, tangentsAt(ptbB[zSlots + 1L], P1), P1)
  Bred0 <- Bred0[apply(Bred0, 1L, function(x) any(x %% P1 != 0)), , drop = FALSE]
  if (nrow(Bred0) == 0L) { if (listAnchors) return(integer(0)); return(fb("forward: empty residual")) }
  rr0 <- .sym_phys_rref(Bred0, physSet, nz, P1)
  anchors <- rr0$piv[rr0$piv %in% physSet]
  if (listAnchors) return(anchors)
  if (!length(anchors)) return(fb("forward: no physical residual anchor"))
  ac <- if (!is.null(acIn)) as.integer(acIn) else anchors[1]
  arow <- match(ac, rr0$piv)
  if (is.na(arow)) return(fb("forward: requested anchor is not a residual pivot"))
  v0 <- as.integer(.sym_mulmod(rr0$R[arow, ], .sym_invmod(rr0$R[arow, ac + 1L] %% P1, P1), P1))
  dirSupport <- sort(unique(stripName(znames[which(v0 %% P1 != 0)])))
  suppNames <- unique(stripName(znames[c(ac, setdiff(which(v0 %% P1 != 0) - 1L, ac)) + 1L]))

  # 2. turnover transversal (avoids the direction's support + the recast coords)
  keepFree <- unique(c(suppNames, grep("^_E_|^_L_", leafNames, value = TRUE)))
  tsv <- setNames(as.list(((seq_along(realStateNames) * 7919L + 11L) %% (P1 - 1L)) + 1L), realStateNames)
  tpv <- setNames(as.list(((seq_along(solveParamNames) * 104729L + 3L) %% (P1 - 1L)) + 1L), solveParamNames)
  ts <- tryCatch(sd$solveForwardModular(models[[1]], realStateNames, solveParamNames, tsv, tpv, P1,
                   forcings = if (length(solveHeld)) solveHeld else NULL,
                   keepFree = as.list(keepFree)), error = function(e) NULL)
  if (is.null(ts) || !isTRUE(ts$ok)) return(fb(paste("forward: no turnover transversal;",
                                                     if (is.null(ts)) "solve error" else ts$why)))
  solveRates <- as.character(ts$solveRates)
  if (rdiag) message(sprintf("[fwd] anchor=%s support=%d transversal(%d)={%s}",
                             znames[ac + 1L], length(suppNames), length(solveRates),
                             paste(solveRates, collapse = ",")))

  # 3. forward base valid at every prime (the linear solve is generically non-singular mod p)
  base0 <- NULL
  for (att in seq_len(80L)) { cand <- draw(npt)
    if (all(vapply(primes, function(pp) !is.null(kxF(cand, pp, ac, solveRates)), logical(1)))) { base0 <- cand; break } }
  if (is.null(base0)) return(fb("forward: no base valid at all primes"))

  # 4. forward relevance over EVERY free leaf (all but the solved turnover rates): params, log
  # leaves _L_, the recast E = base^exp coordinates _E_ (the Hill terms depend on these), and the
  # states. A missed relevant leaf makes an entry base0-specific -- caught by the fresh-POINT
  # verify in step 6, but scanning everything up front is what makes the reported form universal.
  cand <- setdiff(seq_along(leafNames), match(solveRates, leafNames))
  cand <- cand[!is.na(cand)]
  v0f <- kxF(base0, P1, ac, solveRates)
  physSupp <- setdiff(intersect(which(v0f %% P1 != 0) - 1L, physSet), ac)
  relBy <- replicate(length(physSupp), integer(0), simplify = FALSE)
  for (li in cand) { if (.sym_expired(ctrl)) return(fb("forward: timeout in relevance"))
    ch <- rep(FALSE, length(physSupp))
    for (dv in c(3L, 7L, 11L)) { pt <- as.numeric(base0); pt[li] <- pt[li] + dv
      v <- kxF(pt, P1, ac, solveRates); if (is.null(v)) next
      ch <- ch | ((v[physSupp + 1L] - v0f[physSupp + 1L]) %% P1 != 0) }
    for (i in which(ch)) relBy[[i]] <- c(relBy[[i]], li) }
  unionRel <- sort(unique(unlist(relBy)))
  if (!length(unionRel)) return(fb("forward: residual constant in every leaf"))
  maxRel <- max(vapply(relBy, length, integer(1)))
  if (rdiag) message(sprintf("[fwd] physSupp=%d unionRel=%d maxRel=%d", length(physSupp), length(unionRel), maxRel))
  if (nzchar(Sys.getenv("DMOD_SYM_FWDREL"))) {
    for (i in seq_along(physSupp))
      message(sprintf("[fwdrel] %-26s : %s", stripName(znames[physSupp[i] + 1L]),
                      paste(leafNames[relBy[[i]]], collapse = ", ")))
    return(fb("fwdrel diagnostic"))
  }
  # the forward direction legitimately couples more leaves than the eliminated path (it frees the
  # states plus BOTH recast partners E = base^exp and L = log base), so the whole-direction cap is
  # relaxed here; the fit size is still bounded by the per-entry cap relevanceCapSparse (maxRel).
  fwdCapDir <- max(as.integer(ctrl$relevanceCapDir), 48L)
  if (length(unionRel) > fwdCapDir || maxRel > ctrl$relevanceCapSparse)
    return(fb(sprintf("forward: couples %d leaves (entry up to %d)", length(unionRel), maxRel)))

  # 5. SHARED bank, filled LAZILY: forward points are valid at ALL primes (the linear solve has no
  # per-prime degeneracy), so the same points serve every prime. Each per-point kernel pass is the
  # dominant (serial) cost, so we top the bank up per fit-degree on demand instead of pre-filling
  # for the worst-case degree -- a low-degree direction then needs far fewer passes.
  dCap <- max(1L, min(3L, as.integer(ctrl$degreeCap)))
  needOf <- function(k, d) 2L * nrow(.sym_mono_table(k, d)) + 20L
  bankU <- matrix(0L, 0L, length(unionRel)); bankV <- lapply(primes, function(.) matrix(0L, 0L, nz))
  fillBankTo <- function(n) {
    tries <- 0L
    while (nrow(bankU) < n && tries < 12L * n) { if (.sym_expired(ctrl)) break
      tries <- tries + 1L; pt <- as.numeric(base0); pt[unionRel] <- draw(length(unionRel))
      vs <- lapply(primes, function(pp) kxF(pt, pp, ac, solveRates))
      if (any(vapply(vs, is.null, logical(1)))) next
      bankU <<- rbind(bankU, pt[unionRel])
      for (jp in seq_along(primes)) bankV[[jp]] <<- rbind(bankV[[jp]], vs[[jp]]) }
    if (rdiag && nrow(bankU) >= n) message(sprintf("[fwd] bank at %d points (all %d primes)", nrow(bankU), length(primes)))
    nrow(bankU) >= n }

  recEntry <- function(reli, c9) {
    reliCols <- match(reli, unionRel); vars <- leafNames[reli]
    for (d in 0:dCap) {
      mons <- .sym_mono_table(length(reli), d); nMon <- nrow(mons)
      if (!fillBankTo(needOf(length(reli), d))) next
      need <- nrow(bankU)
      refFree <- NULL; coefRes <- matrix(0L, 2L * nMon, length(primes)); ok <- TRUE
      for (jp in seq_along(primes)) { pp <- primes[jp]
        sU <- matrix(as.integer(bankU[seq_len(need), reliCols, drop = FALSE]), need)
        rv <- as.integer(bankV[[jp]][seq_len(need), c9 + 1L])
        fit <- symFitRational(sU, matrix(as.integer(mons), nMon), rv, pp)
        if (!identical(fit$status, "ok")) { ok <- FALSE; break }
        raw <- as.numeric(fit$coeffs); if (is.null(refFree)) refFree <- fit$freeCol
        dn <- raw[refFree + 1L] %% pp; if (dn == 0) { ok <- FALSE; break }
        coefRes[, jp] <- as.integer(.sym_mulmod(raw, .sym_invmod(dn, pp), pp)) }
      if (!ok) next
      rec <- tryCatch(sd$symRatReconBig(coefRes, as.integer(primes)), error = function(e) NULL)
      if (is.null(rec) || any(rec$den == "0")) next
      numI <- seq_len(nMon); denI <- nMon + seq_len(nMon)
      numStr <- .sym_poly_string(rec$num[numI], rec$den[numI], mons, vars)
      denStr <- .sym_poly_string(rec$num[denI], rec$den[denI], mons, vars)
      return(list(expr = if (denStr == "1") numStr else paste0("(", numStr, ")/(", denStr, ")"), reli = vars)) }
    NULL }

  entries <- list(); entryReli <- list()
  entries[[znames[ac + 1L]]] <- "1"; entryReli[[znames[ac + 1L]]] <- character(0)
  for (i in seq_along(physSupp)) { if (.sym_expired(ctrl)) return(fb("forward: timeout in fitting"))
    c9 <- physSupp[i]; nm <- znames[c9 + 1L]; reli <- relBy[[i]]
    if (!length(reli)) {
      if (!fillBankTo(1L)) return(fb("forward: empty bank"))
      resc <- vapply(seq_along(primes), function(jp) bankV[[jp]][1, c9 + 1L], integer(1))
      rc <- tryCatch(sd$symRatReconBig(matrix(as.integer(resc), 1L), as.integer(primes)), error = function(e) NULL)
      if (is.null(rc) || rc$den[1] == "0") return(fb("forward: a constant entry could not be lifted"))
      entries[[nm]] <- if (rc$den[1] == "1") rc$num[1] else paste0(rc$num[1], "/", rc$den[1])
      entryReli[[nm]] <- character(0); next }
    rr <- recEntry(reli, c9)
    if (rdiag) message(sprintf("[fwd]   entry %-24s nvar=%d -> %s", stripName(nm), length(reli),
                               if (is.null(rr)) sprintf("NO FIT (deg>%d)", dCap) else "ok"))
    if (is.null(rr)) return(fb(sprintf("forward: entry %s not a bounded rational (deg>%d)", nm, dCap)))
    entries[[nm]] <- if (is.null(spy)) rr$expr else .sym_simplify(rr$expr, spy); entryReli[[nm]] <- rr$reli }

  # 6. verify at a fresh prime AND a fully FRESH point (every leaf redrawn, not just the
  # relevant ones): evaluating each entry from its relevant leaves alone must still reproduce the
  # kernel residue when the NON-relevant leaves also differ. This certifies the reported form is
  # base-independent (the excluded leaves really are irrelevant), not an artifact of base0's slice.
  ptv <- NULL
  for (t in seq_len(400L)) { pt <- draw(npt)
    if (!is.null(kxF(pt, qv, ac, solveRates))) { ptv <- pt; break } }
  if (is.null(ptv)) return(fb("forward: no fresh verification point"))
  vver <- kxF(ptv, qv, ac, solveRates)
  for (nm in names(entries)) { col <- match(nm, znames); rl <- entryReli[[nm]]
    pred <- tryCatch(sd$evalRationalMod(entries[[nm]], as.list(rl),
              as.list(as.integer(ptv[match(rl, leafNames)])), qv), error = function(e) NULL)
    okv <- !is.null(pred) && ((as.integer(pred) - vver[col]) %% qv == 0)
    if (rdiag) message(sprintf("[fwd]   verify %-24s -> %s", stripName(nm), if (okv) "ok" else "FAIL"))
    if (!okv) return(fb("forward: reconstructed direction failed fresh-point verification")) }
  if (rdiag) message("[fwd] VERIFIED (base-independent) -- closing direction (states are resting-level symbols)")
  list(support = sort(names(entries)), vector = entries, type = "general", closedForm = TRUE)
}


# Reconstruct EVERY residual (non-scaling) direction the coupled/recast case leaves to
# the forward path. One shared anchor set (distinct physical pivots via .sym_phys_rref),
# then the single-direction forward per anchor, so a residual space of dimension > 1 is
# recovered as that many DISTINCT directions instead of the first one repeated. Returns a
# list of directions (each closed or support-only), or NULL if no anchor set forms.
.sym_perprime_forward_multi <- function(residualFree, sc, kcall, kcallFwd, znames, zSlots,
                                        leafNames, nz, scaling, stateColNames, paramNames,
                                        recast, sd, spy, ctrl, physCols, models,
                                        realStateNames, solveParamNames, solveHeld) {
  if (!length(residualFree)) return(list())
  # reconstruct the WHOLE non-scaling residual (all free columns mod the exact
  # scalings): a distinct physical anchor per residual direction. The caller adopts
  # the set only if every one closes, and replaces BOTH the peel and the per-column
  # results, so the peeled/forward split never has to be reconciled.
  fwd1 <- function(acIn, listAnchors)
    .sym_perprime_forward(residualFree[1], sc, kcall, kcallFwd, znames, zSlots, leafNames,
                          nz, scaling, stateColNames, paramNames, recast, sd, spy, ctrl,
                          physCols, models, realStateNames, solveParamNames, solveHeld,
                          acIn = acIn, listAnchors = listAnchors)
  anchors <- fwd1(NULL, TRUE)
  if (!length(anchors)) return(NULL)
  lapply(anchors, function(a) fwd1(a, FALSE))
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


# ---- verification: nullspace membership, fresh-prime re-check, saturation guard ------

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


# STRICT verification for the per-prime (coupled steady-state) reconstruction. There
# a random point rarely admits a modular steady state, so the lenient verifiers
# above accept INCONCLUSIVELY (they cannot reject when the fresh-point solve fails) --
# which would let a spurious low-degree fit through (e.g. a false "constant" that the
# few solvable points happened to agree on). Here retry fresh points until one solves
# pivot-consistently at the verify prime, then REQUIRE the closed form to reproduce
# the nullspace there (a garbage fit does not evaluate to the true residues at a fresh
# point and is rejected). No valid point found within the budget also fails, so the
# per-prime path never reports an unverified closed form. Returns TRUE only on a
# decisive pass at a genuine steady-state point.
.sym_verify_perprime <- function(entry, f, znames, leafNames, point0, NtUsed,
                                 kcall, pivots, pool, poolNext, nz, sd,
                                 relLeaves = NULL, tries = 80L) {
  if (is.null(sd) || is.null(entry$vector)) return(FALSE)
  q <- .symVerifyPrime
  pn <- poolNext
  # perturb only the direction's relevant leaves (as the reconstruction did): a
  # full random point almost never admits an interior modular steady state, so it
  # would never yield a verification point; perturbing the relevant leaves off the
  # base matches the sampling and solves at a workable rate.
  pertIdx <- if (is.null(relLeaves)) seq_along(point0) else relLeaves
  for (t in seq_len(tries)) {
    pt <- as.numeric(point0)
    if (length(pertIdx)) pt[pertIdx] <- pool(pn + seq_along(pertIdx) - 1L)
    pn <- pn + length(pertIdx) + 1L
    rq <- kcall(pt, q, as.integer(NtUsed))
    # an all-constant direction (no relevant leaf) is checked once at the base point
    if (!length(pertIdx) && !isTRUE(rq$ok)) return(FALSE)
    if (!isTRUE(rq$ok) || !identical(as.integer(rq$pivots), as.integer(pivots))) next
    env <- setNames(as.list(pt[seq_along(leafNames)]), leafNames)
    v <- integer(nz); v[f + 1L] <- 1L; bad <- FALSE
    for (nm in names(entry$vector)) {
      col <- match(nm, znames); if (is.na(col)) next
      got <- .sym_eval_modq(entry$vector[[nm]], env, q, sd)
      if (is.na(got)) { bad <- TRUE; break }
      v[col] <- as.integer(got %% q)
    }
    if (bad) next
    freeCols <- setdiff(seq_len(nz) - 1L, as.integer(rq$pivots))
    recon <- numeric(nz)
    for (fcol in freeCols) {
      vf <- v[fcol + 1L] %% q
      if (vf != 0) recon <- (recon + .sym_mulmod(.sym_null_residues(rq, fcol, q), vf, q)) %% q
    }
    return(all((recon - v) %% q == 0))
  }
  FALSE
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

  # plateauNeed consecutive non-growing Lie orders required before declaring the rank
  # saturated. A single plateau stops early on models whose rank has an intermediate
  # plateau then grows again (e.g. a Hill exponent that only becomes observable through a
  # high-order derivative -- the SMAD symbolic rank runs ...227,227,230,230,230,233...),
  # spuriously reporting those params as non-identifiable. Default 3 clears the
  # intermediate plateaus seen in practice; the verify = TRUE cross-check (integer
  # exponents saturate cleanly, no such plateaus) is the backstop for anything a finite
  # plateau would still miss. DMOD_SYM_LIEPLATEAU overrides; DMOD_SYM_LIEDIAG traces.
  plateauNeed <- max(1L, as.integer(Sys.getenv("DMOD_SYM_LIEPLATEAU", "3")))
  lieDiag <- nzchar(Sys.getenv("DMOD_SYM_LIEDIAG"))
  saturateNt <- function(point, Mtot) {
    prev <- -1L; Nt <- 1L; res <- NULL; flat <- 0L; ranks <- integer(0)
    repeat {
      r <- kcall(point, P, Nt, Mtot)
      if (!isTRUE(r$ok)) return(NULL)
      res <- r; ranks <- c(ranks, r$rank)
      flat <- if (r$rank == prev) flat + 1L else 0L
      if (r$rank >= nz || (flat >= plateauNeed && Nt >= 2L) || Nt > nz + 1L) break
      prev <- r$rank; Nt <- Nt + 1L
    }
    if (lieDiag) message("[liediag] Mtot=", Mtot, " ranks by Lie order: ",
                         paste(ranks, collapse = ","), " (nz=", nz, ")")
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


# Nullspace basis (columns) of an arbitrary GF(p) matrix M with nc columns.
.sym_null_basis_of <- function(M, nc, p) {
  rr <- .sym_rref_modp(M, p)
  free <- setdiff(seq_len(nc) - 1L, rr$piv)
  if (!length(free)) return(matrix(0L, nc, 0L))
  vapply(free, function(fc) {
    v <- integer(nc); v[fc + 1L] <- 1L
    for (ri in seq_along(rr$piv))
      v[rr$piv[ri] + 1L] <- as.integer((p - rr$R[ri, fc + 1L]) %% p)
    v
  }, integer(nc))
}


# RREF row-basis (rank x nz) of the intersection of the column spans of A and B
# over GF(p), via the nullspace of [A | -B]: a combination (x, y) there gives a
# vector A x = B y common to both spans. The canonical RREF makes the result
# comparable across primes/points.
.sym_intersect_rref <- function(A, B, nz, p) {
  a <- ncol(A); b <- ncol(B)
  if (a == 0L || b == 0L) return(matrix(0L, 0L, nz))
  Am <- matrix(as.numeric(A) %% p, nz, a)
  M <- cbind(Am, (p - matrix(as.numeric(B) %% p, nz, b)) %% p)   # [A | -B]
  nb <- .sym_null_basis_of(M, a + b, p)
  if (ncol(nb) == 0L) return(matrix(0L, 0L, nz))
  X <- matrix(as.numeric(nb[seq_len(a), , drop = FALSE]), a)
  inter <- .sym_mulmod_matmat(Am, X, p)                          # nz x k, = A x
  rr <- .sym_rref_modp(t(inter), p)
  rr$R[seq_len(rr$rank), , drop = FALSE]
}


# GF(p) matrix product A (m x k) %*% B (k x n), split-multiply to stay exact.
.sym_mulmod_matmat <- function(A, B, p) {
  m <- nrow(A); n <- ncol(B); k <- ncol(A)
  out <- matrix(0, m, n)
  for (j in seq_len(n)) {
    acc <- numeric(m)
    for (l in seq_len(k)) acc <- (acc + .sym_mulmod(A[, l], B[l, j], p)) %% p
    out[, j] <- acc
  }
  out
}


# Exact translation lattice: the CONSTANT tangents that lie in the observability
# nullspace at EVERY point. A translation z_i -> z_i + eps*a_i has a point-
# independent tangent a, so it survives the intersection of the nullspaces over
# several generic points; a scaling's tangent w*z varies with the point and drops
# out. The surviving GF(p) subspace is the translation lattice; each basis vector's
# constant components are lifted to exact rationals by single-prime rational
# reconstruction. Reconstruction-free (no rational fit over the whole coordinate),
# so the directions carry the same exactness as the scalings. Returns the
# translation entries and the Bmat with their tangents appended (excluding them
# from the residual reconstruction), or NULL if generic points could not be
# gathered (caller then skips the peel -- no regression).
.sym_peel_translations <- function(sc, kcall, freeCols, N, nz, P, ctrl, Bmat, znames) {
  if (!ncol(N)) return(list(translations = list(), Bmat = Bmat))
  npts <- if (is.null(ctrl$translPoints)) 3L else as.integer(ctrl$translPoints)
  C <- N
  nL <- length(sc$point0)
  base <- sc$poolNext + 4096L                    # clear of the reconstruction probes
  for (t in seq_len(npts)) {
    if (!ncol(C)) break
    cand <- NULL
    for (att in seq_len(ctrl$probeRetries)) {
      off <- base + ((t - 1L) * ctrl$probeRetries + (att - 1L)) * (nL + 1L)
      pert <- sc$point0
      for (li in seq_len(nL)) pert[li] <- sc$pool(off + li)
      cnd <- tryCatch(kcall(pert, P, sc$NtUsed), error = function(e) NULL)
      if (!is.null(cnd) && isTRUE(cnd$ok) &&
          identical(as.integer(cnd$pivots), as.integer(sc$pivots))) { cand <- cnd; break }
    }
    if (is.null(cand)) return(NULL)              # graceful skip: no regression
    Nk <- .sym_nullspace_basis(cand, freeCols, P)
    rows <- .sym_intersect_rref(C, Nk, nz, P)
    C <- if (nrow(rows)) t(rows) else matrix(0L, nz, 0L)
  }
  if (!ncol(C)) return(list(translations = list(), Bmat = Bmat))
  latt <- .sym_rref_modp(t(C), P)                # canonical RREF row-basis mod P
  rows <- latt$R[seq_len(latt$rank), , drop = FALSE]
  translations <- list()
  for (r in seq_len(nrow(rows))) {
    tang <- as.integer(rows[r, ] %% P)
    if (all(tang == 0) || .sym_in_span(Bmat, tang, nz, P)) next
    rec <- tryCatch(symRatRecon(matrix(tang, ncol = 1L), as.integer(P)),
                    error = function(e) NULL)
    if (is.null(rec) || any(rec$den == "0")) next
    num <- rec$num; den <- rec$den
    comps <- list()
    for (i in seq_len(nz)) {
      if (num[i] == "0") next
      comps[[znames[i]]] <- if (den[i] == "1") num[i] else paste0(num[i], "/", den[i])
    }
    if (!length(comps)) next
    Bmat <- cbind(Bmat, tang)
    translations[[length(translations) + 1L]] <- list(
      support = sort(names(comps)), vector = comps,
      type = "translation", closedForm = TRUE)
  }
  list(translations = translations, Bmat = Bmat)
}


# Classify one closed-form direction into {scaling, translation, affine,
# polynomial, general} on its canonical poly-primitive generator, via the Python
# classifier. A generator is defined only up to a nonzero function h(z); the
# classifier fixes that gauge, so the modular and symbolic engines land on the
# same representative (a "disguised" scaling like xi = (1, -ktl/ktx, m/ktx) is
# recognised as the scaling (ktx, -ktl, m)). A native scaling with an integer
# weight is canonicalised too; one with a symbolic (Hill) weight is left as-is,
# since the classifier reads it as a degree-2 polynomial. Returns the direction
# with $type/$vector/$degree updated (scaling: $vector = integer weights, so the
# weight convention downstream is preserved; every other class: $vector = the
# canonical generator components).
.sym_classify_direction <- function(d, sd, degreeCap) {
  if (is.null(d$vector)) return(d)
  if (isTRUE(d$type == "scaling")) {
    wint <- suppressWarnings(as.integer(as.character(unlist(d$vector))))
    if (anyNA(wint)) return(d)                       # Hill / symbolic weight: keep
    lit <- setNames(as.list(paste0(as.character(unlist(d$vector)), "*",
                                   names(d$vector))), names(d$vector))
  } else {
    lit <- setNames(as.list(as.character(unlist(d$vector))), names(d$vector))
  }
  cls <- tryCatch(sd$classifyDirection(lit, as.integer(degreeCap)),
                  error = function(e) NULL)
  if (is.null(cls) || is.null(cls$type)) return(d)
  d$type <- as.character(cls$type)
  d$degree <- if (!is.null(cls$degree)) as.integer(cls$degree) else NULL
  if (identical(d$type, "scaling")) {
    d$vector <- lapply(cls$weights, function(w) as.character(w))
    d$support <- sort(names(cls$weights))
  } else {
    d$vector <- lapply(cls$components, function(x) as.character(x))
    d$support <- sort(names(cls$components))
  }
  d
}


# Post-process every closed-form direction through the classifier so both the
# modular and symbolic engines report the same canonical generator and symmetry
# class. Support-only directions (no closed form) are left untouched. Runs at the
# top-level observability returns, so it covers both symEngine paths.
.sym_relabel_directions <- function(nonId, sd, degreeCap) {
  if (is.null(sd) || !length(nonId)) return(nonId)
  lapply(nonId, function(d) if (is.null(d$vector)) d
                            else .sym_classify_direction(d, sd, degreeCap))
}


# Optional certificate (reconstControl(certifyPoly = TRUE)): run the polynomial
# Lie-symmetry engine at the relevant degree and flag each affine/polynomial
# direction that is a nonzero constant linear combination of the returned
# generators -- i.e. an exact polynomial Lie POINT symmetry (the strict
# determining-equation notion), not merely an observability non-identifiability.
# Sets $certified (TRUE/FALSE); a FALSE just means the direction is not a strict
# polynomial symmetry, it stays a valid non-identifiability. `modelLines`/`obsLines`
# are the already-serialised f and g (toLines is a local closure of the caller).
.sym_certify_poly <- function(nonId, modelLines, obsLines, forcings, fixed,
                              parameters, ctrl, sd) {
  cand <- Filter(function(d) isTRUE(d$type %in% c("affine", "polynomial")), nonId)
  if (!length(cand) || is.null(sd)) return(nonId)
  deg <- if (!is.null(ctrl$certifyPolyDeg)) as.integer(ctrl$certifyPolyDeg)
         else max(vapply(cand, function(d)
           if (is.null(d$degree)) 1L else as.integer(d$degree), integer(1)))
  gens <- tryCatch({
    reticulate::py_capture_output(
      r <- sd$symmetryDetectiondMod(
        model = modelLines, observation = obsLines,
        ansatz = "multi", pMax = as.integer(max(1L, deg)),
        inputs = if (length(forcings)) forcings else NULL, fixed = fixed,
        allTrafos = TRUE, lieOrder = 0L, exact = TRUE, verify = FALSE,
        backend = "symengine",
        parameters = if (length(parameters)) parameters else NULL, method = "liesym"))
    r
  }, error = function(e) NULL)
  genVecs <- if (is.null(gens)) list()
             else Filter(Negate(is.null), lapply(gens, function(g) g$infinitesimals))
  lapply(nonId, function(d) {
    if (!isTRUE(d$type %in% c("affine", "polynomial"))) return(d)
    d$certified <- length(genVecs) > 0L &&
      isTRUE(tryCatch(sd$certifyInSpan(d$vector, genVecs)$certified,
                      error = function(e) FALSE))
    d
  })
}


# Peel the scaling directions common to all conditions. Each generator is
# projected onto the unknown coordinates z (state weights outside z, e.g. under a
# steady-state constraint, drop out), turned into its tangent w_c * z_c at the
# base point, verified to lie in the nullspace, and kept only if independent of
# the scalings already taken. Returns the scaling entries and their tangent span.
.sym_peel_scalings <- function(scalRes, znames, nz, zval, P, N, sd = NULL) {
  scaling <- list()
  Bmat <- matrix(0L, nz, 0L)
  inSpan <- function(M, x) .sym_in_span(M, x, nz, P)
  # A weight may be symbolic in a free exponent (a Hill scaling xi_kinh = -nhill*kinh).
  # Evaluate it at the base point (each coordinate's finite-field value) for the
  # tangent/validation, but report the symbolic weight verbatim.
  env <- as.list(setNames(as.numeric(zval), znames))
  evalW <- function(w) {
    wi <- suppressWarnings(as.numeric(w))
    if (!is.na(wi) && wi == round(wi)) return(as.numeric(wi %% P))
    if (is.null(sd)) return(NA_real_)
    v <- .sym_eval_modq(w, env, P, sd)
    if (is.null(v) || is.na(v)) NA_real_ else as.numeric(v)
  }
  for (d in scalRes$nonIdentifiable) {
    cols <- match(names(d$vector), znames)
    keep <- !is.na(cols)
    if (!any(keep)) next
    wsym <- as.character(unlist(d$vector))[keep]
    wv <- vapply(wsym, evalW, numeric(1))            # weights evaluated mod P
    if (anyNA(wv)) next
    v <- integer(nz); v[cols[keep]] <- as.integer(wv %% P)
    if (all(v == 0)) next
    tangent <- as.integer((as.numeric(v) * zval) %% P)
    if (all(tangent == 0) || !inSpan(N, tangent) || inSpan(Bmat, tangent)) next
    Bmat <- cbind(Bmat, tangent)
    supp <- names(d$vector)[keep]
    ord <- order(supp)
    scaling[[length(scaling) + 1L]] <- list(
      support = supp[ord],
      vector = setNames(as.list(wsym[ord]), supp[ord]),
      type = "scaling", closedForm = TRUE)
  }
  list(scaling = scaling, Bmat = Bmat)
}


# Expand a multi-condition scaling result onto the wide per-condition joint
# coordinates: a state weight applies identically to each of the state's K
# per-condition columns (they are log-normalised, so one integer weight is shared);
# parameter weights are left as they are. A held-variable pivot state (heldParamOf:
# pivot -> its shared initial-value parameter) is NOT a wide coordinate; its weight is
# the scaling of its resting value, so it is carried by that shared parameter (a single
# column, not per-condition) -- this is how the pool scaling closes on the pivot's
# initial value instead of being dropped.
.sym_joint_expand_scal <- function(scalRes, stateBase, Kc, heldParamOf = character(0)) {
  scalRes$nonIdentifiable <- lapply(scalRes$nonIdentifiable, function(d) {
    vec <- list()
    for (nm in names(d$vector)) {
      w <- d$vector[[nm]]
      if (nm %in% names(heldParamOf))
        vec[[heldParamOf[[nm]]]] <- w
      else if (nm %in% stateBase)
        for (m in seq_len(Kc)) vec[[paste0(nm, "|c", m)]] <- w
      else vec[[nm]] <- w
    }
    d$vector <- vec
    d$support <- names(vec)
    d
  })
  scalRes
}


# ---- parallel steady-state solve pool (PSOCK workers) --------------------------------

# Worker-side steady-state solve for the parallel warm pool. Runs in a fresh
# PSOCK process that has imported the Python engine into `.sym_worker_sd` (a
# worker global set at cluster creation); every argument is plain serialisable
# data, so no R closure or Python handle crosses the process boundary. Returns the
# solve list (valBy, dfJx/dfJt for the joint constraint) or NULL on any failure.
.sym_solve_worker <- function(job, cargs) {
  sd <- tryCatch(get(".sym_worker_sd", envir = globalenv()), error = function(e) NULL)
  if (is.null(sd)) return(NULL)
  tryCatch(sd$solveSteadyStateModular(
    model = job$model, stateNames = cargs$stateNames,
    paramNames = cargs$paramNames, paramVals = job$paramVals, prime = job$p,
    forcings = cargs$forcings, t0events = job$t0events,
    recast = cargs$recast, lVals = job$lVals, jointMode = TRUE,
    heldStates = job$heldVals),
    error = function(e) NULL)
}


# A PSOCK worker pool that each imports the Python symmetryDetection engine, used
# to fill the per-point steady-state solve bank in parallel on every platform
# (mclapply forks only on unix; the closed-form reconstruction of an equilibrated
# model at scale is dominated by these independent coupled solves). The workers run
# Python only -- they never touch the dMod C++ kernel or R closures -- so the pool
# is cheap to stand up and robust. Returns the cluster, or NULL when parallelism is
# off or a worker cannot bring up its interpreter (the caller then stays serial).
.sym_make_solve_cluster <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n <= 1L) return(NULL)
  if (!requireNamespace("parallel", quietly = TRUE)) return(NULL)
  cl <- tryCatch(parallel::makeCluster(n), error = function(e) NULL)
  if (is.null(cl)) return(NULL)
  # bring up each worker's Python interpreter and import the engine into a worker
  # global. clusterCall (not clusterEvalq, which some R builds do not export) runs
  # the initialiser on every node; it is reparented to the global env so it
  # serialises without dragging this frame's cluster handle. Returns TRUE or an
  # error string per worker, so a Python that will not come up degrades to serial.
  codeDir <- system.file("code", package = "dMod")
  pyPath  <- Sys.getenv("RETICULATE_PYTHON")
  initFun <- function(codeDir, pyPath) {
    tryCatch({
      if (!requireNamespace("reticulate", quietly = TRUE)) stop("no reticulate")
      if (nzchar(pyPath)) {
        Sys.setenv(RETICULATE_PYTHON = pyPath)
        suppressWarnings(suppressMessages(
          tryCatch(reticulate::use_python(pyPath, required = FALSE),
                   error = function(e) NULL)))
      }
      sysmod <- reticulate::import("sys", convert = TRUE)
      if (!(codeDir %in% sysmod$path)) sysmod$path <- c(codeDir, sysmod$path)
      assign(".sym_worker_sd",
             reticulate::import("symmetryDetectionVersion2", convert = TRUE),
             envir = globalenv())
      TRUE
    }, error = function(e) conditionMessage(e))
  }
  environment(initFun) <- globalenv()
  res <- tryCatch(parallel::clusterCall(cl, initFun, codeDir, pyPath),
                  error = function(e) list(conditionMessage(e)))
  ok <- all(vapply(res, isTRUE, logical(1)))
  if (!isTRUE(ok)) {
    if (nzchar(Sys.getenv("DMOD_SYM_TIMING")))
      message("[sym] solve pool init failed: ",
              paste(unique(vapply(res, function(r)
                if (isTRUE(r)) "" else as.character(r)[1], character(1))), collapse = "; "))
    tryCatch(parallel::stopCluster(cl), error = function(e) NULL)
    return(NULL)
  }
  cl
}


# ---- observability engine: build [Obs; df], peel scalings, reconstruct residuals -----

# Multi-condition analytic observability. `multi` is the list returned by
# compileObservabilityTapeMulti: per-condition tapes over a shared coordinate
# space. The observability rows of all conditions are stacked and reduced once,
# so the verdict and every reconstructed direction reflect the intersection
# nullspace. Known inputs are baked into each condition and are not leaves, so a
# direction can only involve genuine parameters and free initial values.
.observability_analytic_multi <- function(multi, spy = NULL,
                                          closedForm = FALSE, sd = NULL, cores = 1,
                                          equilZeroStates = character(0),
                                          t0events = list(), nConditions = NULL,
                                          chainOf = NULL, nGaps = 0L,
                                          implicitSteadyState = FALSE,
                                          control = reconstControl(), verify = FALSE,
                                          translationsOnly = FALSE) {
  ctrl <- control
  jointSS <- isTRUE(multi$jointSteadyState) && isTRUE(implicitSteadyState)
  # `cores` drives two NESTED parallelism axes without oversubscribing: the coarse
  # per-point solve fork takes the whole budget (coresGLp), while the inner
  # observability kernel (coresCall, set below) threads over conditions/segments in the
  # non-forked phases and drops to serial only inside that fork. The batch paths (plain
  # ODE, gap seed-batch) are a single OpenMP call and use the full budget there.
  cores <- as.integer(max(1L, cores))
  coresGLp <- cores
  # The joint/gap reconstruction solves each sample point's steady state serially in
  # kcall (the batch kernel only covers the plain ODE path). Those solves are the
  # dominant cost at scale (relevance probe + sample bank) and are independent, so
  # fork the batch over the sample points. Fork is correct with reticulate/sympy on
  # unix (a fresh child per batch, results are plain R values); serial elsewhere. A
  # child crash surfaces as try-error -> NULL, handled like a failed solve downstream.
  parMap <- if (coresGLp > 1L && .Platform$OS.type == "unix")
    function(xs, f) lapply(
      parallel::mclapply(xs, f, mc.cores = coresGLp, mc.preschedule = TRUE),
      function(o) if (inherits(o, "try-error")) NULL else o)
    else function(xs, f) lapply(xs, f)
  # On unix parMap forks the whole per-point kcall (mclapply). Fork is unavailable
  # elsewhere (Windows), so there the dominant per-point cost -- the coupled
  # steady-state solve -- is filled in parallel through a PSOCK pool of Python
  # interpreters instead (warmSolves below populates the solve cache; the kernel and
  # reduce then run on the master over cached solves). The cache dedups either way.
  useSolveCluster <- coresGLp > 1L && .Platform$OS.type != "unix" && isTRUE(closedForm)
  cl <- NULL
  warmSolves <- function(pts, primes) invisible()   # replaced in the jointSS block
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
  # a per-call kernel parallelises over the conditions/segments (chain groups or
  # tapes), so more threads than units cannot help and oversubscribe: each of the
  # many small reconstruction solves then pays the spawn/barrier cost of the idle
  # surplus threads. Cap at the unit count. The one place this NESTS inside the
  # per-point fork (the parMap in kbatch) drops it to 1 there so fork x thread does not
  # oversubscribe; every non-forked phase (saturation, verify, the plain batch) keeps
  # the full unit budget.
  nKernelUnits <- if (hasGaps) length(chainGroups) else length(multi$tapes)
  coresCall <- as.integer(max(1L, min(cores, nKernelUnits)))
  # the recast coordinates E = base^exp, L = log base are appended to the leaf
  # space as extra sampled coordinates; reconstruction fits the rational
  # dependence on them and back-substitutes, recovering closed forms that involve
  # a free exponent. With no recast, the augmented space equals the leaf space.
  recast <- list()
  nAug <- nLeaves
  leafNamesAug <- leafNames
  # leaves that are auxiliary coordinates (per-condition states, recast E/L) rather
  # than physical parameters; excluded from the relevance-cap gate in joint mode
  auxLeaves <- integer(0)
  # transient recast (free power/Hill exponent without equilibrate): E = base^exp and
  # L = log(base) are ordinary free-initial-value leaves tied to (base, exp) by the
  # recast relation rows stacked onto the observability codistribution. Set below.
  recastTransient <- FALSE
  recastAtomNames <- character(0)

  tapeFields <- function(t) {
    out <- list(
      op = as.integer(t$op), a = as.integer(t$a), b = as.integer(t$b),
      cnum = as.character(t$cnum), cden = as.character(t$cden),
      stateSlots = as.integer(t$stateSlots), fOut = as.integer(t$fOut),
      gOut = as.integer(t$gOut), icLeaf = as.integer(t$icLeaf),
      icNum = as.character(t$icNum), icDen = as.character(t$icDen))
    # a segment carries an IC tape seeding the state initial values (and their duals):
    # free/carry values, doses and resets, or -- for a joint equilibrate anchor -- the
    # identity seed R fills with the resting state per sample point.
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

  # batched per-chunk evaluator for the coupled+gap reconstruction loop (set inside
  # the joint block below when hasGaps); NULL keeps the serial parMap fallback.
  kchunk <- NULL

  if (ssConstraint) {
    # equilibrate mode: the states stay free coordinates, seeded on-manifold to the
    # resting state x* per (point, prime), and f = 0 enters as the stacked df tangency
    # rows of the joint determining system below.
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
    # held-variable initial-value parameters carry a state's name; the pivot state is
    # handled through the held-state channel of the solve (a frozen value), so its name
    # must NOT also be passed as a solver parameter -- it would duplicate the state gen
    # in the steady-state system. It stays a tape/report coordinate.
    heldSolveNames <- if (length(multi$heldStateParams))
      as.character(unlist(multi$heldStateParams)) else character(0)
    solveParamNames <- setdiff(c(paramNames, genNames), heldSolveNames)
    # in joint mode the recast coordinates E and L are already free-state leaves
    # (their own z-columns), so they must NOT be appended again as augmented
    # coordinates; the eliminated path augments the leaf space with them instead.
    if (length(recast) && !jointSS) {
      leafNamesAug <- c(leafNames, genNames, lNames)
      nAug <- length(leafNamesAug)
    }
    # The steady state enters implicitly via the joint determining system below:
    # the states stay coordinates and f = 0 is a tangency constraint. The former
    # explicit (eliminated) path (seed x* and its IFT duals into an icSeed) was
    # removed; an explicit steady state is supplied through `trafo` (steadyStates()).
    ssWhy <- NULL              # last steady-state failure reason, for diagnostics
    if (jointSS) {
      # implicit/joint determining system: keep the states as coordinates (seeded
      # on-manifold to x*), build the observability over (x, theta), and stack the
      # resting-manifold tangency rows df_rest = [Jx | Jt]. The nullspace of the
      # combined [Obs ; df] is the joint symmetry space; the scaling peel then
      # recovers the low-degree (often integer-weight) directions the eliminated
      # path cannot. No gap propagation in this first cut.
      slotOfName <- function(nm) { i <- match(nm, leafNames)
        if (is.na(i)) NULL else i }
      # held-variable pivots (equilibrate + reduceCQ = FALSE): one conserved-moiety
      # pivot per moiety keeps its resting value as a shared initial-value PARAMETER
      # (heldParamOf: pivot -> that param). The point solve freezes the pivot to that
      # parameter's residue (read from the parameter leaf, since the pivot is no longer
      # a free-state leaf) and solves the remaining states; the moiety freedom is then
      # reported on the parameter rather than a `total`.
      heldParamOf <- if (length(multi$heldStateParams))
        unlist(multi$heldStateParams) else character(0)   # named: pivot -> p_e
      heldNames <- setdiff(intersect(names(heldParamOf), realStateNames), solveHeld)
      jointPV <- function(point, p) {
        val <- function(nm) { s <- slotOfName(nm)
          if (is.null(s)) 0 else as.numeric(point[s]) %% p }
        pv <- as.list(setNames(vapply(solveParamNames, val, numeric(1)), solveParamNames))
        lv <- if (length(lNames))
          as.list(setNames(vapply(lNames, val, numeric(1)), lNames)) else NULL
        # freeze value of a pivot = residue of ITS initial-value parameter leaf
        hv <- if (length(heldNames))
          as.list(setNames(vapply(heldNames,
            function(nm) val(heldParamOf[[nm]]), numeric(1)), heldNames)) else NULL
        list(paramVals = pv, lVals = lv, heldVals = hv)
      }
      # Per-condition state coordinates. Each equilibrate condition c has its OWN
      # resting steady state x*_c, so its state perturbation xi_x,c is an INDEPENDENT
      # coordinate (only the parameters theta are shared across conditions). A single
      # shared state column would force xi_x,c = x*_c * w for one weight w and miss
      # every non-scaling multi-condition direction (silent over-identifiability). So
      # every state coordinate (the real states plus the recast partners E = base^exp
      # and L = log base) becomes per-condition; only the real params stay shared.
      # States and E are log-normalised per condition (dual x*_c) so a scaling's
      # weight is the same integer in each of its per-condition columns and the peel
      # recovers it; L is a log-shift, so it stays additive (dual 1), not normalised.
      perCondCols <- which(znames %in% as.character(multi$zStateNames))   # states + E + L
      logCols <- which(znames %in% setdiff(as.character(multi$zStateNames), lNames))  # states + E
      jointStateSlot <- zSlots[logCols] + 1L        # 1-based leaf slot to read x*_c
      equilConds <- which(isEquilTape)              # one anchor tape per condition
      Kc <- length(equilConds)
      sharedIdx <- setdiff(seq_len(nz), perCondCols)   # real params only
      nShared <- length(sharedIdx); nSt <- length(perCondCols)
      nzWide <- nShared + Kc * nSt
      # local znames-column (1..nz) -> wide column, for the mi-th equilibrate condition
      wideCols <- lapply(seq_len(Kc), function(mi) {
        wc <- integer(nz)
        wc[sharedIdx] <- seq_len(nShared)
        wc[perCondCols] <- nShared + (mi - 1L) * nSt + seq_len(nSt)
        wc
      })
      # wide coordinate labels/leaf slots: shared columns keep theirs; each state
      # coordinate is duplicated per condition (decorated name, same leaf slot -> the
      # value is overridden to 1 downstream, so all K copies share dual 1)
      stateBase <- znames[perCondCols]
      znamesWide <- c(znames[sharedIdx],
                      unlist(lapply(seq_len(Kc), function(mi)
                        paste0(stateBase, "|c", mi))))
      zSlotsWide <- c(zSlots[sharedIdx], rep(zSlots[perCondCols], times = Kc))
      zStateNamesWide <- if (nSt) znamesWide[(nShared + 1L):nzWide] else character(0)
      # local (per-condition) width/labels/slots the kernel and df builder use; kept
      # separate from the wide nz/znames/zSlots that overwrite the outer names below
      # (R closures would otherwise see the widened values at call time)
      nzL <- nz; znamesL <- znames; zSlotsL <- zSlots
      # Recast relation rows. E = base^exp and L = log(base) are free coordinates in
      # the joint system, but the df tangency only carries their (redundant) time
      # derivatives, so without the algebraic ties E and L are spuriously free and
      # everything looks non-identifiable. Add the linearised relations per condition:
      #   E:  xi_E/E   = exp * xi_base/base + log(base) * xi_exp   (log-weight columns)
      #   L:  xi_L     = xi_base/base
      # exp (the exponent parameter) and log(base) enter as numeric coefficients from
      # the sample point (log(base) is the generic L coordinate). E and base are
      # log-normalised columns, L is additive, exp is an additive param column.
      recastRel <- lapply(recast, function(rc) list(
        baseCol = match(as.character(rc$base), znamesL),
        ECol    = match(as.character(rc$E),    znamesL),
        LCol    = match(as.character(rc$L),    znamesL),
        expCol  = match(as.character(rc$exp),  znamesL),
        expSlot = slotOfName(as.character(rc$exp)),
        LSlot   = slotOfName(as.character(rc$L)),
        # a state base is log-normalised (its column carries xi_base/base); a parameter
        # base -- a Michaelis constant that only appears under the exponent, kept as a
        # coordinate so the pool/Km co-scaling is reported -- is not, so its column
        # carries xi_base and the relation coefficient picks up a 1/base factor.
        baseSlot  = slotOfName(as.character(rc$base)),
        baseParam = as.character(rc$base) %in% paramNames))
      # constraint columns that legitimately map to no z-coordinate: forcings and
      # sink-cluster states (held at 0) and any `fixed` leaf. A df column outside
      # this set that fails to map into znames is a lost constraint (a bug) - flag it.
      dfDroppable <- unique(c(as.character(forcings), as.character(equilZeroStates),
                              setdiff(as.character(leafNames), znames)))
      # df constraint rows for one condition over the nz znames columns
      dfRowsCond <- function(sol, p) {
        Jx <- sol$dfJx; Jt <- sol$dfJt
        scn <- as.character(sol$dfStateCols); pcn <- as.character(sol$dfParamCols)
        rows <- vector("list", length(Jx))
        for (i in seq_along(Jx)) {
          row <- numeric(nzL); xi <- as.numeric(Jx[[i]])
          for (j in seq_along(scn)) { col <- match(scn[j], znamesL)
            if (!is.na(col)) row[col] <- (row[col] + xi[j]) %% p
            else if (!scn[j] %in% dfDroppable && nzchar(Sys.getenv("DMOD_JOINT_DIAG")))
              message("[jointdiag] df drops non-fixed state column: ", scn[j]) }
          for (th in pcn) {
            # a held-variable pivot column is labelled by the pivot state name; it maps
            # to that pivot's shared initial-value parameter (the reduced tangency folds
            # the pivot's df column onto its p_e coordinate)
            thz <- if (th %in% names(heldParamOf)) heldParamOf[[th]] else th
            col <- match(thz, znamesL)
            if (!is.na(col)) row[col] <- (row[col] + as.numeric(Jt[[th]])[i]) %% p
            else if (!th %in% dfDroppable && nzchar(Sys.getenv("DMOD_JOINT_DIAG")))
              message("[jointdiag] df drops non-fixed param column: ", th) }
          rows[[i]] <- row
        }
        if (!length(rows)) matrix(0, 0, nzL) else do.call(rbind, rows)
      }
      # The coupled steady-state solve of one condition depends on the point ONLY
      # through that condition's steady-state parameter values (paramVals), log
      # coordinates (lVals) and any held-variable pivot residues (heldVals) at the
      # prime -- not the observation scales, the other free initial values or any leaf
      # outside f. So the solve is memoised on that subvector: a relevance probe that
      # perturbs a non-solve leaf, or two sample points that agree on the solve
      # parameters, reuse one (expensive, often Groebner) solve. Negatives are cached
      # too so a doomed perturbation is not re-solved. The cache is the shared sink the
      # parallel warm pool fills.
      jointSolveCache <- new.env(parent = emptyenv())
      solveKey <- function(p, pv, ci)
        paste(ci, p, paste0(unlist(pv$paramVals), collapse = ","),
              paste0(unlist(pv$lVals), collapse = ","),
              paste0(unlist(pv$heldVals), collapse = ","), sep = "|")
      # constant (point-independent) arguments of every solve, shipped once to the pool
      solveConst <- list(stateNames = realStateNames, paramNames = solveParamNames,
                         forcings = if (length(solveHeld)) solveHeld else NULL,
                         recast = if (length(recast)) recast else NULL)
      solveRaw <- function(p, pv, ci) {
        evC <- if (ci <= length(t0events) && length(t0events[[ci]]))
          t0events[[ci]] else NULL
        tryCatch(sd$solveSteadyStateModular(
          model = models[[ci]], stateNames = realStateNames,
          paramNames = solveParamNames, paramVals = pv$paramVals, prime = p,
          forcings = if (length(solveHeld)) solveHeld else NULL, t0events = evC,
          recast = if (length(recast)) recast else NULL, lVals = pv$lVals,
          jointMode = TRUE,
          heldStates = if (length(heldNames)) pv$heldVals else NULL),
          error = function(e) NULL)
      }
      # Fill the solve cache for a batch of points in parallel over the PSOCK pool
      # (Windows path). Enumerate the distinct uncached (condition, subvector) jobs,
      # solve them across the Python workers, and store each result (or a negative)
      # under its key. A no-op without a pool: jointSolveCond then solves on demand.
      warmSolves <- function(pts, primes) {
        if (is.null(cl)) return(invisible())
        jobs <- list(); seen <- new.env(parent = emptyenv())
        for (idx in seq_along(pts)) {
          p <- primes[[idx]]; pt <- pts[[idx]]
          for (ci in equilConds) {
            pv <- jointPV(pt, p); key <- solveKey(p, pv, ci)
            if (!is.null(jointSolveCache[[key]]) || !is.null(seen[[key]])) next
            seen[[key]] <- TRUE
            evC <- if (ci <= length(t0events) && length(t0events[[ci]]))
              t0events[[ci]] else NULL
            jobs[[length(jobs) + 1L]] <- list(key = key, ci = ci, p = p,
              model = models[[ci]], paramVals = pv$paramVals, lVals = pv$lVals,
              heldVals = if (length(heldNames)) pv$heldVals else NULL,
              t0events = evC)
          }
        }
        if (!length(jobs)) return(invisible())
        # reparent the worker to the global env so it serialises self-contained (the
        # Python-only workers have no dMod namespace to resolve it against); its body
        # uses only base ops and the worker-global .sym_worker_sd handle
        worker <- .sym_solve_worker
        environment(worker) <- globalenv()
        res <- tryCatch(parallel::parLapply(cl, jobs, worker, cargs = solveConst),
                        error = function(e) vector("list", length(jobs)))
        for (i in seq_along(jobs))
          jointSolveCache[[jobs[[i]]$key]] <-
            if (is.null(res[[i]]) || !isTRUE(res[[i]]$ok)) list(ok = FALSE) else res[[i]]
        invisible()
      }
      # solve one equilibrate condition and return its solution plus the point with
      # this condition's on-manifold x* written into the state leaves (cache-backed)
      jointSolveCond <- function(point, p, ci) {
        pv <- jointPV(point, p)
        key <- solveKey(p, pv, ci)
        sol <- jointSolveCache[[key]]
        if (is.null(sol)) {
          sol <- solveRaw(p, pv, ci)
          jointSolveCache[[key]] <- if (is.null(sol) || !isTRUE(sol$ok)) list(ok = FALSE) else sol
          if (is.null(sol) || !isTRUE(sol$ok)) return(NULL)
        } else if (!isTRUE(sol$ok)) return(NULL)
        ptc <- as.numeric(point[seq_len(nLeaves)])
        # seed the state leaves with the PRE-event resting value (valBy): the IC tape
        # applies the event map E(x_ss) itself, and the df constraint is linearised at
        # the same pre-event x_ss - so obs and df share one operating point
        vb <- sol$valBy
        for (k in seq_along(realStateNames)) {
          s <- slotOfName(realStateNames[k]); v <- vb[[realStateNames[k]]]
          if (!is.null(s) && !is.null(v)) ptc[s] <- as.numeric(v) %% p
        }
        list(sol = sol, ptc = ptc)
      }
      # FORWARD solve variant: the states + free params are CHOSEN (read off the point's
      # leaves) and f = 0 is solved LINEARLY for a turnover subset of rate constants. Unlike
      # the backward Groebner solve this is valid at every prime, so a coordinate the direction
      # depends on (a resting state x*, e.g. C3) is a FREE, SHARED sample coordinate instead of
      # a per-prime-inconsistent x*(theta). `solveRates` is the fixed turnover transversal
      # (chosen once to avoid the direction's support). Returns the same {sol, ptc} shape as the
      # backward solve, with the SOLVED rates written into their leaves and the chosen states.
      jointSolveCondFwd <- function(point, p, ci, solveRates) {
        rd <- function(nm) { s <- slotOfName(nm); if (is.null(s)) 0 else as.numeric(point[s]) %% p }
        sv <- as.list(setNames(vapply(realStateNames, rd, numeric(1)), realStateNames))
        pv <- as.list(setNames(vapply(solveParamNames, rd, numeric(1)), solveParamNames))
        sol <- tryCatch(sd$solveForwardModular(models[[ci]], realStateNames, solveParamNames,
                          sv, pv, p, forcings = if (length(solveHeld)) solveHeld else NULL,
                          solveRates = solveRates), error = function(e) NULL)
        if (is.null(sol) || !isTRUE(sol$ok)) return(NULL)
        ptc <- as.numeric(point[seq_len(nLeaves)])
        for (r in names(sol$rates)) { s <- slotOfName(r)
          if (!is.null(s)) ptc[s] <- as.numeric(sol$rates[[r]]) %% p }
        for (nm in realStateNames) { s <- slotOfName(nm); v <- sol$valBy[[nm]]
          if (!is.null(s) && !is.null(v)) ptc[s] <- as.numeric(v) %% p }
        list(sol = sol, ptc = ptc)
      }
      # each condition is solved and observed at its OWN resting state; the state
      # columns are then log-normalised (multiplied by x*_c) so a scaling weight is
      # condition-independent, and the per-condition [Obs ; df] blocks are stacked.
      #
      # one condition's [Obs ; df] blocks embedded in the wide space, given its kernel
      # result `obs` and solve `sc0`. Shared verbatim by the serial kcall4 below and
      # the batched kchunk, so both paths are byte-identical. NULL = degenerate point.
      oneCondBlocks <- function(mi, obs, sc0, p) {
        if (!isTRUE(obs$ok)) return(NULL)
        oR <- matrix(as.numeric(obs$R), nrow = obs$rank, ncol = nzL)
        dR <- dfRowsCond(sc0$sol, p)
        # a state whose on-manifold value hits 0 mod p (unlucky point/prime, or a
        # genuinely-zero state) makes its log-normalised column degenerate - reject
        # the whole point so the saturator resamples at a fresh point/prime
        xvals <- as.numeric(sc0$ptc[jointStateSlot]) %% p
        if (any(xvals == 0)) {
          if (nzchar(Sys.getenv("DMOD_SYM_FWDDIAG")))
            message("[fwddiag] zero log-normal coord at slots ",
                    paste(jointStateSlot[xvals == 0], collapse = ","))
          return(NULL) }
        for (m in seq_along(logCols)) {
          xv <- xvals[m]
          oR[, logCols[m]] <- .sym_mulmod(oR[, logCols[m]], xv, p)
          if (nrow(dR)) dR[, logCols[m]] <- .sym_mulmod(dR[, logCols[m]], xv, p)
        }
        # embed this condition's [Obs ; df] into its own wide state block (shared
        # param/L columns overlap, state columns go to block mi)
        wc <- wideCols[[mi]]
        bl <- list()
        eO <- matrix(0, nrow(oR), nzWide); eO[, wc] <- oR
        bl[[length(bl) + 1L]] <- eO
        if (nrow(dR)) {
          eD <- matrix(0, nrow(dR), nzWide); eD[, wc] <- dR
          bl[[length(bl) + 1L]] <- eD
        }
        # tie E and L back to the base for this condition (numeric exp/log(base)
        # coefficients from this condition's sample point)
        for (rc in recastRel) {
          expv <- as.numeric(sc0$ptc[rc$expSlot]) %% p
          Lv   <- as.numeric(sc0$ptc[rc$LSlot]) %% p
          # 1/base factor for a non-log-normalised parameter base (state base: 1)
          bScale <- 1
          if (isTRUE(rc$baseParam)) {
            bv <- as.numeric(sc0$ptc[rc$baseSlot]) %% p
            if (bv == 0) return(NULL)                  # degenerate point, resample
            bScale <- .sym_invmod(bv, p)
          }
          r1 <- numeric(nzWide)                       # E = base^exp
          r1[wc[rc$ECol]]    <- 1
          r1[wc[rc$baseCol]] <- (p - .sym_mulmod(expv, bScale, p)) %% p   # -exp/base
          r1[wc[rc$expCol]]  <- (-Lv) %% p
          r2 <- numeric(nzWide)                       # L = log(base)
          r2[wc[rc$LCol]]    <- 1
          r2[wc[rc$baseCol]] <- (p - bScale) %% p                         # -1/base
          bl[[length(bl) + 1L]] <- rbind(r1, r2)
        }
        bl
      }
      # one condition's observability kernel at its seeded resting state: a chain of
      # segments across post-t0 event gaps, or the single equilibrate segment.
      condObs <- function(ci, ptc, p, Nt, Mtot) {
        if (hasGaps)
          symObsNullChain(list(tapes[chainGroups[[chainOf[ci]]]]), nLeaves, nStates,
                          zSlotsL, as.integer(ptc %% p), p, as.integer(Nt),
                          as.integer(Mtot), 1L)
        else
          symObsNullMulti(list(tapes[[ci]]), nLeaves, nStates, zSlotsL,
                          as.integer(ptc %% p), p, as.integer(Nt), 1L)
      }
      kcall4 <- function(point, p, Nt, Mtot = 0L, solveFn = jointSolveCond) {
        blocks <- list()
        for (mi in seq_len(Kc)) {
          ci <- equilConds[mi]
          sc0 <- solveFn(point, p, ci)
          if (is.null(sc0)) { ssWhy <<- "joint solve failed"; return(list(ok = FALSE)) }
          b <- oneCondBlocks(mi, condObs(ci, sc0$ptc, p, Nt, Mtot), sc0, p)
          if (is.null(b)) return(list(ok = FALSE))
          blocks <- c(blocks, b)
        }
        # final stacked reduction over GF(p): the compiled kernel (symRrefMod)
        # runs per accepted sample here; .sym_rref_modp is the identical R fallback
        rr <- symRrefMod(do.call(rbind, blocks), p)
        list(ok = TRUE, R = rr$R,
             pivots = as.integer(rr$piv), rank = as.integer(rr$rank), dim = nzWide)
      }
      # Batched twin of the serial parMap(kcall4) loop for the coupled + gap path: in
      # joint mode symObsNullChain runs once per (point, condition), single-threaded --
      # the dominant reconstruction cost, run serially on Windows (parMap is serial
      # there). Solve every (point, condition) seed (warmSolves fills the cache in
      # parallel), evaluate ALL chain kernels in one OpenMP batch over the seeds, then
      # assemble + reduce each point in R (cheap). Byte-identical to looping kcall4.
      if (hasGaps) {
        chainsList <- lapply(seq_len(Kc), function(mi)
          tapes[chainGroups[[chainOf[equilConds[mi]]]]])
        # CPU-bound batch over evals: one OpenMP call over the pre-solved seeds (no
        # per-point fork here), so it uses the full `cores` budget, capped per call.
        coresChunk <- cores
        kchunk <- function(pointList, primeVec, Nt) {
          warmSolves(pointList, primeVec)
          nP <- length(pointList)
          perCond <- vector("list", nP)   # per point: list of Kc sc0, or NULL if any fails
          seedRows <- list(); evChain <- integer(0); evPrime <- numeric(0)
          for (i in seq_len(nP)) {
            pt <- pointList[[i]]; pp <- primeVec[[i]]
            per <- vector("list", Kc); okAll <- TRUE
            for (mi in seq_len(Kc)) {
              sc0 <- jointSolveCond(pt, pp, equilConds[mi])
              if (is.null(sc0)) { okAll <- FALSE; break }
              per[[mi]] <- sc0
            }
            if (!okAll) next
            perCond[[i]] <- per
            for (mi in seq_len(Kc)) {
              seedRows[[length(seedRows) + 1L]] <- as.integer(per[[mi]]$ptc %% pp)
              evChain <- c(evChain, mi - 1L); evPrime <- c(evPrime, pp)
            }
          }
          kr <- if (length(seedRows))
            symObsNullChainSeedBatch(chainsList, as.integer(evChain),
              do.call(rbind, seedRows), as.numeric(evPrime), nLeaves, nStates,
              zSlotsL, as.integer(Nt), as.integer(MtotUsed),
              as.integer(min(coresChunk, length(evChain))))
            else list()
          out <- vector("list", nP); e <- 0L
          for (i in seq_len(nP)) {
            if (is.null(perCond[[i]])) { out[[i]] <- list(ok = FALSE); next }
            blocks <- list(); bad <- FALSE
            for (mi in seq_len(Kc)) {
              e <- e + 1L
              b <- oneCondBlocks(mi, kr[[e]], perCond[[i]][[mi]], primeVec[[i]])
              if (is.null(b)) bad <- TRUE else blocks <- c(blocks, b)
            }
            out[[i]] <- if (bad) list(ok = FALSE) else {
              rr <- symRrefMod(do.call(rbind, blocks), primeVec[[i]])
              list(ok = TRUE, R = rr$R, pivots = as.integer(rr$piv),
                   rank = as.integer(rr$rank), dim = nzWide) }
          }
          out
        }
      }
      # from here on the analysis runs over the wide per-condition coordinate space;
      # the kernel/df closures above keep the local names (nzL/znamesL/zSlotsL)
      nz <- nzWide; znames <- znamesWide; zSlots <- zSlotsWide
      multi$zStateNames <- zStateNamesWide
      jointStateBase <- stateBase        # original state names, for the scaling peel
      jointKc <- Kc
      # a joint direction is reported in parameter space; its per-condition state and
      # recast leaves are auxiliary. A probe that perturbs such a leaf shifts the
      # pivots and marks it relevant to every entry, inflating the coupling count, so
      # exclude them from the direction/entry relevance gate (the reconstruction still
      # uses every relevant leaf). The dense-fit threshold is lifted by the state count
      # so a small joint direction (few params plus its state/recast leaves) still
      # takes the dense fit rather than the weaker sparse path.
      auxLeaves <- which(leafNamesAug %in% setdiff(leafNames, paramNames))
      ctrl$relevanceCap <- ctrl$relevanceCap + nSt
    }
  } else if (length(multi$powerRecast)) {
    # transient recast: a free power/Hill exponent without equilibrate. E = base^exp
    # and L = log(base) are ordinary free-initial-value leaves (their own z-columns);
    # the observability codistribution is stacked with the linearised recast relations
    #   E:  xi_E - (exp*E/base) xi_base - (E*log base) xi_exp = 0
    #   L:  xi_L - (1/base) xi_base = 0
    # (the log-differential identity, exact at the generic sample point regardless of
    # whether E, L take their tied values there -- base, log(base) and base^exp are
    # algebraically independent), then reduced. The physical report drops E, L.
    recastTransient <- TRUE
    recast <- multi$powerRecast
    for (i in seq_along(recast)) recast[[i]]$inverted <- FALSE
    recastAtomNames <- as.character(multi$recastAtomNames)
    auxLeaves <- which(leafNamesAug %in% recastAtomNames)
    slotOfL <- function(nm) { i <- match(nm, leafNames); if (is.na(i)) NA_integer_ else i }
    colOfN  <- function(nm) { i <- match(nm, znames);    if (is.na(i)) NA_integer_ else i }
    recastRelT <- lapply(recast, function(rc) list(
      Ecol = colOfN(rc$E), baseCol = colOfN(rc$base), expCol = colOfN(rc$exp),
      Lcol = colOfN(rc$L), Eslot = slotOfL(rc$E), baseSlot = slotOfL(rc$base),
      expSlot = slotOfL(rc$exp), Lslot = slotOfL(rc$L)))
    kcall4 <- function(point, p, Nt, Mtot = 0L) {
      pt <- as.integer(point)
      o <- if (hasGaps)
        symObsNullChain(lapply(chainGroups, function(idx) tapes[idx]), nLeaves,
                        nStates, zSlots, pt, p, as.integer(Nt), as.integer(Mtot), coresCall)
        else symObsNullMulti(tapes, nLeaves, nStates, zSlots, pt, p, as.integer(Nt), coresCall)
      if (!isTRUE(o$ok)) return(list(ok = FALSE))
      oR <- matrix(as.numeric(o$R), nrow = o$rank, ncol = nz)
      rel <- list(); seenL <- integer(0)
      for (rc in recastRelT) {
        base0 <- as.numeric(point[rc$baseSlot]) %% p
        if (base0 == 0) return(list(ok = FALSE))          # degenerate point, resample
        invb  <- .sym_invmod(base0, p)
        E0    <- as.numeric(point[rc$Eslot])  %% p
        exp0  <- as.numeric(point[rc$expSlot]) %% p
        L0    <- as.numeric(point[rc$Lslot])  %% p
        r1 <- numeric(nz)                                 # E = base^exp
        r1[rc$Ecol]    <- 1
        r1[rc$baseCol] <- (p - .sym_mulmod(exp0, .sym_mulmod(E0, invb, p), p)) %% p
        r1[rc$expCol]  <- (p - .sym_mulmod(E0, L0, p)) %% p
        rel[[length(rel) + 1L]] <- r1
        if (!(rc$Lcol %in% seenL)) {                      # L = log(base), once per base
          seenL <- c(seenL, rc$Lcol)
          r2 <- numeric(nz)
          r2[rc$Lcol]    <- 1
          r2[rc$baseCol] <- (p - invb) %% p
          rel[[length(rel) + 1L]] <- r2
        }
      }
      rr <- symRrefMod(rbind(oR, do.call(rbind, rel)), p)
      list(ok = TRUE, R = rr$R, pivots = as.integer(rr$piv),
           rank = as.integer(rr$rank), dim = nz)
    }
  } else {
    kcall4 <- function(point, p, Nt, Mtot = 0L) {
      pt <- as.integer(point)
      if (hasGaps)
        symObsNullChain(lapply(chainGroups, function(idx) tapes[idx]), nLeaves,
                        nStates, zSlots, pt, p, as.integer(Nt), as.integer(Mtot), coresCall)
      else
        symObsNullMulti(tapes, nLeaves, nStates, zSlots, pt, p, as.integer(Nt), coresCall)
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
  # joint mode: the state columns are log-normalised in kcall4, so their z-value is
  # the (dimensionless) weight coordinate with unit value; set the base point's
  # state slots to 1 so the peel forms tangent = weight * 1 for a state column and
  # weight * paramvalue for a parameter column. This overwrites the resting-state
  # values in point0, so keep the ORIGINAL solved point (its modular steady state is
  # cached) for the verify guard, which re-evaluates the kernel at higher Lie order.
  point0Solved <- sc$point0
  if (jointSS)
    for (col in which(znames %in% as.character(multi$zStateNames)))
      sc$point0[zSlots[col] + 1L] <- 1
  # the reconstruction samples at the saturated gap order, baked into a 3-arg kcall
  MtotUsed <- if (is.null(sc$MtotUsed)) 0L else as.integer(sc$MtotUsed)
  kcall <- function(point, p, Nt) kcall4(point, p, Nt, MtotUsed)
  # forward-sampling kernel (joint mode only): builds the same [Obs ; df] kernel but with the
  # states/params CHOSEN off the point and a turnover-rate subset solved linearly (valid at
  # every prime), so the gauge-robust reconstruction can sample a direction whose entries
  # depend on the resting state on a SHARED slice across primes. `solveRates` is the fixed
  # transversal the caller determines (once, avoiding the direction's support).
  kcallFwd <- if (jointSS) function(point, p, Nt, solveRates)
      kcall4(point, p, Nt, MtotUsed, function(pt, pp, cc) jointSolveCondFwd(pt, pp, cc, solveRates))
    else NULL
  # batched solve over many (point, prime) pairs for the reconstruction's sample
  # bank. The single-segment, non-constraint path shares one tape set across all
  # points, so the whole batch runs in one OpenMP-over-points kernel call; the
  # constraint and gap paths need per-point seeding, so they fall back to looping
  # the per-call kcall (which still threads cores over conditions/segments).
  # transient recast must loop the per-call kcall too: the batch kernel builds only
  # the observability rows and would omit the stacked recast-relation rows.
  canBatch <- !ssConstraint && !hasGaps && !recastTransient
  kbatch <- function(pointList, primeVec, Nt) {
    if (canBatch && length(pointList) > 0L) {
      M <- do.call(rbind, lapply(pointList, function(pp)
        as.integer(pp[seq_len(nLeaves)])))
      symObsNullBatch(tapes, nLeaves, nStates, zSlots, M, as.numeric(primeVec),
                      as.integer(Nt), cores)
    } else if (!is.null(kchunk) && length(pointList) > 0L &&
               !nzchar(Sys.getenv("DMOD_SYM_NOCHUNK"))) {
      # coupled + gap path: batch every (point, condition) chain-kernel evaluation in
      # one OpenMP call over the pre-solved seeds instead of the serial per-point loop
      # (the dominant cost, previously serial on Windows). kchunk warms the solves.
      # DMOD_SYM_NOCHUNK forces the serial fallback (byte-identical cross-check).
      kchunk(pointList, primeVec, Nt)
    } else {
      # fill the coupled-solve cache for the whole batch in parallel (Windows pool),
      # then run the per-point kernel/reduce over the cached solves. When parMap forks
      # (unix), the fork provides the parallelism, so the inner kernel runs serial here
      # to avoid fork x thread oversubscription (restored on exit for the non-forked
      # phases). On Windows parMap is serial and the master kernel keeps the full budget.
      warmSolves(pointList, primeVec)
      if (coresGLp > 1L && .Platform$OS.type == "unix") {
        ccSaved <- coresCall; coresCall <<- 1L
        on.exit(coresCall <<- ccSaved, add = TRUE)
      }
      parMap(seq_along(pointList),
             function(i) kcall(pointList[[i]], primeVec[[i]], Nt))
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
  # joint mode reports in PARAMETER space (states are auxiliary coordinates): a full
  # wide rank means an empty joint nullspace, hence every parameter identifiable.
  # (The non-identifiable case reprojects dim/rank below after the nullspace is known.)
  if (jointSS) {
    physParams0 <- setdiff(znames, as.character(multi$zStateNames))
    result$dim <- length(physParams0)
    if (sc$rank == nz) result$rank <- length(physParams0)
  }
  # transient recast reports in the physical space (real states + parameters); the
  # recast atoms E = base^exp, L = log(base) are auxiliary. A full augmented rank ties
  # every atom, so the physical space is fully identifiable.
  if (recastTransient) {
    physCoords0 <- setdiff(znames, recastAtomNames)
    result$dim <- length(physCoords0)
    if (sc$rank == nz) result$rank <- length(physCoords0)
  }
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
    # baked (non-coordinate) leaves. In joint mode znames is the WIDE decorated set,
    # so measure against the original coordinate names (shared params/L + state bases)
    fixed <- if (jointSS)
      setdiff(leafNames, c(setdiff(znames, as.character(multi$zStateNames)), jointStateBase))
      else setdiff(leafNames, znames)
    # the recast atoms let the peel impose c_E = exp*c_base and recover the
    # parameter-weighted (Hill) scalings over Q(exp), instead of leaving them to the
    # expensive rational fit (they are exact and need no finite-field sampling)
    scalRes <- tryCatch(sd$scalingSymmetriesMulti(
      perCondModel = modelLines, perCondObs = obsLines,
      inputs = if (length(inputs)) inputs else NULL,
      fixed = if (length(fixed)) fixed else NULL,
      recast = if (length(recast)) recast else NULL), error = function(e) NULL)
    if (!is.null(scalRes)) {
      # a scaling from the (original-name) integer kernel applies its state weight to
      # every per-condition column of that state
      if (jointSS) scalRes <- .sym_joint_expand_scal(scalRes, jointStateBase, jointKc,
                                                     heldParamOf)
      peel <- .sym_peel_scalings(scalRes, znames, nz, sc$point0[zSlots + 1L], P, N,
                                 sd = sd)
      scaling <- peel$scaling; Bmat <- peel$Bmat
    }
  }

  scalCols <- ncol(Bmat)          # scaling tangents only (fixed before translations)

  # exact translation lattice for method = "translation": the constant (additive)
  # directions, extracted before reconstruction. Their tangents join Bmat so the
  # residual reconstruction skips them, like the scalings. A graceful skip (NULL,
  # e.g. no generic point on a coupled steady state) falls back to the ordinary
  # reconstruction, which still finds them.
  translations <- list()
  if (translationsOnly) {
    tp <- tryCatch(.sym_peel_translations(sc, kcall, freeCols, N, nz, P, ctrl, Bmat,
                                          znames),
                   error = function(e) NULL)
    if (!is.null(tp)) {
      translations <- tp$translations; Bmat <- tp$Bmat
      # method = "translation": return only the exact lattice, no rational fit of
      # the residual (general) directions
      if (translationsOnly) { result$nonIdentifiable <- translations; return(result) }
    }
  }

  # free directions not spanned by the scalings/translations are the residual ones
  residualFree <- integer(0)
  for (fc in freeCols) {
    bf <- .sym_null_residues(sc$ref, fc, P)
    if (!.sym_in_span(Bmat, bf, nz, P)) {
      residualFree <- c(residualFree, fc); Bmat <- cbind(Bmat, bf) }
  }
  scalRows <- if (scalCols > 0L) t(Bmat[, seq_len(scalCols), drop = FALSE])
              else matrix(0L, 0L, nz)

  if (isTRUE(closedForm)) {
    .t0 <- Sys.time()
    to <- if (is.null(ctrl$timeout)) Inf else ctrl$timeout
    ctrl$deadline <- if (is.finite(to)) .t0 + to else NULL
    .tlog <- if (nzchar(Sys.getenv("DMOD_SYM_TIMING")))
      function(msg) message(sprintf("[sym %6.1fs] %s",
                                    as.numeric(Sys.time() - .t0, units = "secs"), msg))
      else function(msg) invisible()
    # stand up the parallel solve pool (Windows: fork is unavailable, so the coupled
    # per-point solves that dominate the probe and sample bank are filled across a
    # PSOCK pool of Python interpreters). Only worth it when there are residual
    # directions to reconstruct and jointSS solves them per point.
    if (useSolveCluster && jointSS && length(residualFree)) {
      cl <- .sym_make_solve_cluster(coresGLp)
      if (!is.null(cl)) on.exit(tryCatch(parallel::stopCluster(cl),
                                         error = function(e) NULL), add = TRUE)
      .tlog(sprintf("solve pool: %s", if (is.null(cl)) "unavailable (serial)"
                    else sprintf("%d workers", coresGLp)))
    }
    .tlog(sprintf("start: %d residual direction(s)", length(residualFree)))
    # shared relevance probe: a single solve at a one-leaf perturbation yields the
    # residues of every direction, so the per-leaf scan is done once here rather
    # than once per direction (the dominant cost at scale).
    # a perturbation that shifts the pivot set yields an uncomparable probe, which
    # would mark the leaf relevant to every entry and inflate the fit past its caps;
    # retry with fresh random values until the pivots match before accepting the
    # pessimistic last candidate, so genuinely irrelevant leaves stay out of the fit
    # round-based over retries: each (leaf, attempt) draws a fixed pool index, so a
    # whole round's seeds are pre-solved in parallel (over the sample points) before probing.
    # Leaves that resolve drop out; the final result is gauge-certified, so the
    # changed sampling order does not affect it.
    probeNext <- sc$poolNext
    relProbe <- vector("list", nAug)
    pending <- seq_len(nAug)
    # the relevance probe only feeds the per-direction reconstruction; when there are
    # no residual directions to reconstruct (all identifiability directions came from
    # the scaling pass / base rank) relProbe is never consumed, so skip the probe --
    # its per-leaf solves + kernels otherwise dominate the run at scale
    if (length(residualFree)) for (att in seq_len(ctrl$probeRetries)) {
      if (!length(pending) || .sym_expired(ctrl)) break
      perts <- lapply(pending, function(li) {
        pert <- sc$point0
        pert[li] <- sc$pool(probeNext + (li - 1L) * ctrl$probeRetries + (att - 1L))
        pert
      })
      # the per-leaf solves are independent -> one batched kernel call. On the
      # plain path this is a single OpenMP-over-leaves symObsNullBatch (also on
      # Windows, where parMap is serial); the coupled/gap path falls back inside
      # kbatch to warmSolves (Windows pool) + per-point kcall exactly as before.
      cands <- kbatch(perts, rep(P, length(perts)), sc$NtUsed)
      resolved <- logical(length(pending))
      for (j in seq_along(pending)) {
        li <- pending[j]; pert <- perts[[j]]; cand <- cands[[j]]
        relProbe[[li]] <- list(
          rp = cand, zvals = if (length(zSlots)) pert[zSlots + 1L] else NULL,
          pertval = pert[li])
        if (!is.null(cand) && isTRUE(cand$ok) &&
            identical(as.integer(cand$pivots), as.integer(sc$pivots)))
          resolved[j] <- TRUE
      }
      pending <- pending[!resolved]
    }
    poolNext <- probeNext + nAug * ctrl$probeRetries
    .tlog("relevance probe done")

    # reconstruct one direction in a given gauge (free-column when residueFn is
    # NULL, decoupled-canonical otherwise): interpolate, certify against the
    # nullspace at a fresh prime, and on failure downgrade to support-only. Also
    # certifies that the reconstructed vector actually lies in the nullspace at a
    # fresh point (gauge-independent), which rejects a self-consistent but
    # base-point-contaminated canonical representative.
    zvals0 <- if (length(zSlots)) sc$point0[zSlots + 1L] else NULL
    # physical-parameter support columns (0-based) for the per-prime joint path: the
    # per-condition state / recast columns are auxiliary and reconstructed only
    # implicitly, so the per-prime fit skips them to keep each entry narrow.
    physColsPP <- if (jointSS)
      which(!(znames %in% as.character(multi$zStateNames))) - 1L else NULL
    reconstructOne <- function(fc, rfn, logCoords = FALSE, sharedBank = NULL,
                               fastOnly = FALSE, perPrime = FALSE) {
      if (.sym_expired(ctrl)) {
        v <- if (!is.null(rfn)) rfn(sc$ref, P, zvals0) else .sym_null_residues(sc$ref, fc, P)
        if (is.null(v)) v <- .sym_null_residues(sc$ref, fc, P)
        return(list(support = sort(znames[v != 0]), type = "general",
                    closedForm = FALSE,
                    reason = "reconstruction time budget (reconstControl(timeout=)) exceeded"))
      }
      # the equilibrate/joint (coupled steady-state) path uses per-prime
      # reconstruction: a random perturbation solves at only some primes, so points
      # are collected per prime and the entry fits lifted by CRT (the free-column
      # gauge only; the log/canon rescues keep the shared-bank interpolation)
      dir <- if (isTRUE(perPrime))
        .sym_interpolate_perprime(fc, sc$ref, sc$pivots, znames, zSlots,
                                  leafNamesAug, nAug, sc$point0, sc$pool, poolNext,
                                  sc$NtUsed, kcall, kbatch, spy, relProbe, ctrl,
                                  auxLeaves = auxLeaves, physCols = physColsPP)
        else
        .sym_interpolate_direction(fc, sc$ref, sc$pivots, znames, zSlots,
                                   leafNamesAug, nAug, sc$point0, sc$pool,
                                   poolNext, sc$NtUsed, kcall, spy, relProbe,
                                   rfn, ctrl, kbatch, sharedBank, fastOnly,
                                   auxLeaves = auxLeaves)
      poolNext <<- dir$poolNext
      e <- dir$entry
      # the per-prime reconstructor tags the direction's relevant leaves for the
      # strict verifier (which perturbs only those to find a genuine steady-state
      # point); strip the internal tag so it never leaks into the reported direction
      ppRel <- e$relevantLeaves; e$relevantLeaves <- NULL
      # log-gauge entries are eta = xi / z; turn them back into the tangent xi
      # before verifying, then certify in original coordinates (gauge-independent)
      if (logCoords && isTRUE(e$closedForm))
        e$vector <- .sym_logcoord_backsub(e$vector, spy)
      ok <- isTRUE(e$closedForm) && (
        if (isTRUE(perPrime))
          # coupled path: strict verification at a genuine steady-state point (the
          # lenient verifiers accept inconclusively when the fresh-point solve fails,
          # which would pass a spurious per-prime fit)
          .sym_verify_perprime(e, fc, znames, leafNamesAug, sc$point0, sc$NtUsed,
                               kcall, sc$pivots, sc$pool, poolNext, nz, sd,
                               relLeaves = ppRel)
        else if (logCoords)
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

    allClosed <- function(rec) all(vapply(rec, function(e) isTRUE(e$closedForm),
                                          logical(1)))

    # Peel the minimal-support cocircuits that close on their own via the cheap
    # log-coordinate monomial read-off (a parameter-weighted scaling, e.g. a Hill
    # exponent). Only the directions no narrow cocircuit spans reach the wide
    # free-column fit below, so a single genuinely loop-wide direction no longer
    # forces every simple one into the expensive gauge.
    peeled <- list()
    msPeel <- .sym_minsupport_gauge(residualFree, scalRows, P, nz, sc, freeCols,
                                    candCap = ctrl$minsupportCandCap)
    if (length(msPeel$anchors) && !is.null(msPeel$vectors)) {
      peelVec <- matrix(0L, nz, 0L)
      for (gi in seq_along(msPeel$anchors)) {
        if (is.null(msPeel$residueFns[[gi]])) next
        e <- reconstructOne(msPeel$anchors[gi], msPeel$residueFns[[gi]],
                            logCoords = TRUE, fastOnly = TRUE)
        if (isTRUE(e$closedForm)) {
          peeled[[length(peeled) + 1L]] <- e
          peelVec <- cbind(peelVec, msPeel$vectors[, gi])
        }
      }
      if (length(peeled)) {
        # keep a basis completion of the residual set independent of the peeled span
        keep <- integer(0); span <- peelVec
        for (fc in residualFree) {
          bf <- .sym_null_residues(sc$ref, fc, P)
          if (.sym_in_span(span, bf, nz, P)) next
          keep <- c(keep, fc); span <- cbind(span, bf)
        }
        residualFree <- keep
      }
    }
    .tlog(sprintf("peel done: %d peeled, %d remaining", length(peeled), length(residualFree)))

    # minimal-support gauge run over the remaining residual set: a weighted scaling
    # whose weight is a free parameter (e.g. a Hill exponent) is a sparse circuit
    # of the full nullspace, recovered in log coordinates. `fastOnly` reads it off
    # the base point and the shared probe as a monomial with no kernel sampling.
    ms <- .sym_minsupport_gauge(residualFree, scalRows, P, nz, sc, freeCols,
                                candCap = ctrl$minsupportCandCap)
    msTry <- function(fastOnly) {
      if (length(ms$anchors) != length(residualFree) ||
          all(vapply(ms$residueFns, is.null, logical(1)))) return(NULL)
      rec <- lapply(seq_along(ms$anchors), function(gi)
        reconstructOne(ms$anchors[gi], ms$residueFns[[gi]], logCoords = TRUE,
                       fastOnly = fastOnly))
      if (allClosed(rec)) rec else NULL
    }

    # try the no-sampling monomial read-off first: it closes the parameter-weighted
    # scalings without paying for the sampling gauges below, which would otherwise
    # interpolate their loop-spanning free-column representative (expensive, and
    # doomed for these directions).
    interp <- msTry(TRUE)
    .tlog(if (is.null(interp)) "ms fastOnly: no close, going dense"
          else "ms fastOnly: all closed")
    if (is.null(interp)) {
      # raw free-column gauge; exact for every model whose directions the gauge does
      # not entangle. One shared dense sample bank over the union of all directions'
      # relevant leaves serves every direction: the expensive kernel is evaluated once
      # per point and reused across directions.
      metas <- lapply(residualFree, function(fc)
        .sym_direction_relevance(fc, sc$ref, sc$pivots, zSlots, nAug, sc$point0,
                                 relProbe, NULL))
      needs <- vapply(metas, function(m) .sym_dense_need(m$relByEntry, ctrl)$maxNeed,
                      integer(1))
      denseDir <- which(needs > 0L)
      .tlog(sprintf("relevance/need: %d dense dir(s), needs=[%s]",
                    length(denseDir), paste(needs, collapse = ",")))
      bank <- NULL
      if (length(denseDir)) {
        unionLeaves <- sort(unique(unlist(lapply(metas[denseDir], function(m) m$relevant))))
        bk <- .sym_build_shared_bank(unionLeaves, max(needs[denseDir]), sc$point0,
                                     sc$pool, poolNext, kbatch, sc$pivots, sc$NtUsed,
                                     length(.symPrimes), ctrl)
        poolNext <- bk$poolNext
        if (isTRUE(bk$ok)) bank <- bk
        .tlog(sprintf("shared bank built (%d leaves, ok=%s)",
                      length(unionLeaves), isTRUE(bk$ok)))
      }
      # A coupled steady-state (equilibrate) model rarely fills the shared all-prime
      # bank: a random perturbation admits an interior modular steady state at only
      # some primes, so a point that solves at every prime at once is vanishingly rare
      # (the fast-bail keeps this cheap). Those directions are reconstructed PER PRIME
      # instead -- each prime fills its own solvable points and the entry fits are
      # lifted by CRT. The shared bank still serves the easy equilibrate models, so
      # per-prime is a strict fallback and never regresses them.
      bankMissing <- ssConstraint && length(denseDir) > 0L && is.null(bank)
      interp <- lapply(seq_along(residualFree), function(ii) {
        fc <- residualFree[ii]
        r <- if (bankMissing) reconstructOne(fc, NULL, perPrime = TRUE)
             else reconstructOne(fc, NULL, sharedBank = bank)
        if (!isTRUE(r$closedForm) && ssConstraint && !bankMissing)
          r <- reconstructOne(fc, NULL, perPrime = TRUE)
        .tlog(sprintf("dense dir %d/%d done (closed=%s)", ii, length(residualFree),
                      isTRUE(r$closedForm)))
        r
      })

      # Coupled/recast residual directions the pivot-pinned gauge cannot sample (its
      # relevance probe misses the recast leaf, since perturbing it shifts the pivots, so
      # the standard/per-prime path returns them support-only regardless of whether the
      # shared bank filled): reconstruct the WHOLE residual set with the forward path,
      # each direction on a distinct physical anchor. Adopt only when it closes EVERY
      # residual (each entry is certified base-independently inside), so a partial/failed
      # pass leaves the honest support-only result untouched.
      if (ssConstraint && !is.null(sd) && !is.null(kcallFwd) && length(recast) &&
          !allClosed(interp)) {
        # forward reconstructs the whole non-scaling residual (one anchor per direction).
        # The constant (scaling) directions are already peeled in log coordinates, and the
        # forward path -- built for leaf-DEPENDENT (x*-dependent) entries -- declines them
        # ("residual constant in every leaf"); it closes exactly the leaf-dependent
        # residuals. Adopt those as the per-column results (keeping the peel) when their
        # count matches the residual set, each certified base-independently inside.
        fwd <- tryCatch(
          .sym_perprime_forward_multi(residualFree, sc, kcall, kcallFwd, znames, zSlots,
            leafNamesAug, nz, scaling, as.character(multi$zStateNames),
            as.character(multi$paramNames), recast, sd, spy, ctrl, physColsPP,
            models, realStateNames, solveParamNames, solveHeld),
          error = function(e) NULL)
        fwdClosed <- if (is.null(fwd)) list()
                     else Filter(function(e) isTRUE(e$closedForm), fwd)
        if (length(fwdClosed) == length(residualFree) && length(fwdClosed)) {
          interp <- lapply(fwdClosed, function(e) {
            if (length(recast)) e$vector <- .sym_recast_backsub(e$vector, recast, sd)
            e })
          .tlog(sprintf("forward-multi closed %d residual direction(s)", length(fwdClosed)))
        }
      }

      # rescue: the free-column gauge may entangle a few genuinely independent
      # directions; retry each still-open one modulo the exact scalings, in the
      # decoupled-canonical gauge then in log coordinates (a parameter-weighted scaling
      # stays sparse there, back-substituted after). Per-direction, not all-or-nothing:
      # each closed form is certified against the nullspace at a fresh prime, so a
      # partial adoption is sound.
      rescueEach <- function(gauge, logCoords = FALSE) {
        if (length(gauge$anchors) != length(residualFree)) return(invisible())
        for (ii in which(!vapply(interp, function(e) isTRUE(e$closedForm), logical(1)))) {
          if (is.null(gauge$residueFns[[ii]]) || .sym_expired(ctrl)) next
          rec <- reconstructOne(gauge$anchors[ii], gauge$residueFns[[ii]],
                                logCoords = logCoords)
          if (isTRUE(rec$closedForm)) interp[[ii]] <<- rec
        }
        invisible()
      }
      # the canon/log/ms sampling rescues interpolate through the ALL-PRIME shared
      # bank, which cannot fill for a coupled steady-state model; skip them there (the
      # per-prime path above is the coupled-case route) so they never grind, and run
      # them only for the plain (bank-backed) free-column entanglement they are for.
      if (!bankMissing) {
        if (!allClosed(interp))
          rescueEach(.sym_canon_gauge(residualFree, scalRows, P, nz, sc))
        if (!allClosed(interp))
          rescueEach(.sym_logcoord_gauge(residualFree, scalRows, P, nz, sc, zvals0), TRUE)
        # final rescue: the minimal-support gauge with kernel sampling, for a pinned
        # entry that is a bounded-degree rational rather than a bare monomial.
        if (!allClosed(interp)) {
          .tlog("rescues exhausted, trying ms sampling")
          msrec <- msTry(FALSE)
          if (!is.null(msrec)) interp <- msrec
        }
      }
    }
    .tlog(sprintf("reconstruction done: %d/%d closed",
                  sum(vapply(interp, function(e) isTRUE(e$closedForm), logical(1))),
                  length(interp)))
    result$nonIdentifiable <- c(scaling, translations, peeled, interp)
  } else {
    support <- lapply(residualFree, function(fc) {
      v <- .sym_null_residues(sc$ref, fc, P)
      list(support = sort(znames[v != 0]), type = "general", closedForm = FALSE)
    })
    result$nonIdentifiable <- c(scaling, translations, support)
  }
  # joint mode: drop directions supported only on the auxiliary recast / state
  # coordinates with no physical-parameter content (e.g. the trivial log-coordinate
  # shift L = log(base), which is unobservable by construction). A reported joint
  # direction must move at least one genuine parameter.
  if (jointSS) {
    stateCoords <- as.character(multi$zStateNames)
    physParams <- setdiff(znames, stateCoords)
    keep <- vapply(result$nonIdentifiable, function(d)
      any(d$support %in% physParams), logical(1))
    result$nonIdentifiable <- result$nonIdentifiable[keep]
    # report dim/rank/identifiable in PARAMETER space (the physical question), not
    # the enlarged joint (x, theta) space: project the joint nullspace onto the
    # parameter coordinates and take its rank. This matches the eliminated path's
    # verdict (the state coordinates are determined by the parameters on the manifold).
    paramCols <- which(znames %in% physParams)
    Njoint <- .sym_nullspace_basis(sc$ref, setdiff(0:(nz - 1L), sc$pivots), P)
    paramNull <- if (length(paramCols) && ncol(Njoint))
      .sym_rref_modp(Njoint[paramCols, , drop = FALSE], P)$rank else 0L
    result$dim <- length(physParams)
    result$rank <- as.integer(length(physParams) - paramNull)
    result$identifiable <- (paramNull == 0L)
    # project each joint direction onto the parameters: the state coordinates are
    # determined by the parameters on the resting manifold (xi_x = dx . xi_theta),
    # so the parameter part fully specifies the direction and the states are
    # redundant in the reported closed form. Report the parameter components; keep
    # the (implicit) state components in $stateVector. If a parameter entry's VALUE
    # genuinely references a state, no pure-parameter closed form exists - then the
    # honest closed answer keeps that state symbol (the whole point of the implicit
    # path) and the direction is flagged jointForm = TRUE.
    result$nonIdentifiable <- lapply(result$nonIdentifiable, function(d) {
      d$support <- sort(intersect(d$support, physParams))
      if (!is.null(d$vector)) {
        isState <- names(d$vector) %in% stateCoords
        d$stateVector <- d$vector[isState]
        pv <- d$vector[!isState]
        refsState <- any(vapply(pv, function(e)
          length(intersect(getSymbols(as.character(e)), stateCoords)) > 0L, logical(1)))
        d$vector <- pv
        if (refsState) d$jointForm <- TRUE
      }
      d
    })
  }
  # transient recast: report in the physical space (real states + parameters), with
  # the recast atoms E = base^exp and L = log(base) as auxiliary coordinates. Unlike
  # the joint mode above, the real state initial values ARE physical coordinates (a
  # free-IC symmetry such as FB d/dFB is legitimate), so only the E, L atoms are
  # dropped. Their direction components are determined by the recast relation and
  # redundant; the reconstructed entry VALUES already had E, L back-substituted to
  # base^exp, log(base), so no auxiliary symbol survives in the reported closed form.
  if (recastTransient) {
    physCoords <- setdiff(znames, recastAtomNames)
    keep <- vapply(result$nonIdentifiable, function(d)
      any(d$support %in% physCoords), logical(1))
    result$nonIdentifiable <- result$nonIdentifiable[keep]
    physCols <- which(znames %in% physCoords)
    Nrec <- .sym_nullspace_basis(sc$ref, setdiff(0:(nz - 1L), sc$pivots), P)
    physNull <- if (length(physCols) && ncol(Nrec))
      .sym_rref_modp(Nrec[physCols, , drop = FALSE], P)$rank else 0L
    result$dim <- length(physCoords)
    result$rank <- as.integer(length(physCoords) - physNull)
    result$identifiable <- (physNull == 0L)
    result$nonIdentifiable <- lapply(result$nonIdentifiable, function(d) {
      d$support <- sort(intersect(d$support, physCoords))
      if (!is.null(d$vector)) {
        isAtom <- names(d$vector) %in% recastAtomNames
        d$recastVector <- d$vector[isAtom]
        d$vector <- d$vector[!isAtom]
      }
      d
    })
  }
  # Schwartz-Zippel saturation guard: re-saturate the SAME kernel with a wider plateau
  # (more consecutive non-growing Lie orders required) and check the rank does not climb
  # past the reported value. The plateau rule is the one step that can stop early and
  # silently over-report; a wider re-saturation catches an intermediate plateau it slipped
  # through. This is one saturation pass -- find a modular point, extend the Lie jet -- not
  # a second analysis (no peeling, no reconstruction). Only reached when non-identifiable.
  if (isTRUE(verify))
    result$verification <- tryCatch(
      .sym_sz_saturation_guard(kcall, point0Solved, sc$NtUsed, sc$rank),
      error = function(e) list(ok = NA, method = "saturation guard",
                               reason = conditionMessage(e)))
  result
}


# ---- reconstruction internals: rational-entry interpolation & sampling ---------------

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
  for (pj in seq_len(np)[-1]) {
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
    if (.sym_expired(ctrl)) return(NULL)
    target <- min(if (have == 0L) 16L else 2L * have, maxLen)
    for (pj in seq_len(np)) {
      p <- .symPrimes[pj]
      for (k in (have + 1L):target) {
        if (.sym_expired(ctrl)) return(NULL)
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
    if (.sym_expired(ctrl)) return(NULL)
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
    if (.sym_expired(ctrl)) return(NULL)
    target <- min(if (have == 0L) 16L else 2L * have, maxLen)
    for (pj in seq_len(np)) {
      p <- .symPrimes[pj]
      for (k in (have + 1L):target) {
        if (.sym_expired(ctrl)) return(NULL)
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
  # A pivot-shifted probe is an UNUSABLE sample, not proof the leaf is relevant.
  # Two gauges can read it optimistically (skip) and lean on the fresh-point
  # certification to reject any under-fit: the support-pinned gauge (fixed support)
  # and the free-column gauge (its reconstruction is verified against the nullspace
  # at a fresh prime, so a skipped genuinely-relevant leaf makes the fit inconsistent
  # or fails verification -> support-only, the same verdict a pessimistic reading
  # reaches via the caps, but WITHOUT the over-count a recast-induced pivot shift
  # inflicts on every parameter). A canonical/log residue gauge stays pessimistic.
  optimistic <- isTRUE(attr(residueFn, "pinnedSupport")) || is.null(residueFn)
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
      if (optimistic) next
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
                                   pivots, NtUsed, nP, ctrl = NULL) {
  nu <- length(union)
  U <- matrix(0L, 0L, nu); points <- list(); rps <- list()
  tries <- 0L
  .bd <- nzchar(Sys.getenv("DMOD_SYM_BANKDIAG"))
  .nSolveFail <- 0L; .nPivMismatch <- 0L; .nGood <- 0L
  while (length(points) < needPts && tries < 30L * needPts) {
    # early bail: a direction whose perturbations rarely admit a modular steady state
    # that solves at EVERY prime at once (the coupled equilibrate case) yields no
    # all-prime point. Detect that fast -- 48 candidates with none good means the
    # all-prime rate is < ~2%, so the shared bank cannot fill and the caller should
    # fall back to the per-prime path -- instead of grinding the full 30x budget.
    if (length(points) == 0L && tries >= min(2L * needPts, 48L)) {
      if (.bd) message(sprintf("[bankdiag/bail] union=%d needPts=%d tries=%d good=0 solveFail=%d pivMismatch=%d",
                               nu, needPts, tries, .nSolveFail, .nPivMismatch))
      break
    }
    if (!is.null(ctrl) && .sym_expired(ctrl)) {
      if (.bd) message(sprintf("[bankdiag/expired] union=%d needPts=%d tries=%d good=%d solveFail=%d pivMismatch=%d",
                               nu, needPts, tries, .nGood, .nSolveFail, .nPivMismatch))
      return(list(ok = FALSE, poolNext = poolNext))
    }
    chunk <- max(1L, min(needPts - length(points) + 5L, 30L * needPts - tries))
    # keep the probing chunks small until the first point lands, so the fast-bail
    # above triggers after ~2 small batches rather than one huge coupled-solve batch
    if (length(points) == 0L) chunk <- min(chunk, 24L)
    # under an active deadline, bound the batch so a single pooled seed-solve cannot
    # overrun it: the loop-top expiry check then fires within one small batch
    if (!is.null(ctrl) && !is.null(ctrl$deadline)) chunk <- min(chunk, 96L)
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
        if (!isTRUE(rp$ok)) { good <- FALSE; if (.bd) .nSolveFail <- .nSolveFail + 1L; break }
        if (!identical(as.integer(rp$pivots), as.integer(pivots))) {
          good <- FALSE; if (.bd) .nPivMismatch <- .nPivMismatch + 1L; break
        }
        rpc[[j]] <- rp
      }
      if (!good) next
      if (.bd) .nGood <- .nGood + 1L
      U <- rbind(U, uv[[ci]])
      points[[length(points) + 1L]] <- cand[[ci]]
      rps[[length(rps) + 1L]] <- rpc
    }
  }
  if (.bd) message(sprintf("[bankdiag] union=%d needPts=%d tries=%d good=%d solveFail=%d pivMismatch=%d",
                           nu, needPts, tries, .nGood, .nSolveFail, .nPivMismatch))
  if (length(points) < needPts) return(list(ok = FALSE, poolNext = poolNext))
  list(ok = TRUE, U = U, points = points, rps = rps, union = union,
       poolNext = poolNext)
}


# Read each entry of a support-pinned direction as a Laurent monomial coeff *
# prod(leaf^exp) over the entry's relevant leaves. The integer exponents come from
# the ratio between the base point and each leaf's (already pivot-matched) probe at
# one prime; the rational coefficient is lifted from the base point across primes.
# Returns the closed-form entry list, or NULL if any entry is not such a monomial
# (then the caller interpolates). No fresh kernel sampling beyond the base point.
.sym_pinned_monomials <- function(f, pivots, znames, leafNames, point0,
                                  NtUsed, kcall, nv, base_nv, supportCols,
                                  relByEntry, relProbe, spy, maxExp = 6L) {
  P <- .symPrimes[1]
  powmod <- function(b, e, p) {
    r <- 1; b <- b %% p
    while (e > 0) { if (e %% 2 == 1) r <- .sym_mulmod(r, b, p); e <- e %/% 2
                    if (e > 0) b <- .sym_mulmod(b, b, p) }
    r
  }
  # smallest a in [-maxExp, maxExp] with base^a == target (mod p), or NULL
  findExp <- function(base, target, p) {
    target <- target %% p
    if (target == 1) return(0L)
    acc <- 1
    for (a in seq_len(maxExp)) { acc <- .sym_mulmod(acc, base, p)
      if (acc == target) return(a) }
    acc <- 1; bi <- .sym_invmod(base, p)
    for (a in seq_len(maxExp)) { acc <- .sym_mulmod(acc, bi, p)
      if (acc == target) return(-a) }
    NULL
  }

  entries <- list()
  for (i in seq_along(supportCols)) {
    ci <- supportCols[i]; reli <- relByEntry[[i]]
    e0 <- base_nv[ci + 1L] %% P
    if (e0 == 0) return(NULL)
    exps <- integer(length(reli))
    for (j in seq_along(reli)) {
      l <- reli[j]; pr <- relProbe[[l]]
      if (is.null(pr$rp) || !isTRUE(pr$rp$ok) ||
          !identical(as.integer(pr$rp$pivots), as.integer(pivots))) return(NULL)
      ev <- nv(pr$rp, P, pr$zvals)
      if (is.null(ev)) return(NULL)
      bval <- as.numeric(point0[l]) %% P; pval <- as.numeric(pr$pertval) %% P
      if (bval == 0 || pval == 0) return(NULL)
      a <- findExp(.sym_mulmod(pval, .sym_invmod(bval, P), P),
                   .sym_mulmod(ev[ci + 1L] %% P, .sym_invmod(e0, P), P), P)
      if (is.null(a)) return(NULL)
      exps[j] <- a
    }
    # rational coefficient: entry(base) / prod(leaf^exp) at the base point, lifted
    # from every prime whose base solve stays pivot-consistent (at least the anchor
    # prime, whose residue is already in hand); a small scaling weight needs one.
    monoDenom <- function(pj) {
      d <- 1
      for (j in seq_along(reli)) {
        bv <- as.numeric(point0[reli[j]]) %% pj
        d <- .sym_mulmod(d, if (exps[j] >= 0) powmod(bv, exps[j], pj)
                            else powmod(.sym_invmod(bv, pj), -exps[j], pj), pj)
      }
      d
    }
    residP <- .sym_mulmod(e0, .sym_invmod(monoDenom(P), P), P)
    usePrimes <- P; resid <- residP
    for (pj in .symPrimes[-1]) {
      rp <- kcall(point0, pj, NtUsed)
      if (is.null(rp) || !isTRUE(rp$ok) ||
          !identical(as.integer(rp$pivots), as.integer(pivots))) next
      nvp <- nv(rp, pj)
      if (is.null(nvp)) next
      usePrimes <- c(usePrimes, pj)
      resid <- c(resid, .sym_mulmod(nvp[ci + 1L] %% pj, .sym_invmod(monoDenom(pj), pj), pj))
    }
    rec <- symRatRecon(matrix(as.integer(resid), 1L), as.integer(usePrimes))
    if (rec$den[1] == "0") return(NULL)
    coef <- if (rec$den[1] == "1") rec$num[1] else paste0(rec$num[1], "/(", rec$den[1], ")")
    # assemble coeff * prod(leaf^exp) as a Laurent monomial string
    num <- character(0); den <- character(0)
    for (j in seq_along(reli)) {
      v <- leafNames[reli[j]]; a <- exps[j]
      if (a == 0) next
      e <- abs(a)
      term <- if (e == 1) v else paste0(v, "^", e)
      if (a > 0) num <- c(num, term) else den <- c(den, term)
    }
    numS <- paste(c(if (coef != "1" || !length(num)) coef, num), collapse = "*")
    expr <- if (!length(den)) numS
            else paste0("(", numS, ")/(", paste(den, collapse = "*"), ")")
    entries[[znames[ci + 1L]]] <- if (is.null(spy)) expr else .sym_simplify(expr, spy)
  }
  entries[[znames[f + 1L]]] <- "1"
  list(support = sort(names(entries)), vector = entries, type = "general",
       closedForm = TRUE)
}


# Per-prime reconstruction of one direction, for the equilibrate/joint (coupled
# steady-state) path. A random multi-parameter perturbation admits an interior
# modular steady state at only SOME of the reconstruction primes -- the coupled
# variety is sparse over any single GF(p) -- so a point that solves at every prime
# at once is vanishingly rare and the shared all-prime bank cannot fill (verified:
# perturbing the direction's leaves solves at >= 1 of 4 primes ~70% of the time but
# at all 4 ~0%). Here each prime fills its OWN solvable, pivot-consistent points,
# every entry is fit as a rational function over that prime's points, and the
# per-prime coefficients are lifted to the rationals by Chinese remaindering (the
# fitted function is the same across primes; only the evaluation points differ, so
# the CRT of its per-prime coefficients is exact). Degree is escalated per entry,
# so a low-degree entry needs only a few points. A prime whose coupled solve never
# succeeds is dropped and the CRT runs over the rest. Returns the same shape as
# .sym_interpolate_direction. `f` is the free (anchor) column, fixed to 1.
.sym_interpolate_perprime <- function(f, ref, pivots, znames, zSlots, leafNames,
                                      nLeaves, point0, pool, poolNext, NtUsed,
                                      kcall, kbatch, spy, relProbe, ctrl,
                                      auxLeaves = integer(0), physCols = NULL) {
  nP <- length(.symPrimes)
  zvalsOf <- function(pt) if (is.null(zSlots)) NULL else pt[zSlots + 1L]
  nv <- function(rp, p) .sym_null_residues(rp, f, p)
  rel <- .sym_direction_relevance(f, ref, pivots, zSlots, nLeaves, point0,
                                  relProbe, NULL)
  base_nv <- rel$base_nv; supportCols <- rel$supportCols
  relByEntry <- rel$relByEntry
  # In joint mode reconstruct only the physical-parameter support columns; the wide
  # per-condition state / recast columns are auxiliary (the direction is reported in
  # parameter space and they are projected out downstream), and dropping them keeps
  # every entry as narrow as its own parameters. The anchor stays column `f` (whose
  # value is 1), so a parameter entry is xi_param normalised on that free column.
  if (!is.null(physCols)) {
    keep <- which(supportCols %in% physCols)
    supportCols <- supportCols[keep]; relByEntry <- relByEntry[keep]
  }
  relevant <- sort(unique(unlist(relByEntry)))
  fallback <- function(reason) list(poolNext = poolNext,
    entry = list(support = sort(znames[base_nv != 0]), type = "general",
                 closedForm = FALSE, reason = reason))
  timedOut <- "reconstruction time budget (reconstControl(timeout=)) exceeded"
  if (.sym_expired(ctrl)) return(fallback(timedOut))
  ex <- function(v) if (length(auxLeaves)) setdiff(v, auxLeaves) else v
  maxRelPhys <- max(0L, vapply(relByEntry, function(e) length(ex(e)), integer(1)))
  if (nzchar(Sys.getenv("DMOD_SYM_SPARSEDIAG")))
    message(sprintf("perprime f=%d: support=%d relevant=%d maxRelPhys=%d",
                    f, length(supportCols), length(relevant), maxRelPhys))
  if (length(ex(relevant)) > ctrl$relevanceCapDir || maxRelPhys > ctrl$relevanceCapSparse)
    return(fallback(sprintf("direction couples %d parameters (an entry up to %d)",
                            length(ex(relevant)), maxRelPhys)))

  # a constant entry (no relevant leaf) lifts from the base point -- but the base
  # point, like any point, solves at only SOME primes here, so lift from the primes
  # where it does and CRT over those (a small integer weight needs one; more primes
  # only widen the reconstructible height). Requires perprimeMinPrimes to guard a
  # too-short CRT, unless a single prime already gives an exact small rational.
  constEntry <- function(col) {
    vals <- integer(0); prs <- numeric(0)
    for (pj in .symPrimes) {
      rp <- kcall(point0, pj, NtUsed)
      if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots))) next
      v <- nv(rp, pj); if (is.null(v)) next
      vals <- c(vals, v[col + 1L]); prs <- c(prs, pj)
    }
    if (!length(prs)) return(NULL)
    rec <- symRatRecon(matrix(as.integer(vals), 1L), as.integer(prs))
    if (rec$den[1] == "0") return(NULL)
    if (rec$den[1] == "1") rec$num[1] else paste0(rec$num[1], "/", rec$den[1])
  }
  if (!length(relevant)) {
    entries <- list()
    for (i in seq_along(supportCols)) {
      e <- constEntry(supportCols[i])
      if (is.null(e)) return(fallback("a constant entry could not be lifted"))
      entries[[znames[supportCols[i] + 1L]]] <- e
    }
    entries[[znames[f + 1L]]] <- "1"
    return(list(poolNext = poolNext, entry = list(support = sort(names(entries)),
      vector = entries, type = "general", closedForm = TRUE,
      relevantLeaves = integer(0))))
  }

  # per-prime point banks over the relevant leaves, grown on demand
  nrel <- length(relevant)
  bankU <- lapply(seq_len(nP), function(.) matrix(0L, 0L, nrel))
  bankR <- lapply(seq_len(nP), function(.) matrix(0L, 0L, length(supportCols)))
  triesP <- integer(nP); deadP <- logical(nP)
  ensurePrime <- function(j, need) {
    p <- .symPrimes[j]
    while (nrow(bankU[[j]]) < need && !deadP[j]) {
      if (.sym_expired(ctrl)) break
      have <- nrow(bankU[[j]])
      chunk <- max(1L, min(need - have + 4L, 128L))
      uv <- lapply(seq_len(chunk), function(.) {
        u <- pool(poolNext + seq_len(nrel) - 1L); poolNext <<- poolNext + nrel; u })
      cand <- lapply(uv, function(u) { pt <- point0; pt[relevant] <- u; pt })
      res <- kbatch(cand, rep(p, chunk), NtUsed)
      for (ci in seq_len(chunk)) {
        rp <- res[[ci]]
        if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots))) next
        vv <- nv(rp, p); if (is.null(vv)) next
        bankU[[j]] <<- rbind(bankU[[j]], uv[[ci]])
        bankR[[j]] <<- rbind(bankR[[j]], vv[supportCols + 1L])
      }
      triesP[j] <<- triesP[j] + chunk
      # this prime never yields a solvable/pivot-consistent point: drop it
      if (nrow(bankU[[j]]) == 0L && triesP[j] >= 8L * need + 32L) { deadP[j] <<- TRUE; break }
      if (triesP[j] >= 80L * need + 256L) break
    }
    nrow(bankU[[j]]) >= need
  }

  entries <- list()
  for (i in seq_along(supportCols)) {
    reli <- relByEntry[[i]]
    if (!length(reli)) {
      e <- constEntry(supportCols[i])
      if (is.null(e)) return(fallback("a constant entry could not be lifted"))
      entries[[znames[supportCols[i] + 1L]]] <- e
      next
    }
    cols_i <- match(reli, relevant)
    fitted <- NULL
    for (degree in 0:ctrl$degreeCap) {
      if (.sym_expired(ctrl)) return(fallback(timedOut))
      mons <- .sym_mono_table(length(reli), degree)
      need <- 2L * nrow(mons) - 1L + ctrl$sampleSlack
      if (need > ctrl$perprimeCap) break
      live <- which(vapply(seq_len(nP), function(j) !deadP[j] && ensurePrime(j, need),
                           logical(1)))
      if (length(live) < ctrl$perprimeMinPrimes) next
      # fit num/den over each prime's OWN points. symFitRational normalises to its
      # own free coefficient, which the point set can pick differently per prime; the
      # num/den are only defined up to a common scale, so before the CRT re-normalise
      # every prime to ONE shared monomial (the first fit's free column, which lies in
      # the true support and so is nonzero at every prime). Then the coefficient
      # vectors are one consistent representative and the CRT is exact.
      raw <- vector("list", length(live)); refFree <- NULL; okAll <- TRUE
      for (jj in seq_along(live)) { j <- live[jj]
        su <- matrix(as.integer(bankU[[j]][seq_len(need), cols_i, drop = FALSE]), need)
        rv <- as.integer(bankR[[j]][seq_len(need), i])
        fit <- symFitRational(su, matrix(as.integer(mons), nrow(mons)), rv, .symPrimes[j])
        if (!identical(fit$status, "ok")) { okAll <- FALSE; break }
        raw[[jj]] <- as.numeric(fit$coeffs)
        if (is.null(refFree)) refFree <- fit$freeCol
      }
      if (!okAll) next
      coefRes <- matrix(0L, 2L * nrow(mons), length(live))
      for (jj in seq_along(live)) { j <- .symPrimes[live[jj]]
        d <- raw[[jj]][refFree + 1L] %% j
        if (d == 0) { okAll <- FALSE; break }              # refFree not in support here
        coefRes[, jj] <- as.integer(.sym_mulmod(raw[[jj]], .sym_invmod(d, j), j))
      }
      if (!okAll) next
      rec <- symRatRecon(coefRes, as.integer(.symPrimes[live]))
      if (any(rec$den == "0")) next
      nMon <- nrow(mons); numI <- seq_len(nMon); denI <- nMon + seq_len(nMon)
      fitted <- list(numN = rec$num[numI], numD = rec$den[numI],
                     denN = rec$num[denI], denD = rec$den[denI], mons = mons)
      if (nzchar(Sys.getenv("DMOD_SYM_SPARSEDIAG")))
        message(sprintf("  perprime entry %d/%d col=%s reli=%d: fit at degree %d (%d primes, %d pts/prime)",
                        i, length(supportCols), znames[supportCols[i] + 1L],
                        length(reli), degree, length(live), need))
      break
    }
    if (is.null(fitted)) {
      if (nzchar(Sys.getenv("DMOD_SYM_SPARSEDIAG")))
        message(sprintf("  perprime entry %d/%d col=%s reli=%d: FAILED (degreeCap/perprimeCap)",
                        i, length(supportCols), znames[supportCols[i] + 1L], length(reli)))
      return(fallback(paste("a joint entry could not be fit per-prime within",
                            "degreeCap/perprimeCap; raise them or the direction is",
                            "not a bounded-degree rational")))
    }
    vars <- leafNames[reli]
    numStr <- .sym_poly_string(fitted$numN, fitted$numD, fitted$mons, vars)
    denStr <- .sym_poly_string(fitted$denN, fitted$denD, fitted$mons, vars)
    expr <- if (denStr == "1") numStr else paste0("(", numStr, ")/(", denStr, ")")
    if (!is.null(spy)) expr <- .sym_simplify(expr, spy)
    entries[[znames[supportCols[i] + 1L]]] <- expr
  }
  entries[[znames[f + 1L]]] <- "1"
  list(poolNext = poolNext, entry = list(support = sort(names(entries)),
    vector = entries, type = "general", closedForm = TRUE,
    relevantLeaves = relevant))
}


.sym_interpolate_direction <- function(f, ref, pivots, znames, zSlots, leafNames,
                                       nLeaves, point0, pool, poolNext, NtUsed,
                                       kcall, spy, relProbe = NULL,
                                       residueFn = NULL, ctrl = reconstControl(),
                                       kbatch = NULL, sharedBank = NULL,
                                       fastOnly = FALSE, auxLeaves = integer(0)) {
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
  timedOut <- "reconstruction time budget (reconstControl(timeout=)) exceeded"
  if (.sym_expired(ctrl)) return(fallback(timedOut))

  # a constant entry (no relevant leaf) is reconstructed from the base point. The
  # base-point kernel is identical across support columns and primes, so evaluate
  # it once (batched over primes -> one OpenMP call on the plain path) and index
  # the cached residues per column.
  constBase <- NULL
  constEntry <- function(col) {
    if (is.null(constBase)) {
      kr <- kbatch(rep(list(point0), length(.symPrimes)), .symPrimes, NtUsed)
      constBase <<- lapply(seq_along(.symPrimes), function(k) {
        rp <- kr[[k]]
        if (!isTRUE(rp$ok) || !identical(as.integer(rp$pivots), as.integer(pivots)))
          return(NULL)
        nv(rp, .symPrimes[k])
      })
    }
    res <- vapply(seq_along(.symPrimes), function(k) {
      nvp <- constBase[[k]]
      if (is.null(nvp)) NA_real_ else nvp[col + 1L]
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

  # pinned-support fast path: a weighted scaling's entry is a Laurent monomial in
  # its relevant leaves (e.g. eta = -nhill), so read it off the base point and the
  # already-matched per-leaf probes instead of sampling the kernel afresh. This
  # avoids resampling a recast Hill exponent, whose generic modular values almost
  # never keep the pivots consistent. Falls through to interpolation if any entry
  # is not a bounded-degree monomial.
  if (isTRUE(attr(residueFn, "pinnedSupport")) && !is.null(relProbe)) {
    mono <- .sym_pinned_monomials(f, pivots, znames, leafNames, point0,
                                  NtUsed, kcall, nv, base_nv, supportCols,
                                  relByEntry, relProbe, spy)
    if (!is.null(mono)) return(list(poolNext = poolNext, entry = mono))
    # the cheap monomial read-off failed; under fastOnly do not fall through to
    # kernel sampling (the caller will try the sampling gauges instead)
    if (isTRUE(fastOnly))
      return(fallback("pinned entry is not a bounded-degree monomial"))
  }

  maxRel <- max(0L, vapply(relByEntry, length, integer(1)))
  # the relevance caps guard against expensive fits over many PARAMETERS. In the
  # joint/implicit mode the per-condition state and recast leaves (auxLeaves) inflate
  # the count although the direction is reported in parameter space, so exclude them
  # from the gate (the reconstruction itself still uses every relevant leaf).
  ex <- function(v) if (length(auxLeaves)) setdiff(v, auxLeaves) else v
  relPhys <- ex(relevant)
  maxRelPhys <- max(0L, vapply(relByEntry, function(e) length(ex(e)), integer(1)))
  if (nzchar(Sys.getenv("DMOD_SYM_SPARSEDIAG")))
    message(sprintf("interp f=%d: relevant=%d (phys %d) maxRel=%d (phys %d) supportCols=%d (capDir=%d capSparse=%d)",
                    f, length(relevant), length(relPhys), maxRel, maxRelPhys,
                    length(supportCols), ctrl$relevanceCapDir, ctrl$relevanceCapSparse))
  if (length(relPhys) > ctrl$relevanceCapDir || maxRelPhys > ctrl$relevanceCapSparse)
    return(fallback(sprintf(
      paste("direction couples %d parameters (an entry up to %d), above",
            "relevanceCapDir=%d / relevanceCapSparse=%d; raise the relevant cap"),
      length(relPhys), maxRelPhys, ctrl$relevanceCapDir, ctrl$relevanceCapSparse)))
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
        if (.sym_expired(ctrl)) return(fallback(timedOut))
        # early bail: relevant leaves that leave NO valid steady state under
        # perturbation (a transcendental exponent confound, where the held recast E is
        # inconsistent with the perturbed parameters) yield zero usable points. Give up
        # after a small budget rather than grinding the full 30x, so a doomed rescue
        # gauge falls to support-only fast instead of dominating the runtime.
        if (nrow(sampleU) == 0L && tries >= 2L * maxNeed)
          return(fallback("no valid steady state under perturbation of the relevant leaves"))
        chunk <- max(1L, min(maxNeed - nrow(sampleU) + ctrl$sampleSlack,
                             30L * maxNeed - tries))
        if (!is.null(ctrl$deadline)) chunk <- min(chunk, 96L)
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
      if (is.null(e) && .sym_expired(ctrl)) return(fallback(timedOut))
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
                                    equilibrate = FALSE,
                                    condSubs = NULL, condInitial = NULL) {
  # a per-condition trafo list (condSubs / condInitial) sets the condition count when
  # no grid is given, and must align with the grid rows when one is
  nGrid <- if (is.null(conditions)) 0L else nrow(as.data.frame(conditions))
  Kcond <- max(1L, length(condSubs), length(condInitial), nGrid)
  if (nGrid && length(condSubs) && length(condSubs) != nGrid)
    stop("symmetryDetection(): the per-condition `trafo` list length (",
         length(condSubs), ") must match the condition grid rows (", nGrid, ").",
         call. = FALSE)
  if (is.null(conditions)) conditions <- data.frame(row.names = as.character(seq_len(Kcond)))
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
    v <- if (val %in% cols) cell(ci, val) else val
    # also apply this condition's trafo-list substitutions, so an event value the
    # per-condition `trafo` bakes -- a known dose (init_TGFb -> 1), a switch level,
    # a knockdown rename -- is resolved here rather than left as a free symbol. This
    # is the event-value counterpart of the grid cell lookup; without it a dose the
    # trafo pins to a constant would surface as a spurious non-identifiability.
    if (length(condSubs) && ci <= length(condSubs) && length(condSubs[[ci]])) {
      sub <- condSubs[[ci]]
      v <- replaceSymbols(names(sub), unlist(sub), v)
    }
    v
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
    # a per-condition trafo augments/overrides the grid substitutions for this
    # condition (a parameter baked to a number or an expression, exactly like a cell)
    if (length(condSubs) && !is.null(condSubs[[k]]))
      for (nm in names(condSubs[[k]])) subsBase[[nm]] <- condSubs[[k]][[nm]]
    # this condition's initial conditions: the per-condition trafo list wins, else
    # the single shared `initial`
    initK <- if (length(condInitial)) condInitial[[k]] else initial
    initK <- if (is.null(initK)) NULL else as.eqnvec(initK)
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
          hasInit <- !is.null(initK) && X %in% names(initK)
          ops <- leIdx[evVar[leIdx] == X]
          if (!isForcing && !hasInit && !length(ops)) next
          cur <- if (isForcing) "0" else if (hasInit) as.character(initK[[X]]) else X
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


# the partial-derivative operator d/d<var>, printed with the glyph U+2202 built
# from an escape so this source file stays pure ASCII (a literal glyph would trip
# the R CMD check non-ASCII warning, which fails CI).
.sym_partial <- function(var) paste0("\u2202/\u2202", var)

# ---- result display: the shared print()/summary() renderer ---------------------------

# the generator component xi_i of a scaling from its weight w and coordinate:
# w = 1 -> "var", w = -1 -> "-var", else "w*var"
.sym_scaling_component <- function(w, var) {
  w <- as.character(w)
  if (w == "1") var
  else if (w == "-1") paste0("-", var)
  else paste0(w, "*", var)
}

# one signed generator term "<coef> d/d<var>": pull a leading '-' out for the
# join, elide a unit coefficient, and parenthesise a multi-term (sum) coefficient.
.sym_signed_term <- function(xi, var) {
  xi <- trimws(xi)
  neg <- startsWith(xi, "-")
  mag <- if (neg) trimws(sub("^-", "", xi)) else xi
  bare <- gsub("\\*\\*", "^", mag)                    # exponents are not sums
  needParen <- grepl("[+]", bare) || grepl(".[-]", bare)
  coef <- if (mag == "1") ""
          else paste0(if (needParen) paste0("(", mag, ")") else mag, " ")
  list(neg = neg, text = paste0(coef, .sym_partial(var)))
}

# join signed terms into "a - b + c"
.sym_join_generator <- function(signed) {
  if (!length(signed)) return("0")
  out <- paste0(if (signed[[1]]$neg) "-" else "", signed[[1]]$text)
  for (k in seq_along(signed)[-1])
    out <- paste0(out, if (signed[[k]]$neg) " - " else " + ", signed[[k]]$text)
  out
}

# the symmetry-class tag, carrying the degree for a polynomial direction
.sym_class_tag <- function(d) {
  if (isTRUE(d$type == "polynomial") && !is.null(d$degree))
    sprintf("polynomial (%d)", d$degree)
  else d$type
}

# one formatted line for a symmetry generator: the differential generator
# sum_i xi_i d/dz_i followed by its class tag. $generator always holds the
# components xi_i directly (a scaling's integer weights are expanded to
# xi_i = w_i z_i at the finalisation boundary), so the printer is class-agnostic.
.sym_direction_line <- function(d) {
  if (is.null(d$generator))
    return(paste0("[", d$type, ", support only] involves: ",
                  paste(d$support, collapse = ", ")))
  signed <- lapply(names(d$generator),
                   function(k) .sym_signed_term(as.character(d$generator[[k]]), k))
  tag <- .sym_class_tag(d)
  if (isTRUE(d$certified))  tag <- paste0(tag, ", certified")
  if (isFALSE(d$verified))  tag <- paste0(tag, ", unverified")
  paste0(.sym_join_generator(signed), "   [", tag, "]")
}

# the plural group heading for a symmetry type
.sym_type_label <- function(t) switch(t,
  scaling = "Scalings", affine = "Affine", translation = "Translations",
  polynomial = "Polynomial", general = "General",
  paste0(toupper(substring(t, 1, 1)), substring(t, 2)))

# a non-negative integer as Unicode subscript digits (U+2080..U+2089), built from
# escapes so this source stays pure ASCII. Each digit is a single display column,
# so nchar() (which the label padding relies on) stays correct.
.sym_subscript <- function(n) {
  sub <- c("\u2080", "\u2081", "\u2082", "\u2083", "\u2084",
           "\u2085", "\u2086", "\u2087", "\u2088", "\u2089")
  paste(sub[as.integer(strsplit(as.character(n), "")[[1]]) + 1L], collapse = "")
}

# left-justify to a display width. R's sprintf/formatC pad by BYTES, which the
# multibyte d/d glyph (U+2202, 3 bytes) breaks; nchar() counts display columns.
.sym_ljust <- function(x, w) paste0(x, strrep(" ", pmax(0L, w - nchar(x))))

# Render one generator sum_i xi_i d/dz_i as aligned text lines: `perRow` terms per
# line (so the coefficients stay readable), the d/d operators padded into columns
# so they sit directly under each other down the block. Returns a character vector
# of lines; the first carries `prefix` ("Xk = "), the rest its blank indent.
.sym_format_generator <- function(generator, prefix, perRow = 3L) {
  nm <- names(generator)
  cells <- character(0); ops <- character(0)
  for (i in seq_along(nm)) {
    comp <- trimws(as.character(generator[[i]]))
    neg  <- startsWith(comp, "-")
    mag  <- if (neg) trimws(sub("^-", "", comp)) else comp
    bare <- gsub("\\*\\*", "^", mag)                    # exponents are not sums
    coef <- if (mag == "1") ""
            else if (grepl("[+]", bare) || grepl(".[-]", bare)) paste0("(", mag, ")")
            else mag
    conn <- if (i == 1L) { if (neg) "-" else "" } else if (neg) "- " else "+ "
    cells <- c(cells, paste0(conn, coef))
    ops   <- c(ops, .sym_partial(nm[i]))
  }
  n <- length(cells)
  if (!n) return(paste0(prefix, "0"))
  cw <- ow <- integer(perRow)                           # per-column cell / op widths
  for (c in seq_len(perRow)) {
    idx <- if (c <= n) seq.int(c, n, by = perRow) else integer(0)
    if (length(idx)) { cw[c] <- max(nchar(cells[idx])); ow[c] <- max(nchar(ops[idx])) }
  }
  cont <- strrep(" ", nchar(prefix))
  vapply(seq_len(ceiling(n / perRow)), function(r) {
    parts <- character(0)
    for (c in seq_len(perRow)) {
      i <- (r - 1L) * perRow + c
      if (i > n) break
      parts <- c(parts, paste0(.sym_ljust(cells[i], cw[c]), " ", .sym_ljust(ops[i], ow[c])))
    }
    paste0(if (r == 1L) prefix else cont, trimws(paste(parts, collapse = "  "), "right"))
  }, character(1))
}

# Print the symmetries grouped by type (Scalings / Affine / ...), each generator
# labelled X1, X2, ... contiguously, with a trailing flag block for a degree /
# certificate / verification note. In verbose mode the finite transformation of a
# polynomial generator and the reason a closed form was missed are shown below it.
.sym_cat_generators <- function(object, verbose = FALSE) {
  syms <- object$symmetries
  if (!length(syms)) return(invisible())
  types <- vapply(syms, function(d) as.character(d$type), character(1))
  ord <- order(match(types, c("scaling", "affine", "translation", "polynomial", "general")),
               seq_along(syms))                          # contiguous numbering per group
  syms <- syms[ord]; types <- types[ord]
  labels <- vapply(seq_along(syms), function(i) paste0("X", .sym_subscript(i)),
                   character(1))
  lw <- max(nchar(labels))
  cur <- ""
  for (i in seq_along(syms)) {
    d <- syms[[i]]
    if (types[i] != cur) { cur <- types[i]; cat(.sym_type_label(cur), ":\n", sep = "") }
    prefix <- paste0("  ", .sym_ljust(labels[i], lw), " = ")
    if (is.null(d$generator)) {
      cat(prefix, paste(d$support, collapse = ", "), "   [support only]\n", sep = "")
      if (isTRUE(verbose) && !is.null(d$reason))
        cat(strrep(" ", nchar(prefix)), "reason: ", d$reason, "\n", sep = "")
      next
    }
    lines <- .sym_format_generator(d$generator, prefix)
    flag <- character(0)
    if (isTRUE(d$type == "polynomial") && !is.null(d$degree)) flag <- c(flag, sprintf("deg %d", d$degree))
    if (isTRUE(d$certified)) flag <- c(flag, "certified")
    if (isFALSE(d$verified)) flag <- c(flag, "unverified")
    if (length(flag)) lines[1] <- paste0(lines[1], "   [", paste(flag, collapse = ", "), "]")
    cat(lines, sep = "\n"); cat("\n")
    if (isTRUE(verbose) && !is.null(d$transformation)) {
      tr <- d$transformation
      cat(strrep(" ", nchar(prefix)), "transformation: ",
          paste(vapply(names(tr), function(k) paste0(k, " -> ", tr[[k]]), character(1)),
                collapse = ",  "), "\n", sep = "")
    }
  }
}

# "<n> <singular|plural>", picking the number-agreeing noun explicitly (nicer than
# a "(s)" suffix)
.sym_plural <- function(n, one, many) sprintf("%d %s", n, if (n == 1L) one else many)

# The result section (shared by print() and summary()): the one-line verdict, then
# the grouped generators. This is all print() shows -- summary() prints the header
# and computation block above it.
.sym_cat_result <- function(object, verbose = FALSE) {
  m <- object$method; isObs <- m %in% c("observability", "translation")
  n <- length(object$symmetries)
  if (isObs && isTRUE(object$identifiable)) {
    cat(sprintf("Result:  structurally locally identifiable (rank %d / %d)\n",
                object$rank, object$dim))
    return(invisible())
  }
  if (isObs) {
    dirs <- if (m == "translation") .sym_plural(n, "translation direction", "translation directions")
            else .sym_plural(n, "non-identifiable direction", "non-identifiable directions")
    pre <- if (m == "translation") "translation lattice, " else ""
    cat(sprintf("Result:  %srank %d / %d  --  %s\n\n", pre, object$rank, object$dim, dirs))
  } else if (m == "scaling") {
    cat(sprintf("Result:  %s (exact integer kernel)\n\n",
                .sym_plural(n, "scaling symmetry", "scaling symmetries")))
  } else if (n == 0L) {
    cat(sprintf("Result:  no Lie-symmetry generator found up to pMax=%s\n",
                format(object$info$settings$pMax)),
        "         (a non-exhaustive search; not a proof of identifiability)\n", sep = "")
    return(invisible())
  } else {
    cat(sprintf("Result:  %s\n\n",
                .sym_plural(n, "polynomial Lie-symmetry generator",
                            "polynomial Lie-symmetry generators")))
  }
  .sym_cat_generators(object, verbose)
}

# ---- finalisation: normalise every engine's raw result into the public object -----
# Turn one raw engine direction/generator into the public $symmetries element.
# The raw shapes differ (observability/scaling directions carry $vector; the
# polynomial engine carries $infinitesimals); both map to a single schema with
# $generator (the xi_i components) and, for a scaling, the integer $weights.
.sym_public_symmetry <- function(d) {
  gen <- if (!is.null(d$vector)) d$vector else d$infinitesimals
  weights <- NULL
  # observability/scaling engines encode a scaling as integer weights in $vector
  # (xi_i = w_i z_i); expand to explicit components and keep the weights too. The
  # polynomial engine already delivers xi_i directly ($infinitesimals), no weights.
  if (isTRUE(d$type == "scaling") && !is.null(d$vector) && is.null(d$infinitesimals)) {
    weights <- d$vector
    gen <- setNames(as.list(vapply(names(d$vector),
             function(k) .sym_scaling_component(d$vector[[k]], k), character(1))),
             names(d$vector))
  }
  structure(list(
    type           = sub("^Type: ", "", as.character(d$type)),
    generator      = gen,
    weights        = weights,
    degree         = d$degree,
    support        = if (!is.null(d$support)) d$support else names(gen),
    # a direction is explicit iff it carries a full generator (all components),
    # whatever produced it -- a modular reconstruction, an exact scaling, or a
    # symbolic solve; a support-only direction has none. Derived from the generator
    # itself, not the internal reconstruction flag, so it is engine-agnostic.
    explicit       = !is.null(gen),
    reason         = d$reason,
    certified      = isTRUE(d$certified),
    transformation = d$transformation,
    verified       = if (is.null(d$verified)) NA else isTRUE(d$verified)
  ), class = "symmetryGenerator")
}

# Assemble the class-"symmetryDetection" object every engine return funnels
# through. Top level carries only the verdict (method / identifiable / rank /
# dim / symmetries); everything about *how* it was computed lives in $info.
# `identifiable` is a real rank verdict only for observability/translation; the
# scaling and polynomial engines are non-exhaustive (they find only scalings, or
# only generators up to pMax), so there "nothing found" is no proof -> NA.
.sym_finalize <- function(raw, method, settings, call, elapsed = NA_real_) {
  isObs   <- method %in% c("observability", "translation")
  rawSyms <- if (isObs || method == "scaling") raw$nonIdentifiable else raw
  if (is.null(rawSyms)) rawSyms <- list()
  syms    <- lapply(rawSyms, .sym_public_symmetry)

  rank <- if (!is.null(raw$rank)) as.integer(raw$rank) else NA_integer_
  dim  <- if (!is.null(raw$dim))  as.integer(raw$dim)  else NA_integer_
  identifiable <- if (isObs) isTRUE(rank == dim) else NA

  engine <- switch(method,
    scaling    = "integer-kernel",
    polynomial = "lie-ansatz",
    if (identical(raw$engine, "symbolic")) "symbolic" else "modular")

  info <- list(
    engine       = engine,
    lieOrderUsed = raw$lieOrderUsed,
    gapOrderUsed = raw$gapOrderUsed,
    conditions   = raw$conditions,
    segments     = raw$segments,
    settings     = settings,
    elapsed      = elapsed,
    verification = raw$verification)

  structure(list(
    method       = method,
    identifiable = identifiable,
    rank         = rank,
    dim          = dim,
    symmetries   = syms,
    info         = info,
    call         = call), class = "symmetryDetection")
}

# ---- the single report renderer shared by print() and summary() -----------------
# a compact "k=v" join of the settings relevant to `method` (more when verbose)
.sym_settings_line <- function(s, method, verbose) {
  keys <- switch(method,
    observability = , translation =
      if (verbose) c("reduceCQ", "equilibrate", "reconstruct", "verify",
                     "symEngine", "certifyPoly", "degreeCap")
      else c("reduceCQ", "equilibrate", "reconstruct"),
    scaling    = if (verbose) "reduceCQ" else character(0),
    polynomial = if (verbose) c("ansatz", "pMax", "polyBackend")
                 else c("ansatz", "pMax"))
  keys <- keys[vapply(keys, function(k) !is.null(s[[k]]), logical(1))]
  if (!length(keys)) return(NULL)
  paste(vapply(keys, function(k) paste0(k, "=", format(s[[k]])), character(1)),
        collapse = ", ")
}

# the Schwartz-Zippel saturation-guard verdict, one line (+ detail when verbose)
.sym_guard_lines <- function(v, verbose) {
  if (!is.list(v)) return(NULL)
  head <- if (isTRUE(v$ok))
      sprintf("saturation guard: PASSED (%s)", v$reason)
    else if (isFALSE(v$ok))
      sprintf("saturation guard: FAILED -- %s (directions may be over-reported)", v$reason)
    else
      sprintf("saturation guard: inconclusive (%s)",
              if (!is.null(v$reason)) v$reason else "unavailable")
  out <- head
  if (isTRUE(verbose) && !is.null(v$kernelRank))
    out <- c(out, sprintf("  kernel rank %d (extended %d), checked to Lie order %d%s",
                          v$kernelRank, v$kernelRankExtended, v$ordersChecked,
                          if (!is.na(v$growAt)) sprintf(", grew at %d", v$growAt) else ""))
  out
}

.sym_report <- function(object, verbose = FALSE) {
  m    <- object$method
  info <- object$info
  s    <- info$settings
  isObs <- m %in% c("observability", "translation")
  bar  <- strrep("-", 60)

  engLabel <- switch(as.character(info$engine),
    modular          = "modular (GF(p) + CRT)",
    symbolic         = "symbolic (pure sympy)",
    `integer-kernel` = "integer kernel (exact)",
    `lie-ansatz`     = paste0("Lie ansatz (", if (!is.null(s$polyBackend)) s$polyBackend
                                              else "symengine", ")"),
    as.character(info$engine))

  cat(bar, "\n", sep = "")
  cat(sprintf("symmetryDetection  |  method: %s   engine: %s\n", m, engLabel))
  cat(bar, "\n", sep = "")

  # ---- computation report ----
  comp <- character(0)
  if (isObs && !is.null(info$lieOrderUsed))
    comp <- c(comp, sprintf("Lie order %d (gap order %d)",
                            info$lieOrderUsed, if (!is.null(info$gapOrderUsed)) info$gapOrderUsed else 0L))
  if (!is.null(info$conditions) && info$conditions > 1L)
    comp <- c(comp, paste0(.sym_plural(info$conditions, "condition", "conditions"), ", ",
                           .sym_plural(info$segments, "segment", "segments")))
  setline <- .sym_settings_line(s, m, verbose)
  if (!is.null(setline)) comp <- c(comp, paste0("settings: ", setline))
  comp <- c(comp, .sym_guard_lines(info$verification, verbose))
  if (!is.null(info$elapsed) && is.finite(info$elapsed) && info$elapsed >= 0.05)
    comp <- c(comp, sprintf("elapsed: %.1fs", info$elapsed))
  if (length(comp)) {
    cat("Computation:\n")
    for (l in comp) cat("  ", l, "\n", sep = "")
  }
  cat("\n")
  .sym_cat_result(object, verbose)
  invisible(object)
}

#' @export
# print() is deliberately terse: just the verdict and the grouped generators
# (the "Result" section), no header or computation block -- that is what
# summary() adds.
print.symmetryDetection <- function(x, ...) { .sym_cat_result(x); invisible(x) }

#' @export
summary.symmetryDetection <- function(object, verbose = FALSE, ...) .sym_report(object, verbose)
