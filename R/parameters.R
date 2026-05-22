## Functions to generate parameter transformation ----

#' Generate a parameter transformation function
#' 
#' @description
#' This function provides a unified interface for generating condition-specific
#' parameter transformations, as commonly required in ODE-based modeling workflows.
#' 
#' `P()` can operate in three modes:
#' 
#' - **Explicit mode** (`method = "explicit"`, see [Pexpl]):
#'   Inner parameters are directly computed from symbolic expressions.
#'
#' - **Implicit mode** (`method = "implicit"`, see [Pimpl]):
#'   Steady states found via nonlinear root finding.
#'
#' - **Equilibrate mode** (`method = "equilibrate"`, see [Pequil]):
#'   Steady states found via ODE integration to equilibrium.
#'   This is used automatically for [eqnlist] entries when no method is specified.
#'
#' @param trafo object of class [eqnvec], named character, [eqnlist], or list thereof.
#' @param parameters character vector
#' @param condition character, the condition for which the transformation is generated
#' @param compile logical, compile the function (see [CppODE::funCpp])
#' @param modelname character, see [CppODE::funCpp]
#' @param method character, one of \code{"explicit"}, \code{"implicit"}, or \code{"equilibrate"}.
#'   If \code{NULL} (default), auto-selects \code{"equilibrate"} for [eqnlist] entries 
#'   and \code{"explicit"} for all others.
#' @param cores Number of cores for the per-condition `mclapply()`.
#'   `NULL` (default) auto-detects via [detectFreeCores]; forced to 1 on
#'   Windows.
#' @param verbose Print out information during compilation
#' @param deriv Logical, attach first-order parameter sensitivities to the
#'   result. Default `TRUE`.
#' @param deriv2 Logical, additionally attach the second-order sensitivity
#'   `attr(., "deriv2")` of shape `[innerPar, outerPar, outerPar]`. Requires
#'   `deriv = TRUE`. Default `FALSE`.
#' @param ... Additional arguments passed to the underlying transformation function
#'   ([Pexpl], [Pimpl], or [Pequil]).
#'
#' @return
#' An object of class [parfn], representing the parameter transformation.
#'
#' @seealso
#' [Pexpl], [Pimpl], [Pequil], [parfn]
#'
#' @export
P <- function(trafo = NULL, parameters = NULL, condition = NULL,
              compile = FALSE, modelname = NULL, method = NULL,
              cores = NULL, verbose = FALSE,
              deriv = TRUE, deriv2 = FALSE, ...) {

  if (is.null(trafo)) return()
  if (isTRUE(deriv2) && !isTRUE(deriv))
    stop("P(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)

  # Wrap single trafo in named list
  if (!is.list(trafo) || inherits(trafo, "eqnlist") || inherits(trafo, "eqnvec")) {
    trafo_list <- list(trafo)
    names(trafo_list) <- condition
  } else { trafo_list <- trafo }

  # Resolve `cores`. mclapply forks aren't available on Windows, so we cap at
  # 1 there unconditionally. On POSIX, defer to detectFreeCores() unless the
  # caller passed an explicit value — and call detectFreeCores() at most once,
  # since it has visible side effects (warnings on Windows, SSH on remotes).
  on_windows <- Sys.info()[['sysname']] == "Windows"
  if (on_windows) {
    cores <- 1L
  } else if (is.null(cores)) {
    cores <- detectFreeCores()
  } else {
    cores <- min(detectFreeCores(), cores)
  }

  # Always do codegen-only inside mclapply; actual compile + dyn.load must
  # happen in the parent process below, because dyn.load in a forked worker
  # is lost when the fork exits.
  result <- Reduce("+", mclapply(seq_along(trafo_list), function(i) {

    tr <- trafo_list[[i]]
    cond <- names(trafo_list[i])

    # Auto-select method per entry
    m <- if (!is.null(method)) match.arg(method, c("explicit", "implicit", "equilibrate"))
    else if (inherits(tr, "eqnlist")) "equilibrate"
    else "explicit"

    switch(m,
           explicit    = Pexpl(as.eqnvec(tr), parameters = parameters, condition = cond, compile = FALSE,
                               modelname = modelname, verbose = verbose,
                               deriv = deriv, deriv2 = deriv2, ...),
           implicit    = Pimpl(trafo = as.eqnvec(tr), parameters = parameters, condition = cond,
                               compile = FALSE, modelname = modelname, verbose = verbose,
                               deriv = deriv, deriv2 = deriv2, ...),
           equilibrate = Pequil(trafo = tr, parameters = parameters, condition = cond,
                                compile = FALSE, modelname = modelname, verbose = verbose,
                                deriv = deriv, deriv2 = deriv2, ...)
    )

  }, mc.cores = cores))

  if (compile) compile(result, cores = cores, output = modelname, verbose = verbose)

  result
}


#' Parameter transformation (explicit, algebraic)
#'
#' Builds `p_inner = f(p_outer)` from symbolic expressions and returns a
#' [parfn] whose evaluation attaches the Jacobian (and optionally the
#' Hessian). Backed by [CppODE::funCpp] in either forward-mode AD or
#' SymPy mode.
#'
#' @param trafo Named character or [eqnvec]. Names are inner parameters,
#'   values are expressions in the outer parameters.
#' @param parameters Outer-parameter names. Defaults to `getSymbols(trafo)`.
#' @param attach.input If `TRUE`, append the outer inputs to the output.
#' @param condition Condition label.
#' @param compile,modelname,verbose Forwarded to [CppODE::funCpp].
#' @param deriv If `TRUE` (default), attach the Jacobian `attr(., "deriv")`
#'   of shape `[p, theta]`.
#' @param deriv2 If `TRUE`, attach the Hessian `attr(., "deriv2")` of
#'   shape `[p, theta, theta]`. Requires `deriv = TRUE`.
#' @param derivMode `"dual"` (forward-mode AD, needs `compile = TRUE`) or
#'   `"symbolic"` (SymPy + analytic chain rule).
#'
#' @return A [parfn].
#' @seealso [Pimpl], [Pequil], [P].
#' @importFrom CppODE funCpp
#' @export
Pexpl <- function(trafo, parameters = NULL, attach.input = FALSE, condition = NULL,
                  compile = FALSE, modelname = NULL, verbose = FALSE,
                  deriv = TRUE, deriv2 = FALSE, derivMode = c("dual", "symbolic")) {

  derivMode <- match.arg(derivMode)
  emit_d1   <- isTRUE(deriv)
  emit_d2   <- isTRUE(deriv2)
  if (emit_d2 && !emit_d1)
    stop("Pexpl(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)

  if (is.null(parameters)) {
    parameters <- getSymbols(trafo)
  } else {
    identity <- setNames(parameters[!(parameters %in% names(trafo))],
                         parameters[!(parameters %in% names(trafo))])
    trafo <- c(trafo, identity)
    parameters <- getSymbols(trafo)
  }

  if (is.null(modelname)) modelname <- "expl_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

  PEval <- suppressWarnings(CppODE::funCpp(
    unclass(trafo), variables = NULL, parameters = parameters, fixed = NULL,
    compile = compile, modelname = modelname, outdir = getwd(),
    verbose = verbose, convenient = FALSE, derivMode = derivMode,
    deriv = emit_d1, deriv2 = emit_d2))

  fun <- PEval$func; jac <- PEval$jac; hess <- PEval$hess; evaluate <- PEval$evaluate
  use_ad     <- derivMode == "dual"
  ad_symbol  <- paste0(modelname, "_eval_ad")
  ad2_symbol <- paste0(modelname, "_eval_ad2")

  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {

    if (deriv2 && !emit_d2)
      stop("Pexpl was built with deriv2 = FALSE; rebuild with deriv2 = TRUE.", call. = FALSE)
    if (!emit_d1) deriv <- FALSE  # constructor-level gate: no first-order
    if (deriv2 && !deriv) deriv <- TRUE

    p <- c(pars, fixed)
    ad_ok  <- use_ad && !is.null(evaluate) && is.loaded(ad_symbol)
    ad_ok2 <- ad_ok && emit_d2 && is.loaded(ad2_symbol)
    if (deriv2 && !ad_ok2 && use_ad)
      stop("Pexpl(deriv2 = TRUE) needs the compiled AD2 entry; rebuild with compile = TRUE.", call. = FALSE)

    Jac <- NULL; Hess <- NULL

    if (ad_ok && deriv) {
      ## AD path. The dual-mode entry reads `params` and `dP` positionally
      ## against the codegen order, so both must be restricted to (and
      ## reordered to) `parameters`.
      dP  <- attr(pars, "deriv")
      dP2 <- if (deriv2) attr(pars, "deriv2") else NULL
      if (is.null(dP)) {
        active <- setdiff(parameters, names(fixed))
        dP <- diag(length(active)); dimnames(dP) <- list(active, active)
      }
      out <- evaluate(NULL, p[parameters], dX = NULL, dP = dP, dX2 = NULL, dP2 = dP2,
                      deriv2 = deriv2, attach.input = attach.input,
                      fixed = intersect(names(fixed), parameters))
      pinnerVal <- out$y[1, ]
      if (!is.null(out$dy))
        Jac <- matrix(out$dy, dim(out$dy)[2], dim(out$dy)[3],
                      dimnames = list(dimnames(out$dy)[[2]], dimnames(out$dy)[[3]]))
      if (deriv2 && !is.null(out$d2y))
        Hess <- array(out$d2y, dim(out$d2y)[2:4], dimnames = dimnames(out$d2y)[2:4])
    } else {
      pinnerVal <- fun(NULL, p, attach.input = attach.input, fixed = names(fixed))[, ]
      if (deriv && !is.null(jac)) {
        Jac <- as.matrix(jac(NULL, p, attach.input = attach.input, fixed = names(fixed))[1, , ])
        dP  <- attr(pars, "deriv")
        if (!is.null(dP)) {
          Jac <- Jac %*% dP[colnames(Jac), , drop = FALSE]
          dimnames(Jac) <- list(names(pinnerVal), colnames(dP))
        }
      }
      if (deriv2) {
        if (is.null(hess))
          stop("Pexpl(deriv2 = TRUE) requires hess(); rebuild with deriv2 = TRUE.", call. = FALSE)
        H4 <- hess(NULL, p, dX = NULL, dP = attr(pars, "deriv"),
                   dX2 = NULL, dP2 = attr(pars, "deriv2"),
                   attach.input = attach.input, fixed = names(fixed))
        Hess <- array(H4, dim(H4)[2:4], dimnames = dimnames(H4)[2:4])
      }
    }

    if (any(is.nan(pinnerVal)))
      stop("Inner parameter(s) evaluate to NaN:\n\t",
           paste(names(pinnerVal)[is.nan(pinnerVal)], collapse = "\n\t"),
           ".\nLikely cause: division by zero or missing inputs.", call. = FALSE)

    Jac_keep  <- if (deriv  && !is.null(Jac))  Jac[rowSums(Jac != 0) > 0, , drop = FALSE] else FALSE
    Hess_keep <- if (deriv2 && !is.null(Hess))
      (if (is.matrix(Jac_keep)) Hess[rownames(Jac_keep), , , drop = FALSE] else Hess) else FALSE
    pinner <- as.parvec(pinnerVal, deriv = Jac_keep, deriv2 = Hess_keep)

    if (attach.input && !all(names(pars) %in% names(pinnerVal)))
      pinner <- c(pinner, as.parvec(pars[setdiff(names(pars), names(pinnerVal))],
                                    deriv  = if (deriv)  NULL else FALSE,
                                    deriv2 = if (deriv2) NULL else FALSE))
    pinner
  }

  attr(p2p, "equations")   <- as.eqnvec(trafo)
  attr(p2p, "parameters")  <- parameters
  attr(p2p, "modelname")   <- modelname
  attr(p2p, "compileInfo") <- collectCompileInfo(fun, jac, hess, evaluate)
  parfn(p2p, parameters, condition)
}


#' Substitute conserved quantities into ODE rates
#'
#' Each independent conservation law `sum_k c_k * x_k = total_i` removes one
#' degree of freedom from the system: exactly one species `x_e` is eliminated
#' and substituted as `x_e = total_i - (rest_of_CQ)` throughout `f`. Species
#' already eliminated by a previous CQ are skipped (a CQ with no remaining
#' eliminatable species is redundant and dropped). When the user lists
#' species in `parameters` they are kept as pass-through inputs and the
#' eliminator picks a different candidate from the same CQ when possible.
#'
#' @param smatrix Stoichiometric matrix, or `NULL` (then only warns if
#'   parameter-states still occur unsubstituted).
#' @param f,states,parameters Equation set, state names, user parameters.
#' @return `list(f, parameters, cq_info, elim_states)`.
#' @keywords internal
.detect_and_substitute_cq <- function(smatrix, f, states, parameters) {
  cq_info <- list(); elim_states <- character(0)

  if (!is.null(smatrix)) {
    cq <- conservedQuantities(smatrix)
    if (!is.null(cq) && nrow(cq) > 0) {
      if (is.null(parameters)) parameters <- character(0)
      ## Per-CQ species sets, so we can look ahead to species that future
      ## CQs still need to keep free.
      cq_sets <- lapply(seq_len(nrow(cq)), function(i)
        intersect(getSymbols(as.character(cq[i, 1])), states))
      substitutions <- list()
      for (i in seq_len(nrow(cq))) {
        cq_species <- cq_sets[[i]]
        cq_expr    <- as.character(cq[i, 1])
        cand       <- setdiff(cq_species, elim_states)
        if (length(cand) < 2L) next  # CQ already constrained or redundant
        ## Avoid picking a pivot that another CQ further down will also
        ## need to eliminate; that produces mutually-recursive
        ## substitutions whose textual expansion never terminates. After
        ## that filter, prefer species the user did NOT list as
        ## pass-throughs; alphabetic order is the last tie-breaker.
        future <- unique(unlist(cq_sets[(i + 1L):length(cq_sets)],
                                use.names = FALSE))
        safe <- setdiff(cand, future)
        pool <- if (length(safe)) safe else cand
        not_user <- setdiff(pool, parameters)
        e <- if (length(not_user)) sort(not_user)[1L] else sort(pool)[1L]
        total_name <- paste0("total_", i)
        recon_expr <- paste0("(", total_name, " - (", replaceSymbols(e, "0", cq_expr), "))")
        substitutions[[e]] <- recon_expr
        elim_states <- c(elim_states, e)
        parameters  <- union(parameters, c(e, total_name))
        cq_info[[length(cq_info) + 1L]] <- list(
          total_name = total_name, elim_state = e,
          recon_expr = setNames(recon_expr, e))
      }
      if (length(substitutions))
        f <- replaceSymbols(names(substitutions),
                            unname(unlist(substitutions)), f)
    }
  } else if (!is.null(parameters)) {
    ## Without an smatrix we cannot auto-detect CQs; warn if the user nominated
    ## a state as parameter while leaving it in the dependent equations.
    param_states <- intersect(parameters, states)
    if (length(param_states)) {
      remaining <- f[setdiff(states, param_states)]
      still <- if (length(remaining))
        param_states[vapply(param_states, function(ps) any(ps %in% getSymbols(remaining)), logical(1))]
        else character(0)
      if (length(still))
        warning("States in 'parameters' still appear in dependent equations: ",
                paste(still, collapse = ", "),
                ". Provide an eqnlist for automatic CQ substitution.", call. = FALSE)
    }
  }
  list(f = f, parameters = parameters, cq_info = cq_info, elim_states = elim_states)
}


#' Detect states that are structurally zero in steady state
#'
#' Three-layer structural test on the stoichiometric matrix of an `eqnlist`,
#' iterated until stable (matches the a-priori-zero detection in AlyssaPetit
#' v1.2):
#'
#' \enumerate{
#'   \item *Neg-only column*: the state's stoichiometric column has only
#'     non-positive entries and at least one strictly negative entry. The
#'     state has only outflux, so it must be zero at steady state.
#'   \item *Pos-only column with single-state feeder*: the state's column has
#'     only non-negative entries and at least one strictly positive entry.
#'     Each feeding flux must be zero in steady state; if any feeding flux's
#'     rate expression involves exactly one state, that state must be zero.
#'   \item *Sink cluster (LP)*: subsets of states whose combined mass leaks
#'     monotonically (mass-balance LP via `lpSolve::lp`, only when the
#'     `lpSolve` package is installed; reduces to layer 1 for singletons).
#'     Catches cases like `{TGFb, R1_TGFb, R1_TGFb_int, ...}` where every
#'     individual state's column looks balanced but the cluster as a whole
#'     degrades.
#' }
#'
#' For every state found to be zero, the column is removed from the
#' stoichiometric matrix and the state symbol is substituted by `"0"` in the
#' remaining rate expressions; reactions whose rate becomes structurally
#' zero after substitution are dropped entirely. Re-detection then runs on
#' the reduced system, since each removal can expose new zero-states.
#'
#' @param eqnlist_obj An object of class `eqnlist`.
#' @return List with elements:
#'   \describe{
#'     \item{`zero_states`}{Character vector of state names found to be zero.}
#'     \item{`eqnlist`}{Reduced `eqnlist` with zero-state columns and their
#'       trivial reactions removed and rates simplified.}
#'   }
#' @keywords internal
.zeroStatesFromSmatrix <- function(eqnlist_obj) {
  S0 <- eqnlist_obj$smatrix
  if (is.null(S0) || ncol(S0) == 0L || nrow(S0) == 0L)
    return(list(zero_states = character(0), eqnlist = eqnlist_obj))

  S <- suppressWarnings(matrix(as.numeric(S0), nrow = nrow(S0), ncol = ncol(S0),
                               dimnames = dimnames(S0)))
  S[is.na(S)] <- 0
  rates       <- eqnlist_obj$rates
  description <- eqnlist_obj$description
  reactionCompartment <- eqnlist_obj$reactionCompartment

  zero_states <- character(0)

  .is_struct_zero <- function(expr) {
    v <- tryCatch({
      syms <- getSymbols(expr)
      if (!length(syms)) eval(parse(text = expr))
      else eval(parse(text = replaceSymbols(syms, rep("1", length(syms)), expr)))
    }, error = function(e) NA_real_)
    isTRUE(v == 0)
  }

  ## Drop reactions whose rate is structurally zero a priori (e.g.
  ## `addReaction("", "X", "0")` patterns users add to lock in a
  ## compartment for X). Their stoichiometric +1 entries otherwise mask
  ## the structural sink-cluster LP downstream.
  drop_rate0 <- vapply(rates, .is_struct_zero, logical(1))
  if (any(drop_rate0)) {
    keep <- !drop_rate0
    S     <- S[keep, , drop = FALSE]
    rates <- rates[keep]
    description <- description[keep]
    if (!is.null(reactionCompartment))
      reactionCompartment <- reactionCompartment[keep]
  }

  .neg_col <- function(M) {
    ## A column with no strictly positive entry has no influx; the state
    ## can only stay constant or decrease. In steady state it must be zero.
    ## Includes all-zero columns (state cascaded out by a prior reduction).
    for (j in seq_len(ncol(M))) {
      col <- M[, j]
      if (!any(col > 0)) return(j)
    }
    NA_integer_
  }

  .pos_col_zero_state <- function(M, rates_chr) {
    cn <- colnames(M)
    for (j in seq_len(ncol(M))) {
      col <- M[, j]
      if (any(col < 0) || !any(col > 0)) next
      for (k in which(col > 0)) {
        in_rate <- intersect(getSymbols(rates_chr[k]), cn)
        if (length(in_rate) == 1L) return(in_rate)
      }
    }
    NA_character_
  }

  .sink_cluster <- function(M, eps = 1e-8, Mbig = 1e4) {
    if (!requireNamespace("lpSolve", quietly = TRUE)) return(integer(0))
    nF <- nrow(M); nS <- ncol(M)
    if (nF == 0L || nS == 0L) return(integer(0))
    c_obj <- colSums(M)
    id    <- diag(nS)
    for (i in seq_len(nS)) {
      lb <- rep(0, nS); ub <- rep(Mbig, nS); lb[i] <- 1; ub[i] <- 1
      res <- tryCatch(
        lpSolve::lp("min", c_obj,
                    rbind(M, id, id),
                    c(rep("<=", nF), rep(">=", nS), rep("<=", nS)),
                    c(rep(0,  nF), lb, ub)),
        error = function(e) NULL)
      if (!is.null(res) && res$status == 0 && res$objval < -eps)
        return(which(res$solution > eps))
    }
    integer(0)
  }

  .zero_out <- function(j) {
    state_name <- colnames(S)[j]
    zero_states <<- c(zero_states, state_name)
    keep_rxn <- rep(TRUE, nrow(S))
    for (k in seq_len(nrow(S))) {
      ## Two ways the state can enter a reaction: stoichiometrically (column
      ## entry != 0, e.g. an educt or product) or kinetically (appears in
      ## the rate expression, e.g. an enzyme/modifier). Both cases must be
      ## substituted; matching AlyssaPetit's F.subs(state, 0) check.
      in_stoich <- S[k, j] != 0
      in_rate   <- state_name %in% getSymbols(rates[k])
      if (!in_stoich && !in_rate) next
      new_rate <- if (in_rate) replaceSymbols(state_name, "0", rates[k]) else rates[k]
      if (.is_struct_zero(new_rate)) keep_rxn[k] <- FALSE
      else rates[k] <<- new_rate
    }
    S    <<- S[keep_rxn, -j, drop = FALSE]
    rates       <<- rates[keep_rxn]
    description <<- description[keep_rxn]
    if (!is.null(reactionCompartment))
      reactionCompartment <<- reactionCompartment[keep_rxn]
  }

  repeat {
    progressed <- FALSE
    repeat {
      step <- FALSE
      jneg <- .neg_col(S)
      if (!is.na(jneg)) { .zero_out(jneg); step <- TRUE; progressed <- TRUE; next }
      pzs <- .pos_col_zero_state(S, rates)
      if (!is.na(pzs)) {
        j <- match(pzs, colnames(S))
        if (!is.na(j)) { .zero_out(j); step <- TRUE; progressed <- TRUE; next }
      }
      if (!step) break
    }
    sink <- .sink_cluster(S)
    if (!length(sink)) break
    for (j in sort(sink, decreasing = TRUE)) { .zero_out(j); progressed <- TRUE }
    if (!progressed) break
  }

  if (!length(zero_states))
    return(list(zero_states = character(0), eqnlist = eqnlist_obj))

  S_out <- S
  S_out[S_out == 0] <- NA
  storage.mode(S_out) <- storage.mode(S0)

  new_obj <- eqnlist_obj
  new_obj$smatrix     <- S_out
  new_obj$states      <- colnames(S_out)
  new_obj$rates       <- rates
  new_obj$description <- description
  if (!is.null(reactionCompartment))
    new_obj$reactionCompartment <- reactionCompartment
  if (!is.null(new_obj$compartmentOf))
    new_obj$compartmentOf <- new_obj$compartmentOf[colnames(S_out)]

  list(zero_states = zero_states, eqnlist = new_obj)
}


#' Reset warm-start caches in `Pequil`/`Pimpl` parameter transformations
#'
#' Both [Pequil] and [Pimpl] cache the previous root (and its
#' sensitivities) as a warm start for the next call, gated by
#' `keep.root = TRUE`. The cache is keyed to the closure of a single
#' `parfn` and persists across calls; for workflows that cross basins of
#' attraction (`mstrust`, parameter-grid sweeps, repeated profile
#' likelihoods after a structural change) a stale cache can pin the
#' solver in the wrong region.
#'
#' `resetWarmStarts()` clears the cache on the supplied function and on
#' every `Pequil`/`Pimpl` parfn reachable through composition (e.g.
#' `obj` constructed via `normL2(data, g * x * p)` — calling
#' `resetWarmStarts(obj)` walks the closure environments and resets the
#' cache on `p` as well).
#'
#' @param fn A `parfn`, `prdfn`, `obsfn`, `objfn`, or composed `fn`.
#' @param verbose Print one-line summary of cleared caches. Default
#'   `TRUE`.
#' @return Invisibly, a character vector of labels for the cleared
#'   caches (one entry per Pequil/Pimpl parfn touched). Length zero if
#'   nothing was found.
#'
#' @export
resetWarmStarts <- function(fn, verbose = TRUE) {
  if (!is.function(fn))
    stop("`fn` must be a function (parfn / prdfn / obsfn / objfn / composed fn).",
         call. = FALSE)

  visited_envs   <- new.env(parent = emptyenv())
  invoked_resets <- new.env(parent = emptyenv())
  labels         <- character(0)

  call_reset <- function(r) {
    key <- format(environment(r))
    if (key %in% names(invoked_resets)) return()
    assign(key, TRUE, envir = invoked_resets)
    new_labels <- tryCatch(r(),
                           error = function(e) sprintf("<reset error: %s>",
                                                       conditionMessage(e)))
    labels <<- c(labels, as.character(new_labels))
  }

  walk <- function(x) {
    env <- NULL
    if (is.function(x)) {
      r <- attr(x, "resetWarmStart")
      if (!is.null(r) && is.function(r)) call_reset(r)
      env <- environment(x)
    } else if (is.environment(x)) {
      env <- x
    } else if (is.list(x)) {
      for (v in x) walk(v)
      return()
    } else {
      return()
    }
    if (is.null(env)) return()
    key <- format(env)
    if (key %in% names(visited_envs)) return()
    assign(key, TRUE, envir = visited_envs)
    for (nm in ls(env, all.names = TRUE)) {
      val <- tryCatch(get(nm, envir = env, inherits = FALSE),
                      error = function(e) NULL)
      if (is.function(val) || is.environment(val) || is.list(val)) walk(val)
    }
  }
  walk(fn)

  if (isTRUE(verbose)) {
    if (length(labels))
      message("resetWarmStarts: cleared ", length(labels),
              " warm-start cache(s):\n  ",
              paste(labels, collapse = "\n  "))
    else
      message("resetWarmStarts: no warm-start caches found.")
  }
  invisible(labels)
}


#' Solve `A %*% X = B` with SVD pseudoinverse fallback
#'
#' Used by [Pimpl] to invert `df/dx` at the steady state. Well-conditioned
#' systems take a standard LU solve; ill-conditioned or rank-deficient
#' systems fall back to the Moore-Penrose pseudoinverse (minimum-norm
#' least-squares solution, equal to the IFT sensitivity restricted to the
#' constraint manifold), with a warning listing the null-space
#' direction(s) in `row_names` coordinates so the missing conserved
#' quantity or redundant equation can be identified.
#'
#' @keywords internal
.pimpl_solve_dfdx <- function(A, B, row_names = rownames(A), warn_rcond = 1e-10) {
  if (nrow(A) == 0L) return(B)
  sv  <- svd(A); d <- sv$d
  tol <- max(dim(A)) * d[1L] * .Machine$double.eps
  rnk <- sum(d > tol)
  rc  <- if (d[1L] > 0) d[length(d)] / d[1L] else 0

  if (rnk == nrow(A) && (is.finite(rc) && rc >= warn_rcond))
    return(solve(A, B))

  ## Pseudoinverse path. Warn once, listing the rank-deficient direction(s).
  if (rnk < nrow(A)) {
    nd <- sv$v[, (rnk + 1L):ncol(sv$v), drop = FALSE]
    rownames(nd) <- row_names
    warning(.pimpl_format_singularity(d, rnk, nd), call. = FALSE)
  } else {
    warning(sprintf("df/dx is ill-conditioned (rcond = %.2e); using SVD pseudoinverse.", rc),
            call. = FALSE)
  }
  inv_d <- ifelse(d > tol, 1 / d, 0)
  X <- sv$v %*% (inv_d * crossprod(sv$u, B))
  dimnames(X) <- list(row_names, colnames(B))
  X
}

#' @keywords internal
.pimpl_format_singularity <- function(d, rnk, null_dirs) {
  lines <- vapply(seq_len(ncol(null_dirs)), function(j) {
    v   <- null_dirs[, j]
    sig <- abs(v) > 0.05 * max(abs(v)); if (!any(sig)) sig <- abs(v) == max(abs(v))
    paste0("  null #", j, ": ",
           paste(sprintf("%+.3f*%s", v[sig], rownames(null_dirs)[sig]), collapse = " "))
  }, character(1))
  paste0(
    sprintf("df/dx is rank-deficient (rank %d of %d; smallest sv %.2e of largest %.2e); ",
            rnk, length(d), d[length(d)], d[1L]),
    "using SVD pseudoinverse (minimum-norm sensitivity on the constraint manifold).\n",
    "Null-space direction(s) in dependent-state coordinates:\n",
    paste(lines, collapse = "\n"),
    "\nLikely cause: an unmodelled conserved quantity, redundant equations, or a ",
    "continuum of steady states. Pass an `eqnlist` so `Pimpl` can auto-detect CQs ",
    "for cleaner sensitivities."
  )
}


#' Parameter transformation (implicit, root-finding)
#'
#' Solves `f(x, p) = 0` for the dependent states `x` via [nleqslv::nleqslv]
#' (multistart on failure), then returns a [parfn] over the outer inputs
#' carrying the IFT-derived Jacobian (and Hessian, with `deriv2 = TRUE`).
#' For [eqnlist] inputs, conserved quantities are detected and eliminated
#' species are replaced by `total_*` parameters; their values are
#' reconstructed from the solved root.
#'
#' @param trafo Named character / [eqnvec] / [eqnlist].
#' @param parameters Outer parameters. Auto-extended with state names that
#'   need to be eliminated for CQ substitution.
#' @param forcings Forcing names; replaced by 0 and removed from the system.
#' @param condition Condition label.
#' @param keep.root Cache the root as warm-start for the next call.
#' @param positive Solve in log-space to enforce positivity.
#' @param compile,modelname,verbose Forwarded to [CppODE::funCpp].
#' @param deriv If `TRUE` (default), attach first-order parameter
#'   sensitivities `attr(., "deriv")` via the implicit function theorem.
#' @param deriv2 Emit `attr(., "deriv2")` via the implicit function theorem
#'   (closed form: one extra `df/dx`-solve per parameter pair). Requires
#'   `funCpp` to expose `hess()`. Same shape as in [Pexpl] / [Pequil].
#' @param controlsNleqslv Solver-control overrides (merged into the
#'   defaults below).
#'
#' @return A [parfn].
#' @seealso [Pexpl], [Pequil], [P].
#' @export
#' @import nleqslv
#' @importFrom digest digest
Pimpl <- function(trafo, parameters = NULL, forcings = NULL, condition = NULL,
                  keep.root = TRUE, positive = TRUE, compile = FALSE,
                  modelname = NULL, verbose = FALSE, deriv = TRUE, deriv2 = FALSE,
                  controlsNleqslv = list()) {

  emit_d1 <- isTRUE(deriv)
  emit_d2 <- isTRUE(deriv2)
  if (emit_d2 && !emit_d1)
    stop("Pimpl(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)
  ## Newton iteration needs the analytical Jacobian regardless of
  ## emit_d1; the construct-time deriv flag only gates the *output*
  ## IFT chain-rule (so a deriv = FALSE parfn skips that work).

  smatrix <- NULL
  zero_states <- character(0)
  if (inherits(trafo, "eqnlist")) {
    if (!is.null(forcings))
      trafo$rates <- replaceSymbols(forcings, rep("0", length(forcings)), trafo$rates)
    zs <- .zeroStatesFromSmatrix(trafo)
    zero_states <- zs$zero_states
    trafo <- zs$eqnlist
    smatrix <- trafo$smatrix
  }
  trafo  <- as.eqnvec(trafo)
  states <- names(trafo)

  if (!is.null(forcings)) {
    if (is.null(smatrix))
      trafo <- replaceSymbols(forcings, rep("0", length(forcings)), trafo)
    trafo  <- trafo[setdiff(names(trafo), forcings)]
    states <- names(trafo)
  }

  ## States with rhs == "0" become parameters (no equation to solve).
  const_states <- states[vapply(unclass(trafo), function(x)
    tryCatch(identical(eval(parse(text = x)), 0), error = function(e) FALSE), logical(1))]
  if (length(const_states)) parameters <- union(parameters %||% character(0), const_states)

  cq <- .detect_and_substitute_cq(smatrix, trafo, states, parameters)
  trafo <- cq$f; parameters <- cq$parameters
  cq_info <- cq$cq_info; elim_states <- cq$elim_states

  dependent <- setdiff(states, parameters)
  if (!length(dependent))
    stop("No dependent states to solve for. All states appear in 'parameters'.", call. = FALSE)

  parameters <- if (is.null(parameters))
    getSymbols(trafo, exclude = dependent)
  else union(getSymbols(trafo, exclude = dependent), parameters)
  parms_all  <- setdiff(parameters, dependent)
  n_dep      <- length(dependent)

  if (is.null(modelname)) modelname <- "impl_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

  ## Combined evaluator: dep residuals + recon expressions for eliminated states.
  all_exprs <- unclass(trafo[dependent])
  if (length(cq_info))
    all_exprs <- c(all_exprs, setNames(
      vapply(cq_info, function(ci) ci$recon_expr, character(1)),
      vapply(cq_info, function(ci) ci$elim_state, character(1))))
  n_all <- length(all_exprs)

  ## Eliminated states still live in `parameters` (so `total_*` ride along)
  ## but no longer occur in any expression; CppODE's symbolic deriv2 path
  ## errors on unused parameters, so we strip them here.
  parms_all <- intersect(parms_all, getSymbols(all_exprs))

  PEval <- suppressWarnings(CppODE::funCpp(
    all_exprs, variables = dependent, parameters = parms_all, fixed = NULL,
    compile = compile, modelname = modelname, outdir = getwd(),
    verbose = verbose, convenient = FALSE,
    deriv = TRUE, deriv2 = emit_d2, derivMode = "symbolic"))

  jac_cols <- c(dependent, parms_all)

  ## funCpp returns [1, n_all, ...] arrays; we squeeze the leading obs axis.
  X <- function(x) matrix(x[dependent], 1, dimnames = list(NULL, dependent))
  eval_f <- function(x, p) {
    F <- PEval$func(X(x), p[parms_all])
    setNames(as.numeric(F[1, ]), dimnames(F)[[2]])[seq_len(n_dep)]
  }
  eval_J_full <- function(x, p) {
    J <- PEval$jac(X(x), p[parms_all])
    matrix(c(J), dim(J)[2], dim(J)[3], dimnames = list(dimnames(J)[[2]], dimnames(J)[[3]]))
  }
  eval_J <- function(x, p) eval_J_full(x, p)[seq_len(n_dep), , drop = FALSE]
  eval_H <- function(x, p) {
    if (is.null(PEval$hess)) return(NULL)
    H4 <- PEval$hess(X(x), p[parms_all])
    array(c(H4), dim(H4)[2:4], dimnames = dimnames(H4)[2:4])
  }

  cache <- new.env(parent = emptyenv())
  cache$guess <- NULL; cache$failed_pv_hash <- NULL; cache$failed_result <- NULL

  controls <- modifyList(list(
    keep.root = keep.root, positive = positive, nstarts = 100,
    ## Biological steady-state magnitudes routinely span ~10 orders of
    ## magnitude when rate constants vary by 2-3 orders themselves. Default
    ## to log-uniform sampling over [1e-5, 1e5] when positive = TRUE so
    ## Newton can land in the basin even for stiff systems.
    lower = if (positive) 1e-5 else 0,
    upper = if (positive) 1e5  else 100,
    debugPlot = FALSE,
    method = "Newton", global = "dbldog", xscalm = "fixed",
    xtol = 1e-4, ftol = 1e-2, btol = 1e-3, cndtol = 1e-12,
    maxit = 200L, allowSingular = TRUE),
    controlsNleqslv)
  pimpl_keys   <- c("keep.root","positive","nstarts","lower","upper","debugPlot")
  nleqslv_args <- c("method","global","xscalm")

  expand_bounds <- function(b, dep, default_val) {
    if (is.null(names(b)) || length(b) == 1L)
      return(setNames(rep(b[1L], length(dep)), dep))
    out <- setNames(rep(default_val, length(dep)), dep)
    nm  <- intersect(names(b), dep); out[nm] <- b[nm]; out
  }

  ## One nleqslv attempt from x0 (in log-space iff `positive`).
  solve_once <- function(x0, pv, positive, top, ctrl) {
    if (positive) {
      s0 <- x0; s0[s0 <= 0] <- 1
      sol <- nleqslv::nleqslv(
        x   = log(s0),
        fn  = function(lx) { names(lx) <- dependent; eval_f(exp(lx), pv) },
        jac = function(lx) {
          names(lx) <- dependent; xv <- exp(lx)
          eval_J(xv, pv)[, dependent, drop = FALSE] %*% diag(xv, n_dep)
        },
        method = top$method, global = top$global, xscalm = top$xscalm, control = ctrl)
      root <- setNames(exp(sol$x), dependent)
    } else {
      sol <- nleqslv::nleqslv(
        x   = x0,
        fn  = function(x) { names(x) <- dependent; eval_f(x, pv) },
        jac = function(x) { names(x) <- dependent; eval_J(x, pv)[, dependent, drop = FALSE] },
        method = top$method, global = top$global, xscalm = top$xscalm, control = ctrl)
      root <- setNames(sol$x, dependent)
    }
    res <- tryCatch(eval_f(root, pv), error = function(e) setNames(rep(Inf, n_dep), dependent))
    list(root = root, sol = sol, res = res,
         maxres = max(abs(res)), termcd = sol$termcd, iter = sol$iter)
  }

  ## Diagnostic plot for multistart termination codes (debugPlot = TRUE).
  waterfall_plot <- function(log_df, ftol) {
    lbl <- c("1" = "converged", "2" = "xtol (f may be large)", "3" = "stalled",
             "4" = "maxit exceeded", "5" = "ill-conditioned",
             "6" = "singular", "7" = "unusable Jacobian")
    col <- c("converged" = "#2ca02c", "xtol (f may be large)" = "#ff7f0e",
             "stalled" = "#d62728", "maxit exceeded" = "#9467bd",
             "ill-conditioned" = "#8c564b", "singular" = "#e377c2",
             "unusable Jacobian" = "#7f7f7f")
    log_df$termcd_label <- factor(lbl[as.character(log_df$termcd)], levels = lbl)
    print(ggplot2::ggplot(log_df, ggplot2::aes(x = factor(rank), y = maxres, fill = termcd_label)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::geom_hline(yintercept = ftol, linetype = "dashed", color = "steelblue", linewidth = 0.6) +
      ggplot2::annotate("text", x = nrow(log_df), y = ftol, label = paste0("ftol = ", ftol),
                        hjust = 1, vjust = -0.5, color = "steelblue", size = 3) +
      ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-12),
                                  breaks = c(0, 10^(-10:4)),
                                  labels = function(x) ifelse(x == 0, "0", scales::scientific(x))) +
      ggplot2::scale_fill_manual(values = col, name = "termination", drop = TRUE) +
      ggplot2::labs(x = "index (sorted by residual)", y = "max |f(x)|",
                    title = "Pimpl multistart diagnostics") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom",
                     axis.text.x = if (nrow(log_df) > 30) ggplot2::element_blank()
                                   else ggplot2::element_text(size = 7)))
  }

  ## IFT-based Jacobian and (optionally) Hessian over `input_cols`. The
  ## chain-rule contractions are:
  ##   dx*/dp_a       = -(df/dx)^{-1}  df/dp_a
  ##   d2 x*/dp_a dpb = -(df/dx)^{-1} [ f_xx(dx_a, dx_b)
  ##                                  + f_xp_b dx_a + f_xp_a dx_b
  ##                                  + f_pp_{ab} ]
  ## Eliminated species recon: x_e = R(x*, p), so dx_e = R_x dx + R_p, and
  ## d2 x_e correspondingly with an extra R_x . d2 x term. Incoming dP, dP2
  ## chain-rule are applied last. Returns rows trimmed to non-zero entries.
  build_derivs <- function(root, pv, out, p, emptypars, fixed, dP, dP2,
                           deriv, want_d2) {
    Jfull <- eval_J(root, pv)
    dfdx  <- Jfull[, dependent, drop = FALSE]
    dfdp  <- Jfull[, parms_all, drop = FALSE]

    dxdp <- if (length(parms_all))
      .pimpl_solve_dfdx(dfdx, -dfdp, row_names = dependent)
    else matrix(numeric(0), n_dep, 0, dimnames = list(dependent, character(0)))

    input_cols <- setdiff(names(p), c(dependent, names(fixed)))
    n_in <- length(input_cols); n_par <- length(parms_all)
    par_input <- intersect(parms_all, input_cols)

    jacobian <- matrix(0, length(out), n_in, dimnames = list(names(out), input_cols))
    ep <- intersect(emptypars, input_cols)
    if (length(ep)) jacobian[cbind(ep, ep)] <- 1
    cd <- intersect(colnames(dxdp), input_cols)
    if (length(cd)) jacobian[rownames(dxdp), cd] <- dxdp[, cd, drop = FALSE]

    recon_jac <- NULL
    if (length(elim_states)) {
      Jall      <- eval_J_full(root, pv)
      recon_jac <- Jall[(n_dep + 1L):n_all, , drop = FALSE]
      rownames(recon_jac) <- elim_states
      dep_jac  <- jacobian[dependent, input_cols, drop = FALSE]
      elim_jac <- recon_jac[, dependent, drop = FALSE] %*% dep_jac
      if (length(par_input))
        elim_jac[, par_input] <- elim_jac[, par_input, drop = FALSE] +
          recon_jac[, parms_all, drop = FALSE][, par_input, drop = FALSE]
      jacobian[elim_states, ] <- elim_jac
    }

    hessian <- NULL
    if (want_d2) {
      H_all <- eval_H(root, pv)
      if (is.null(H_all))
        stop("Pimpl(deriv2 = TRUE) requires hess(); rebuild with deriv2 = TRUE.", call. = FALSE)

      f_xx <- H_all[seq_len(n_dep), dependent, dependent, drop = FALSE]
      f_xp <- H_all[seq_len(n_dep), dependent, parms_all, drop = FALSE]
      f_pp <- H_all[seq_len(n_dep), parms_all, parms_all, drop = FALSE]

      ## IFT RHS: T1 (f_xx sandwich by dxdp) + sym(T2 from f_xp) + f_pp,
      ## all batched over the residual equation index.
      T1 <- if (n_par > 0L) t(dxdp) %bmm% f_xx %bmm% dxdp
            else array(0, c(n_dep, 0L, 0L))
      T2 <- if (n_par > 0L) t(dxdp) %bmm% f_xp
            else array(0, c(n_dep, 0L, 0L))
      RHS_k <- T1 + T2 + aperm(T2, c(1L, 3L, 2L)) + f_pp
      dimnames(RHS_k) <- list(dependent, parms_all, parms_all)

      d2xdp2 <- if (n_par > 0L)
        array(.pimpl_solve_dfdx(dfdx, -matrix(RHS_k, n_dep, n_par * n_par),
                                row_names = dependent),
              c(n_dep, n_par, n_par),
              dimnames = list(dependent, parms_all, parms_all))
        else array(0, c(n_dep, 0L, 0L))

      hess_arr <- array(0, c(length(out), n_in, n_in),
                        dimnames = list(names(out), input_cols, input_cols))
      if (length(par_input))
        hess_arr[dependent, par_input, par_input] <-
          d2xdp2[dependent, par_input, par_input, drop = FALSE]

      if (length(elim_states) && n_par > 0L) {
        r_xx <- H_all[(n_dep + 1L):n_all, dependent, dependent, drop = FALSE]
        r_xp <- H_all[(n_dep + 1L):n_all, dependent, parms_all, drop = FALSE]
        r_pp <- H_all[(n_dep + 1L):n_all, parms_all, parms_all, drop = FALSE]
        rj_x <- recon_jac[, dependent, drop = FALSE]

        A1 <- t(dxdp) %bmm% r_xx %bmm% dxdp
        A2 <- t(dxdp) %bmm% r_xp
        ## rj_x %*% d2xdp2 contracted over the dependent axis (= batch of d2xdp2)
        contract <- array(rj_x %*% matrix(d2xdp2, n_dep, n_par * n_par),
                          c(length(elim_states), n_par, n_par))
        d2elim <- A1 + A2 + aperm(A2, c(1L, 3L, 2L)) + r_pp + contract
        dimnames(d2elim) <- list(elim_states, parms_all, parms_all)
        if (length(par_input))
          hess_arr[elim_states, par_input, par_input] <-
            d2elim[, par_input, par_input, drop = FALSE]
      }
      hessian <- hess_arr
    }

    if (!is.null(dP)) {
      dPsub <- submatrix(dP, rows = colnames(jacobian))
      th    <- colnames(dPsub); n_th <- length(th)
      if (want_d2 && !is.null(hessian)) {
        ## Chain rule: H_new[i] = dP^T H[i] dP, batched -> two %bmm% calls.
        new_hess <- t(dPsub) %bmm% hessian %bmm% dPsub
        dimnames(new_hess) <- list(names(out), th, th)
        if (!is.null(dP2)) {
          dP2sub <- dP2[input_cols, th, th, drop = FALSE]
          new_hess <- new_hess + array(
            jacobian %*% matrix(dP2sub, n_in, n_th * n_th),
            c(length(out), n_th, n_th),
            dimnames = list(names(out), th, th))
        }
        hessian <- new_hess
      }
      jacobian <- jacobian %*% dPsub
    }

    keep <- rowSums(jacobian != 0) > 0
    jacobian <- jacobian[keep, , drop = FALSE]
    if (!is.null(hessian)) hessian <- hessian[keep, , , drop = FALSE]
    list(jacobian = jacobian, hessian = hessian)
  }


  ## Reconstruct eliminated species at `root`. Returns named numeric.
  reconstruct <- function(root, pv) {
    if (!length(elim_states)) return(numeric(0))
    Fv <- PEval$func(X(root), pv[parms_all])
    setNames(as.numeric(Fv[1, ])[(n_dep + 1L):n_all], elim_states)
  }

  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
    if (deriv2 && !emit_d2)
      stop("Pimpl was built with deriv2 = FALSE; rebuild with deriv2 = TRUE.", call. = FALSE)
    if (!emit_d1) deriv <- FALSE
    if (deriv2 && !deriv) deriv <- TRUE

    p   <- pars
    dP  <- attr(p, "deriv")
    dP2 <- if (deriv2) attr(p, "deriv2") else NULL

    top  <- controls[intersect(names(controls), nleqslv_args)]
    ctrl <- controls[setdiff(names(controls), c(pimpl_keys, nleqslv_args))]
    keep.root <- controls$keep.root; positive <- controls$positive
    nstarts   <- controls$nstarts;   ftol     <- controls$ftol
    debugPlot <- controls$debugPlot

    if (!is.null(fixed)) {
      p <- p[!names(p) %in% names(fixed)]
      p <- c(p, fixed)
    }
    emptypars <- setdiff(names(p), c(dependent, names(fixed)))
    miss <- setdiff(dependent, names(p)); if (length(miss)) p[miss] <- 1
    pv <- p[parms_all]
    pv_hash <- digest::digest(pv, algo = "xxhash64")

    zero_vec <- if (length(zero_states))
      setNames(rep(0, length(zero_states)), zero_states) else NULL

    ## Fast path: previous multistart already exhausted on this pv.
    if (!is.null(cache$failed_pv_hash) && identical(pv_hash, cache$failed_pv_hash)) {
      warning("Multistart previously failed for these parameters. ",
              "Returning cached (inaccurate) solution.", call. = FALSE)
      root <- cache$failed_result$root
      out  <- c(root, reconstruct(root, pv), zero_vec,
                p[setdiff(names(p), c(names(root), elim_states, zero_states))])
      d    <- tryCatch(build_derivs(root, pv, out, p, emptypars, fixed, dP, dP2, deriv, deriv2),
                       error = function(e) NULL)
      return(as.parvec(out,
                       deriv  = if (deriv  && !is.null(d)) d$jacobian else NULL,
                       deriv2 = if (deriv2 && !is.null(d)) d$hessian else if (deriv2) NULL else FALSE))
    }

    ## Multistart: cache -> user init -> random sweep, possibly cycling
    ## through alternative nleqslv methods.
    log <- list(); best <- NULL
    record <- function(r, label) {
      if (is.null(r)) return()
      log[[length(log) + 1L]] <<- list(label = label, maxres = r$maxres,
                                       termcd = r$termcd, iter = r$iter)
      if (is.null(best) || r$maxres < best$maxres) best <<- r
    }
    try_solve <- function(x0, label, override = NULL)
      record(tryCatch(solve_once(x0, pv, positive,
                                 if (is.null(override)) top else override, ctrl),
                      error = function(e) NULL), label)

    if (!is.null(cache$guess)) {
      x0 <- p[dependent]
      cd <- intersect(dependent, names(cache$guess))
      if (length(cd)) x0[cd] <- cache$guess[cd]
      try_solve(x0, "cache")
    }
    if (is.null(best) || best$maxres > ftol) try_solve(p[dependent], "user")
    if ((is.null(best) || best$maxres > ftol) && nstarts > 1L) {
      ## When the system is solved in log-space (positive = TRUE), draw the
      ## random starts log-uniform over [lower, upper] so they span several
      ## orders of magnitude — orders matter much more than ranges in
      ## biological networks. Otherwise stay with linear runif.
      lo <- expand_bounds(controls$lower, dependent, 0)
      hi <- expand_bounds(controls$upper, dependent, 10)
      if (positive) {
        lo_log <- log(pmax(lo, .Machine$double.eps))
        hi_log <- log(pmax(hi, .Machine$double.eps * 10))
      }
      ## A few alternative globalizations to cycle through when nleqslv
      ## stalls on stiff biological systems.
      altMethods <- list(top,
                         modifyList(top, list(global = "pwldog")),
                         modifyList(top, list(method = "Broyden", global = "dbldog")))
      for (i in seq_len(nstarts)) {
        x0_rand <- if (positive)
          setNames(exp(runif(n_dep, lo_log, hi_log)), dependent)
        else
          setNames(runif(n_dep, lo, hi), dependent)
        try_solve(x0_rand, paste0("random_", i),
                  override = altMethods[[((i - 1L) %% length(altMethods)) + 1L]])
        if (!is.null(best) && best$maxres <= ftol) break
      }
    }

    if (debugPlot && nstarts > 1L && length(log)) {
      log_df <- do.call(rbind, lapply(log, as.data.frame, stringsAsFactors = FALSE))
      log_df <- log_df[order(log_df$maxres), ]; log_df$rank <- seq_len(nrow(log_df))
      waterfall_plot(log_df, ftol)
    }

    if (is.null(best)) {
      warning("All solve attempts failed. Returning input values.", call. = FALSE)
      out <- c(p, zero_vec)
      return(as.parvec(out, deriv = NULL))
    }

    sol <- best$sol; root <- best$root
    if (sol$termcd != 1L && best$maxres > ftol)
      warning("nleqslv did not converge (code ", sol$termcd, "): ", sol$message, call. = FALSE)
    if (best$maxres > ftol) {
      ord  <- order(abs(best$res), decreasing = TRUE)
      topn <- min(10L, length(best$res))
      tbl  <- paste0(sprintf("  %-25s  |f| = %s",
                             names(best$res)[ord[1:topn]],
                             formatC(abs(best$res[ord[1:topn]]), format = "e", digits = 2)),
                     collapse = "\n")
      warning("Best residual norm ", formatC(best$maxres, format = "e", digits = 2),
              " (ftol = ", formatC(ftol, format = "e", digits = 2),
              "). Solution may be inaccurate.\nLargest residuals (top ", topn, "):\n", tbl,
              call. = FALSE)
    }

    out <- c(root, reconstruct(root, pv), zero_vec,
             p[setdiff(names(p), c(names(root), elim_states, zero_states))])

    if (keep.root) cache$guess <- out
    if (best$maxres > ftol && nstarts > 1L) {
      cache$failed_pv_hash <- pv_hash; cache$failed_result <- best
    } else { cache$failed_pv_hash <- NULL; cache$failed_result <- NULL }

    d <- tryCatch(
      build_derivs(root, pv, out, p, emptypars, fixed, dP, dP2, deriv, deriv2),
      error = function(e) {
        warning("Pimpl: IFT-based sensitivities unavailable at the current ",
                "root (", conditionMessage(e),
                "). Returning value only.", call. = FALSE)
        NULL
      })
    as.parvec(out,
              deriv  = if (deriv  && !is.null(d)) d$jacobian else NULL,
              deriv2 = if (deriv2 && !is.null(d)) d$hessian  else if (deriv2) NULL else FALSE)
  }

  attr(p2p, "equations")   <- as.eqnvec(trafo)
  attr(p2p, "parameters")  <- parameters
  attr(p2p, "modelname")   <- modelname
  attr(p2p, "compileInfo") <- collectCompileInfo(PEval$func, PEval$jac, PEval$hess)
  attr(p2p, "resetWarmStart") <- local({
    cache_ref <- cache; mn <- modelname; cond <- condition
    function() {
      cache_ref$guess <- NULL
      cache_ref$failed_pv_hash <- NULL
      cache_ref$failed_result <- NULL
      paste0("Pimpl(", mn, if (!is.null(cond)) paste0(":", cond) else "", ")")
    }
  })
  parfn(p2p, parameters, condition)
}


#' Parameter transformation (steady states via pre-equilibration)
#'
#' Returns a [parfn] that maps outer parameters to the steady state found
#' by integrating the ODE from `start.time` to `end.time`. Conserved
#' quantities are not eliminated; pick the basin of attraction by choosing
#' the dependent-state initial values in the input parvec. The Jacobian
#' (and optional Hessian) come from CppODE's analytical sensitivity
#' integration; the chain rule with incoming `dP`/`dP2` is applied here.
#'
#' @param trafo Named character / [eqnvec] / [eqnlist].
#' @param parameters Outer parameter names. States listed here are not
#'   integrated; they act as initial conditions and pass through.
#' @param forcings Forcing names; replaced by 0.
#' @param condition Condition label.
#' @param attach.input Append pass-through inputs to the output.
#' @param keep.root Warm-start subsequent calls from the cached steady
#'   state, its sensitivities, and re-use the previous result if the
#'   inputs are unchanged.
#' @param controlsODE Overrides for the ODE solver controls.
#' @param start.time,end.time Integration window; the root event fires at
#'   `end.time` if the steady state has not been reached.
#' @param compile,modelname,verbose Forwarded to [CppODE::CppODE].
#' @param deriv If `TRUE` (default), attach first-order parameter
#'   sensitivities `attr(., "deriv")` of the steady state.
#' @param deriv2 Emit second-order sensitivities `attr(., "deriv2")`;
#'   requires the model to be compiled with deriv2 support and
#'   `deriv = TRUE`.
#' @param ... Forwarded to [CppODE::CppODE].
#'
#' @return A [parfn].
#' @seealso [Pexpl], [Pimpl], [P].
#' @import CppODE
#' @importFrom digest digest
#' @export
Pequil <- function(trafo, parameters = NULL, forcings = NULL, condition = NULL,
                   attach.input = TRUE, start.time = -1e7, end.time = 0,
                   keep.root = TRUE, controlsODE = list(),
                   compile = FALSE, modelname = NULL, verbose = FALSE,
                   deriv = TRUE, deriv2 = FALSE, ...) {

  emit_d1 <- isTRUE(deriv)
  emit_d2 <- isTRUE(deriv2)
  if (emit_d2 && !emit_d1)
    stop("Pequil(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)

  zero_states <- character(0)
  if (inherits(trafo, "eqnlist")) {
    if (!is.null(forcings)) trafo$rates <- replaceSymbols(forcings, "0", trafo$rates)
    zs <- .zeroStatesFromSmatrix(trafo)
    zero_states <- zs$zero_states
    trafo <- zs$eqnlist
    f <- as.eqnvec(trafo)
  } else if (inherits(trafo, "eqnvec") || is.character(trafo)) {
    f <- as.eqnvec(replaceSymbols(forcings, "0", trafo))
  } else {
    stop("'trafo' must be an eqnlist, eqnvec or character vector", call. = FALSE)
  }

  states <- names(f)
  const_states <- states[vapply(unclass(f), function(x)
    tryCatch(identical(eval(parse(text = x)), 0), error = function(e) FALSE), logical(1))]
  if (length(const_states)) parameters <- union(parameters %||% character(0), const_states)

  dependent <- setdiff(states, parameters)
  if (!length(dependent)) stop("No dependent states left to equilibrate.", call. = FALSE)
  n_dep <- length(dependent)
  f_red <- f[dependent]

  parameters <- if (is.null(parameters)) getSymbols(f_red, exclude = dependent)
                else union(getSymbols(f_red, exclude = dependent), parameters)
  parms_all  <- setdiff(parameters, dependent)

  if (is.null(modelname)) modelname <- "equil_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

  ## Three CppODE models so callers can pay for the deriv-order they actually
  ## use at evaluation time.
  dotArgs <- list(...); dotArgs[["deriv2"]] <- NULL
  base <- c(list(rhs = unclass(f_red), rootfunc = "equilibrate", compile = compile,
                 outdir = getwd(), useDenseOutput = FALSE, verbose = verbose), dotArgs)
  model    <- do.call(CppODE::CppODE, c(base, list(deriv = FALSE, deriv2 = FALSE,
                                                   modelname = modelname)))
  model_s  <- if (emit_d1)
    do.call(CppODE::CppODE, c(base, list(deriv = TRUE,  deriv2 = FALSE,
                                         modelname = paste0(modelname, "_s"),
                                         fixed = names(f)))) else NULL
  model_s2 <- if (emit_d2)
    do.call(CppODE::CppODE, c(base, list(deriv = TRUE, deriv2 = TRUE,
                                         modelname = paste0(modelname, "_s2"),
                                         fixed = names(f)))) else NULL
  all_sens <- if (emit_d1) attr(model_s, "dimNames")$sens else character(0)

  ode_ctrl <- modifyList(list(abstol = 1e-6, reltol = 1e-6, maxsteps = 1e6L,
                              maxprogress = 100L, hini = 0, roottol = 1e-6, maxroot = 1L),
                         controlsODE)

  cache <- new.env(parent = emptyenv())
  cache$yini <- NULL; cache$sensini <- NULL; cache$sens2ini <- NULL
  cache$last_hash <- NULL; cache$last_result <- NULL

  default_sens <- matrix(0, n_dep, length(all_sens), dimnames = list(dependent, all_sens))
  diag_vars <- intersect(dependent, all_sens)
  if (length(diag_vars)) default_sens[cbind(diag_vars, diag_vars)] <- 1
  default_sens2 <- if (emit_d2)
    array(0, c(n_dep, length(all_sens), length(all_sens)),
          dimnames = list(dependent, all_sens, all_sens)) else NULL

  controls <- c(list(keep.root = keep.root, attach.input = attach.input,
                     start.time = start.time, end.time = end.time), ode_ctrl)

  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
    if (deriv2 && !emit_d2)
      stop("Pequil(deriv2 = TRUE) requires the model to be built with deriv2 = TRUE.",
           call. = FALSE)
    if (!emit_d1) deriv <- FALSE
    if (deriv2 && !deriv) deriv <- TRUE

    p   <- pars
    dP  <- attr(p, "deriv")
    dP2 <- if (deriv2) attr(p, "deriv2") else NULL
    keep.root    <- controls$keep.root
    attach.input <- controls$attach.input

    if (!is.null(fixed)) {
      p <- p[!names(p) %in% names(fixed)]
      p <- c(p, fixed)
    }
    emptypars <- setdiff(names(p), c(dependent, names(fixed)))
    miss <- setdiff(dependent, names(p)); if (length(miss)) p[miss] <- 1

    pv_hash <- NULL
    if (keep.root) {
      pv_hash <- digest::digest(list(p[dependent], p[parms_all], fixed, deriv, deriv2),
                                algo = "xxhash64")
      if (!is.null(cache$last_hash) && identical(pv_hash, cache$last_hash) &&
          !is.null(cache$last_result))
        return(cache$last_result)
      if (!is.null(cache$yini)) p[dependent] <- cache$yini
    }

    fixed_char  <- if (!is.null(fixed)) intersect(names(fixed), all_sens) else NULL
    active_sens <- if (length(fixed_char)) all_sens[-match(fixed_char, all_sens)] else all_sens
    n_active    <- length(active_sens)

    s1ini <- if (deriv  && keep.root && !is.null(cache$sensini))
               cache$sensini[, active_sens, drop = FALSE]
             else if (deriv)
               default_sens[, active_sens, drop = FALSE]
    s2ini <- if (deriv2 && keep.root && !is.null(cache$sens2ini))
               cache$sens2ini[, active_sens, active_sens, drop = FALSE]
             else if (deriv2)
               default_sens2[, active_sens, active_sens, drop = FALSE]

    sens_model <- if (deriv2) model_s2 else if (deriv) model_s else model
    res <- tryCatch(
      withCallingHandlers(
        CppODE::solveODE(
          sens_model,
          times = c(controls$start.time, controls$end.time),
          parms = c(p[dependent], p[parms_all]),
          sens1ini = s1ini, sens2ini = s2ini,
          fixed = if (deriv || deriv2) fixed_char,
          roottol = controls$roottol, abstol = controls$abstol, reltol = controls$reltol,
          maxsteps = as.integer(controls$maxsteps),
          maxprogress = as.integer(controls$maxprogress),
          hini = controls$hini, maxroot = as.integer(controls$maxroot)),
        warning = function(w) { warning(w$message, call. = FALSE); invokeRestart("muffleWarning") }),
      error = function(e) { warning("ODE integration failed: ", e$message, call. = FALSE); NULL })

    zero_vec <- if (length(zero_states))
      setNames(rep(0, length(zero_states)), zero_states) else NULL

    if (is.null(res)) {
      out <- if (attach.input)
               c(p[dependent], zero_vec, p[setdiff(names(p), dependent)])
             else c(p[dependent], zero_vec)
      return(as.parvec(out, deriv = NULL, deriv2 = NULL))
    }

    last <- length(res$time)
    if (last >= 1L && res$time[last] == controls$start.time)
      stop("Pequil: ODE solver made no progress from start.time = ", controls$start.time,
           ". The initial conditions and rate parameters likely produce a pathologically ",
           "stiff or unbounded system. Try realistic parameter magnitudes, loosen ",
           "`controlsODE` (`abstol`, `reltol`, `maxprogress`), or check for missing ",
           "forcings.", call. = FALSE)
    if (res$time[last] == end.time)
      warning("Steady state not reached within integration time.", call. = FALSE)

    digits <- floor(-log10(controls$roottol)) + 1L
    root <- setNames(round(res$variable[last, ], digits), dependent)
    out  <- if (attach.input)
              c(root, zero_vec, p[setdiff(names(p), dependent)])
            else c(root, zero_vec)

    if (keep.root) {
      cache$yini <- root
      cache$sensini <- if (!is.null(res$sens1)) {
        s <- default_sens; s[, active_sens] <- res$sens1[last, , ]; s
      } else NULL
      cache$sens2ini <- if (deriv2 && !is.null(res$sens2)) {
        s <- default_sens2; s[, active_sens, active_sens] <- res$sens2[last, , , ]; s
      } else NULL
    }

    if (!deriv || is.null(res$sens1)) {
      result <- as.parvec(out, deriv = NULL, deriv2 = NULL)
    } else {
      ## Don't round the sensitivities: rounding to roottol precision
      ## flushes legitimate small sensitivities (e.g. for states whose
      ## steady-state value is below threshold but still has a non-zero
      ## chain-rule contribution) to exact zero, which then makes the
      ## state appear "fixed" in as.parvec and creates a
      ## parameter-dependent `fixed` set.
      sens_final <- matrix(res$sens1[last, , ], n_dep, n_active,
                           dimnames = list(dependent, active_sens))
      input_cols <- setdiff(names(p), c(dependent, names(fixed)))
      jacobian <- matrix(0, length(out), length(input_cols),
                         dimnames = list(names(out), input_cols))
      if (attach.input) {
        idx <- intersect(emptypars, input_cols)
        if (length(idx)) jacobian[cbind(idx, idx)] <- 1
      }
      sr <- intersect(dependent, rownames(sens_final))
      sc <- intersect(input_cols, colnames(sens_final))
      if (length(sr) && length(sc)) jacobian[sr, sc] <- sens_final[sr, sc, drop = FALSE]

      hess_attr <- NULL
      if (deriv2 && !is.null(res$sens2)) {
        sens2_final <- array(res$sens2[last, , , ],
                             c(n_dep, n_active, n_active),
                             dimnames = list(dependent, active_sens, active_sens))
        hess_arr <- array(0, c(length(out), length(input_cols), length(input_cols)),
                          dimnames = list(names(out), input_cols, input_cols))
        if (length(sr) && length(sc))
          hess_arr[sr, sc, sc] <- sens2_final[sr, sc, sc, drop = FALSE]
        if (!is.null(dP)) {
          dPsub <- submatrix(dP, rows = input_cols)
          th    <- colnames(dPsub); n_th <- length(th)
          new_hess <- t(dPsub) %bmm% hess_arr %bmm% dPsub
          dimnames(new_hess) <- list(names(out), th, th)
          if (!is.null(dP2)) {
            dP2sub <- dP2[input_cols, th, th, drop = FALSE]
            new_hess <- new_hess + array(
              jacobian %*% matrix(dP2sub, length(input_cols), n_th * n_th),
              c(length(out), n_th, n_th),
              dimnames = list(names(out), th, th))
          }
          jacobian  <- jacobian %*% dPsub
          hess_attr <- new_hess
        } else hess_attr <- hess_arr
      } else if (!is.null(dP)) {
        jacobian <- jacobian %*% submatrix(dP, rows = input_cols)
      }

      keep <- rowSums(jacobian != 0) > 0
      hess_keep <- if (!is.null(hess_attr)) hess_attr[keep, , , drop = FALSE] else FALSE
      result <- as.parvec(out, deriv = jacobian[keep, , drop = FALSE], deriv2 = hess_keep)
    }

    if (keep.root) { cache$last_hash <- pv_hash; cache$last_result <- result }
    result
  }

  attr(p2p, "equations")   <- as.eqnvec(f)
  attr(p2p, "parameters")  <- parameters
  attr(p2p, "modelname")   <- modelname
  attr(p2p, "compileInfo") <- collectCompileInfo(model, model_s, model_s2)
  attr(p2p, "resetWarmStart") <- local({
    cache_ref <- cache; mn <- modelname; cond <- condition
    function() {
      cache_ref$yini <- NULL
      cache_ref$sensini <- NULL
      cache_ref$sens2ini <- NULL
      cache_ref$last_hash <- NULL
      cache_ref$last_result <- NULL
      paste0("Pequil(", mn, if (!is.null(cond)) paste0(":", cond) else "", ")")
    }
  })
  parfn(p2p, parameters, condition)
}


#' Construct and modify parameter transformations
#'
#' Symbolic helpers used by [P()] and [Xs()] to build, substitute, and
#' branch transformation rules. The condition table from [branch()] is
#' stored on `attr(., "tree")`.
#'
#' - `define` resets the LHS of `expr` to its RHS.
#' - `insert` substitutes the RHS for the LHS wherever it occurs.
#' - `branch` duplicates a trafo across conditions, optionally applying
#'   per-condition substitutions taken from `table`.
#'
#' @param trafo Named character / [eqnvec], or a list thereof.
#' @param expr `"lhs ~ rhs"` formula string.
#' @param table Condition table (row per condition, column per parameter)
#'   carried as the `tree` attribute when branching.
#' @param conditions Condition names; default `rownames(table)`.
#' @param apply One of `"nothing"`, `"insert"`, `"define"`; how the
#'   `table` entries are folded into each branch.
#' @param conditionMatch Regex on condition names; restricts the operation.
#' @param ... Named values to substitute into `expr` symbols.
#'
#' @return Same shape as `trafo` (or a per-condition list if branched).
#' @export
#' @example inst/examples/define.R
define <- function(trafo, expr, ..., conditionMatch = NULL) {
  if (missing(trafo)) trafo <- NULL
  tree <- attr(trafo, "tree")
  if (is.list(trafo) && is.null(names(trafo)))
    stop("If trafo is a list, elements must be named.", call. = FALSE)
  if (is.list(trafo) && !all(names(trafo) %in% rownames(tree)))
    stop("List names must be a subset of rownames(attr(trafo, 'tree')).", call. = FALSE)
  mytrafo <- if (is.list(trafo)) trafo else list(trafo)

  dots <- substitute(alist(...))
  out  <- lapply(seq_along(mytrafo), function(i) {
    row <- if (is.list(trafo)) tree[names(mytrafo)[i], , drop = FALSE]
           else tree[1, , drop = FALSE]
    if (!is.null(conditionMatch) && !str_detect(rownames(row), conditionMatch))
      return(mytrafo[[i]])
    with(row, do.call(repar,
      c(list(expr = expr, trafo = mytrafo[[i]], reset = TRUE), eval(dots))))
  })
  names(out) <- names(mytrafo)
  if (!is.list(trafo)) out <- out[[1]]
  attr(out, "tree") <- tree
  out
}


#' @export
#' @rdname define
insert <- function(trafo, expr, ..., conditionMatch = NULL) {
  if (missing(trafo)) trafo <- NULL
  tree <- attr(trafo, "tree")
  if (is.list(trafo) && is.null(names(trafo)))
    stop("If trafo is a list, elements must be named.", call. = FALSE)
  if (is.list(trafo) && !all(names(trafo) %in% rownames(tree)))
    stop("List names must be a subset of rownames(attr(trafo, 'tree')).", call. = FALSE)
  mytrafo <- if (is.list(trafo)) trafo else list(trafo)

  dots <- substitute(alist(...))
  out  <- lapply(seq_along(mytrafo), function(i) {
    cur <- mytrafo[[i]]
    row <- if (is.list(trafo)) tree[names(mytrafo)[i], , drop = FALSE]
           else tree[1, , drop = FALSE]
    if (!is.null(conditionMatch) && !str_detect(rownames(row), conditionMatch))
      return(cur)
    with(row, {
      ## Caller may pass logical dots to gate substitution per condition,
      ## and non-logical dots to substitute symbols in `expr`. Logical dots
      ## are stripped before forwarding to `repar`.
      .apply <- function() {
        d <- eval(dots)
        if (!length(d)) return(do.call(repar, list(expr = expr, trafo = cur)))
        d_eval  <- lapply(d, function(x) eval.parent(x, 3))
        is_log  <- vapply(d_eval, is.logical, logical(1))
        gate    <- do.call(c, d[is_log])
        if (!is.null(gate) && any(!gate)) return(cur)
        do.call(repar, c(list(expr = expr, trafo = cur), d_eval[!is_log]))
      }
      .apply()
    })
  })
  names(out) <- names(mytrafo)
  if (!is.list(trafo)) out <- out[[1]]
  attr(out, "tree") <- tree
  out
}


#' @export
#' @rdname define
branch <- function(trafo, table = NULL,
                   conditions = rownames(table),
                   apply = c("nothing", "insert", "define")) {
  apply <- match.arg(apply)
  if (is.null(table) && is.null(conditions)) return(trafo)
  if (is.null(conditions)) conditions <- paste0("C", seq_len(nrow(table)))
  if (is.null(table))      table      <- data.frame(condition = conditions,
                                                    row.names = conditions)
  rownames(table) <- conditions

  out <- setNames(lapply(conditions, function(x) trafo), conditions)
  attr(out, "tree") <- table
  if (apply == "nothing") return(out)

  for (cn in conditions) {
    row <- table[cn, !colnames(table) %in% c("condition", "conditions"), drop = FALSE]
    for (par in colnames(row)) {
      val <- row[[par]]; if (is.na(val)) next
      single <- out[[cn]]; attr(single, "tree") <- row
      out[[cn]] <- if (apply == "insert") insert(single, paste0(par, " ~ ", val))
                   else                    define(single, paste0(par, " ~ ", val))
      attr(out[[cn]], "tree") <- NULL
    }
  }
  out
}



#' Reparameterization
#'
#' Replaces symbols on either side of `"lhs ~ rhs"`. With `reset = TRUE`
#' the LHS entry of `trafo` is overwritten by the RHS (per row when `...`
#' supplies vector replacements); otherwise the LHS is substituted into
#' `trafo` wherever it occurs. Symbols separated by `_` are recognised as
#' compound identifiers, e.g. `Delta_x` -> symbols `"Delta"` and `"x"`.
#'
#' @param expr `"lhs ~ rhs"` string (or a formula).
#' @param trafo Character / [eqnvec] / list. `NULL` builds a fresh trafo
#'   from the LHS.
#' @param ... Named character/numeric vectors; each row supplies one
#'   substitution.
#' @param reset Overwrite (`TRUE`) or substitute into (`FALSE`).
#' @return Same shape as `trafo`.
#' @export
#' @importFrom stats as.formula
#' @examples
#' innerpars   <- letters[1:3]
#' constraints <- c(a = "b + c")
#' mycondition <- "cond1"
#' trafo <- repar("x ~ x",        x = innerpars)
#' trafo <- repar("x ~ y",        trafo, x = names(constraints), y = constraints)
#' trafo <- repar("x ~ exp(x)",   trafo, x = innerpars)
#' trafo <- repar("x ~ x + Delta_x_condition",
#'                trafo, x = innerpars, condition = mycondition)
repar <- function(expr, trafo = NULL, ..., reset = FALSE) {
  if (inherits(expr, "formula")) expr <- deparse(expr)
  parsed <- as.character(stats::as.formula(gsub("_", ":", expr, fixed = TRUE)))
  lhs <- parsed[2]; rhs <- parsed[3]

  args <- lapply(list(...), as.character)
  if (length(args)) {
    reps <- as.data.frame(args, stringsAsFactors = FALSE)
    apply_repl <- function(side) vapply(seq_len(nrow(reps)), function(i)
      gsub(":", "_", replaceSymbols(colnames(reps), reps[i, ], side), fixed = TRUE),
      character(1))
    lhs <- apply_repl(lhs); rhs <- apply_repl(rhs)
  } else {
    lhs <- gsub(":", "_", lhs, fixed = TRUE)
    rhs <- gsub(":", "_", rhs, fixed = TRUE)
  }

  if (is.null(trafo))                           as.eqnvec(structure(lhs, names = lhs))
  else if (is.list(trafo)      && !reset)       lapply(trafo, function(t) replaceSymbols(lhs, rhs, t))
  else if (is.character(trafo) && !reset)       replaceSymbols(lhs, rhs, trafo)
  else if (is.list(trafo)      &&  reset)       lapply(trafo, function(t) { t[lhs] <- rhs; t })
  else { trafo[lhs] <- rhs; trafo }
}

paste_ <- function(...) paste(..., sep = "_")