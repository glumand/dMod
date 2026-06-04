## Functions to generate parameter transformation ----

#' Generate a parameter transformation function
#'
#' Unified entry to the three backends: explicit ([Pexpl], algebraic),
#' implicit ([Pimpl], root-finding) and equilibrate ([Pequil], ODE
#' pre-integration). `method = NULL` picks `"equilibrate"` for
#' [eqnlist] entries, `"explicit"` otherwise.
#'
#' @param trafo An [eqnvec], named character, [eqnlist], or list thereof.
#' @param parameters Outer-parameter names.
#' @param condition Condition label.
#' @param compile,modelname,verbose Forwarded to [CppODE::funCpp].
#' @param method One of `"explicit"`, `"implicit"`, `"equilibrate"`, or `NULL`.
#' @param cores Per-condition `mclapply()` cores. `NULL` auto-detects via
#'   [detectFreeCores]; capped at 1 on Windows.
#' @param deriv,deriv2 Attach first/second-order sensitivities. `deriv2`
#'   requires `deriv = TRUE`.
#' @param ... Forwarded to the chosen backend.
#'
#' @return A [parfn].
#' @seealso [Pexpl], [Pimpl], [Pequil], [parfn]
#' @export
P <- function(trafo = NULL, parameters = NULL, condition = NULL,
              compile = FALSE, modelname = NULL, method = NULL,
              cores = NULL, verbose = FALSE,
              deriv = TRUE, deriv2 = FALSE, ...) {

  if (is.null(trafo)) return()
  if (isTRUE(deriv2) && !isTRUE(deriv))
    stop("P(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)

  if (!is.list(trafo) || inherits(trafo, "eqnlist") || inherits(trafo, "eqnvec")) {
    trafo_list <- list(trafo); names(trafo_list) <- condition
  } else trafo_list <- trafo

  ## detectFreeCores() has side effects (warnings, SSH on remotes); call once.
  cores <- if (Sys.info()[['sysname']] == "Windows") 1L
           else if (is.null(cores)) detectFreeCores()
           else min(detectFreeCores(), cores)

  ## Codegen-only inside mclapply: dyn.load in a forked worker is lost when
  ## the fork exits, so the actual compile happens in the parent below.
  result <- Reduce("+", mclapply(seq_along(trafo_list), function(i) {
    tr   <- trafo_list[[i]]
    cond <- names(trafo_list[i])
    m <- if (!is.null(method)) match.arg(method, c("explicit", "implicit", "equilibrate"))
         else if (inherits(tr, "eqnlist")) "equilibrate" else "explicit"
    switch(m,
      explicit    = Pexpl(as.eqnvec(tr), parameters = parameters, condition = cond,
                          compile = FALSE, modelname = modelname, verbose = verbose,
                          deriv = deriv, deriv2 = deriv2, ...),
      implicit    = Pimpl(trafo = tr, parameters = parameters, condition = cond,
                          compile = FALSE, modelname = modelname, verbose = verbose,
                          deriv = deriv, deriv2 = deriv2, ...),
      equilibrate = Pequil(trafo = tr, parameters = parameters, condition = cond,
                           compile = FALSE, modelname = modelname, verbose = verbose,
                           deriv = deriv, deriv2 = deriv2, ...))
  }, mc.cores = cores))

  if (compile) compile(result, cores = cores, output = modelname, verbose = verbose)
  result
}


#' Parameter transformation (explicit, algebraic)
#'
#' Builds `p_inner = f(p_outer)` from symbolic expressions via
#' [CppODE::funCpp], in forward-mode AD or SymPy mode. The returned
#' [parfn] attaches the Jacobian and, optionally, the Hessian.
#'
#' @param trafo Named character / [eqnvec]; names are inner parameters,
#'   values are expressions in the outer parameters.
#' @param parameters Outer parameters; defaults to `getSymbols(trafo)`.
#' @param attach.input Append outer inputs to the output.
#' @param condition Condition label.
#' @param compile,modelname,verbose Forwarded to [CppODE::funCpp].
#' @param deriv,deriv2 Attach `attr(., "deriv")` `[p, theta]` and/or
#'   `attr(., "deriv2")` `[p, theta, theta]`. `deriv2` needs `deriv = TRUE`.
#' @param derivMode `"dual"` (AD, needs `compile = TRUE`) or `"symbolic"`.
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
    if (!emit_d1) deriv <- FALSE
    if (deriv2 && !deriv) deriv <- TRUE

    p <- c(pars, fixed)
    ad_ok  <- use_ad && !is.null(evaluate) && is.loaded(ad_symbol)
    ad_ok2 <- ad_ok && emit_d2 && is.loaded(ad2_symbol)
    if (deriv2 && !ad_ok2 && use_ad)
      stop("Pexpl(deriv2 = TRUE) needs the compiled AD2 entry; rebuild with compile = TRUE.", call. = FALSE)

    Jac <- NULL; Hess <- NULL

    if (ad_ok && deriv) {
      ## Dual-mode entry reads `params` and `dP` positionally against the
      ## codegen order, so reorder both to `parameters`.
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


#' Conserved-quantity coefficient matrix and pivot choice
#'
#' Builds the linear CQ coefficient matrix `C` (rows = totals, columns =
#' participating species) and picks one pivot species per total by Gaussian
#' elimination, so that `C[, pivots]` is invertible even for overlapping
#' totals. Tie-breaks prefer species not needed by later totals, then
#' non-`parameters`, then alphabetical order.
#'
#' @param totals Named list of CQ expressions (from [getTotals()]).
#' @param states State names participating in the model.
#' @param parameters Outer parameters, used only for the pivot tie-break.
#' @return `list(C_mat, pivots, cq_sets)`. `pivots` is aligned to `totals`
#'   with `NA` where a total has fewer than two free species.
#' @keywords internal
.cq_pivot_decomposition <- function(totals, states, parameters = character(0)) {
  cq_sets <- lapply(totals, function(expr) intersect(getSymbols(expr), states))
  all_cq_species <- unique(unlist(cq_sets, use.names = FALSE))
  n_cq <- length(totals)
  C_mat <- matrix(0, n_cq, length(all_cq_species),
                  dimnames = list(names(totals), all_cq_species))
  for (i in seq_len(n_cq)) {
    expr <- parse(text = totals[[i]])
    e0 <- setNames(as.list(rep(0, length(all_cq_species))), all_cq_species)
    base <- eval(expr, envir = e0)
    for (s in cq_sets[[i]]) {
      e1 <- e0; e1[[s]] <- 1
      C_mat[i, s] <- eval(expr, envir = e1) - base
    }
  }

  C_red <- C_mat
  pivots <- rep(NA_character_, n_cq)
  for (i in seq_len(n_cq)) {
    cand <- setdiff(cq_sets[[i]], pivots[!is.na(pivots)])
    if (length(cand) < 2L) next
    nz <- cand[abs(C_red[i, cand]) > 1e-12]
    if (!length(nz)) next
    future <- if (i < n_cq)
      unique(unlist(cq_sets[(i + 1L):n_cq], use.names = FALSE)) else character(0)
    pool <- if (length(safe <- setdiff(nz, future))) safe else nz
    nu   <- setdiff(pool, parameters)
    e    <- if (length(nu)) sort(nu)[1L] else sort(pool)[1L]
    if (i < n_cq) {
      piv <- C_red[i, e]
      for (j in (i + 1L):n_cq)
        if (abs(C_red[j, e]) > 1e-12)
          C_red[j, ] <- C_red[j, ] - (C_red[j, e] / piv) * C_red[i, ]
    }
    pivots[i] <- e
  }
  list(C_mat = C_mat, pivots = pivots, cq_sets = cq_sets)
}

#' Substitute conserved quantities into ODE rates
#'
#' For each conserved quantity one pivot species is eliminated.
#' `expressInTotals = TRUE` substitutes `x_e -> total_i - rest` and adds a
#' `total_i` parameter; `FALSE` promotes the pivot to a pass-through
#' parameter so its redundant rate equation drops.
#'
#' @param totals Named list of CQ expressions (from [getTotals()]).
#'   Empty list disables CQ handling.
#' @param has_smatrix Whether the caller had structural info available;
#'   gates the diagnostic warning when `totals` is empty.
#' @param f,states,parameters Equation set, state names, user parameters.
#' @param expressInTotals See above.
#' @return `list(f, parameters, cq_info, elim_states)`. `cq_info`/
#'   `elim_states` are non-empty only in `TRUE` mode.
#' @keywords internal
.detect_and_substitute_cq <- function(totals, has_smatrix, f, states, parameters,
                                      expressInTotals = TRUE) {
  cq_info <- list(); elim_states <- character(0)

  if (length(totals)) {
    if (is.null(parameters)) parameters <- character(0)
    dec <- .cq_pivot_decomposition(totals, states, parameters)
    C_mat <- dec$C_mat
    substitutions <- list()
    for (i in seq_along(totals)) {
      e <- dec$pivots[i]
      if (is.na(e)) next
      if (expressInTotals) {
        total_name <- names(totals)[i]
        coef_e <- C_mat[i, e]
        rest <- paste0("(", replaceSymbols(e, "0", totals[[i]]), ")")
        recon_expr <- if (isTRUE(all.equal(coef_e, 1)))
          paste0("(", total_name, " - ", rest, ")")
        else
          paste0("((", total_name, " - ", rest, ") / (", coef_e, "))")
        substitutions[[e]] <- recon_expr
        elim_states <- c(elim_states, e)
        parameters  <- union(parameters, total_name)
        cq_info[[length(cq_info) + 1L]] <- list(
          total_name = total_name, elim_state = e,
          recon_expr = setNames(recon_expr, e))
      } else {
        parameters <- union(parameters, e)
      }
    }
    if (length(substitutions)) {
      keys <- names(substitutions)
      vals <- unname(unlist(substitutions))
      for (.iter in seq_along(substitutions)) {
        fNew <- replaceSymbols(keys, vals, f)
        if (identical(fNew, f)) break
        f <- fNew
      }
      for (i in seq_along(cq_info)) {
        rec <- cq_info[[i]]$recon_expr
        for (.iter in seq_along(substitutions)) {
          recNew <- replaceSymbols(keys, vals, rec)
          if (identical(recNew, rec)) break
          rec <- recNew
        }
        cq_info[[i]]$recon_expr <- setNames(rec, cq_info[[i]]$elim_state)
      }
    }
  } else if (!has_smatrix && !is.null(parameters)) {
    param_states <- intersect(parameters, states)
    if (length(param_states)) {
      remaining <- f[setdiff(states, param_states)]
      still <- if (length(remaining))
        param_states[vapply(param_states, function(ps)
          any(ps %in% getSymbols(remaining)), logical(1))] else character(0)
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
#' Iterated three-layer test on the stoichiometric matrix of an `eqnlist`
#' (matches AlyssaPetit v1.2):
#' \enumerate{
#'   \item *Neg-only column*: only outflux, must be zero at SS.
#'   \item *Pos-only column with single-state feeder*: each feeding flux
#'     vanishes at SS; if a feeder's rate involves exactly one state, that
#'     state must be zero.
#'   \item *Sink cluster (LP)*: subsets whose combined mass leaks
#'     monotonically (mass-balance LP via `lpSolve::lp`, only when installed).
#' }
#' For each zero-state the column is dropped, the state symbol is substituted
#' by `"0"` in remaining rates, structurally-zero reactions are removed, and
#' detection re-runs (removals can expose new zero-states).
#'
#' @param eqnlist_obj An [eqnlist].
#' @return `list(zero_states, eqnlist)` with the reduced system.
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

  ## Drop a-priori-zero rates: their +1 stoichiometry otherwise masks the
  ## sink-cluster LP downstream.
  drop_rate0 <- vapply(rates, .is_struct_zero, logical(1))
  if (any(drop_rate0)) {
    keep <- !drop_rate0
    S     <- S[keep, , drop = FALSE]
    rates <- rates[keep]
    description <- description[keep]
    if (!is.null(reactionCompartment))
      reactionCompartment <- reactionCompartment[keep]
  }

  ## No influx (column has no strictly positive entry) => zero at SS;
  ## also catches all-zero columns from prior reductions.
  .neg_col <- function(M) {
    for (j in seq_len(ncol(M))) if (!any(M[, j] > 0)) return(j)
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

  ## Substitute a state -> 0 stoichiometrically AND kinetically; drop any
  ## reaction whose rate is structurally zero after substitution.
  .zero_out <- function(j) {
    state_name <- colnames(S)[j]
    zero_states <<- c(zero_states, state_name)
    keep_rxn <- rep(TRUE, nrow(S))
    for (k in seq_len(nrow(S))) {
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
#' Walks `fn` and its closure environments and clears the warm-start
#' cache on every reachable [Pequil] / [Pimpl] parfn. Use before workflows
#' that cross basins of attraction (multistart, profile after a structural
#' change) where a stale root pins the solver in the wrong region.
#'
#' @param fn A `parfn`, `prdfn`, `obsfn`, `objfn`, or composed `fn`.
#' @param verbose Print one-line summary of cleared caches.
#' @return Invisibly, labels of the cleared caches (empty if none found).
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


## Solve A %*% X = B; LU when well-conditioned, Moore-Penrose pseudoinverse
## (minimum-norm IFT sensitivity on the constraint manifold) otherwise. The
## warning lists null-space directions in `row_names` coordinates so the
## missing CQ or redundant equation is identifiable.
#' @keywords internal
.pimpl_solve_dfdx <- function(A, B, row_names = rownames(A), warn_rcond = 1e-10) {
  if (nrow(A) == 0L) return(B)
  sv  <- svd(A); d <- sv$d
  tol <- max(dim(A)) * d[1L] * .Machine$double.eps
  rnk <- sum(d > tol)
  rc  <- if (d[1L] > 0) d[length(d)] / d[1L] else 0

  if (rnk == nrow(A) && is.finite(rc) && rc >= warn_rcond)
    return(solve(A, B))

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

## Shared SS preamble: coerce to eqnvec, zero+drop forcings (a state with
## rhs == 0 stays at its initial value so it can be cut), run sink-state
## detection, promote constant-rhs states to parameters, and eliminate
## CQs according to `expressInTotals` (see `.detect_and_substitute_cq`).
## Returns the normalised record used by Pimpl / Pequil.
#' Normalise steady-state inputs
#'
#' Coerces `trafo` to an `eqnvec`, zeroes forcings, drops structurally-zero
#' states and promotes constant-rhs states to parameters. In the default
#' `fullsystem = FALSE` mode it eliminates conserved quantities via
#' [.detect_and_substitute_cq] (used by [Pimpl] and `Pequil` without
#' `expressInTotals`). In `fullsystem = TRUE` mode it leaves `f` intact and
#' returns the CQ pivot decomposition so the caller can integrate the full
#' system and inject the totals through initial conditions.
#'
#' @param trafo An [eqnlist], [eqnvec] or named character vector.
#' @param parameters Outer parameters (may be `NULL`).
#' @param forcings Forcing names; zeroed and removed.
#' @param expressInTotals Passed to [.detect_and_substitute_cq] when
#'   `fullsystem = FALSE`.
#' @param fullsystem If `TRUE`, skip substitution and return the full `f`
#'   plus `C_mat`, `pivots`, `moiety_species` and `totals`.
#' @return A named list with the normalised `trafo`, `states`, `dependent`,
#'   `parameters`, `parms_all`, `zero_states` and (mode-specific) CQ fields.
#' @keywords internal
.normalize_ss_inputs <- function(trafo, parameters, forcings, expressInTotals = TRUE,
                                 fullsystem = FALSE) {
  smatrix <- NULL; zero_states <- character(0)
  original_params <- character(0)
  totals <- list()

  if (inherits(trafo, "eqnlist")) {
    original_params <- setdiff(getParameters(trafo), trafo$states)
    if (!is.null(forcings))
      trafo$rates <- replaceSymbols(forcings, rep("0", length(forcings)), trafo$rates)
    zs <- .zeroStatesFromSmatrix(trafo)
    zero_states <- zs$zero_states
    trafo       <- zs$eqnlist
    smatrix     <- trafo$smatrix
    if (!length(trafo$states))
      stop("All states are structurally zero in steady state; no dynamical ",
           "state remains to solve for. The network likely has irreversible ",
           "drains without matching influx (an open system with trivial ",
           "all-zero equilibrium).", call. = FALSE)
    totals      <- getTotals(trafo)
  } else if (inherits(trafo, "eqnvec") || is.character(trafo)) {
    if (!is.null(forcings))
      trafo <- replaceSymbols(forcings, rep("0", length(forcings)), trafo)
  } else stop("'trafo' must be an eqnlist, eqnvec or character vector", call. = FALSE)
  trafo <- as.eqnvec(trafo)
  if (!is.null(forcings)) trafo <- trafo[setdiff(names(trafo), forcings)]

  states <- names(trafo)
  const_states <- states[vapply(unclass(trafo), function(x)
    tryCatch(identical(eval(parse(text = x)), 0),
             error = function(e) FALSE), logical(1))]
  if (length(const_states))
    parameters <- union(parameters %||% character(0), const_states)

  if (fullsystem && length(totals)) {
    dec    <- .cq_pivot_decomposition(totals, states, parameters %||% character(0))
    pivots <- dec$pivots[!is.na(dec$pivots)]
    dependent <- setdiff(states, parameters %||% character(0))
    if (!length(dependent))
      stop("No dynamical states to integrate. All states appear in 'parameters'.",
           call. = FALSE)
    parameters <- Reduce(union, list(getSymbols(trafo[dependent], exclude = dependent),
                                     parameters %||% character(0),
                                     original_params, names(totals)))
    return(list(trafo = trafo, states = states, zero_states = zero_states,
                dependent = dependent, parameters = parameters,
                parms_all = setdiff(parameters, dependent), smatrix = smatrix,
                totals = totals, C_mat = dec$C_mat, pivots = pivots,
                moiety_species = colnames(dec$C_mat)))
  }

  cq <- .detect_and_substitute_cq(totals, !is.null(smatrix), trafo, states,
                                  parameters, expressInTotals)
  trafo       <- cq$f
  parameters  <- cq$parameters
  cq_info     <- cq$cq_info
  elim_states <- cq$elim_states

  dependent <- setdiff(states, c(parameters, elim_states))
  if (!length(dependent))
    stop("No dependent states to solve for. All states appear in 'parameters'.",
         call. = FALSE)

  parameters <- Reduce(union, list(getSymbols(trafo, exclude = dependent),
                                   parameters %||% character(0),
                                   original_params))
  parms_all  <- setdiff(parameters, dependent)

  list(trafo = trafo, states = states, zero_states = zero_states,
       dependent = dependent, parameters = parameters, parms_all = parms_all,
       smatrix = smatrix, cq_info = cq_info, elim_states = elim_states)
}

#' @keywords internal
.expand_bounds <- function(b, dep, default_val) {
  if (is.null(names(b)) || length(b) == 1L)
    return(setNames(rep(b[1L], length(dep)), dep))
  out <- setNames(rep(default_val, length(dep)), dep)
  nm  <- intersect(names(b), dep); out[nm] <- b[nm]; out
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
#' Returns a [parfn] over the outer inputs. On call, the parfn solves
#' `f(x, p) = 0` for the dependent states via [nleqslv::nleqslv], warm
#' starting from the cached root when available and falling back to
#' multistart otherwise. The IFT-derived Jacobian (and Hessian when
#' `deriv2 = TRUE`) are attached to the result. For [eqnlist] inputs with
#' conserved moieties, see `expressInTotals`.
#'
#' @param trafo Named character / [eqnvec] / [eqnlist].
#' @param parameters Outer parameters; auto-extended with the `total_*`
#'   conserved-quantity parameters.
#' @param forcings Forcing names; zeroed and removed.
#' @param condition Condition label.
#' @param keep.root Cache the root as warm-start for the next call.
#' @param expressInTotals If `TRUE` (default), every conserved moiety stays a
#'   solve variable and one redundant rate equation per conserved quantity is
#'   replaced by the algebraic conservation constraint `sum(c_k x_k) = total_X`,
#'   adding `total_X` as a parameter. The moiety is solved in log space, so all
#'   species stay positive and none is reconstructed by subtraction; the
#'   conservation then holds to the solver tolerance (`controlsNleqslv$ftol`).
#'   If `FALSE`, the pivot species per conserved quantity becomes a pass-through
#'   parameter and its redundant equation is dropped.
#' @param compile,modelname,verbose Forwarded to [CppODE::funCpp].
#' @param deriv,deriv2 Attach first/second-order IFT sensitivities.
#'   `deriv2` requires `funCpp` to expose `hess()`.
#' @param controlsMS Multistart controls. Recognised keys: `nStarts`
#'   (default `100L`; `1L` disables multistart), `positive` (default
#'   `TRUE`; selects nleqslv's log-space transform and log-uniform
#'   random starts), `lower`/`upper` (scalar or named vector of bounds
#'   for the random sweep), `debugPlot` (default `FALSE`; emit a
#'   waterfall plot of multistart termination codes).
#' @param controlsNleqslv nleqslv tuning: `method`, `global`, `xscalm`,
#'   `xtol`, `ftol`, `btol`, `cndtol`, `maxit`, `allowSingular`. Residuals
#'   are scaled per equation by their turnover `sum_k |df_i/dx_k| |x_k|`, so
#'   `ftol` is a relative criterion that converges multi-scale systems
#'   uniformly.
#'
#' @return A [parfn].
#' @seealso [Pexpl], [Pequil], [P].
#' @export
#' @import nleqslv
#' @importFrom digest digest
Pimpl <- function(trafo, parameters = NULL, forcings = NULL, condition = NULL,
                  keep.root = TRUE, expressInTotals = TRUE, compile = FALSE,
                  modelname = NULL, verbose = FALSE, deriv = TRUE, deriv2 = FALSE,
                  controlsMS = list(), controlsNleqslv = list()) {

  emit_d1 <- isTRUE(deriv)
  emit_d2 <- isTRUE(deriv2)
  if (emit_d2 && !emit_d1)
    stop("Pimpl(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)

  norm <- .normalize_ss_inputs(trafo, parameters, forcings,
                               expressInTotals = expressInTotals,
                               fullsystem = isTRUE(expressInTotals))
  states      <- norm$states
  zero_states <- norm$zero_states
  dependent   <- norm$dependent
  parameters  <- norm$parameters
  parms_all   <- norm$parms_all

  if (is.null(modelname)) modelname <- "impl_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

  if (!is.null(norm$pivots) && length(norm$pivots)) {
    tn   <- names(norm$totals)
    cons <- setNames(vapply(seq_along(norm$totals), function(i)
      paste0("((", norm$totals[[i]], ")/(", tn[i], ") - 1)"), character(1)), tn)
    all_exprs <- c(unclass(norm$trafo[setdiff(dependent, norm$pivots)]), cons)
  } else {
    all_exprs <- unclass(norm$trafo[dependent])
  }
  n_dep <- length(dependent)
  parms_all <- intersect(parms_all, getSymbols(all_exprs))

  PEval <- suppressWarnings(CppODE::funCpp(
    all_exprs, variables = dependent, parameters = parms_all, fixed = NULL,
    compile = compile, modelname = modelname, outdir = getwd(),
    verbose = verbose, convenient = FALSE,
    deriv = TRUE, deriv2 = emit_d2, derivMode = "symbolic"))

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
  cache$guess <- NULL

  ## Biological SS magnitudes span ~10 orders when rate constants span 2-3:
  ## default log-uniform sampling over [1e-5, 1e5] when positive = TRUE.
  ms <- modifyList(list(nStarts = 100L, positive = TRUE,
                        lower = 1e-5, upper = 1e5, debugPlot = FALSE),
                   controlsMS)
  if (!ms$positive) {
    if (identical(ms$lower, 1e-5)) ms$lower <- 0
    if (identical(ms$upper, 1e5))  ms$upper <- 100
  }
  nleq <- modifyList(list(method = "Newton", global = "dbldog", xscalm = "fixed",
                          xtol = 1e-4, ftol = 1e-2, btol = 1e-3,
                          cndtol = 1e-12, maxit = 200L, allowSingular = TRUE),
                     controlsNleqslv)
  nleqslv_top <- c("method", "global", "xscalm")

  turnover <- function(x, pv)
    pmax(as.numeric(abs(eval_J(x, pv)[, dependent, drop = FALSE]) %*% abs(x)),
         .Machine$double.eps)

  solve_once <- function(x0, pv, positive, top, ctrl) {
    nl <- function(start, scale) {
      if (positive) {
        s0 <- start; s0[s0 <= 0] <- 1
        sol <- nleqslv::nleqslv(
          x   = log(s0),
          fn  = function(lx) { x <- exp(lx); names(x) <- dependent; eval_f(x, pv) / scale },
          jac = function(lx) { x <- exp(lx); names(x) <- dependent
                               sweep(eval_J(x, pv)[, dependent, drop = FALSE], 2L, x, `*`) / scale },
          method = top$method, global = top$global, xscalm = top$xscalm, control = ctrl)
        list(root = setNames(exp(sol$x), dependent), sol = sol)
      } else {
        sol <- nleqslv::nleqslv(
          x   = start,
          fn  = function(x) { names(x) <- dependent; eval_f(x, pv) / scale },
          jac = function(x) { names(x) <- dependent; eval_J(x, pv)[, dependent, drop = FALSE] / scale },
          method = top$method, global = top$global, xscalm = top$xscalm, control = ctrl)
        list(root = setNames(sol$x, dependent), sol = sol)
      }
    }
    r1  <- nl(x0, rep(1, n_dep))
    sc  <- turnover(r1$root, pv)
    res <- tryCatch(eval_f(r1$root, pv) / sc, error = function(e) setNames(rep(Inf, n_dep), dependent))
    if (max(abs(res)) <= ctrl$ftol)
      return(list(root = r1$root, sol = r1$sol, res = res,
                  maxres = max(abs(res)), termcd = r1$sol$termcd, iter = r1$sol$iter))
    r2  <- nl(r1$root, sc)
    res <- tryCatch(eval_f(r2$root, pv) / sc, error = function(e) setNames(rep(Inf, n_dep), dependent))
    list(root = r2$root, sol = r2$sol, res = res,
         maxres = max(abs(res)), termcd = r2$sol$termcd, iter = r1$sol$iter + r2$sol$iter)
  }

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

    hessian <- NULL
    if (want_d2) {
      H_all <- eval_H(root, pv)
      if (is.null(H_all))
        stop("Pimpl(deriv2 = TRUE) requires hess(); rebuild with deriv2 = TRUE.", call. = FALSE)

      f_xx <- H_all[seq_len(n_dep), dependent, dependent, drop = FALSE]
      f_xp <- H_all[seq_len(n_dep), dependent, parms_all, drop = FALSE]
      f_pp <- H_all[seq_len(n_dep), parms_all, parms_all, drop = FALSE]

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
      hessian <- hess_arr
    }

    if (!is.null(dP)) {
      dPsub <- submatrix(dP, rows = colnames(jacobian))
      th    <- colnames(dPsub); n_th <- length(th)
      if (want_d2 && !is.null(hessian)) {
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


  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
    if (deriv2 && !emit_d2)
      stop("Pimpl was built with deriv2 = FALSE; rebuild with deriv2 = TRUE.", call. = FALSE)
    if (!emit_d1) deriv <- FALSE
    if (deriv2 && !deriv) deriv <- TRUE

    p   <- pars
    dP  <- attr(p, "deriv")
    dP2 <- if (deriv2) attr(p, "deriv2") else NULL

    top  <- nleq[intersect(names(nleq), nleqslv_top)]
    ctrl <- nleq[setdiff(names(nleq), nleqslv_top)]
    positive  <- ms$positive
    nStarts   <- ms$nStarts
    debugPlot <- ms$debugPlot
    ftol      <- nleq$ftol

    if (!is.null(fixed)) {
      p <- p[!names(p) %in% names(fixed)]
      p <- c(p, fixed)
    }
    emptypars <- setdiff(names(p), c(dependent, names(fixed)))
    miss <- setdiff(dependent, names(p)); if (length(miss)) p[miss] <- 1
    pv <- p[parms_all]

    zero_vec <- if (length(zero_states))
      setNames(rep(0, length(zero_states)), zero_states) else NULL

    ## Multistart: cache -> user init -> random sweep, cycling through
    ## alternative nleqslv globalizations.
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
    if ((is.null(best) || best$maxres > ftol) && nStarts > 1L) {
      lo <- .expand_bounds(ms$lower, dependent, 0)
      hi <- .expand_bounds(ms$upper, dependent, 10)
      if (positive) {
        lo_log <- log(pmax(lo, .Machine$double.eps))
        hi_log <- log(pmax(hi, .Machine$double.eps * 10))
      }
      altMethods <- list(top,
                         modifyList(top, list(global = "pwldog")),
                         modifyList(top, list(method = "Broyden", global = "dbldog")))
      for (i in seq_len(nStarts)) {
        x0_rand <- if (positive)
          setNames(exp(runif(n_dep, lo_log, hi_log)), dependent)
        else
          setNames(runif(n_dep, lo, hi), dependent)
        try_solve(x0_rand, paste0("random_", i),
                  override = altMethods[[((i - 1L) %% length(altMethods)) + 1L]])
        if (!is.null(best) && best$maxres <= ftol) break
      }
    }

    if (debugPlot && nStarts > 1L && length(log)) {
      log_df <- do.call(rbind, lapply(log, as.data.frame, stringsAsFactors = FALSE))
      log_df <- log_df[order(log_df$maxres), ]; log_df$rank <- seq_len(nrow(log_df))
      waterfall_plot(log_df, ftol)
    }

    if (is.null(best))
      stop("Pimpl: all ", length(log), " solve attempt(s) failed (no usable ",
           "nleqslv result). The residual system may be degenerate at the ",
           "current parameter values, or the basin may lie outside ",
           "[", format(min(ms$lower)), ", ", format(max(ms$upper)),
           "]. Increase `controlsMS$nStarts`, widen ",
           "`controlsMS$lower`/`upper`, or check the model.",
           call. = FALSE)

    if (best$maxres > ftol) {
      ord  <- order(abs(best$res), decreasing = TRUE)
      topn <- min(10L, length(best$res))
      tbl  <- paste0(sprintf("  %-25s  rel|f| = %s",
                             names(best$res)[ord[1:topn]],
                             formatC(abs(best$res[ord[1:topn]]), format = "e", digits = 2)),
                     collapse = "\n")
      stop("Pimpl: best relative residual ", formatC(best$maxres, format = "e", digits = 2),
           " exceeds ftol = ", formatC(ftol, format = "e", digits = 2),
           " after ", length(log), " attempt(s) (nleqslv termcd ",
           best$sol$termcd, "). No steady state reached.\n",
           "Largest relative residuals (top ", topn, "):\n", tbl, call. = FALSE)
    }

    sol <- best$sol; root <- best$root

    out <- c(root, zero_vec,
             p[setdiff(names(p), c(names(root), zero_states))])

    if (keep.root) cache$guess <- out

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

  attr(p2p, "equations")   <- as.eqnvec(all_exprs)
  attr(p2p, "parameters")  <- parameters
  attr(p2p, "modelname")   <- modelname
  attr(p2p, "compileInfo") <- collectCompileInfo(PEval$func, PEval$jac, PEval$hess)
  attr(p2p, "resetWarmStart") <- local({
    cache_ref <- cache; mn <- modelname; cond <- condition
    function() {
      cache_ref$guess <- NULL
      paste0("Pimpl(", mn, if (!is.null(cond)) paste0(":", cond) else "", ")")
    }
  })
  parfn(p2p, parameters, condition)
}


#' Steady-state transformation with conserved moieties expressed as totals
#'
#' Builds the [Pequil] parfn for the `expressInTotals = TRUE` case. The full
#' (uneliminated) system is integrated to its fixed point, so mass-action
#' positivity and the conservation laws hold automatically. Each integration
#' start is a random non-negative point on the conservation manifold fixed by
#' the totals (the mass is distributed across the moiety, not placed in a
#' single species), anchored at a private species and spread through the null
#' space of the CQ matrix. Total sensitivities are obtained from the pivot
#' species' initial-condition sensitivities via the constant map
#' `C[, pivots]^{-1}`, independent of the starting distribution.
#'
#' @param norm Record from [.normalize_ss_inputs] with `fullsystem = TRUE`.
#' @param ms Resolved multistart controls.
#' @param emit_d1,emit_d2 Whether first/second-order sensitivities are built.
#' @param attach.input,keep.root,controlsODE,compile,modelname,condition,verbose,start.time,end.time
#'   As in [Pequil].
#' @param dotArgs Extra arguments forwarded to [CppODE::CppODE].
#' @return A [parfn].
#' @keywords internal
.Pequil_totals <- function(norm, ms, emit_d1, emit_d2, attach.input, keep.root,
                           controlsODE, compile, modelname, condition, verbose,
                           start.time, end.time, dotArgs) {
  f           <- norm$trafo
  states      <- norm$states
  zero_states <- norm$zero_states
  dependent   <- norm$dependent
  parameters  <- norm$parameters
  totals      <- norm$totals
  C_mat       <- norm$C_mat
  pivots      <- norm$pivots
  moiety      <- norm$moiety_species
  n_dep       <- length(dependent)
  total_names <- names(totals)
  nonmoiety   <- setdiff(dependent, moiety)

  model_params <- setdiff(getSymbols(unclass(f[dependent])), dependent)
  Cp_inv <- solve(C_mat[, pivots, drop = FALSE])
  dimnames(Cp_inv) <- list(pivots, total_names)

  shared_count <- colSums(abs(C_mat) > 1e-12)
  private <- vapply(seq_along(totals), function(i) {
    cand <- colnames(C_mat)[abs(C_mat[i, ]) > 1e-12 & shared_count == 1L]
    if (length(cand)) sort(cand)[1L] else NA_character_
  }, character(1))
  if (anyNA(private))
    stop("Pequil(expressInTotals = TRUE): conserved quantity '",
         total_names[which(is.na(private))[1L]], "' has no private species ",
         "(a form occurring only in that total) to seed its initial condition. ",
         "Fully-shared moiety systems are unsupported in totals mode; use ",
         "expressInTotals = FALSE.", call. = FALSE)
  private_coef <- vapply(seq_along(totals), function(i) C_mat[i, private[i]], numeric(1))

  sv <- svd(C_mat, nu = 0L, nv = ncol(C_mat))
  rk <- sum(sv$d > max(dim(C_mat)) * .Machine$double.eps * sv$d[1L])
  Cnull <- if (rk < ncol(C_mat)) sv$v[, (rk + 1L):ncol(C_mat), drop = FALSE] else
    matrix(0, ncol(C_mat), 0L)
  rownames(Cnull) <- colnames(C_mat)

  moiety_ic <- function(tot) {
    x <- setNames(rep(0, length(moiety)), moiety)
    x[private] <- as.numeric(tot[total_names]) / private_coef
    if (ncol(Cnull)) {
      dir  <- setNames(as.numeric(Cnull %*% runif(ncol(Cnull), -1, 1)), rownames(Cnull))[names(x)]
      neg  <- dir < -1e-12
      amax <- if (any(neg)) min(-x[neg] / dir[neg]) else max(x)
      x <- pmax(x + runif(1L, 0, amax) * dir, 0)
    }
    x
  }

  if (is.null(modelname)) modelname <- "equil_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

  base <- c(list(rhs = unclass(f[dependent]), rootfunc = "equilibrate", compile = compile,
                 outdir = getwd(), useDenseOutput = FALSE, verbose = verbose), dotArgs)
  fixed_states <- setdiff(dependent, pivots)
  model    <- do.call(CppODE::CppODE, c(base, list(deriv = FALSE, deriv2 = FALSE,
                                                   modelname = modelname)))
  model_s  <- if (emit_d1)
    do.call(CppODE::CppODE, c(base, list(deriv = TRUE, deriv2 = FALSE,
                                         modelname = paste0(modelname, "_s"),
                                         fixed = fixed_states))) else NULL
  model_s2 <- if (emit_d2)
    do.call(CppODE::CppODE, c(base, list(deriv = TRUE, deriv2 = TRUE,
                                         modelname = paste0(modelname, "_s2"),
                                         fixed = fixed_states))) else NULL
  all_sens <- if (emit_d1) attr(model_s, "dimNames")$sens else character(0)

  kin_sens   <- intersect(model_params, all_sens)
  outer_sens <- c(total_names, kin_sens)
  Tmat <- matrix(0, length(all_sens), length(outer_sens),
                 dimnames = list(all_sens, outer_sens))
  Tmat[pivots, total_names] <- Cp_inv
  if (length(kin_sens)) Tmat[cbind(kin_sens, kin_sens)] <- 1

  ode_ctrl <- modifyList(list(abstol = 1e-6, reltol = 1e-6, maxsteps = 1e6L,
                              maxprogress = 100L, hini = 0, roottol = 1e-6, maxroot = 1L),
                         controlsODE)
  controls <- c(list(keep.root = keep.root, attach.input = attach.input,
                     start.time = start.time, end.time = end.time), ode_ctrl)

  cache <- new.env(parent = emptyenv())
  cache$yini <- cache$last_hash <- cache$last_result <- NULL

  default_sens <- matrix(0, n_dep, length(all_sens), dimnames = list(dependent, all_sens))
  if (length(pivots)) default_sens[cbind(pivots, pivots)] <- 1
  default_sens2 <- if (emit_d2)
    array(0, c(n_dep, length(all_sens), length(all_sens)),
          dimnames = list(dependent, all_sens, all_sens)) else NULL

  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
    if (deriv2 && !emit_d2)
      stop("Pequil(deriv2 = TRUE) requires the model to be built with deriv2 = TRUE.",
           call. = FALSE)
    if (!emit_d1) deriv <- FALSE
    if (deriv2 && !deriv) deriv <- TRUE
    p   <- pars
    dP  <- attr(p, "deriv")
    dP2 <- if (deriv2) attr(p, "deriv2") else NULL
    if (!is.null(fixed)) { p <- p[!names(p) %in% names(fixed)]; p <- c(p, fixed) }

    tot       <- p[total_names]
    emptypars <- setdiff(names(p), c(dependent, names(fixed)))

    pv_hash <- NULL
    if (keep.root) {
      pv_hash <- digest::digest(list(tot, p[model_params], fixed, deriv, deriv2),
                                algo = "xxhash64")
      if (!is.null(cache$last_hash) && identical(pv_hash, cache$last_hash) &&
          !is.null(cache$last_result))
        return(cache$last_result)
    }

    sens_model <- if (deriv2) model_s2 else if (deriv) model_s else model
    run_attempt <- function(y0) {
      tryCatch(
        CppODE::solveODE(
          sens_model, times = c(controls$start.time, controls$end.time),
          parms = c(y0, p[model_params]),
          sens1ini = if (deriv) default_sens else NULL,
          sens2ini = if (deriv2) default_sens2 else NULL,
          roottol = controls$roottol, abstol = controls$abstol, reltol = controls$reltol,
          maxsteps = as.integer(controls$maxsteps),
          maxprogress = as.integer(controls$maxprogress),
          hini = controls$hini, maxroot = as.integer(controls$maxroot),
          onFailure = "silent"),
        error = function(e) NULL)
    }
    is_success <- function(r) {
      if (is.null(r) || is.null(r$diagnostics)) return(FALSE)
      rc <- r$diagnostics$return_code
      if (!is.null(rc) && rc < 0L) return(FALSE)
      length(r$time) >= 1L &&
        r$time[length(r$time)] < controls$end.time - .Machine$double.eps
    }

    y0 <- setNames(rep(1, n_dep), dependent)
    y0[moiety] <- moiety_ic(tot)
    if (keep.root && !is.null(cache$yini)) y0[nonmoiety] <- cache$yini[nonmoiety]
    res <- run_attempt(y0)

    if (!is_success(res) && ms$nStarts > 1L) {
      lo <- .expand_bounds(ms$lower, nonmoiety, 0)
      hi <- .expand_bounds(ms$upper, nonmoiety, 10)
      if (ms$positive) {
        lo_log <- log(pmax(lo, .Machine$double.eps))
        hi_log <- log(pmax(hi, .Machine$double.eps * 10))
      }
      for (i in seq_len(ms$nStarts - 1L)) {
        y0 <- setNames(rep(0, n_dep), dependent)
        y0[moiety] <- moiety_ic(tot)
        y0[nonmoiety] <- if (ms$positive)
          exp(runif(length(nonmoiety), lo_log, hi_log))
        else runif(length(nonmoiety), lo, hi)
        res <- run_attempt(y0)
        if (is_success(res)) break
      }
    }

    zero_vec <- if (length(zero_states))
      setNames(rep(0, length(zero_states)), zero_states) else NULL
    if (!is_success(res)) {
      rc <- if (!is.null(res) && !is.null(res$diagnostics))
              as.character(res$diagnostics$return_code) else "exception"
      stop("Pequil: no steady state reached after ", ms$nStarts,
           " integration attempt(s) (last return_code: ", rc, "). Either no stable ",
           "fixed point exists in this regime, the totals admit no non-negative ",
           "steady state, or the ODE is too stiff. Increase `controlsMS$nStarts`, ",
           "widen `controlsMS$lower`/`upper`, or relax `controlsODE`.", call. = FALSE)
    }

    last <- length(res$time)
    digits <- floor(-log10(controls$roottol)) + 1L
    root <- setNames(round(res$variable[last, ], digits), dependent)
    out  <- if (attach.input)
              c(root, zero_vec, p[setdiff(names(p), c(dependent, zero_states))])
            else c(root, zero_vec)
    if (keep.root) cache$yini <- root

    if (!deriv || is.null(res$sens1)) {
      result <- as.parvec(out, deriv = NULL, deriv2 = NULL)
    } else {
      sens_outer <- matrix(res$sens1[last, , ], n_dep, length(all_sens),
                           dimnames = list(dependent, all_sens)) %*% Tmat
      input_cols <- setdiff(names(p), c(dependent, names(fixed)))
      jacobian <- matrix(0, length(out), length(input_cols),
                         dimnames = list(names(out), input_cols))
      if (attach.input) {
        idx <- intersect(emptypars, input_cols)
        if (length(idx)) jacobian[cbind(idx, idx)] <- 1
      }
      sc <- intersect(input_cols, colnames(sens_outer))
      if (length(sc)) jacobian[dependent, sc] <- sens_outer[dependent, sc, drop = FALSE]

      hess_attr <- NULL
      if (deriv2 && !is.null(res$sens2)) {
        ns <- length(all_sens)
        sens2 <- array(res$sens2[last, , , ], c(n_dep, ns, ns),
                       dimnames = list(dependent, all_sens, all_sens))
        hess_arr <- array(0, c(length(out), length(input_cols), length(input_cols)),
                          dimnames = list(names(out), input_cols, input_cols))
        oc <- colnames(Tmat)
        hess_arr[dependent, oc, oc] <- t(Tmat) %bmm% (sens2 %bmm% Tmat)
        if (!is.null(dP)) {
          dPsub <- submatrix(dP, rows = input_cols)
          th    <- colnames(dPsub); n_th <- length(th)
          new_hess <- t(dPsub) %bmm% hess_arr %bmm% dPsub
          dimnames(new_hess) <- list(names(out), th, th)
          if (!is.null(dP2)) {
            dP2sub <- dP2[input_cols, th, th, drop = FALSE]
            new_hess <- new_hess + array(
              jacobian %*% matrix(dP2sub, length(input_cols), n_th * n_th),
              c(length(out), n_th, n_th), dimnames = list(names(out), th, th))
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

  attr(p2p, "equations")   <- as.eqnvec(f[dependent])
  attr(p2p, "parameters")  <- parameters
  attr(p2p, "modelname")   <- modelname
  attr(p2p, "compileInfo") <- collectCompileInfo(model, model_s, model_s2)
  attr(p2p, "resetWarmStart") <- local({
    cache_ref <- cache; mn <- modelname; cond <- condition
    function() {
      cache_ref$yini <- cache_ref$last_hash <- cache_ref$last_result <- NULL
      paste0("Pequil(", mn, if (!is.null(cond)) paste0(":", cond) else "", ")")
    }
  })
  parfn(p2p, parameters, condition)
}


#' Parameter transformation (steady states via pre-equilibration)
#'
#' Returns a [parfn] over the outer inputs. On call, the parfn integrates
#' the ODE from `start.time` to `end.time`, warm starting from the cached
#' root when available and falling back to multistart on the initial
#' conditions otherwise. The Jacobian (and Hessian when `deriv2 = TRUE`)
#' come from CppODE's analytical sensitivity integration. Conserved
#' quantities are detected; how they are parametrised is controlled by
#' `expressInTotals`.
#'
#' @param trafo Named character / [eqnvec] / [eqnlist].
#' @param parameters Outer parameters; listed states pass through as
#'   initial conditions instead of being integrated.
#' @param forcings Forcing names; zeroed and removed.
#' @param condition Condition label.
#' @param attach.input Append pass-through inputs to the output.
#' @param keep.root Warm-start subsequent calls and re-use the cached
#'   result when inputs are unchanged.
#' @param expressInTotals If `FALSE` (default), the eliminated species per
#'   conserved quantity becomes a pass-through parameter held constant during
#'   integration. If `TRUE`, the full (uneliminated) system is integrated from
#'   initial conditions that distribute each conserved total across its moiety
#'   on the conservation manifold, so mass-action positivity and the
#'   conservation laws hold automatically (no reconstructed species can turn
#'   negative). The totals become outer parameters; their sensitivities are
#'   obtained from the pivot species' initial-condition sensitivities. Requires
#'   every conserved quantity to have a private species (a form occurring only
#'   in that total) to anchor a feasible initial condition.
#' @param controlsODE Overrides for the ODE solver controls.
#' @param start.time,end.time Integration window; the equilibrate root
#'   event fires before `end.time` on success.
#' @param controlsMS Multistart controls. Recognised keys: `nStarts`
#'   (default `10L`; `1L` disables multistart), `positive` (default
#'   `TRUE`; draws log-uniform random initial conditions over
#'   `[lower, upper]`), `lower`/`upper` (scalar or named vector of
#'   bounds for the random sweep). Under `expressInTotals = TRUE` the sweep
#'   covers the non-conserved states only; moiety species are restarted on
#'   the conservation manifold fixed by the totals.
#' @param compile,modelname,verbose Forwarded to [CppODE::CppODE].
#' @param deriv,deriv2 Attach first/second-order sensitivities; `deriv2`
#'   requires the model built with `deriv2 = TRUE`.
#' @param ... Forwarded to [CppODE::CppODE].
#'
#' @return A [parfn].
#' @seealso [Pexpl], [Pimpl], [P].
#' @import CppODE
#' @importFrom digest digest
#' @export
Pequil <- function(trafo, parameters = NULL, forcings = NULL, condition = NULL,
                   attach.input = TRUE, start.time = 0, end.time = 1e10,
                   keep.root = TRUE, expressInTotals = FALSE,
                   controlsODE = list(), controlsMS = list(),
                   compile = FALSE, modelname = NULL, verbose = FALSE,
                   deriv = TRUE, deriv2 = FALSE, ...) {

  ms <- modifyList(list(nStarts = 10L, positive = TRUE,
                        lower = 1e-5, upper = 1e5), controlsMS)

  emit_d1 <- isTRUE(deriv)
  emit_d2 <- isTRUE(deriv2)
  if (emit_d2 && !emit_d1)
    stop("Pequil(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)

  norm <- .normalize_ss_inputs(trafo, parameters, forcings,
                               expressInTotals = expressInTotals,
                               fullsystem = isTRUE(expressInTotals))

  if (!is.null(norm$pivots) && length(norm$pivots))
    return(.Pequil_totals(norm, ms, emit_d1, emit_d2, attach.input, keep.root,
                          controlsODE, compile, modelname, condition, verbose,
                          start.time, end.time, list(...)))

  f           <- norm$trafo
  states      <- norm$states
  zero_states <- norm$zero_states
  dependent   <- norm$dependent
  parameters  <- norm$parameters
  parms_all   <- intersect(norm$parms_all, getSymbols(norm$trafo[norm$dependent]))
  n_dep       <- length(dependent)
  f_red       <- f[dependent]

  if (is.null(modelname)) modelname <- "equil_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

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
  cache$yini <- cache$sensini <- cache$sens2ini <-
    cache$last_hash <- cache$last_result <- NULL

  default_sens <- matrix(0, n_dep, length(all_sens),
                         dimnames = list(dependent, all_sens))
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

    sens_model <- if (deriv2) model_s2 else if (deriv) model_s else model

    run_attempt <- function(y0_dep, use_cache_sens) {
      s1 <- if (deriv  && keep.root && use_cache_sens && !is.null(cache$sensini))
              cache$sensini[, active_sens, drop = FALSE]
            else if (deriv)
              default_sens[, active_sens, drop = FALSE]
      s2 <- if (deriv2 && keep.root && use_cache_sens && !is.null(cache$sens2ini))
              cache$sens2ini[, active_sens, active_sens, drop = FALSE]
            else if (deriv2)
              default_sens2[, active_sens, active_sens, drop = FALSE]
      tryCatch(
        CppODE::solveODE(
          sens_model,
          times = c(controls$start.time, controls$end.time),
          parms = c(y0_dep, p[parms_all]),
          sens1ini = s1, sens2ini = s2,
          fixed = if (deriv || deriv2) fixed_char,
          roottol = controls$roottol, abstol = controls$abstol, reltol = controls$reltol,
          maxsteps = as.integer(controls$maxsteps),
          maxprogress = as.integer(controls$maxprogress),
          hini = controls$hini, maxroot = as.integer(controls$maxroot),
          onFailure = "silent"),
        error = function(e) NULL)
    }
    is_success <- function(r) {
      if (is.null(r) || is.null(r$diagnostics)) return(FALSE)
      rc <- r$diagnostics$return_code
      if (!is.null(rc) && rc < 0L) return(FALSE)
      length(r$time) >= 1L &&
        r$time[length(r$time)] < controls$end.time - .Machine$double.eps
    }

    res <- run_attempt(p[dependent], use_cache_sens = TRUE)

    if (!is_success(res) && ms$nStarts > 1L) {
      lo <- .expand_bounds(ms$lower, dependent, 0)
      hi <- .expand_bounds(ms$upper, dependent, 10)
      if (ms$positive) {
        lo_log <- log(pmax(lo, .Machine$double.eps))
        hi_log <- log(pmax(hi, .Machine$double.eps * 10))
      }
      for (i in seq_len(ms$nStarts - 1L)) {
        y0_rand <- if (ms$positive)
          setNames(exp(runif(n_dep, lo_log, hi_log)), dependent)
        else
          setNames(runif(n_dep, lo, hi), dependent)
        res <- run_attempt(y0_rand, use_cache_sens = FALSE)
        if (is_success(res)) break
      }
    }

    zero_vec <- if (length(zero_states))
      setNames(rep(0, length(zero_states)), zero_states) else NULL

    if (!is_success(res)) {
      rc <- if (!is.null(res) && !is.null(res$diagnostics))
              as.character(res$diagnostics$return_code) else "exception"
      stop("Pequil: no steady state reached after ", ms$nStarts,
           " integration attempt(s) (last return_code: ", rc, "). ",
           "Either no stable fixed point exists in this parameter regime, ",
           "the basin lies outside [", format(min(ms$lower)), ", ",
           format(max(ms$upper)), "], or the ODE is too stiff for the current ",
           "`controlsODE`. Increase `controlsMS$nStarts`, widen ",
           "`controlsMS$lower`/`upper`, or relax ",
           "`controlsODE$abstol`/`reltol`/`maxprogress`.", call. = FALSE)
    }

    last <- length(res$time)

    digits <- floor(-log10(controls$roottol)) + 1L
    root <- setNames(round(res$variable[last, ], digits), dependent)
    out  <- if (attach.input)
              c(root, zero_vec, p[setdiff(names(p), c(dependent, zero_states))])
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

  attr(p2p, "equations")   <- as.eqnvec(f_red)
  attr(p2p, "parameters")  <- parameters
  attr(p2p, "modelname")   <- modelname
  attr(p2p, "compileInfo") <- collectCompileInfo(model, model_s, model_s2)
  attr(p2p, "resetWarmStart") <- local({
    cache_ref <- cache; mn <- modelname; cond <- condition
    function() {
      cache_ref$yini <- cache_ref$sensini <- cache_ref$sens2ini <-
        cache_ref$last_hash <- cache_ref$last_result <- NULL
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
    .currentTrafo   <- mytrafo[[i]]
    .currentSymbols <- if (is.null(.currentTrafo)) NULL else getSymbols(.currentTrafo)
    row <- if (is.list(trafo)) tree[names(mytrafo)[i], , drop = FALSE]
           else tree[1, , drop = FALSE]
    if (!is.null(conditionMatch) && !str_detect(rownames(row), conditionMatch))
      return(.currentTrafo)
    with(row, do.call(repar,
      c(list(expr = expr, trafo = .currentTrafo, reset = TRUE), eval(dots))))
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
    .currentTrafo   <- mytrafo[[i]]
    .currentSymbols <- if (is.null(.currentTrafo)) NULL else getSymbols(.currentTrafo)
    row <- if (is.list(trafo)) tree[names(mytrafo)[i], , drop = FALSE]
           else tree[1, , drop = FALSE]
    if (!is.null(conditionMatch) && !str_detect(rownames(row), conditionMatch))
      return(.currentTrafo)
    with(row, {
      ## Caller may pass logical dots to gate substitution per condition,
      ## and non-logical dots to substitute symbols in `expr`. Logical dots
      ## are stripped before forwarding to `repar`.
      .apply <- function() {
        d <- eval(dots)
        if (!length(d)) return(do.call(repar, list(expr = expr, trafo = .currentTrafo)))
        d_eval  <- lapply(d, function(x) eval.parent(x, 3))
        is_log  <- vapply(d_eval, is.logical, logical(1))
        gate    <- do.call(c, d[is_log])
        if (!is.null(gate) && any(!gate)) return(.currentTrafo)
        do.call(repar, c(list(expr = expr, trafo = .currentTrafo), d_eval[!is_log]))
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