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
#' @param cores Number of cores for parallel method call per condition
#' @param verbose Print out information during compilation
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
              cores = detectFreeCores(), verbose = FALSE, ...) {
  
  if (is.null(trafo)) return()
  
  # Wrap single trafo in named list
  if (!is.list(trafo) || inherits(trafo, "eqnlist") || inherits(trafo, "eqnvec")) {
    trafo_list <- list(trafo)
    names(trafo_list) <- condition
  } else { trafo_list <- trafo }
  
  if (Sys.info()[['sysname']] == "Windows") cores <- 1
  
  Reduce("+", mclapply(seq_along(trafo_list), function(i) {
    
    tr <- trafo_list[[i]]
    cond <- names(trafo_list[i])
    
    # Auto-select method per entry
    m <- if (!is.null(method)) match.arg(method, c("explicit", "implicit", "equilibrate"))
    else if (inherits(tr, "eqnlist")) "equilibrate" 
    else "explicit"
    
    switch(m,
           explicit    = Pexpl(as.eqnvec(tr), parameters = parameters, condition = cond, compile = compile, 
                               modelname = modelname, verbose = verbose, ...),
           implicit    = Pimpl(trafo = as.eqnvec(tr), parameters = parameters, condition = cond, 
                               compile = compile, modelname = modelname, verbose = verbose, ...),
           equilibrate = Pequil(trafo = tr, parameters = parameters, condition = cond, 
                                compile = compile, modelname = modelname, verbose = verbose, ...)
    )
    
  }, mc.cores = min(detectFreeCores(), cores)))
}


#' Parameter transformation (explicit)
#'
#' Constructs a parameter transformation function that maps **outer parameters**
#' \eqn{p_{\text{outer}}} to **inner parameters** \eqn{p_{\text{inner}}}
#' according to symbolic expressions.
#'
#' @description
#' The explicit parameter transformation defines a direct, algebraic mapping
#'
#' \deqn{p_{\text{inner}} = \mathrm{parfn}(p_{\text{outer}}),}
#'
#' where \eqn{\mathrm{parfn}} is a vector-valued function composed from symbolic
#' expressions. Each element of `trafo` defines one component of
#' \eqn{p_{\text{inner}}}.
#'
#' The **Jacobian** is obtained by **symbolic differentiation**.
#' It is attached as attribute `"deriv"`
#' to the resulting function output and automatically composed when
#' transformations are combined via the [parfn] interface.
#'
#' @param trafo `eqnvec` or named character vector.
#' Names correspond to **inner parameters**; each element defines how it depends
#' on **outer parameters**.
#' @param parameters Character vector of outer parameter names. If omitted,
#' all symbols in `trafo` are used.
#' @param attach.input Logical. If `TRUE`, include unchanged input parameters
#' in the output vector (identity mapping).
#' @param condition Character label for which the transformation is generated.
#' @param compile Logical. If `TRUE`, compile the transformation via [funCpp]
#' for faster evaluation.
#' @param modelname Base name for generated C++ code if `compile = TRUE`.
#' @param verbose Logical. Print compiler messages.
#' @param derivMode Character. Selects the derivative backend used by [funCpp]
#'   to evaluate the transformation Jacobian. One of `"symbolic"` (default,
#'   classical SymPy Jacobian — appropriate for the typically small parameter
#'   transformations), `"ad"` (forward-mode automatic differentiation via
#'   `jac_chain`; requires `compile = TRUE`), or `"none"` (no derivatives).
#'
#' @return
#' A function of class [parfn].
#'
#' @seealso
#' [Pimpl] for implicit (steady-state) parameter transformations,
#' [P] for automatic mode selection.
#'
#' @importFrom CppODE funCpp
#' @export
Pexpl <- function(trafo, parameters = NULL, attach.input = FALSE, condition = NULL,
                  compile = FALSE, modelname = NULL, verbose = FALSE,
                  derivMode = c("symbolic", "ad", "none")) {

  derivMode <- match.arg(derivMode)

  # Determine parameter sets
  if (is.null(parameters)) {
    parameters <- getSymbols(trafo)
  } else {
    identity <- parameters[!(parameters %in% names(trafo))]
    names(identity) <- identity
    trafo <- c(trafo, identity)
    parameters <- getSymbols(trafo)
  }

  # Model name with condition label
  if (is.null(modelname)) modelname <- "expl_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

  # Build compiled (or fallback R) evaluator for transformation
  PEval <- suppressWarnings(
    CppODE::funCpp(
      unclass(trafo),
      variables  = NULL,
      parameters = parameters,
      fixed      = NULL,
      compile    = compile,
      modelname  = modelname,
      outdir     = getwd(),
      verbose    = verbose,
      convenient = FALSE,
      derivMode  = derivMode
    )
  )

  fun       <- PEval$func
  jac       <- PEval$jac
  jac_chain <- PEval$jac_chain
  use_ad    <- derivMode == "ad"

  # Define returned parameter transformation function
  p2p <- function(pars, fixed = NULL, deriv = TRUE) {

    # Prepare pars
    p <- c(pars, fixed)

    Jac <- NULL
    if (use_ad && deriv && !is.null(jac_chain)) {
      # AD path: jac_chain returns y and dy already chain-ruled w.r.t. theta.
      out <- jac_chain(NULL, p, dX = NULL, dP = attr(pars, "deriv"),
                       attach.input = attach.input, fixed = names(fixed))
      pinnerVal <- out$y[, 1]
      if (!is.null(out$dy)) {
        Jac <- matrix(out$dy, dim(out$dy)[1], dim(out$dy)[2],
                      dimnames = list(dimnames(out$dy)[[1]], dimnames(out$dy)[[2]]))
      }
    } else {
      # Symbolic path (also serves "both" and "none" via NULL jac).
      pinnerVal <- fun(NULL, p, attach.input = attach.input, fixed = names(fixed))[,]
      if (deriv && !is.null(jac)) {
        Jac <- as.matrix(jac(NULL, p, attach.input = attach.input, fixed = names(fixed))[,,1])
        dP <- attr(pars, "deriv")
        if (!is.null(dP)) {
          Jac <- Jac %*% dP[colnames(Jac), , drop = FALSE]
          dimnames(Jac) <- list(names(pinnerVal), colnames(dP))
        }
      }
    }

    if (any(is.nan(pinnerVal))) {
      stop(
        paste0(
          "The following inner parameter(s) evaluate to NaN:\n\t",
          paste0(names(pinnerVal)[is.nan(pinnerVal)], collapse = "\n\t"),
          ".\nLikely cause: division by zero or missing inputs."
        )
      )
    }

    # Assemble result
    pinner <- as.parvec(pinnerVal, deriv = if (deriv && !is.null(Jac)) Jac[rowSums(Jac != 0) > 0, , drop = FALSE] else FALSE)

    if (attach.input && !all(names(pars) %in% names(pinnerVal))) {
      pinner <- c(pinner,
                  as.parvec(pars[setdiff(names(pars), names(pinnerVal))],
                            deriv = if (deriv) NULL else FALSE))
    }

    pinner
  }

  attr(p2p, "equations")  <- as.eqnvec(trafo)
  attr(p2p, "parameters") <- parameters
  attr(p2p, "modelname")  <- modelname
  attr(p2p, "compileInfo") <- collectCompileInfo(fun, jac, jac_chain)

  parfn(p2p, parameters, condition)
}


#' Detect conserved quantities and substitute eliminated species
#'
#' For each conserved quantity (e.g. A + B = const), one species is eliminated
#' from the equations by substitution. The eliminated species is replaced by
#' \code{(total_i - rest)}, where \code{total_i} becomes a new parameter.
#'
#' @param smatrix Stoichiometric matrix (or NULL).
#' @param f Named character vector / eqnvec of equations.
#' @param states Character vector of state names.
#' @param parameters Character vector of user-specified parameters (modified in place).
#' @return A list with components:
#'   \item{f}{Modified equations with substitutions applied.}
#'   \item{parameters}{Updated parameter vector (eliminated species added).}
#'   \item{cq_info}{List of per-CQ info: \code{total_name}, \code{elim_state},
#'     \code{recon_expr} (the reconstruction expression as eqnvec for funCpp),
#'     \code{recon_vars} (dependent species appearing in reconstruction),
#'     \code{recon_parms} (parameters appearing in reconstruction, including total).}
#'   \item{elim_states}{Character vector of all eliminated state names.}
#' @keywords internal
.detect_and_substitute_cq <- function(smatrix, f, states, parameters) {
  
  cq_info <- list()
  elim_states <- character(0)
  
  if (!is.null(smatrix)) {
    cq <- conservedQuantities(smatrix)
    if (!is.null(cq) && nrow(cq) > 0) {
      if (is.null(parameters)) parameters <- character(0)
      
      substitutions <- list()  # elim_state -> replacement expression
      
      for (i in seq_len(nrow(cq))) {
        cq_expr <- as.character(cq[i, 1])
        cq_species <- intersect(getSymbols(cq_expr), states)
        n_need <- length(cq_species) - 1
        already <- intersect(cq_species, parameters)
        
        # Auto-fill: ensure n-1 species from this CQ are in parameters
        if (length(already) < n_need) {
          fill <- setdiff(cq_species, already)[seq_len(n_need - length(already))]
          parameters <- union(parameters, fill)
        }
        
        elim <- intersect(cq_species, parameters)
        
        for (e in elim) {
          total_name <- if (length(elim) > 1) paste0("total_", i, "_", e) else paste0("total_", i)
          
          # Build replacement: e = total_name - (cq_expr with e set to 0)
          rest_expr <- replaceSymbols(e, "0", cq_expr)
          recon_expr <- paste0("(", total_name, " - (", rest_expr, "))")
          
          substitutions[[e]] <- recon_expr
          elim_states <- c(elim_states, e)
          
          # Identify which dependent and which parameters appear in recon_expr
          recon_syms <- getSymbols(recon_expr)
          
          cq_info[[length(cq_info) + 1]] <- list(
            total_name  = total_name,
            elim_state  = e,
            recon_expr  = setNames(recon_expr, e)
          )
        }
      }
      
      # Apply all substitutions to equations at once
      if (length(substitutions) > 0) {
        f <- replaceSymbols(names(substitutions),
                            unname(unlist(substitutions)), f)
      }
    }
  } else if (!is.null(parameters)) {
    # eqnvec/character: warn if eliminated species still appear unsubstituted
    param_states <- intersect(parameters, states)
    if (length(param_states)) {
      remaining_eqs <- f[setdiff(states, param_states)]
      if (length(remaining_eqs) > 0) {
        still_present <- param_states[vapply(param_states, function(ps) {
          any(ps %in% getSymbols(remaining_eqs))
        }, logical(1))]
        if (length(still_present))
          warning("States in 'parameters' still appear in dependent equations: ",
                  paste(still_present, collapse = ", "),
                  ". Without an eqnlist, conserved quantities cannot be substituted ",
                  "automatically. Please substitute manually or provide an eqnlist.",
                  call. = FALSE)
      }
    }
  }
  
  list(f = f, parameters = parameters, cq_info = cq_info, elim_states = elim_states)
}

#' Reconstruct eliminated species via compiled function
#' @keywords internal
.recon_elim <- function(recon_eval, root, pv, dependent, parms_all, elim_states) {
  if (is.null(recon_eval) || length(elim_states) == 0) return(numeric(0))
  vals <- recon_eval$func(
    matrix(root[dependent], 1, dimnames = list(NULL, dependent)),
    pv[parms_all]
  )[, 1]
  names(vals) <- elim_states
  vals
}



#' Parameter transformation (implicit, root-finding)
#'
#' @param trafo Named character vector, \code{\link{eqnvec}}, or
#'   \code{\link{eqnlist}} defining the equations to be set to zero.
#'   For \code{eqnlist} inputs, conserved quantities are detected and
#'   eliminated species are substituted automatically. The user then
#'   provides \code{total_*} parameters instead of the eliminated species.
#' @param parameters Character vector, the independent variables.
#' @param forcings Character vector of forcing/dummy state names.
#' @param condition Character, the condition label.
#' @param keep.root Logical, cache the root for warm-starting.
#' @param positive Logical, solve in log-space for positivity.
#' @param compile Logical, compile C++ code.
#' @param modelname Character, base filename for C++ code.
#' @param verbose Logical, print compiler output.
#' @param controlsNleqslv Named list of solver controls.
#'
#' @return A function of class \code{\link{parfn}}.
#' @seealso [Pexpl], [Pequil], [P]
#' @export
#' @import nleqslv
#' @importFrom digest digest
Pimpl <- function(trafo, parameters = NULL, forcings = NULL, condition = NULL, keep.root = TRUE,
                  positive = TRUE, compile = FALSE, modelname = NULL, verbose = FALSE,
                  controlsNleqslv = list()) {
  
  # ---- Accept eqnlist, eqnvec, or named character ----
  smatrix <- NULL
  if (inherits(trafo, "eqnlist")) {
    if (!is.null(forcings))
      trafo$rates <- replaceSymbols(forcings, rep("0", length(forcings)), trafo$rates)
    smatrix <- trafo$smatrix
    trafo <- as.eqnvec(trafo)
  } else {
    trafo <- as.eqnvec(trafo)
  }
  
  states <- names(trafo)
  
  # Replace forcing symbols by zero (for eqnvec/character input)
  if (!is.null(forcings) && is.null(smatrix))
    trafo <- replaceSymbols(forcings, rep("0", length(forcings)), trafo)
  if (!is.null(forcings)) {
    trafo <- trafo[setdiff(names(trafo), forcings)]
    states <- names(trafo)
  }
  
  # Auto-detect constant states (rhs == "0")
  const_states <- states[vapply(unclass(trafo), function(x) {
    tryCatch(identical(eval(parse(text = x)), 0), error = function(e) FALSE)
  }, logical(1))]
  if (length(const_states)) parameters <- union(parameters %||% character(0), const_states)
  
  # ---- Conserved quantities: detect, substitute ----
  cq <- .detect_and_substitute_cq(smatrix, trafo, states, parameters)
  trafo      <- cq$f
  parameters <- cq$parameters
  cq_info    <- cq$cq_info
  elim_states <- cq$elim_states
  
  dependent <- setdiff(states, parameters)
  
  # ---- Input validation ----
  if (length(dependent) == 0)
    stop("No dependent states to solve for. All states appear in 'parameters'.")
  
  # ---- Determine required input parameters ----
  if (is.null(parameters)) {
    parameters <- getSymbols(trafo, exclude = dependent)
  } else {
    parameters <- union(getSymbols(trafo, exclude = dependent), parameters)
  }
  
  parms_all <- setdiff(parameters, dependent)
  n_dep     <- length(dependent)
  
  if (is.null(modelname)) modelname <- "impl_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  # ---- Build combined compiled evaluator for f and recon ----
  # trafo[dependent] already has CQ substitutions applied, so funCpp
  # sees total_* as normal parameters and differentiates correctly.
  # If conserved quantities exist, append reconstruction expressions
  # so that a single funCpp call covers both f(x) and recon(x).
  all_exprs <- unclass(trafo[dependent])
  if (length(cq_info) > 0) {
    recon_exprs <- setNames(
      vapply(cq_info, function(ci) ci$recon_expr, character(1)),
      vapply(cq_info, function(ci) ci$elim_state, character(1))
    )
    all_exprs <- c(all_exprs, recon_exprs)
  }
  n_all <- length(all_exprs)
  all_names <- names(all_exprs)
  
  PEval <- suppressWarnings(CppODE::funCpp(
    all_exprs, variables = dependent, parameters = parms_all,
    fixed = NULL, compile = compile, modelname = modelname,
    outdir = getwd(), verbose = verbose, convenient = FALSE, deriv = TRUE
  ))
  
  jac_cols <- c(dependent, parms_all)
  
  eval_f <- function(x, p) {
    vals <- PEval$func(matrix(x[dependent], 1, dimnames = list(NULL, dependent)), p[parms_all])[, 1]
    vals[seq_len(n_dep)]
  }
  eval_J <- function(x, p) {
    J <- PEval$jac(matrix(x[dependent], 1, dimnames = list(NULL, dependent)), p[parms_all])
    matrix(J[seq_len(n_dep), , 1], n_dep, length(jac_cols), dimnames = list(dependent, jac_cols))
  }
  
  # ---- Cache ----
  cache <- new.env(parent = emptyenv())
  cache$guess <- NULL
  cache$failed_pv_hash <- NULL
  cache$failed_result  <- NULL
  
  # ---- Controls ----
  controls <- list(
    keep.root  = keep.root,
    positive   = positive,
    nstarts    = 100,
    lower      = 0,
    upper      = 100,
    debugPlot  = FALSE,
    method     = "Newton",
    global     = "dbldog",
    xscalm     = "fixed",
    xtol       = 1e-4,
    ftol       = 1e-2,
    btol       = 1e-3,
    cndtol     = 1e-12,
    maxit      = 200L,
    allowSingular = TRUE
  )
  controls <- modifyList(controls, controlsNleqslv)
  
  pimpl_keys   <- c("keep.root", "positive", "nstarts", "lower", "upper", "debugPlot")
  nleqslv_args <- c("method", "global", "xscalm")
  
  # ---- Helper: expand scalar bound ----
  expand_bounds <- function(b, dep, default_val) {
    if (is.null(names(b)) || length(b) == 1)
      return(setNames(rep(b[1], length(dep)), dep))
    out <- setNames(rep(default_val, length(dep)), dep)
    out[intersect(names(b), dep)] <- b[intersect(names(b), dep)]
    out
  }
  
  # ---- Helper: single nleqslv solve ----
  solve_once <- function(x0, pv, positive, nleqslv_top, nleqslv_ctrl) {
    method <- nleqslv_top$method
    global <- nleqslv_top$global
    xscalm <- nleqslv_top$xscalm
    
    if (positive) {
      s0 <- x0; s0[s0 <= 0] <- 1
      sol <- nleqslv::nleqslv(
        x   = log(s0),
        fn  = function(lx) { names(lx) <- dependent; eval_f(exp(lx), pv) },
        jac = function(lx) {
          names(lx) <- dependent
          xv <- exp(lx)
          eval_J(xv, pv)[, dependent, drop = FALSE] %*% diag(xv, n_dep)
        },
        method = method, global = global, xscalm = xscalm,
        control = nleqslv_ctrl)
      root <- setNames(exp(sol$x), dependent)
    } else {
      sol <- nleqslv::nleqslv(
        x   = x0,
        fn  = function(x) { names(x) <- dependent; eval_f(x, pv) },
        jac = function(x) { names(x) <- dependent; eval_J(x, pv)[, dependent, drop = FALSE] },
        method = method, global = global, xscalm = xscalm,
        control = nleqslv_ctrl)
      root <- setNames(sol$x, dependent)
    }
    res <- tryCatch(eval_f(root, pv), error = function(e) setNames(rep(Inf, n_dep), dependent))
    list(root = root, sol = sol, res = res, maxres = max(abs(res)), termcd = sol$termcd, iter = sol$iter)
  }
  
  # ---- Helper: waterfall plot ----
  waterfall_plot <- function(log_df, ftol) {
    termcd_labels <- c(
      "1" = "converged", "2" = "xtol (f may be large)", "3" = "stalled",
      "4" = "maxit exceeded", "5" = "ill-conditioned",
      "6" = "singular", "7" = "unusable Jacobian"
    )
    termcd_colors <- c(
      "converged" = "#2ca02c", "xtol (f may be large)" = "#ff7f0e",
      "stalled" = "#d62728", "maxit exceeded" = "#9467bd",
      "ill-conditioned" = "#8c564b", "singular" = "#e377c2",
      "unusable Jacobian" = "#7f7f7f"
    )
    log_df$termcd_label <- factor(termcd_labels[as.character(log_df$termcd)], levels = termcd_labels)
    P <- ggplot2::ggplot(log_df, ggplot2::aes(x = factor(rank), y = maxres, fill = termcd_label)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::geom_hline(yintercept = ftol, linetype = "dashed", color = "steelblue", linewidth = 0.6) +
      ggplot2::annotate("text", x = nrow(log_df), y = ftol, label = paste0("ftol = ", ftol),
                        hjust = 1, vjust = -0.5, color = "steelblue", size = 3) +
      ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-12),
                                  breaks = c(0, 10^(-10:4)),
                                  labels = function(x) ifelse(x == 0, "0", scales::scientific(x))) +
      ggplot2::scale_fill_manual(values = termcd_colors, name = "termination", drop = TRUE) +
      ggplot2::labs(x = "index (sorted by residual)", y = "max |f(x)|",
                    title = "Pimpl multistart diagnostics") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom",
                     axis.text.x = if (nrow(log_df) > 30) ggplot2::element_blank()
                     else ggplot2::element_text(size = 7))
    print(P)
  }
  
  # ---- Helper: build Jacobian including eliminated species ----
  build_jacobian <- function(root, pv, out, p, emptypars, fixed, dP, deriv) {
    Jfull <- eval_J(root, pv)
    dfdx  <- Jfull[, dependent, drop = FALSE]
    dfdp  <- Jfull[, parms_all, drop = FALSE]
    dxdp  <- solve(dfdx, -dfdp)
    
    input_cols <- setdiff(names(p), c(dependent, names(fixed)))
    jacobian <- matrix(0, length(out), length(input_cols),
                       dimnames = list(names(out), input_cols))
    
    for (ep in intersect(emptypars, input_cols)) jacobian[ep, ep] <- 1
    
    jacobian[rownames(dxdp), intersect(colnames(dxdp), input_cols)] <-
      dxdp[, intersect(colnames(dxdp), input_cols), drop = FALSE]
    
    # Jacobian of eliminated species from combined PEval
    # Rows (n_dep+1):n_all of the full Jacobian are the recon derivatives
    # d(elim)/d(input) = d(elim)/d(dep) * d(dep)/d(input) + d(elim)/d(parms) * I
    if (length(elim_states) > 0) {
      Jall <- PEval$jac(matrix(root[dependent], 1, dimnames = list(NULL, dependent)), pv[parms_all])
      recon_idx <- (n_dep + 1):n_all
      recon_jac <- matrix(Jall[recon_idx, , 1], length(elim_states), length(jac_cols),
                          dimnames = list(elim_states, jac_cols))
      
      drecon_ddep   <- recon_jac[, dependent, drop = FALSE]
      drecon_dparms <- recon_jac[, parms_all, drop = FALSE]
      
      dep_jac <- jacobian[dependent, input_cols, drop = FALSE]
      elim_jac <- drecon_ddep %*% dep_jac
      parm_input <- intersect(parms_all, input_cols)
      if (length(parm_input))
        elim_jac[, parm_input] <- elim_jac[, parm_input, drop = FALSE] +
        drecon_dparms[, parm_input, drop = FALSE]
      jacobian[elim_states, ] <- elim_jac
    }
    
    if (!is.null(dP))
      jacobian <- jacobian %*% submatrix(dP, rows = colnames(jacobian))
    
    nonzero <- rowSums(jacobian != 0) > 0
    jacobian[nonzero, , drop = FALSE]
  }
  
  # ---- Main p2p function ----
  p2p <- function(pars, fixed = NULL, deriv = TRUE) {
    
    p <- pars; dP <- attr(p, "deriv")
    keep.root  <- controls$keep.root
    positive   <- controls$positive
    nstarts    <- controls$nstarts
    ftol       <- controls$ftol
    debugPlot  <- controls$debugPlot
    
    nleqslv_top  <- controls[intersect(names(controls), nleqslv_args)]
    nleqslv_ctrl <- controls[setdiff(names(controls), c(pimpl_keys, nleqslv_args))]
    
    if (!is.null(fixed)) {
      is.fixed <- which(names(p) %in% names(fixed))
      if (length(is.fixed) > 0) p <- p[-is.fixed]
      p <- c(p, fixed)
    }
    
    emptypars <- names(p)[!names(p) %in% c(dependent, names(fixed))]
    if (!all(dependent %in% names(p))) p[setdiff(dependent, names(p))] <- 1
    
    pv <- p[parms_all]
    
    # ---- Failure tracking ----
    pv_hash <- digest::digest(pv, algo = "xxhash64")
    
    if (!is.null(cache$failed_pv_hash) && identical(pv_hash, cache$failed_pv_hash)) {
      warning("Multistart previously failed for these parameters. ",
              "Returning cached (inaccurate) solution.", call. = FALSE)
      fr <- cache$failed_result
      root <- fr$root
      
      # Reconstruct eliminated species from combined PEval
      elim_values <- if (length(elim_states) > 0) {
        vals <- PEval$func(matrix(root[dependent], 1, dimnames = list(NULL, dependent)), pv[parms_all])[, 1]
        setNames(vals[(n_dep + 1):n_all], elim_states)
      } else numeric(0)
      out <- c(root, elim_values, p[setdiff(names(p), c(names(root), names(elim_values)))])
      
      jac <- tryCatch(
        build_jacobian(root, pv, out, p, emptypars, fixed, dP, deriv),
        error = function(e) NULL)
      return(as.parvec(out, deriv = if (deriv && !is.null(jac)) jac else NULL))
    }
    
    # ---- Solve with multistart ----
    all_results <- list()
    best <- NULL
    
    record <- function(result, label) {
      if (is.null(result)) return()
      idx <- length(all_results) + 1L
      all_results[[idx]] <<- list(label = label, maxres = result$maxres,
                                  termcd = result$termcd, iter = result$iter)
      if (is.null(best) || result$maxres < best$maxres) best <<- result
    }
    
    if (!is.null(cache$guess)) {
      x0 <- p[dependent]
      cached_dep <- intersect(dependent, names(cache$guess))
      if (length(cached_dep)) x0[cached_dep] <- cache$guess[cached_dep]
      record(tryCatch(solve_once(x0, pv, positive, nleqslv_top, nleqslv_ctrl), error = function(e) NULL), "cache")
    }
    
    if (is.null(best) || best$maxres > ftol)
      record(tryCatch(solve_once(p[dependent], pv, positive, nleqslv_top, nleqslv_ctrl), error = function(e) NULL), "user")
    
    if ((is.null(best) || best$maxres > ftol) && nstarts > 1) {
      lo <- expand_bounds(controls$lower, dependent, 0)
      hi <- expand_bounds(controls$upper, dependent, 10)
      for (i in seq_len(nstarts)) {
        x0_rand <- setNames(runif(n_dep, lo, hi), dependent)
        record(tryCatch(solve_once(x0_rand, pv, positive, nleqslv_top, nleqslv_ctrl), error = function(e) NULL),
               paste0("random_", i))
        if (!is.null(best) && best$maxres <= ftol) break
      }
    }
    
    if (debugPlot && nstarts > 1 && length(all_results) > 0) {
      log_df <- data.frame(
        index = seq_along(all_results),
        label = vapply(all_results, `[[`, "", "label"),
        maxres = vapply(all_results, `[[`, 0, "maxres"),
        termcd = vapply(all_results, `[[`, 0L, "termcd"),
        iter = vapply(all_results, `[[`, 0L, "iter"),
        stringsAsFactors = FALSE)
      log_df <- log_df[order(log_df$maxres), ]
      log_df$rank <- seq_len(nrow(log_df))
      waterfall_plot(log_df, ftol)
    }
    
    if (is.null(best)) {
      warning("All solve attempts failed. Returning input values.", call. = FALSE)
      return(as.parvec(p, deriv = NULL))
    }
    
    root <- best$root
    sol  <- best$sol
    
    if (sol$termcd != 1 && best$maxres > ftol)
      warning("nleqslv did not converge (code ", sol$termcd, "): ", sol$message, call. = FALSE)
    
    if (best$maxres > ftol) {
      res_vec <- best$res
      res_order <- order(abs(res_vec), decreasing = TRUE)
      top_n <- min(10, length(res_vec))
      res_table <- paste0(
        sprintf("  %-25s  |f| = %s",
                names(res_vec)[res_order[1:top_n]],
                formatC(abs(res_vec[res_order[1:top_n]]), format = "e", digits = 2)),
        collapse = "\n")
      warning("Best residual norm is ",
              formatC(best$maxres, format = "e", digits = 2),
              " (ftol = ", formatC(ftol, format = "e", digits = 2),
              "). Solution may be inaccurate.\n",
              "Largest residuals (top ", top_n, "):\n", res_table, call. = FALSE)
    }
    
    # Reconstruct eliminated species from combined PEval
    elim_values <- if (length(elim_states) > 0) {
      vals <- PEval$func(matrix(root[dependent], 1, dimnames = list(NULL, dependent)), pv[parms_all])[, 1]
      setNames(vals[(n_dep + 1):n_all], elim_states)
    } else numeric(0)
    out <- c(root, elim_values, p[setdiff(names(p), c(names(root), names(elim_values)))])
    
    if (keep.root) cache$guess <- out
    if (best$maxres > ftol && nstarts > 1) {
      cache$failed_pv_hash <- pv_hash; cache$failed_result <- best
    } else {
      cache$failed_pv_hash <- NULL; cache$failed_result <- NULL
    }
    
    jac <- build_jacobian(root, pv, out, p, emptypars, fixed, dP, deriv)
    as.parvec(out, deriv = if (deriv) jac else NULL)
  }
  
  attr(p2p, "equations")  <- as.eqnvec(trafo)
  attr(p2p, "parameters") <- parameters
  attr(p2p, "modelname")  <- modelname
  attr(p2p, "compileInfo") <- collectCompileInfo(PEval$func, PEval$jac)
  parfn(p2p, parameters, condition)
}


#' Parameter transformation (steady states via pre-equilibration)
#'
#' @description
#' Finds steady states by integrating the ODE to equilibrium.
#' All states are integrated directly without conserved-quantity
#' elimination. Initial values for each state are taken from the
#' input parameter vector.
#'
#' @param trafo An \code{\link{eqnlist}}, \code{\link{eqnvec}},
#'   or named character vector.
#' @param parameters Character vector of independent parameters.
#' @param forcings Character vector of forcing names replaced by zero.
#' @param condition Character label for the condition.
#' @param attach.input Logical. Append pass-through parameters to output.
#' @param keep.root Logical. Cache steady state for warm-starting.
#' @param controlsODE Named list of ODE solver options, accessible via
#'   \code{\link{controls}()}.
#' @param compile Logical. Compile the C++ ODE model.
#' @param modelname Character, base name for C++ code files.
#' @param verbose Logical. Print compiler output.
#' @param ... Additional arguments passed to \code{\link[CppODE]{CppODE}}.
#'
#' @return A function of class \code{\link{parfn}}.
#' @seealso [Pexpl], [Pimpl], [P]
#' @import CppODE
#' @importFrom digest digest
#' @export
Pequil <- function(trafo, parameters = NULL, forcings = NULL, condition = NULL, attach.input = TRUE, 
                   start.time = -1e7, end.time = 0, keep.root = TRUE, controlsODE = list(),
                   compile = FALSE, modelname = NULL, verbose = FALSE, ...) {
  
  # ---- Accept eqnlist, eqnvec, or named character ----
  if (inherits(trafo, "eqnlist")) {
    if (!is.null(forcings))
      trafo$rates <- replaceSymbols(forcings, "0", trafo$rates)
    f <- as.eqnvec(trafo)
  } else if (inherits(trafo, "eqnvec") || is.character(trafo)) {
    f <- as.eqnvec(replaceSymbols(forcings, "0", trafo))
  } else {
    stop("'trafo' must be an eqnlist, eqnvec or character vector")
  }
  
  states <- names(f)
  
  # Auto-detect constant states (rhs == "0") -> treated as parameters
  const_states <- states[vapply(unclass(f), function(x) {
    tryCatch(identical(eval(parse(text = x)), 0), error = function(e) FALSE)
  }, logical(1))]
  if (length(const_states)) parameters <- union(parameters %||% character(0), const_states)
  
  # States in parameters are NOT integrated
  dependent <- setdiff(states, parameters)
  if (length(dependent) == 0) stop("No dependent states left to equilibrate.")
  n_dep <- length(dependent)
  f_red <- f[dependent]
  
  if (is.null(parameters)) {
    parameters <- getSymbols(f_red, exclude = dependent)
  } else {
    parameters <- union(getSymbols(f_red, exclude = dependent), parameters)
  }
  parms_all <- setdiff(parameters, dependent)
  
  if (is.null(modelname)) modelname <- "equil_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  # ---- Build CppODE models ----
  dotArgs <- list(...)
  .args <- c(list(rhs = unclass(f_red), rootfunc = "equilibrate", deriv2 = FALSE,
                  compile = compile, outdir = getwd(), useDenseOutput = FALSE, verbose = verbose),
             dotArgs)
  model   <- do.call(CppODE::CppODE, c(.args, list(deriv = FALSE, modelname = modelname)))
  model_s <- do.call(CppODE::CppODE, c(.args, list(deriv = TRUE,  modelname = paste0(modelname, "_s"),
                                                   fixed = names(f))))
  dims    <- attr(model_s, "dim_names")
  all_sens <- dims$sens
  
  # ---- ODE controls ----
  ode_defaults <- list(abstol = 1e-6, reltol = 1e-6, maxsteps = 1e6L,
                       maxprogress = 100L, hini = 0, roottol = 1e-6, maxroot = 1L)
  ode_ctrl <- modifyList(ode_defaults, controlsODE)
  
  # ---- Cache ----
  cache <- new.env(parent = emptyenv())
  cache$yini <- NULL; cache$sensini <- NULL
  cache$last_hash <- NULL; cache$last_result <- NULL
  
  # Default sensitivities
  default_sens <- matrix(0, n_dep, length(all_sens), dimnames = list(dependent, all_sens))
  diag_vars <- intersect(dependent, all_sens)
  if (length(diag_vars)) default_sens[cbind(diag_vars, diag_vars)] <- 1
  
  # ---- Controls ----
  controls <- c(list(keep.root = keep.root, attach.input = attach.input,
                     start.time = start.time, end.time = end.time), ode_ctrl)
  
  # ---- p2p ----
  p2p <- function(pars, fixed = NULL, deriv = TRUE) {
    p <- pars; dP <- attr(p, "deriv")
    keep.root    <- controls$keep.root
    attach.input <- controls$attach.input
    
    if (!is.null(fixed)) {
      is.fixed <- names(p) %in% names(fixed)
      if (any(is.fixed)) p <- p[!is.fixed]
      p <- c(p, fixed)
    }
    emptypars <- names(p)[!names(p) %in% c(dependent, names(fixed))]
    missing_dep <- setdiff(dependent, names(p))
    if (length(missing_dep)) p[missing_dep] <- 1
    
    # Fast path: skip if parameters unchanged
    if (keep.root) {
      pv_hash <- digest::digest(list(p[dependent], p[parms_all], fixed, deriv), algo = "xxhash64")
      if (!is.null(cache$last_hash) && identical(pv_hash, cache$last_hash) &&
          !is.null(cache$last_result))
        return(cache$last_result)
    }
    
    # Warm-start
    if (keep.root && !is.null(cache$yini)) p[dependent] <- cache$yini
    
    # Active sensitivities
    fixed_char <- if (!is.null(fixed)) intersect(names(fixed), all_sens) else NULL
    if (length(fixed_char)) {
      active_sens <- all_sens[-match(fixed_char, all_sens)]
    } else { fixed_char <- NULL; active_sens <- all_sens }
    n_active <- length(active_sens)
    
    s1ini <- NULL
    if (deriv && keep.root && !is.null(cache$sensini))
      s1ini <- cache$sensini[, active_sens, drop = FALSE]
    else if (deriv)
      s1ini <- default_sens[, active_sens, drop = FALSE]
    
    # Integrate
    res <- tryCatch(
      withCallingHandlers(
        CppODE::solveODE(
          if (deriv) model_s else model,
          times = c(controls$start.time, controls$end.time), parms = c(p[dependent], p[parms_all]),
          sens1ini = s1ini, fixed = fixed_char,
          roottol = controls$roottol, abstol = controls$abstol, reltol = controls$reltol,
          maxsteps = as.integer(controls$maxsteps), maxprogress = as.integer(controls$maxprogress),
          hini = controls$hini, maxroot = as.integer(controls$maxroot)),
        warning = function(w) { warning(w$message, call. = FALSE); invokeRestart("muffleWarning") }),
      error = function(e) { warning("ODE integration failed: ", e$message, call. = FALSE); NULL })
    
    if (is.null(res)) {
      out <- if (attach.input) c(p[dependent], p[setdiff(names(p), dependent)]) else p[dependent]
      return(as.parvec(out, deriv = NULL))
    }
    
    last <- length(res$time)
    if (res$time[last] == end.time)
      warning("Steady state not reached within integration time.", call. = FALSE)
    
    root <- round(res$variable[, last], floor(-log10(controls$roottol))+1L)
    names(root) <- dependent
    
    if (attach.input) {
      out <- c(root, p[setdiff(names(p), dependent)])
    } else {
      out <- root
    }
    
    # Update warm-start cache
    if (keep.root) {
      cache$yini <- root
      if (!is.null(res$sens1)) {
        full_new <- default_sens
        full_new[, active_sens] <- res$sens1[, , last]
        cache$sensini <- full_new
      } else cache$sensini <- NULL
    }
    
    if (deriv && !is.null(res$sens1)) {
      sens_final <- round(matrix(res$sens1[, , last], nrow = n_dep, ncol = n_active,
                                 dimnames = list(dependent, active_sens)),
                          floor(-log10(controls$roottol))+1L)
      
      input_cols <- setdiff(names(p), c(dependent, names(fixed)))
      jacobian <- matrix(0, length(out), length(input_cols),
                         dimnames = list(names(out), input_cols))
      
      # Pass-through identity
      if (attach.input) {
        diag_idx <- intersect(emptypars, input_cols)
        if (length(diag_idx)) jacobian[cbind(diag_idx, diag_idx)] <- 1
      }
      
      # Dependent state sensitivities (direct from ODE)
      sr <- intersect(dependent, rownames(sens_final))
      sc <- intersect(input_cols, colnames(sens_final))
      if (length(sr) && length(sc))
        jacobian[sr, sc] <- sens_final[sr, sc, drop = FALSE]
      
      if (!is.null(dP)) jacobian <- jacobian %*% submatrix(dP, rows = input_cols)
      nonzero <- rowSums(jacobian != 0) > 0
      result <- as.parvec(out, deriv = jacobian[nonzero, , drop = FALSE])
    } else {
      result <- as.parvec(out, deriv = NULL)
    }
    
    if (keep.root) { cache$last_hash <- pv_hash; cache$last_result <- result }
    result
  }
  
  attr(p2p, "equations")  <- as.eqnvec(f)
  attr(p2p, "parameters") <- parameters
  attr(p2p, "modelname")  <- modelname
  attr(p2p, "compileInfo") <- collectCompileInfo(model, model_s)
  parfn(p2p, parameters, condition)
}


#' Construct parameter transformations
#'
#' Helper functions to construct and modify symbolic parameter transformations
#' used by prediction functions such as [P()] and [Xs()].
#'
#' The functions [define()], [insert()] and [branch()] operate exclusively on
#' the symbolic level. They are used to build transformation objects that
#' describe how *outer parameters* are expressed in terms of *inner parameters*
#' or constants.
#'
#' No model evaluation, sensitivity calculation or parameter checking is
#' performed by these functions. The resulting transformations are interpreted
#' later when prediction or objective functions are constructed.
#'
#' \describe{
#'   \item{define}{
#'     Reset or redefine a transformation rule by explicitly specifying a new
#'     right-hand side.
#'   }
#'   \item{insert}{
#'     Insert symbolic substitutions into existing transformation rules without
#'     resetting them.
#'   }
#'   \item{branch}{
#'     Duplicate a transformation for multiple conditions and optionally apply
#'     condition-specific substitutions.
#'   }
#' }
#'
#' When transformations are branched, a condition table is stored as metadata
#' (attribute \code{"tree"}) and may be used to restrict subsequent calls to
#' [define()] or [insert()] to specific conditions.
#'
#' @param trafo
#'   A named character vector, an object of class \code{eqnvec}, or a list
#'   thereof representing parameter transformations.
#'
#' @param expr
#'   Character string of the form \code{"lhs ~ rhs"} defining a symbolic
#'   transformation or substitution.
#'
#' @param table
#'   Optional data frame specifying condition-specific substitutions. Rownames
#'   identify conditions; columns correspond to parameter names.
#'
#' @param conditions
#'   Character vector of condition names. If supplied, overrides
#'   \code{rownames(table)}.
#'
#' @param apply
#'   Character string specifying whether and how entries of \code{table} are
#'   applied when branching:
#'   \describe{
#'     \item{"nothing"}{Only duplicate the transformation (default).}
#'     \item{"insert"}{Apply entries via [insert()].}
#'     \item{"define"}{Apply entries via [define()].}
#'   }
#'
#' @param conditionMatch
#'   Optional character string (regular expression). If provided, the operation
#'   is applied only to conditions whose names match this expression.
#'
#' @param ...
#'   Named values used to substitute symbols occurring in \code{expr}.
#'
#' @return
#' An object of the same type as \code{trafo}, possibly expanded to a list if
#' branching has been applied.
#'
#' @export
#' @example inst/examples/define.R
define <- function(trafo, expr, ..., conditionMatch = NULL) {
  
  if (missing(trafo)) trafo <- NULL
  lookuptable <- attr(trafo, "tree")
  
  
  if (is.list(trafo) & is.null(names(trafo)))
    stop("If trafo is a list, elements must be named.")
  
  if (is.list(trafo) & !all(names(trafo) %in% rownames(lookuptable)))
    stop("If trafo is a list and contains a lookuptable (is branched from a tree), the list names must be contained in the rownames of the tree.")
  
  if (!is.list(trafo)) {
    mytrafo <- list(trafo)
  } else {
    mytrafo <- trafo
  }
  
  dots <- substitute(alist(...))
  out <- lapply(1:length(mytrafo), function(i) {
    
    .currentTrafo <- mytrafo[[i]]
    .currentSymbols <- NULL
    if (!is.null(.currentTrafo))
      .currentSymbols <- getSymbols(.currentTrafo)
    
    if (is.list(trafo)) {
      mytable <- lookuptable[names(mytrafo)[i], , drop = FALSE]
    } else {
      mytable <- lookuptable[1, , drop = FALSE]
    }
    
    if ((!is.null(conditionMatch)))
      if((!str_detect(rownames(mytable), conditionMatch))) 
        return(mytrafo[[i]])
    
    with(mytable, {
      args <- c(list(expr = expr, trafo = mytrafo[[i]], reset = TRUE), eval(dots))
      do.call(repar, args)
    })
    
    
  })
  names(out) <- names(mytrafo)
  if (!is.list(trafo)) out <- out[[1]]
  attr(out, "tree") <- lookuptable
  
  return(out)
  
}


#' @export
#' @rdname define
insert <- function(trafo, expr, ..., conditionMatch = NULL) {
  
  
  if (missing(trafo)) trafo <- NULL
  lookuptable <- attr(trafo, "tree")
  
  
  if (is.list(trafo) & is.null(names(trafo)))
    stop("If trafo is a list, elements must be named.")
  
  if (is.list(trafo) & !all(names(trafo) %in% rownames(lookuptable)))
    stop("If trafo is a list and contains a lookuptable (is branched from a tree), the list names must be contained in the rownames of the tree.")
  
  if (!is.list(trafo)) {
    mytrafo <- list(trafo)
  } else {
    mytrafo <- trafo
  }
  
  dots <- substitute(alist(...))
  out <- lapply(1:length(mytrafo), function(i) {
    
    .currentTrafo <- mytrafo[[i]]
    .currentSymbols <- NULL
    if (!is.null(.currentTrafo))
      .currentSymbols <- getSymbols(.currentTrafo)
    
    if (is.list(trafo)) {
      mytable <- lookuptable[names(mytrafo)[i], , drop = FALSE]
    } else {
      mytable <- lookuptable[1, , drop = FALSE]
    }
    
    if ((!is.null(conditionMatch)))
      if((!str_detect(rownames(mytable), conditionMatch))) 
        return(.currentTrafo)
    
    
    
    with(mytable, {
      .fun <- function() {
        # subset conditions by logicals expressions supplied by the dots
        dots_eval <- eval(dots)                                            # convert from substituted to language
        if (length(dots_eval) == 0) 
          return(do.call(repar, list(expr = expr, trafo = .currentTrafo)))
        
        dots_eval_eval <- lapply(dots_eval, function(i) eval.parent(i, 3)) # evaluate the language in the "mytable" frame. parent1: lapply, parent2: .fun, parent3: with
        which_logical <- vapply(dots_eval_eval, function(i) {is.logical(i)}, FUN.VALUE = vector("logical", 1)) # which of the dots are logical
        logical_dots <- dots_eval[which_logical]                           # subset to matching conditions
        matching <- do.call(c, logical_dots)
        if(!is.null(matching)) { # null means no logical dots were supplied
          if (any(!matching)) {  
            return(.currentTrafo) }}
        
        args <- c(list(expr = expr, trafo = .currentTrafo), dots_eval_eval[!which_logical]) # feed the rest of the eval'd dots into repar
        return(do.call(repar, args))
      }
      .fun()
    })
    
    
  })
  names(out) <- names(mytrafo)
  if (!is.list(trafo)) out <- out[[1]]
  attr(out, "tree") <- lookuptable
  
  return(out)
  
}


#' @export
#' @rdname define
branch <- function(
    trafo,
    table = NULL,
    conditions = rownames(table),
    apply = c("nothing", "insert", "define")) {
  
  
  apply <- match.arg(apply)
  
  # --- trivial case ----------------------------------------------------------
  if (is.null(table) && is.null(conditions))
    return(trafo)
  
  # --- normalize inputs ------------------------------------------------------
  if (is.null(conditions))
    conditions <- paste0("C", seq_len(nrow(table)))
  
  if (is.null(table))
    table <- data.frame(condition = conditions, row.names = conditions)
  
  rownames(table) <- conditions
  
  # --- branch trafo -----------------------------------------------------------
  out <- setNames(lapply(conditions, function(x) trafo), conditions)
  attr(out, "tree") <- table
  
  # --- optional application of table -----------------------------------------
  if (apply == "nothing")
    return(out)
  
  for (cn in conditions) {
    row <- table[cn, , drop = FALSE]
    row <- row[, !colnames(row) %in% c("condition", "conditions"), drop = FALSE]
    
    for (par in colnames(row)) {
      val <- row[[par]]
      if (is.na(val)) next
      
      expr <- paste0(par, " ~ ", val)
      
      # Einzelelement mit einzeiliger Dummy-Tree versehen
      single_trafo <- out[[cn]]
      attr(single_trafo, "tree") <- row
      
      if (apply == "insert") {
        out[[cn]] <- insert(single_trafo, expr)
      }
      if (apply == "define") {
        out[[cn]] <- define(single_trafo, expr)
      }
      attr(out[[cn]], "tree") <- NULL
    }
  }
  
  out
}



#' Reparameterization
#' 
#' @param expr character of the form `"lhs ~ rhs"` where `rhs`
#' reparameterizes `lhs`. Both `lhs` and `rhs`
#' can contain a number of symbols whose values need to be passed by the `...` argument.
#' @param trafo character or equation vector or list thereof. The object where the replacement takes place in
#' @param ... pass symbols as named arguments
#' @param reset logical. If true, the trafo element corresponding to lhs is reset according to rhs. 
#' If false, lhs wherever it occurs in the rhs of trafo is replaced by rhs of the formula.
#' @return an equation vector with the reparameterization.
#' @details Left and right-hand side of `expr` are searched for symbols. If separated by
#' "_", symbols are recognized as such, e.g. in `Delta_x` where the symbols are 
#' "Delta" and "x". Each symbol for which values (character or numbers) are passed by the
#' `...` argument is replaced.
#' @export
#' @importFrom stats as.formula
#' @examples
#' innerpars <- letters[1:3]
#' constraints <- c(a = "b + c")
#' mycondition <- "cond1"
#' 
#' trafo <- repar("x ~ x", x = innerpars)
#' trafo <- repar("x ~ y", trafo, x = names(constraints), y = constraints)
#' trafo <- repar("x ~ exp(x)", trafo, x = innerpars)
#' trafo <- repar("x ~ x + Delta_x_condition", trafo, x = innerpars, condition = mycondition)
repar <- function(expr, trafo = NULL, ..., reset = FALSE) {
  
  if (inherits(expr, "formula")) expr <- deparse(expr)
  
  parsed.expr <- as.character(stats::as.formula(gsub("_", ":", expr, fixed = TRUE)))
  lhs <- parsed.expr[2]
  lhs.symbols <- getSymbols(lhs)
  rhs <- parsed.expr[3]
  rhs.symbols <- getSymbols(rhs)
  
  # Make sure that arguments are characters
  args <- lapply(list(...), as.character)
  
  replacements <- as.data.frame(args, stringsAsFactors = FALSE)
  
  lhs <- sapply(1:nrow(replacements), function(i) {
    out <- replaceSymbols(colnames(replacements), replacements[i, ], lhs)
    gsub(":", "_", out, fixed = TRUE)
  })
  
  rhs <- sapply(1:nrow(replacements), function(i) {
    out <- replaceSymbols(colnames(replacements), replacements[i, ], rhs)
    gsub(":", "_", out, fixed = TRUE)
  })
  
  if (is.null(trafo)) {
    trafo <- as.eqnvec(structure(lhs, names = lhs))
  } else if (is.list(trafo) & !reset) {
    trafo <- lapply(trafo, function(t) replaceSymbols(lhs, rhs, t))
  } else if (is.character(trafo) & !reset) {
    trafo <- replaceSymbols(lhs, rhs, trafo)
  } else if (is.list(trafo) & reset) {
    trafo <- lapply(trafo, function(t) {t[lhs] <- rhs; return(t)})
  } else if (is.character(trafo) & reset) {
    trafo[lhs] <- rhs
  }
  
  return(trafo)
  
  
}

paste_ <- function(...) paste(..., sep = "_")