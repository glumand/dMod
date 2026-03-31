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
                  compile = FALSE, modelname = NULL, verbose = FALSE) {
  
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
      convenient = FALSE
    )
  )
  
  fun <- PEval$func
  jac <- PEval$jac
  
  
  # Define returned parameter transformation function
  p2p <- function(pars, fixed = NULL, deriv = TRUE) {
    
    # Prepare pars
    p <- c(pars, fixed)
    
    # Evaluate inner parameters
    pinnerVal <- fun(NULL, p, attach.input = attach.input, fixed = names(fixed))[,]
    
    if (any(is.nan(pinnerVal))) {
      stop(
        paste0(
          "The following inner parameter(s) evaluate to NaN:\n\t",
          paste0(names(pinnerVal)[is.nan(pinnerVal)], collapse = "\n\t"),
          ".\nLikely cause: division by zero or missing inputs."
        )
      )
    }
    
    # Apply chain rule for derivatives
    Jac <- NULL
    if (deriv && !is.null(jac)) {
      Jac <- as.matrix(jac(NULL, p, attach.input = attach.input, fixed = names(fixed))[,,1])
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        Jac <- Jac %*% dP[colnames(Jac), , drop = FALSE]
        dimnames(Jac) <- list(names(pinnerVal), colnames(dP))
      }
    }
    
    # Assemble result
    pinner <- as.parvec(pinnerVal, deriv = if (deriv) Jac[rowSums(Jac != 0) > 0, , drop = FALSE] else FALSE)
    
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
  
  parfn(p2p, parameters, condition)
}


#' Parameter transformation (implicit)
#' 
#' @param trafo Named character vector defining the equations to be set to zero. 
#' Names correspond to dependent variables.
#' @param parameters Character vector, the independent variables.  
#' @param condition character, the condition for which the transformation is generated
#' @param compile Logical, compile the C++ code containing the function and jacobian (see \link{funCpp})
#' @param keep.root logical, The root of the last evaluation of the parameter transformation 
#' function is saved as guess for the next evaluation.
#' @param positive logical, if \code{TRUE} the implicit variables are solved in log-space to 
#' enforce strict positivity of the solution. 
#' @param modelname Character, sets a fixed filename for the C++ file.
#' @param verbose Print compiler output to R command line.
#' @return a function \code{p2p(p, fixed = NULL, deriv = TRUE)} representing the parameter 
#' transformation. Here, \code{p} is a named numeric vector with the values of the outer parameters,
#' \code{fixed} is a named numeric vector with values of the outer parameters being considered
#' as fixed (no derivatives returned) and \code{deriv} is a logical determining whether the Jacobian
#' of the parameter transformation is returned as attribute "deriv".
#' @details Usually, the equations contain the dependent variables, the independent variables and 
#' other parameters. The argument \code{p} of \code{p2p} must provide values for the independent
#' variables and the parameters but ALSO FOR THE DEPENDENT VARIABLES. Those serve as initial guess
#' for the dependent variables. The dependent variables are then numerically computed by 
#' \link[nleqslv]{nleqslv}. The Jacobian of the solution with respect to dependent variables
#' and parameters is computed by the implicit function theorem. The function \code{p2p} returns
#' all parameters as they are with corresponding 1-entries in the Jacobian.
#' 
#' When \code{positive = TRUE}, the solver works in log-space internally: the substitution 
#' \code{x = exp(logx)} is applied, so \code{f(exp(logx), p) = 0} is solved for \code{logx}. 
#' This guarantees \code{x > 0} by construction. The analytic Jacobian for the nonlinear solver 
#' is transformed via the chain rule: \code{df/d(logx) = df/dx * diag(x)}.
#' 
#' @seealso \link{Pexpl} for explicit parameter transformations
#' @examples
#' ########################################################################
#' ## Example 1: Steady-state trafo
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pimpl(f, "A")
#' 
#' p.outerValues <- c(k1 = 1, k2 = 0.1, A = 10, B = 1)
#' P.steadyState(p.outerValues)
#' 
#' ########################################################################
#' ## Example 2: Steady-state trafo combined with log-transform
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pimpl(f, "A")
#' 
#' logtrafo <- c(k1 = "exp(logk1)", k2 = "exp(logk2)", A = "exp(logA)", B = "exp(logB)")
#' P.log <- P(logtrafo)
#' 
#' p.outerValue <- c(logk1 = 1, logk2 = -1, logA = 0, logB = 0)
#' (P.log)(p.outerValue)
#' (P.steadyState * P.log)(p.outerValue)
#' 
#' ########################################################################
#' ## Example 3: Steady-states with conserved quantitites
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B", B = "k1*A - k2*B")
#' replacement <- c(B = "A + B - total")
#' f[names(replacement)] <- replacement
#' 
#' pSS <- Pimpl(f, "total")
#' pSS(c(k1 = 1, k2 = 2, A = 5, B = 5, total = 3))
#' @export
#' @import nleqslv
Pimpl <- function(trafo, parameters = NULL, condition = NULL, keep.root = TRUE, 
                  positive = TRUE, compile = FALSE, modelname = NULL, verbose = FALSE) {
  
  states    <- names(trafo)
  dependent <- setdiff(states, parameters)
  
  # Determine required input parameters (analogous to Pexpl):
  # All symbols in the equations that are not dependent variables.
  # If the user passed `parameters`, those not appearing as states are 
  # added as identity pass-through.
  if (is.null(parameters)) {
    parameters <- getSymbols(trafo, exclude = dependent)
  } else {
    parameters <- union(getSymbols(trafo, exclude = dependent), parameters)
  }
  
  # For funCpp: the dependent variables are "variables", everything else is "parameters"
  parms_all <- setdiff(parameters, dependent)
  n_dep     <- length(dependent)
  
  if (is.null(modelname)) modelname <- "impl_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  # Single funCpp call: Jacobian columns = c(dependent, parms_all)
  PEval <- suppressWarnings(CppODE::funCpp(
    unclass(trafo[dependent]), variables = dependent, parameters = parms_all,
    fixed = NULL, compile = compile, modelname = modelname,
    outdir = getwd(), verbose = verbose, convenient = FALSE, deriv = TRUE
  ))
  
  jac_cols <- c(dependent, parms_all)
  
  eval_f <- function(x, p) {
    PEval$func(matrix(x[dependent], 1, dimnames = list(NULL, dependent)), p[parms_all])[, 1]
  }
  eval_J <- function(x, p) {
    J <- PEval$jac(matrix(x[dependent], 1, dimnames = list(NULL, dependent)), p[parms_all])
    matrix(J[, , 1], n_dep, length(jac_cols), dimnames = list(dependent, jac_cols))
  }
  
  guess    <- NULL
  controls <- list(keep.root = keep.root, positive = positive)
  
  p2p <- function(pars, fixed = NULL, deriv = TRUE) {
    
    p <- pars; dP <- attr(p, "deriv")
    keep.root <- controls$keep.root; positive <- controls$positive
    
    if (!is.null(fixed)) {
      is.fixed <- which(names(p) %in% names(fixed))
      if (length(is.fixed) > 0) p <- p[-is.fixed]
      p <- c(p, fixed)
    }
    
    emptypars <- names(p)[!names(p) %in% c(dependent, names(fixed))]
    if (!is.null(guess)) p[intersect(dependent, names(guess))] <- guess[intersect(dependent, names(guess))]
    if (!all(dependent %in% names(p))) p[setdiff(dependent, names(p))] <- 1
    
    pv <- p[parms_all]
    
    if (positive) {
      s0 <- p[dependent]; s0[s0 <= 0] <- 1
      sol <- nleqslv::nleqslv(
        x   = log(s0),
        fn  = function(lx) { names(lx) <- dependent; eval_f(exp(lx), pv) },
        jac = function(lx) { names(lx) <- dependent; x <- exp(lx); eval_J(x, pv)[, dependent] %*% diag(x, n_dep) },
        method = "Newton", control = list(maxit = 100L, ftol = 1e-12))
      root <- setNames(exp(sol$x), dependent)
    } else {
      sol <- nleqslv::nleqslv(
        x   = p[dependent],
        fn  = function(x) { names(x) <- dependent; eval_f(x, pv) },
        jac = function(x) { names(x) <- dependent; eval_J(x, pv)[, dependent] },
        method = "Newton", control = list(maxit = 100L, ftol = 1e-12))
      root <- setNames(sol$x, dependent)
    }
    if (sol$termcd > 2) warning("nleqslv did not converge (code ", sol$termcd, "): ", sol$message)
    
    out <- c(root, p[setdiff(names(p), names(root))])
    if (keep.root) guess <<- out
    
    # Implicit function theorem: dx/dp = -dfdx^{-1} dfdp
    Jfull <- eval_J(root, pv)
    dxdp  <- solve(Jfull[, dependent], -Jfull[, parms_all])
    
    # Jacobian columns: only actual input parameters (exclude dependent and fixed)
    input_cols <- setdiff(names(p), c(dependent, names(fixed)))
    jacobian <- matrix(0, length(out), length(input_cols), dimnames = list(names(out), input_cols))
    for (ep in intersect(emptypars, input_cols)) jacobian[ep, ep] <- 1
    jacobian[rownames(dxdp), intersect(colnames(dxdp), input_cols)] <- dxdp[, intersect(colnames(dxdp), input_cols)]
    if (!is.null(dP)) jacobian <- jacobian %*% submatrix(dP, rows = colnames(jacobian))
    
    as.parvec(out, deriv = if (deriv) jacobian else NULL)
  }
  
  attr(p2p, "equations")  <- as.eqnvec(trafo)
  attr(p2p, "parameters") <- parameters
  attr(p2p, "modelname")  <- modelname
  parfn(p2p, parameters, condition)
}


#' Parameter transformation (steady states via pre-equilibration)
#'
#' Constructs a parameter transformation that finds steady states by
#' integrating the ODE system until the norm of the right-hand side falls below
#' a tolerance, using \code{CppODE} with \code{rootfunc = "equilibrate"}.
#'
#' @description
#' When \code{trafo} is an \code{\link{eqnlist}}, conserved quantities are
#' automatically detected. For each conserved quantity involving \code{n}
#' species, exactly \code{n-1} are placed into \code{parameters} (either
#' user-supplied or auto-filled). States whose right-hand side evaluates to
#' zero are also automatically treated as parameters.
#'
#' The reduced ODE (only dependent states) is integrated to steady state.
#' Sensitivities are obtained from the ODE sensitivity integration and
#' composed via the chain rule when upstream derivatives are present.
#'
#' When \code{trafo} is an \code{\link{eqnvec}} or named character vector,
#' equilibration proceeds without conserved-quantity checks and a warning
#' is issued.
#'
#' @param trafo An \code{\link{eqnlist}} (recommended) or \code{\link{eqnvec}} /
#'   named character vector. For \code{eqnlist} inputs, conserved quantities
#'   and constant states are detected automatically. For \code{eqnvec} inputs,
#'   a warning is issued since independence of equations cannot be verified.
#' @param parameters Character vector of independent parameters. Species listed
#'   here are excluded from integration and treated as fixed ODE parameters.
#'   For systems with conserved quantities, exactly \code{n-1} of the \code{n}
#'   involved species must appear here (auto-filled if missing).
#' @param forcings Character vector of forcing names. These symbols are replaced
#'   by zero in the reaction rates before equilibration. Only used when
#'   \code{trafo} is an \code{\link{eqnlist}}.
#' @param condition Character label for which the transformation is generated.
#' @param attach.input Logical. If \code{TRUE} (default), pass-through
#'   parameters are appended to the output with identity derivatives.
#'   If \code{FALSE}, only the equilibrated dependent states are returned.
#' @param keep.root Logical. If \code{TRUE} (default), the steady state and
#'   sensitivities are cached for warm-starting subsequent evaluations.
#'   Improves performance during iterative fitting.
#' @param controlsODE Named list of ODE solver options passed to
#'   \code{\link[CppODE]{solveODE}}. Supported entries:
#'   \code{abstol} (default 1e-6), \code{reltol} (default 1e-6),
#'   \code{maxsteps} (default 1e7), \code{maxprogress} (default 100),
#'   \code{hini} (default 0), \code{roottol} (default 1e-6),
#'   \code{maxroot} (default 1).
#' @param compile Logical. If \code{TRUE}, compile the C++ ODE model for
#'   faster evaluation.
#' @param modelname Character, base name for generated C++ code files.
#' @param verbose Logical. If \code{TRUE}, print compiler output.
#'
#' @return A function of class \code{\link{parfn}}.
#'
#' @seealso [Pexpl] for explicit transformations, [Pimpl] for implicit
#'   (root-finding) transformations, [P] for automatic method selection.
#'
#' @import CppODE
#' @export
Pequil <- function(trafo, parameters = NULL, forcings = NULL, condition = NULL,
                   attach.input = TRUE, keep.root = TRUE, controlsODE = list(),
                   compile = FALSE, modelname = NULL, verbose = FALSE) {
  
  if (inherits(trafo, "eqnlist")) {
    # Exclude forcings
    if (!is.null(forcings))
      trafo$rates <- replaceSymbols(forcings, rep("0", length(forcings)), trafo$rates)
    smatrix <- trafo$smatrix
    f <- as.eqnvec(trafo)
  } else if (inherits(trafo, "eqnvec") || is.character(trafo)) {
    warning("'trafo' is not an eqnlist. Conserved quantities cannot be checked. ",
            "Consider using an eqnlist for automatic independence validation.",
            call. = FALSE)
    smatrix <- NULL
    f <- as.eqnvec(trafo)
  } else {
    stop("'trafo' must be an eqnlist or eqnvec")
  }
  
  states <- names(f)
  
  # Auto-detect constant states (rhs == "0") and add to parameters
  const_states <- states[vapply(unclass(f), function(x) {
    tryCatch(identical(eval(parse(text = x)), 0), error = function(e) FALSE)
  }, logical(1))]
  if (length(const_states)) parameters <- union(parameters %||% character(0), const_states)
  
  # Auto-fill parameters from conserved quantities (only if smatrix available)
  if (!is.null(smatrix)) {
    cq <- conservedQuantities(smatrix)
    if (!is.null(cq) && nrow(cq) > 0) {
      if (is.null(parameters)) parameters <- character(0)
      for (i in seq_len(nrow(cq))) {
        cq_species <- intersect(getSymbols(as.character(cq[i, 1])), states)
        n_need <- length(cq_species) - 1
        already <- intersect(cq_species, parameters)
        if (length(already) < n_need) {
          fill <- setdiff(cq_species, already)[seq_len(n_need - length(already))]
          parameters <- union(parameters, fill)
        }
      }
    }
  }
  
  # States in parameters are NOT integrated — they become ODE parameters
  dependent <- setdiff(states, parameters)
  if (length(dependent) == 0) stop("No dependent states left to equilibrate.")
  n_dep <- length(dependent)
  f_red <- f[dependent]
  
  if (is.null(parameters)) { parameters <- getSymbols(f_red, exclude = dependent)
  } else { parameters <- union(getSymbols(f_red, exclude = dependent), parameters) }
  parms_all <- setdiff(parameters, dependent)
  
  if (is.null(modelname)) modelname <- "equil_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  # Build two models: with and without sensitivities
  .args <- list(rhs = unclass(f_red), rootfunc = "equilibrate", deriv2 = FALSE, 
                compile = compile, outdir = getwd(), useDenseOutput = FALSE, verbose = verbose)
  model   <- do.call(CppODE::CppODE, c(.args, list(deriv = FALSE, modelname = modelname)))
  model_s <- do.call(CppODE::CppODE, c(.args, list(deriv = TRUE,  modelname = paste0(modelname, "_s"))))
  dims    <- attr(model_s, "dim_names")
  all_sens <- dims$sens
  
  # Default ODE controls
  ode_ctrl <- modifyList(list(abstol = 1e-6, reltol = 1e-6, maxsteps = 1e6L, 
                              maxprogress = 100L, hini = 0, roottol = 1e-5, maxroot = 1L), controlsODE)
  
  # Cache for warm-starting
  cache <- new.env(parent = emptyenv())
  cache$yini    <- NULL
  cache$sensini <- NULL
  
  # Default sensitivities [n_dep x all_sens]: dx_i/dx_j = delta_ij, rest 0
  default_sens <- matrix(0, n_dep, length(all_sens), dimnames = list(dependent, all_sens))
  diag_vars <- intersect(dependent, all_sens)
  if (length(diag_vars)) default_sens[cbind(diag_vars, diag_vars)] <- 1
  
  controls <- list(keep.root = keep.root, attach.input = attach.input)
  
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
    
    # Warm-start from cache
    if (keep.root && !is.null(cache$yini)) p[dependent] <- cache$yini
    
    # Compute active_sens exactly as solveODE does
    fixed_char <- if (!is.null(fixed)) intersect(names(fixed), all_sens) else NULL
    if (length(fixed_char)) {
      active_sens <- all_sens[-match(fixed_char, all_sens)]
    } else { fixed_char <- NULL; active_sens <- all_sens }
    n_active <- length(active_sens)
    
    # Prepare sens1ini [n_dep x n_active] from cache
    s1ini <- NULL
    if (deriv && keep.root && !is.null(cache$sensini)) {
      s1ini <- cache$sensini[, active_sens, drop = FALSE]
    } else if (deriv) {
      s1ini <- default_sens[, active_sens, drop = FALSE]
    }
    
    # Integrate to steady state
    res <- tryCatch(
      withCallingHandlers(
        CppODE::solveODE(
          if (deriv) model_s else model,
          times = c(0, 1e100), parms = c(p[dependent], p[parms_all]),
          sens1ini = s1ini, fixed = fixed_char,
          roottol = ode_ctrl$roottol, abstol = ode_ctrl$abstol, reltol = ode_ctrl$reltol,
          maxsteps = as.integer(ode_ctrl$maxsteps), maxprogress = as.integer(ode_ctrl$maxprogress),
          hini = ode_ctrl$hini, maxroot = as.integer(ode_ctrl$maxroot)),
        warning = function(w) { warning(w$message, call. = FALSE); invokeRestart("muffleWarning") }),
      error = function(e) { warning("ODE integration failed: ", e$message, call. = FALSE); NULL })
    
    if (is.null(res)) {
      out <- if (attach.input) c(p[dependent], p[setdiff(names(p), dependent)]) else p[dependent]
      return(as.parvec(out, deriv = NULL))
    }
    
    last <- length(res$time)
    if (res$time[last] == 1e15)
      warning("Steady state not reached within integration time.", call. = FALSE)
    
    root <- res$variable[, last]
    names(root) <- dependent
    
    if (attach.input) {
      out <- c(root, p[setdiff(names(p), dependent)])
    } else { out <- root }
    
    # Update cache
    if (keep.root) {
      cache$yini <- root
      if (!is.null(res$sens1)) {
        full_new <- default_sens
        full_new[, active_sens] <- res$sens1[, , last]
        cache$sensini <- full_new
      } else { cache$sensini <- NULL }
    }
    
    if (deriv && !is.null(res$sens1)) {
      sens_final <- matrix(res$sens1[, , last], nrow = n_dep, ncol = n_active,
                           dimnames = list(dependent, active_sens))
      input_cols <- setdiff(names(p), c(dependent, names(fixed)))
      jacobian <- matrix(0, length(out), length(input_cols),
                         dimnames = list(names(out), input_cols))
      if (attach.input) {
        diag_idx <- intersect(emptypars, input_cols)
        if (length(diag_idx)) jacobian[cbind(diag_idx, diag_idx)] <- 1
      }
      sr <- intersect(dependent, rownames(sens_final))
      sc <- intersect(input_cols, colnames(sens_final))
      if (length(sr) && length(sc))
        jacobian[sr, sc] <- sens_final[sr, sc, drop = FALSE]
      if (!is.null(dP)) jacobian <- jacobian %*% submatrix(dP, rows = input_cols)
      nonzero <- rowSums(jacobian != 0) > 0
      as.parvec(out, deriv = jacobian[nonzero, , drop = FALSE])
    } else { as.parvec(out, deriv = NULL) }
  }
  
  attr(p2p, "equations")  <- as.eqnvec(f)
  attr(p2p, "parameters") <- parameters
  attr(p2p, "modelname")  <- modelname
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
  attr(out, "tree") <- table
  
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