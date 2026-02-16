## Functions to generate parameter transformation ----


#' Generate a parameter transformation function
#' 
#' @description
#' This function provides a unified interface for generating condition-specific
#' parameter transformations, as commonly required in ODE-based modeling workflows.
#' 
#' `P()` can operate in two modes:
#' 
#' - **Explicit mode** (`method = "explicit"`, see [Pexpl]):
#'   Inner parameters are directly computed from symbolic expressions,
#'   for example
#'   \deqn{p_{\text{inner}} = \mathrm{parfn}(p_{\text{outer}})}
#'   A common application of the explicit mode is the log-transformation
#'   \deqn{p_{\text{outer}} \mapsto \exp(p_{\text{outer}})}
#'   which ensures positive parameters.
#'
#' - **Implicit mode** (`method = "implicit"`, see [Pimpl]):
#'   Typically used to infer initial values \eqn{p_{\text{ini}}} satisfying the
#'   steady-state condition
#'   \eqn{f(p_{\text{ini}}, p_{\text{dyn}}) = 0}.
#'   This yields an overall **partially implicit mapping**
#'   \deqn{p_{\text{dyn}} \mapsto (p_{\text{ini}}, p_{\text{dyn}})}
#'   where \eqn{f} usually represents the right-hand side (RHS) of an ODE model.
#' 
#' Both transformation types can be combined with other mappings via arithmetic
#' operators (`+` and `*`) thanks to the [parfn] interface.
#'
#' @param trafo object of class [eqnvec] or named character or list thereof. In case,
#' trafo is a list, [P()] is called on each element and conditions are assumed to be
#' the list names.
#' @param parameters character vector
#' @param condition character, the condition for which the transformation is generated
#' @param attach.input attach those incoming parameters to output which are not overwritten by
#' the parameter transformation.
#' @param keep.root logical, applies for \code{method = "implicit"}. The root of the last
#' evaluation of the parameter transformation function is saved as guess for the next 
#' evaluation.
#' @param compile logical, compile the function (see [CppODE::funCpp])
#' @param modelname character, see (see [CppODE::funCpp])
#' @param method character, either \code{"explicit"} or \code{"implicit"}
#' @param cores Number of cores for parallel method call of [Pexpl] or [Pimpl] per condtion
#' @param verbose Print out information during compilation
#'
#' @return
#' An object of class [parfn], representing the parameter transformation.
#'
#' @seealso
#' [Pexpl], [Pimpl], [parfn]
#'
#' @export
P <- function(trafo = NULL, parameters=NULL, condition = NULL, attach.input = FALSE,  keep.root = TRUE, compile = FALSE, 
              modelname = NULL, method = c("explicit", "implicit"), cores = detectFreeCores(), verbose = FALSE) {
  
  if (is.null(trafo)) return()
  if (!is.list(trafo)) {
    trafo <- list(trafo)
    names(trafo) <- condition
  }
  
  method <- match.arg(method)
  if (Sys.info()[['sysname']] == "Windows") cores <- 1
  Reduce("+", mclapply(1:length(trafo), function(i) {
    
    switch(method, 
           explicit = Pexpl(trafo = as.eqnvec(trafo[[i]]), parameters = parameters, attach.input = attach.input, condition = names(trafo[i]), compile = compile, modelname = modelname, verbose = verbose),
           implicit = Pimpl(trafo = as.eqnvec(trafo[[i]]), parameters = parameters, keep.root = keep.root, condition = names(trafo[i]), compile = compile, modelname = modelname, verbose = verbose))
    
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
      
      if (is.na(val))
        next
      
      expr <- paste0(par, " ~ ", val)
      
      if (apply == "insert") {
        out <- insert(out, expr, conditionMatch = cn)
      }
      
      if (apply == "define") {
        out <- define(out, expr, conditionMatch = cn)
      }
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