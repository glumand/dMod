## Methods of class odemodel

#' @export
print.odemodel <- function(x, ...) {

  func      <- x$func
  extended  <- x$extended
  extended2 <- x$extended2

  isCVODE <- inherits(x, "CppODE") && identical(attr(func, "backend"), "cvode")
  isCppODE <- inherits(x, "CppODE") && !isCVODE
  solver <- if (inherits(x, "deSolve")) "deSolve (cOde)"
            else if (isCVODE)           "Sundials (CVODE)"
            else if (isCppODE)          "CppODE"
            else                        "unknown"
  stepper <- attr(func, "method")

  suppressWarnings({

  cat("dMod odemodel\n", sep = "")
  cat("  Solver:  ", solver,
      if (!is.null(stepper)) paste0(" [", stepper, "]") else "", "\n", sep = "")
  cat("  Model:   ", as.character(func), "\n", sep = "")
  if (!is.null(extended)) {
    cat("  Sens1:   ", as.character(extended), "\n", sep = "")
  } else {
    cat("  Sens1:   not compiled (deriv = FALSE)\n", sep = "")
  }
  if (isCppODE) {
    if (!is.null(extended2)) {
      cat("  Sens2:   ", as.character(extended2), "\n", sep = "")
    } else {
      cat("  Sens2:   not compiled (deriv2 = FALSE)\n", sep = "")
    }
  }

  cat("\nEquations:\n", sep = "")
  print(as.eqnvec(attr(func, "equations")))
  cat("\nStates:\n", sep = "")
  print(sort(attr(func, "variables")))
  cat("\nParameters:\n", sep = "")
  print(sort(attr(func, "parameters")))
  forcs <- attr(func, "forcings")
  if (length(forcs) > 0) {
    cat("\nForcings:\n", sep = "")
    print(sort(forcs))
  }

  })

  invisible(x)
}


#' Generate model objects for use in Xs (models with sensitivities)
#'
#' Creates and compiles model objects for systems of ordinary differential equations (ODEs)
#' with optional first- and second-order sensitivities. Depending on the selected solver,
#' the function interfaces either to [cOde::funC()] (for `solver = "deSolve"`)
#' or to [CppODE::CppODE()] (for `solver = "CppODE"` or `solver = "Sundials"`).
#'
#' @param f Something that can be converted to [eqnvec], e.g. a named character vector
#'   specifying the right-hand sides of the ODE system.
#' @param deriv Logical. If `TRUE`, generate first-order sensitivities.
#'   Defaults to `TRUE`.
#' @param deriv2 Logical. If `TRUE`, also generate second-order sensitivities
#'   (requires `solver = "CppODE"`). Implies `deriv = TRUE`. Defaults to
#'   `FALSE`. The CVODE backend (`solver = "Sundials"`) does not support
#'   second-order sensitivities.
#'
#'   When `deriv2 = TRUE`, three native artefacts are emitted:
#'   `<modelname>` (value only), `<modelname>_s` (1st-order sensitivities,
#'   AD path) and `<modelname>_s2` (1st + 2nd-order sensitivities, AD2 path).
#'   Downstream prediction functions then dispatch on the call's `deriv2`
#'   flag, so callers can keep simulating / fitting at the cheaper
#'   1st-order cost via `deriv = TRUE, deriv2 = FALSE` even when the model
#'   was built with `deriv2 = TRUE`.
#' @param forcings Character vector with the names of external forcings.
#' @param events An [eventlist] (or `data.frame` coercible via [as.eventlist]).
#'   Must be defined here — not on [Xs()] — so that the sensitivity equations
#'   are extended consistently.
#' @param outputs Named character vector for additional output variables.
#' @param fixed Character vector with the names of parameters (initial values and dynamic)
#'   for which no sensitivities are required (this speeds up integration).
#' @param estimate Character vector specifying parameters (initial values and dynamic)
#'   for which sensitivities are returned. If specified, `estimate` overwrites `fixed`.
#' @param modelname Character. The base name of the generated C/C++ file.
#' @param solver Character string specifying the solver backend.
#'   One of `"deSolve"`, `"CppODE"` or `"Sundials"`.
#' @param gridpoints Integer specifying the minimum number of internal time points
#'   where the ODE is evaluated.
#' @param verbose Logical. If `TRUE`, print compiler output to the R console.
#' @param ... Additional arguments passed to [cOde::funC()] or
#'   [CppODE::CppODE()].
#'
#' @return list with \code{func} (ODE object) and \code{extended} (ODE+Sensitivities object).
#'   Carries a \code{"compileInfo"} attribute listing source files and per-file
#'   compile/link flags collected from \code{func} and \code{extended}. This is
#'   consumed by [compile()] when the model is later compiled via a prediction
#'   function, so solver-specific linker requirements (e.g. Sundials libraries
#'   for \code{solver = "Sundials"}) are applied to the right files only.
#'
#' @seealso [cOde::funC()], [CppODE::CppODE()]
#'
#' @example inst/examples/odemodel.R
#' @export
odemodel <- function(f, deriv = TRUE, deriv2 = FALSE, forcings=NULL, events = NULL, outputs = NULL,
                     fixed = NULL, estimate = NULL, modelname = "odemodel", solver = c("deSolve", "CppODE", "Sundials"),
                     gridpoints = NULL, verbose = FALSE, ...) {

  f <- as.eqnvec(f)
  solver <- match.arg(solver)

  if (deriv2 && !deriv) {
    warning("`deriv2 = TRUE` implies `deriv = TRUE`. Setting deriv = TRUE.",
            call. = FALSE)
    deriv <- TRUE
  }

  dots <- list(...)
  if (solver == "deSolve" && "nStack" %in% names(dots)) {
    warning("'nStack' is only supported for solver = 'CppODE' or 'Sundials' and will be ignored.",
            call. = FALSE)
    dots <- dots[setdiff(names(dots), "nStack")]
  }

  if (deriv2 && solver == "deSolve")
    stop("Second-order sensitivities require solver = 'CppODE'.")
  if (deriv2 && solver == "Sundials")
    stop("Second-order sensitivities are not available with CVODE; use solver = 'CppODE'.")

  if (solver == "deSolve") {

    if (is.null(gridpoints)) gridpoints <- 2
    func <- do.call(cOde::funC,
                    c(list(f, forcings = forcings, events = events, outputs = outputs,
                           fixed = fixed, modelname = modelname, solver = solver,
                           nGridpoints = gridpoints),
                      dots))
    extended <- NULL
    if (deriv) {
      modelname_s <- paste0(modelname, "_s")
      mystates <- attr(func, "variables")
      myparameters <- attr(func, "parameters")

      if (is.null(estimate) & !is.null(fixed)) {
        mystates <- setdiff(mystates, fixed)
        myparameters <- setdiff(myparameters, fixed)
      }

      if (!is.null(estimate)) {
        mystates <- intersect(mystates, estimate)
        myparameters <- intersect(myparameters, estimate)
      }

      s <- sensitivitiesSymb(f,
                             states = mystates,
                             parameters = myparameters,
                             inputs = attr(func, "forcings"),
                             events = attr(func, "events"),
                             reduce = TRUE)
      fs <- c(f, s)
      outputs <- c(attr(s, "outputs"), attr(func, "outputs"))

      events.sens <- attr(s, "events")
      events.func <- attr(func, "events")
      events <- NULL
      if (!is.null(events.func)) {
        if (is.data.frame(events.sens)) {
          events <- rbind(
            as.eventlist(events.sens),
            as.eventlist(events.func),
            stringsAsFactors = FALSE)
        } else {
          events <- do.call(rbind, lapply(1:nrow(events.func), function(i) {
            rbind(
              as.eventlist(events.sens[[i]]),
              as.eventlist(events.func[i,]),
              stringsAsFactors = FALSE)
          }))
        }

      }

      extended <- do.call(cOde::funC,
                          c(list(fs, forcings = forcings, modelname = modelname_s,
                                 solver = solver, nGridpoints = gridpoints,
                                 events = events, outputs = outputs),
                            dots))
    }
    out <- list(func = func, extended = extended)
    class(out) <- c("deSolve", "odemodel")
  }
  else {
    unsupported_args <- list(
      outputs = outputs,
      estimate = estimate,
      gridpoints = gridpoints
    )

    unsupported <- names(unsupported_args)[
      sapply(unsupported_args, function(arg) !is.null(arg) && !(is.logical(arg) && arg == FALSE))
    ]

    if (length(unsupported) > 0) {
      warning(sprintf("The following arguments are not (yet) supported by CppODE and will be ignored: %s", paste(unsupported, collapse = ", ")), call. = FALSE)
    }
    # nStack is a compile-time setting on the sensitivity system only. CppODE
    # rejects a finite nStack when deriv = FALSE, so strip it from the ...
    # list used for the func compile and keep it on the extended compile.
    dots_func <- dots[setdiff(names(dots), "nStack")]
    if (solver == "CppODE") {
      func <- do.call(CppODE::CppODE,
                      c(list(f, events = events, fixed = fixed, modelname = modelname,
                             outdir = getwd(), deriv = FALSE, verbose = verbose),
                        dots_func))
      extended <- NULL
      extended2 <- NULL
      if (deriv) {
        extended <- do.call(CppODE::CppODE,
                            c(list(f, events = events, fixed = fixed, forcings = forcings,
                                   modelname = paste0(modelname, "_s"), outdir = getwd(),
                                   deriv = TRUE, deriv2 = FALSE, verbose = verbose),
                              dots))
        if (deriv2) {
          extended2 <- do.call(CppODE::CppODE,
                               c(list(f, events = events, fixed = fixed, forcings = forcings,
                                      modelname = paste0(modelname, "_s2"), outdir = getwd(),
                                      deriv = TRUE, deriv2 = TRUE, verbose = verbose),
                                 dots))
        }
      }
      out <- list(func = func, extended = extended, extended2 = extended2)
      class(out) <- c("CppODE", "odemodel")
    } else if (solver == "Sundials") {
      func <- do.call(CppODE::CVODE,
                      c(list(f, events = events, fixed = fixed, modelname = modelname,
                             outdir = getwd(), deriv = FALSE, verbose = verbose),
                        dots_func))
      extended <- NULL
      if (deriv) {
        extended <- do.call(CppODE::CVODE,
                            c(list(f, events = events, fixed = fixed, forcings = forcings,
                                   modelname = paste0(modelname, "_s"), outdir = getwd(),
                                   deriv = TRUE, verbose = verbose),
                              dots))
      }
      out <- list(func = func, extended = extended)
      class(out) <- c("CppODE", "odemodel")
      }
  }
  attr(out, "compileInfo") <- collectCompileInfo(out$func, out$extended, out$extended2)
  return(out)
}
