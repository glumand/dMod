

#' Model prediction function for ODE models. 
#' @description Interface to combine an ODE and its sensitivity equations
#' into one model function `x(times, pars, deriv = TRUE)` returning ODE output and sensitivities.
#' @param odemodel object of class 'odemodel' or 'odemodel++', see [odemodel]
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric).
#' The ODE forcings. Forcing support for the CppODE / Sundials backends depends
#' on the chosen solver method; see [CppODE::CppODE()].
#' @param events An [eventlist] (or `data.frame` coercible via [as.eventlist]).
#' Applied to the forward simulation only — sensitivities are not corrected.
#' Define events on [odemodel()] unless the prediction is used purely for
#' forward simulation.
#' @param names character vector with the states to be returned. If NULL, all states are returned.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @param optionsSens list with arguments to be passed to odeC() for integration of the extended system
#' @param fcontrol list with additional fine-tuning arguments for the forcing interpolation. 
#' See [approxfun][stats::approxfun] for possible arguments.
#' @param ... Additional arguments passed to methods.
#' @return Object of class [prdfn]. When called with transformed parameters
#'   (see [P]), the chain rule is applied automatically; the result is
#'   stored in `attr(., "deriv")` (which then differs from `"sensitivities"`).
#' @export
Xs <- function(odemodel, ...) {
  UseMethod("Xs", odemodel)
}

#' @export
#' @rdname Xs
Xs.deSolve <- function(odemodel, forcings = NULL, events = NULL, names = NULL, condition = NULL,
                       optionsOde = list(method = "lsoda"), optionsSens = list(method = "lsodes"),
                       fcontrol = NULL, ...) {
  
  func <- odemodel$func
  extended <- odemodel$extended
  if (is.null(extended)) warning("Element 'extended' empty. ODE model does not contain sensitivities.")
  
  myforcings <- forcings
  myevents <- events
  myfcontrol <- fcontrol
  
  if (!is.null(attr(func, "events")) & !is.null(myevents))
    warning("Events already defined in odemodel. Additional events in Xs() will be ignored.")
  if (is.null(attr(func, "events")) & !is.null(myevents))
    message("Events should be defined in odemodel(). If defined in Xs(), events will be applied, but sensitivities will not be reset accordingly.")
  
  # Variable and parameter names
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  forcnames <- attr(func, "forcings")
  
  # Variable and parameter names of sensitivities
  sensvar <- attr(extended, "variables")[!attr(extended, "variables") %in% variables]
  senssplit <- strsplit(sensvar, ".", fixed = TRUE)
  senssplit.1 <- unlist(lapply(senssplit, function(v) v[1]))
  senssplit.2 <- unlist(lapply(senssplit, function(v) paste(v[-1], collapse = ".")))
  svariables <- intersect(senssplit.2, variables)
  sparameters <- setdiff(senssplit.2, variables)
  senspars <- c(svariables, sparameters)
  
  # Initial values for sensitivities
  yiniSens <- as.numeric(senssplit.1 == senssplit.2)
  names(yiniSens) <- sensvar
  
  # Only a subset of all variables/forcings is returned
  if (is.null(names)) names <- c(variables, forcnames)
  
  # Controls to be modified from outside
  controls <- list(
    forcings = myforcings,
    events = myevents,
    names = names,
    optionsOde = optionsOde,
    optionsSens = optionsSens,
    fcontrol = myfcontrol
  )
  
  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {

    if (deriv2)
      stop("Xs.deSolve: second-order sensitivities require solver = 'CppODE'.")

    fixedNames <- names(fixed)
    params <- c(unclass(pars), unclass(fixed))
    yini <- params[variables]
    mypars <- params[parameters]
    
    forcings <- controls$forcings
    events <- controls$events
    optionsOde <- controls$optionsOde
    optionsSens <- controls$optionsSens
    fcontrol <- controls$fcontrol
    names <- controls$names
    
    # Add event time points (required by integrator) 
    times <- sort(union(unique(events$time), times))
    
    # Sort event time points
    if (!is.null(events)) events <- events[order(events$time), ]
    
    # Filter sensGrid by fixed parameters
    if (length(fixedNames)) senspars <- setdiff(senspars, fixedNames)
    
    myderivs <- NULL
    if (!deriv) {
      
      # Evaluate model without sensitivities
      if (!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
      out <- suppressWarnings(do.call(odeC, c(list(y = unclass(yini), times = times, func = func, 
                                                   parms = mypars, forcings = forc, 
                                                   events = list(data = events), fcontrol = fcontrol), optionsOde)))
      out <- submatrix(out, cols = c("time", names))
      
    } else {
      
      # Evaluate extended model
      if (!is.null(forcings)) forc <- setForcings(extended, forcings) else forc <- NULL
      outSens <- suppressWarnings(do.call(odeC, c(list(y = c(unclass(yini), yiniSens), times = times, 
                                                       func = extended, parms = mypars, 
                                                       forcings = forc, fcontrol = fcontrol,
                                                       events = list(data = events)), optionsSens)))
      
      out <- submatrix(outSens, cols = c("time", names))

      # Forcings live in `names` so the value matrix can return them, but the
      # extended ODE only carries sensitivities for state variables. Restrict
      # the deriv axis to states that are actually requested.
      svars <- intersect(names, variables)

      # Apply parameter transformation to the derivatives.
      # deSolve's outSens columns are laid out so that array() fills into
      # [time, variable, sensitivity] naturally under column-major.
      sensNames <- as.vector(outer(svars, senspars, paste, sep = "."))
      mysensitivities <- array(outSens[, sensNames],
                               dim = c(nrow(outSens), length(svars), length(senspars)))

      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        dPsub <- dP[senspars, , drop = FALSE]
        if(any(rownames(dP) %in% senspars)) {
          myderivs <- mysensitivities %bmm% dPsub
          dimnames(myderivs) <- list(NULL, svars, colnames(dPsub))
        } else {
          myderivs <- NULL
        }
      } else {
        myderivs <- mysensitivities
        dimnames(myderivs) <- list(NULL, svars, senspars)
      }
      
    }
    
    prdframe(out, deriv = myderivs, parameters = c(pars, fixed))
    
  }
  
  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  attr(P2X, "compileInfo") <- attr(odemodel, "compileInfo")

  prdfn(P2X, c(variables, parameters), condition)

}

#' @export
#' @rdname Xs
Xs.CppODE <- function(odemodel, forcings = NULL, events = NULL, names = NULL, condition = NULL,
                      optionsOde = list(), optionsSens = list(), ...) {
  
  if (!is.null(forcings)) {
    if (!inherits(forcings, "data.frame")) {
      stop("'forcings' must be a data.frame, data.table, or tibble")
    }
    
    if (!all(c("name", "time", "value") %in% names(forcings))) {
      stop("'forcings' must contain columns: ", paste(c("name", "time", "value"), collapse = ", "))
    }
    
    if (!is.character(forcings$name)) stop("'name' must be a character")
    if (!is.numeric(forcings$time)) stop("'time' must be numeric")
    if (!is.numeric(forcings$value)) stop("'value' must be numeric")
    
    required_cols <- c("name", "time", "value")
    if (anyNA(forcings[, required_cols])) stop("'forcings' contains NA values")
    
    forcs <- split(forcings[, c("time", "value")], forcings$name)
  } else {
    forcs <- NULL
  }
  
  if (!is.null(events)) {
    stop("Events should be passed to odemodel() when using solver = 'boost'")
  }
  
  optionsDefault <- list(atol = 1e-6, rtol = 1e-6, maxWithoutProgress = 20L, maxsteps = 1e6L,
                         hini = 0, roottol = 1e-6, maxroot = 1L,
                         usePID = "none", onFailure = "stop", traceFile = NULL)
  
  # Warn about unknown options
  warn_unknown <- function(user, defaults, label) {
    bad <- setdiff(names(user), names(defaults))
    if (length(bad) > 0)
      warning(sprintf("%s: Ignoring unknown option(s): %s", label, paste(bad, collapse = ", ")))
  }
  warn_unknown(optionsOde, optionsDefault, "optionsOde")
  warn_unknown(optionsSens, optionsDefault, "optionsSens")
  
  optionsOde <- modifyList(optionsDefault, optionsOde)
  optionsSens <- modifyList(optionsDefault, optionsSens)
  
  func <- odemodel$func
  extended <- odemodel$extended
  extended2 <- odemodel$extended2
  if (is.null(extended)) warning("Element 'extended' empty. ODE model does not contain sensitivities.")

  # Extract metadata
  paramNames <- c(attr(func, "variables"), attr(func, "parameters"))
  dim_names <- attr(func, "dimNames")
  dim_names_sens <- attr(extended, "dimNames")
  inner_names <- c(attr(func, "variables"), attr(func, "parameters"))

  # Only a subset of all variables is returned
  if (is.null(names)) names <- dim_names$variable else names <- intersect(dim_names$variable, names)
  if (length(names) == 0) stop(paste("Valid names are:", paste(dim_names$variable, collapse = ", ")))

  # Controls to be modified from outside
  controls <- list(
    forcings = forcs,
    names = names,
    optionsOde = optionsOde,
    optionsSens = optionsSens,
    sensnames = dim_names_sens$sens,
    inner_names = inner_names
  )

  has_deriv2 <- !is.null(extended2)

  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {

    if (deriv2 && !has_deriv2)
      stop("Xs.CppODE: model was compiled without deriv2; rebuild via odemodel(..., deriv2 = TRUE).")
    if (deriv2 && !deriv) deriv <- TRUE
    # Pick the cheapest extension that satisfies the requested derivative order.
    sens_model <- if (deriv2) extended2 else extended

    params <- c(unclass(pars), unclass(fixed))
    forcings <- controls$forcings
    names <- controls$names
    optionsOde <- controls$optionsOde
    optionsSens <- controls$optionsSens

    dX <- NULL
    dX2 <- NULL
    if (!deriv) {

      out <- CppODE::solveODE(func, times, params,
                              sens1ini = NULL, sens2ini = NULL, fixed = NULL,
                              forcings = forcings,
                              abstol = optionsOde$atol, reltol = optionsOde$rtol,
                              maxprogress = optionsOde$maxWithoutProgress,
                              maxsteps = optionsOde$maxsteps,
                              hini = optionsOde$hini,
                              roottol = optionsOde$roottol,
                              maxroot = optionsOde$maxroot,
                              usePID = optionsOde$usePID,
                              onFailure = optionsOde$onFailure,
                              traceFile = optionsOde$traceFile)

      out <- cbind(out$time, submatrix(out$variable, cols = names))
      colnames(out)[1] <- "time"

    } else {
      phi_rows <- inner_names
      deriv_in <- attr(pars, "deriv")
      sens1ini <- NULL
      sens2ini <- NULL
      if (!is.null(deriv_in)) {
        present  <- intersect(rownames(deriv_in), phi_rows)
        sens1ini <- matrix(0, length(phi_rows), ncol(deriv_in),
                           dimnames = list(phi_rows, colnames(deriv_in)))
        sens1ini[present, ] <- deriv_in[present, , drop = FALSE]
        if (deriv2) {
          d2 <- attr(pars, "deriv2")
          if (!is.null(d2)) {
            sens2ini <- array(0, c(length(phi_rows), dim(d2)[2], dim(d2)[3]),
                              dimnames = c(list(phi_rows), dimnames(d2)[2:3]))
            sens2ini[present, , ] <- d2[present, , , drop = FALSE]
          }
        }
      }

      outSens <- CppODE::solveODE(sens_model, times, params,
                                  sens1ini = sens1ini,
                                  sens2ini = sens2ini,
                                  fixed = NULL,
                                  forcings = forcings,
                                  abstol = optionsSens$atol, reltol = optionsSens$rtol,
                                  maxprogress = optionsSens$maxWithoutProgress,
                                  maxsteps = optionsSens$maxsteps,
                                  hini = optionsSens$hini,
                                  roottol = optionsSens$roottol,
                                  maxroot = optionsSens$maxroot,
                                  usePID = optionsSens$usePID,
                                  onFailure = optionsSens$onFailure,
                                  traceFile = optionsSens$traceFile)

      out <- cbind(outSens$time, submatrix(outSens$variable, cols = names))
      colnames(out)[1] <- "time"

      # outSens$sens1 is already [time, variable, theta] with theta colnames.
      dX <- outSens$sens1[, names, , drop = FALSE]
      if (deriv2 && !is.null(outSens$sens2))
        dX2 <- outSens$sens2[, names, , , drop = FALSE]

    }

    prdframe(out, deriv = dX, deriv2 = dX2, parameters = c(pars,fixed))

  }
  
  attr(P2X, "parameters") <- paramNames
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  attr(P2X, "compileInfo") <- attr(odemodel, "compileInfo")

  prdfn(P2X, paramNames, condition)
}


#' Model prediction function for ODE models without sensitivities.
#' @description Reduced version of [Xs] that returns the ODE output without
#' first- or second-order sensitivities. Dispatches on the [odemodel] class:
#' the `deSolve` method drives the cOde backend, the `CppODE` method drives
#' the CppODE / Sundials backend.
#' @param odemodel Object of class [odemodel].
#' @param forcings see [Xs].
#' @param events see [Xs].
#' @param condition either NULL (generic prediction for any condition) or a
#' character denoting the condition for which the function makes a prediction.
#' @param optionsOde list with arguments passed to the ODE integrator (deSolve
#' or [CppODE::solveODE]).
#' @param fcontrol list with additional fine-tuning arguments for the forcing
#' interpolation (cOde backend only). See [approxfun][stats::approxfun].
#' @param ... not used.
#' @details Can be used to integrate additional quantities, e.g. fluxes, by
#' adding them to `f`. All quantities not initialised by `pars` are initialised
#' to 0. For more details and the return value see [Xs].
#' @export
Xf <- function(odemodel, ...) {
  UseMethod("Xf", odemodel)
}

#' @export
#' @rdname Xf
Xf.deSolve <- function(odemodel, forcings = NULL, events = NULL, condition = NULL,
                       optionsOde = list(method = "lsoda"), fcontrol = NULL, ...) {

  func <- odemodel$func

  myforcings <- forcings
  myevents <- events
  myfcontrol <- fcontrol

  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  yini <- rep(0, length(variables))
  names(yini) <- variables

  controls <- list(
    forcings = myforcings,
    events = myevents,
    optionsOde = optionsOde,
    fcontrol = myfcontrol
  )

  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {

    if (deriv2)
      stop("Xf: second-order sensitivities are not implemented (use Xs() for deriv2).")

    events <- controls$events
    forcings <- controls$forcings
    optionsOde <- controls$optionsOde
    fcontrol <- controls$fcontrol

    # Xf carries no sensitivities, so fixed/free collapses to a single pars vector.
    pars <- c(unclass(pars), unclass(fixed))

    times <- sort(union(unique(events$time), times))

    yini[names(pars[names(pars) %in% variables])] <- pars[names(pars) %in% variables]
    mypars <- pars[parameters]

    if (!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
    out <- suppressWarnings(do.call(odeC, c(list(y = yini, times = times, func = func,
                                                 parms = mypars, forcings = forc,
                                                 events = list(data = events),
                                                 fcontrol = fcontrol), optionsOde)))

    prdframe(out, deriv = NULL, parameters = pars)
  }

  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  attr(P2X, "compileInfo") <- attr(odemodel, "compileInfo")

  prdfn(P2X, c(variables, parameters), condition)
}

#' @export
#' @rdname Xf
Xf.CppODE <- function(odemodel, forcings = NULL, events = NULL, condition = NULL,
                      optionsOde = list(), ...) {

  if (!is.null(forcings)) {
    if (!inherits(forcings, "data.frame"))
      stop("'forcings' must be a data.frame, data.table, or tibble")
    if (!all(c("name", "time", "value") %in% names(forcings)))
      stop("'forcings' must contain columns: name, time, value")
    if (!is.character(forcings$name)) stop("'name' must be a character")
    if (!is.numeric(forcings$time))   stop("'time' must be numeric")
    if (!is.numeric(forcings$value))  stop("'value' must be numeric")
    if (anyNA(forcings[, c("name", "time", "value")])) stop("'forcings' contains NA values")
    forcs <- split(forcings[, c("time", "value")], forcings$name)
  } else {
    forcs <- NULL
  }

  if (!is.null(events))
    stop("Events must be passed to odemodel() for solver = 'CppODE' / 'Sundials'.")

  optionsDefault <- list(atol = 1e-6, rtol = 1e-6, maxWithoutProgress = 20L, maxsteps = 1e6L,
                         hini = 0, roottol = 1e-6, maxroot = 1L,
                         usePID = "none", onFailure = "stop", traceFile = NULL)
  bad <- setdiff(names(optionsOde), names(optionsDefault))
  if (length(bad))
    warning(sprintf("optionsOde: Ignoring unknown option(s): %s", paste(bad, collapse = ", ")))
  optionsOde <- modifyList(optionsDefault, optionsOde)

  func <- odemodel$func
  paramNames <- c(attr(func, "variables"), attr(func, "parameters"))

  controls <- list(forcings = forcs, optionsOde = optionsOde)

  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {

    if (deriv2)
      stop("Xf: second-order sensitivities are not implemented (use Xs() for deriv2).")

    params <- c(unclass(pars), unclass(fixed))
    forcings <- controls$forcings
    optionsOde <- controls$optionsOde

    out <- CppODE::solveODE(func, times, params,
                            sens1ini = NULL, sens2ini = NULL, fixed = NULL,
                            forcings = forcings,
                            abstol = optionsOde$atol, reltol = optionsOde$rtol,
                            maxprogress = optionsOde$maxWithoutProgress,
                            maxsteps = optionsOde$maxsteps,
                            hini = optionsOde$hini,
                            roottol = optionsOde$roottol,
                            maxroot = optionsOde$maxroot,
                            usePID = optionsOde$usePID,
                            onFailure = optionsOde$onFailure,
                            traceFile = optionsOde$traceFile)

    out <- cbind(out$time, out$variable)
    colnames(out)[1] <- "time"

    prdframe(out, deriv = NULL, parameters = c(pars, fixed))
  }

  attr(P2X, "parameters") <- paramNames
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  attr(P2X, "compileInfo") <- attr(odemodel, "compileInfo")

  prdfn(P2X, paramNames, condition)
}


#' Model prediction function from data.frame
#' 
#' @param data data.frame with columns "name", "time", and row names that 
#' are taken as parameter names. The data frame can contain a column "value"
#' to initialize the parameters.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @return Object of class [prdfn], i.e. 
#' a function `x(times pars, deriv = TRUE, conditions = NULL)`, 
#' see also [Xs]. Attributes are "parameters", the parameter names (row names of
#' the data frame), and possibly "pouter", a named numeric vector which is generated
#' from `data$value`.
#' @examples
#' \dontrun{
#' # Generate a data.frame and corresponding prediction function
#' timesD <- seq(0, 2*pi, 0.5)
#' mydata <- data.frame(name = "A", time = timesD, value = sin(timesD), 
#'                      row.names = paste0("par", 1:length(timesD)))
#' x <- Xd(mydata)
#' 
#' # Evaluate the prediction function at different time points
#' times <- seq(0, 2*pi, 0.01)
#' pouter <- structure(mydata$value, names = rownames(mydata))
#' prediction <- x(times, pouter)
#' plot(prediction)
#' 
#' }
#' @export
Xd <- function(data, condition = NULL) {
  
  states <- unique(as.character(data$name))
  
  
  # List of prediction functions with sensitivities
  predL <- lapply(states, function(s) {
    subdata <- subset(data, as.character(name) == s)
    
    M <- diag(1, nrow(subdata), nrow(subdata))
    parameters.specific <- rownames(subdata)
    if(is.null(parameters.specific)) parameters.specific <- paste("par", s, 1:nrow(subdata), sep = "_")
    sensnames <- paste(s, parameters.specific, sep = ".")
    
    # return function
    out <- function(times, pars) {
      value <- approx(x = subdata$time, y = pars[parameters.specific], xout = times, rule = 2)$y
      grad <- do.call(cbind, lapply(1:nrow(subdata), function(i) {
        approx(x = subdata$time, y = M[, i], xout = times, rule = 2)$y
      }))
      colnames(grad) <- sensnames
      attr(value, "sensitivities") <- grad
      attr(value, "sensnames") <- sensnames
      return(value)
    }
    
    attr(out, "parameters") <- parameters.specific
    
    return(out)
    
  }); names(predL) <- states
  
  # Collect parameters
  parameters <- unlist(lapply(predL, function(p) attr(p, "parameters")))
  
  # Initialize parameters if available
  pouter <- NULL
  if(any(colnames(data) == "value")) 
    pouter <- structure(data$value[match(parameters, rownames(data))], names = parameters)
  
  sensGrid <- expand.grid(states, parameters, stringsAsFactors=FALSE)
  sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")  
  
  
  controls <- list()  
  
  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE){

    if (deriv2)
      stop("Xd: second-order sensitivities are not implemented for data-driven prediction.")
    # `fixed` is accepted for prdfn-wrapper symmetry; Xd is purely
    # data-grid-driven, so any fixed parameters are merged into `pars`
    # for the lookup.
    if (!is.null(fixed)) pars <- c(unclass(pars), unclass(fixed))


    predictions <- lapply(states, function(s) predL[[s]](times, pars)); names(predictions) <- states
    
    out <- cbind(times, do.call(cbind, predictions))
    colnames(out) <- c("time", states)
    
    myderivs <- NULL
    if (deriv) {

      # Fill in sensitivities — column layout is state-fastest, matching the
      # expand.grid(states, parameters) ordering used to build sensNames.
      outSens <- matrix(0, nrow = length(times), ncol = length(sensNames),
                        dimnames = list(NULL, sensNames))
      for (s in states) {
        mysens   <- attr(predictions[[s]], "sensitivities")
        mynames  <- attr(predictions[[s]], "sensnames")
        outSens[, mynames] <- mysens
      }

      # Reshape to 3D [time, state, param] (batch-first), matching Xs.
      myderivs <- array(outSens,
                        dim = c(length(times), length(states), length(parameters)),
                        dimnames = list(NULL, states, parameters))

      # Chain rule via upstream parameter transformation.
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        dPsub <- dP[parameters, , drop = FALSE]
        myderivs <- myderivs %bmm% dPsub
        dimnames(myderivs) <- list(NULL, states, colnames(dPsub))
      }
    }

    prdframe(out, deriv = myderivs, parameters = pars)
    
  }
  
  attr(P2X, "parameters") <- structure(parameters, names = NULL)
  attr(P2X, "pouter") <- pouter
  
  prdfn(P2X, attr(P2X, "parameters"), condition)
  
}


#' Observation functions
#'
#' @description 
#' Creates an object of type [obsfn] that evaluates an observation function
#' and, if requested, its first and second derivatives based on the output of a model 
#' prediction function, see [prdfn], as e.g. produced by [Xs].
#' 
#' @param g Named character vector or [eqnvec] defining the observation function.
#' @param f Named character vector of equations or an object that can be converted 
#' to [eqnvec], or an object of class 'fn'. If `f` is provided, states and parameters 
#' are automatically inferred from `f`.
#' @param states Character vector, alternative definition of state variables, usually 
#' the names of `f`. If both `f` and `states` are provided, the `states` argument 
#' overrides those derived from `f`.
#' @param parameters Character vector, alternative definition of parameters, usually 
#' the symbols contained in `g` and `f` except for `states` and the keyword `time`. 
#' If both `f` and `parameters` are provided, the `parameters` argument overrides those 
#' derived from `f` and `g`.
#' @param condition Either `NULL` (generic prediction for any condition) or a character 
#' string specifying the condition for which the function generates predictions.
#' @param attach.input Logical, indicating whether the original model input
#'   should be included in the output.
#' @param compile Logical, if `TRUE`, the function is compiled (see
#'   [CppODE::funCpp]).
#' @param modelname Character, used if `compile = TRUE`, specifies a fixed
#'   filename for the generated C file.
#' @param verbose Logical, print compiler output to the R console.
#' @param derivMode Character. Jacobian backend: `"dual"` (default,
#'   forward-mode AD; faster for many parameters; requires compiled native
#'   code) or `"symbolic"` (SymPy Jacobian + chain rule against upstream
#'   `dX`/`dP`; pure R).
#' @param deriv Logical. If `TRUE` (default), attach the first-order
#'   sensitivity `attr(., "deriv")` of shape `[time, observable, theta]`.
#' @param deriv2 Logical. If `TRUE`, attach a second-order derivative
#'   `attr(., "deriv2")` array of shape `[time, observable, theta, theta]`.
#'   Requires `deriv = TRUE`. Default `FALSE`.
#'
#' @return
#' An object of class [obsfn], i.e. a function  `g(..., fixed = NULL, deriv = TRUE, condition = NULL, env = NULL)`
#' returning predictions for observables and its derivatives.
#' 
#' @example inst/examples/prediction.R
#' 
#' @importFrom CppODE funCpp
#' @importFrom abind abind
#' @export
Y <- function(g, f = NULL, states = NULL, parameters = NULL,
              condition = NULL, attach.input = TRUE,
              compile = FALSE, modelname = NULL, verbose = FALSE,
              deriv = TRUE, deriv2 = FALSE,
              derivMode = c("dual", "symbolic")) {

  derivMode <- match.arg(derivMode)
  emit_d1 <- isTRUE(deriv)
  emit_d2 <- isTRUE(deriv2)
  if (emit_d2 && !emit_d1)
    stop("Y(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)

  if (is.null(f) && is.null(states) && is.null(parameters))
    stop("Not all three arguments f, states and parameters can be NULL")

  # Define model name with condition suffix
  if (is.null(modelname)) modelname <- "obsfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

  # Identify symbols in g
  symbols <- getSymbols(unclass(g))

  # Infer states and parameters
  if (is.null(f)) {
    states <- union(states, "time")
    parameters <- union(parameters, setdiff(symbols, states))
  } else if (inherits(f, "fn")) {
    myforcings <- Reduce(union, lapply(lapply(attr(f, "mappings"),
                                              function(m) attr(m, "forcings")),
                                       function(ff) as.character(ff$name)))
    mystates <- unique(c(do.call(c, lapply(getEquations(f), names)), "time"))
    if (length(intersect(myforcings, mystates)) > 0)
      stop("Forcings and states overlap in different conditions.")

    mystates <- c(mystates, myforcings)
    myparameters <- setdiff(union(getParameters(f), getSymbols(unclass(g))),
                            c(mystates, myforcings))
    states <- union(mystates, states)
    parameters <- union(myparameters, parameters)
  } else {
    f <- as.eqnvec(f)
    mystates <- union(names(f), "time")
    myparameters <- getSymbols(c(unclass(g), unclass(f)), exclude = mystates)
    states <- union(mystates, states)
    parameters <- union(myparameters, parameters)
  }

  observables <- names(g)
  obsParams <- intersect(symbols, parameters)
  obsStates <- setdiff(symbols, parameters)

  # Compile evaluator for g (value, Jacobian, Hessian, AD chain)
  gEval <- suppressWarnings(
    CppODE::funCpp(
      unclass(g),
      variables  = obsStates,
      parameters = obsParams,
      compile    = compile,
      modelname  = modelname,
      outdir     = getwd(),
      verbose    = verbose,
      convenient = FALSE,
      derivMode  = derivMode,
      deriv      = emit_d1,
      deriv2     = emit_d2
    )
  )

  gfun       <- gEval$func
  gjac       <- gEval$jac
  ghess      <- gEval$hess
  gevaluate  <- gEval$evaluate
  use_ad     <- derivMode == "dual"

  controls <- list(attach.input = attach.input)

  # Core observation mapping function
  X2Y <- function(out, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {

    if (deriv2 && !emit_d2)
      stop("Y() was built with deriv2 = FALSE; rebuild Y() with deriv2 = TRUE.", call. = FALSE)
    if (!emit_d1) deriv <- FALSE
    if (deriv2 && !deriv) deriv <- TRUE

    attach.input <- controls$attach.input
    fixedObsParams <- intersect(union(attr(pars, "fixed"), names(fixed)), obsParams)
    params <- c(unclass(pars), unclass(fixed))

    if (use_ad && deriv && !is.null(gevaluate)) {
      # AD path: evaluate() returns y and dy already chain-ruled via dX/dP seeds.
      dX_full <- attr(out, "deriv")
      dX2_full <- attr(out, "deriv2")
      # Gate dX the same way the symbolic path gates activeS: if no obsStates
      # appear in dX's state dim, it carries no upstream state sensitivity for
      # this observation function. Suppress to avoid spurious theta mismatches
      # against dP (e.g. Xt() returns a deriv array with unrelated layout).
      dX <- dX_full
      if (!is.null(dX) && length(intersect(obsStates, dimnames(dX)[[2]])) == 0)
        dX <- NULL
      dX2 <- dX2_full
      if (!is.null(dX2) && length(intersect(obsStates, dimnames(dX2)[[2]])) == 0)
        dX2 <- NULL
      ad_out <- gevaluate(out[, obsStates, drop = FALSE], params[obsParams],
                          dX = dX, dP = attr(pars, "deriv"),
                          dX2 = dX2, dP2 = attr(pars, "deriv2"),
                          deriv2 = deriv2,
                          attach.input = attach.input, fixed = fixedObsParams)
      # Values: evaluate() returns observables (and pass-through extras when
      # attach.input = TRUE) under attach.input semantics matching gfun.
      gAll <- ad_out$y
      gVal <- gAll[, observables, drop = FALSE]
      if (any(is.nan(gVal)))
        stop("Observable(s) evaluate to NaN: ",
             paste(observables[colSums(is.nan(gVal)) > 0], collapse = ", "),
             "\nLikely cause: division by zero or missing inputs.")
      values <- cbind(time = out[, "time"], gVal)
      if (attach.input) values <- cbind(values, submatrix(out, cols = -1))
      myderivs <- ad_out$dy
      myderivs2 <- if (deriv2) ad_out$d2y else NULL
      # Append pass-through state sensitivities for states that are attached
      # but not consumed by the observables; the AD path only emits sensitivities
      # for obsStates and would otherwise leave those rows missing.
      if (attach.input && !is.null(myderivs) && !is.null(dX_full)) {
        theta <- dimnames(myderivs)[[3]]
        outer_theta <- theta %||% dimnames(dX_full)[[3]]
        missing <- setdiff(outer_theta, dimnames(dX_full)[[3]])
        if (length(missing)) dX_full <- abind::abind(dX_full, array(0, c(dim(dX_full)[1], dim(dX_full)[2], length(missing)), dimnames = list(NULL, NULL, missing)), along = 3)
        already <- intersect(dimnames(myderivs)[[2]], dimnames(dX_full)[[2]])
        add_states <- setdiff(dimnames(dX_full)[[2]], already)
        if (length(add_states))
          myderivs <- abind::abind(myderivs, dX_full[, add_states, outer_theta, drop = FALSE], along = 2)
      }
      if (attach.input && !is.null(myderivs2) && !is.null(dX2_full)) {
        theta <- dimnames(myderivs2)[[3]]
        outer_theta <- theta %||% dimnames(dX2_full)[[3]]
        missing <- setdiff(outer_theta, dimnames(dX2_full)[[3]])
        if (length(missing)) {
          # Pad dim 3 (theta1) with zero blocks; keep dim 4 matching existing dX2_full.
          dX2_full <- abind::abind(dX2_full,
                                   array(0, c(dim(dX2_full)[1], dim(dX2_full)[2], length(missing), dim(dX2_full)[4]),
                                         dimnames = list(NULL, NULL, missing, dimnames(dX2_full)[[4]])),
                                   along = 3)
          # Pad dim 4 (theta2) with zero blocks; dim 3 now includes the missing entries.
          dX2_full <- abind::abind(dX2_full,
                                   array(0, c(dim(dX2_full)[1], dim(dX2_full)[2], dim(dX2_full)[3], length(missing)),
                                         dimnames = list(NULL, NULL, dimnames(dX2_full)[[3]], missing)),
                                   along = 4)
        }
        already <- intersect(dimnames(myderivs2)[[2]], dimnames(dX2_full)[[2]])
        add_states <- setdiff(dimnames(dX2_full)[[2]], already)
        if (length(add_states))
          myderivs2 <- abind::abind(myderivs2, dX2_full[, add_states, outer_theta, outer_theta, drop = FALSE], along = 2)
      }
    } else {
      # Symbolic path (also serves the !compile fallback for AD modes).
      gVal <- gfun(out[, obsStates, drop = FALSE], params[obsParams], attach.input, fixedObsParams)[, observables, drop = FALSE]

      if (any(is.nan(gVal)))
        stop("Observable(s) evaluate to NaN: ",
             paste(observables[colSums(is.nan(gVal)) > 0], collapse = ", "),
             "\nLikely cause: division by zero or missing inputs.")

      values <- cbind(time = out[, "time"], gVal)
      if (attach.input) values <- cbind(values, submatrix(out, cols = -1))

      myderivs <- NULL
      myderivs2 <- NULL
      if (deriv && !is.null(gjac)) {

        dX <- attr(out, "deriv")  # [time, states, theta] state sensitivities
        dP <- attr(pars, "deriv") # [p, theta] parameter transformation Jacobian
        dG <- gjac(out[, obsStates, drop = FALSE], params[obsParams]) # [time, obs, states+params]

        activeP <- setdiff(obsParams, fixedObsParams)
        activeS <- if (!is.null(dX)) intersect(obsStates, dimnames(dX)[[2]]) else character()
        theta <- if (!is.null(dP)) colnames(dP) else if (!is.null(dX)) dimnames(dX)[[3]] else NULL

        # Chain rule: dY/dtheta = dG/dX * dX/dtheta + dG/dP * dP/dtheta
        t1 <- if (length(activeS)) dG[,,activeS,drop=F] %bmm% dX[,activeS,,drop=F] else NULL
        t2 <- if (!is.null(dP) && length(activeP)) dG[,,activeP,drop=F] %bmm% dP[activeP,,drop=F] else NULL

        # Align by theta names before addition
        if (!is.null(t1)) dimnames(t1)[[3]] <- dimnames(dX)[[3]]
        if (!is.null(t2)) dimnames(t2)[[3]] <- colnames(dP)
        myderivs <- if (!is.null(t1) && !is.null(t2)) t1[,,theta,drop=F] + t2[,,theta,drop=F] else t1 %||% t2

        # Fallback: no upstream derivs, return dG/dp directly

        if (is.null(myderivs) && length(activeP)) myderivs <- dG[,,activeP,drop=F]
        if (!is.null(myderivs)) dimnames(myderivs) <- list(NULL, observables, theta)

        # Append original state sensitivities if attach.input
        if (attach.input && !is.null(myderivs) && !is.null(dX)) {
          outer_theta <- theta %||% dimnames(dX)[[3]]
          missing <- setdiff(outer_theta, dimnames(dX)[[3]])
          if (length(missing)) dX <- abind::abind(dX, array(0, c(dim(dX)[1], dim(dX)[2], length(missing)), dimnames = list(NULL, NULL, missing)), along = 3)
          myderivs <- abind::abind(myderivs, dX[, , outer_theta, drop = FALSE], along = 2)
        }

        # Second-order: delegate the full sandwich to ghess(), which applies
        # chain_hess_sym internally. We pass upstream seeds (dX, dP, dX2, dP2)
        # and let CppODE produce d2y already aligned to theta.
        if (deriv2) {
          if (is.null(ghess))
            stop("Y(deriv2 = TRUE) requires hess(); rebuild Y with deriv2 = TRUE.")
          dX2_in <- attr(out, "deriv2")
          dP2_in <- attr(pars, "deriv2")
          gH <- ghess(out[, obsStates, drop = FALSE], params[obsParams],
                      dX = dX, dP = dP, dX2 = dX2_in, dP2 = dP2_in)
          # gH: [time, obs, theta, theta]. Restrict columns to declared observables
          # (gjac path output above uses [, , observables, theta]).
          if (!is.null(gH)) {
            obs_cols <- intersect(observables, dimnames(gH)[[2]])
            if (length(obs_cols)) gH <- gH[, obs_cols, , , drop = FALSE]
            myderivs2 <- gH

            # Pass-through state Hessians for attach.input
            if (attach.input && !is.null(dX2_in)) {
              outer_theta <- dimnames(myderivs2)[[3]] %||% dimnames(dX2_in)[[3]]
              missing <- setdiff(outer_theta, dimnames(dX2_in)[[3]])
              if (length(missing)) dX2_in <- abind::abind(dX2_in,
                                                          array(0, c(dim(dX2_in)[1], dim(dX2_in)[2], length(missing), length(missing)),
                                                                dimnames = list(NULL, NULL, missing, missing)),
                                                          along = 3)
              already <- intersect(dimnames(myderivs2)[[2]], dimnames(dX2_in)[[2]])
              add_states <- setdiff(dimnames(dX2_in)[[2]], already)
              if (length(add_states))
                myderivs2 <- abind::abind(myderivs2, dX2_in[, add_states, outer_theta, outer_theta, drop = FALSE], along = 2)
            }
          }
        }
      }
    }

    prdframe(prediction = values, deriv = myderivs, deriv2 = myderivs2, parameters = c(pars, fixed))
  }

  attr(X2Y, "equations")  <- as.eqnvec(g)
  attr(X2Y, "parameters") <- parameters
  attr(X2Y, "states")     <- states
  attr(X2Y, "modelname")  <- modelname
  attr(X2Y, "compileInfo") <- collectCompileInfo(gfun, gjac, ghess, gevaluate)

  obsfn(X2Y, parameters, condition)
}

 
#' Generate a prediction function that returns times
#'
#' Function to deal with non-ODE models within the framework of dMod. See example.
#'
#' @param condition  either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @return Object of class [prdfn].
#' @examples
#' \dontrun{
#' x <- Xt()
#' g <- Y(c(y = "a*time^2+b"), f = NULL, parameters = c("a", "b"),
#'        compile = TRUE, modelname = "Xt_example_obs")
#'
#' times <- seq(-1, 1, by = .05)
#' pars <- c(a = .1, b = 1)
#'
#' plot((g*x)(times, pars))
#' }
#' @export
Xt <- function(condition = NULL) {
  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
    n_times <- length(times)
    par_names <- names(pars)
    n_pars <- length(par_names)

    out <- matrix(times, ncol = 1, dimnames = list(NULL, "time"))

    # time has no parameter dependence — both sens1 and sens2 are zero arrays
    # in batch-first [time, observable, ...] layout matching Xs.
    sens  <- array(0, dim = c(n_times, 1, n_pars),
                   dimnames = list(NULL, "time", par_names))
    sens2 <- if (deriv2)
      array(0, dim = c(n_times, 1, n_pars, n_pars),
            dimnames = list(NULL, "time", par_names, par_names))
    else NULL

    prdframe(out, deriv = sens, deriv2 = sens2, parameters = pars)
  }
  attr(P2X, "parameters") <- NULL
  attr(P2X, "equations") <- NULL
  attr(P2X, "forcings") <- NULL
  attr(P2X, "events") <- NULL
  prdfn(P2X, NULL, condition)
}


#' An identity function which vanishes upon concatenation of fns
#'
#' @return fn of class idfn
#' @export
#'
#' @examples
#' x <- Xt()
#' id <- Id()
#'
#' (id*x)(1:10, pars = c(a = 1))
#' (x*id)(1:10, pars = c(a = 1))
#' str(id*x)
#' str(x*id)
Id <- function() {
  outfn <- function() return(NULL)
  class(outfn) <- c("idfn", "fn")
  return(outfn)
}
