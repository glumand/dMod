

#' Model prediction function for ODE models. 
#' @description Interface to combine an ODE and its sensitivity equations
#' into one model function `x(times, pars, deriv = TRUE)` returning ODE output and sensitivities.
#' @param odemodel object of class 'odemodel' or 'odemodel++', see [odemodel]
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric).
#' The ODE forcings. Not (yet) implemented for boost::odeint::rosenbrock4
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (numeric, time point), "value" (numeric, value), "method" (character, either
#' "replace", "add" or "multiply"). See [events][deSolve::events].
#' ATTENTION: Sensitivities for event states will only be correctly computed if defined within
#' [odemodel()]. Specify events within `Xs()` only for forward simulation.
#' @param names character vector with the states to be returned. If NULL, all states are returned.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @param optionsSens list with arguments to be passed to odeC() for integration of the extended system
#' @param fcontrol list with additional fine-tuning arguments for the forcing interpolation. 
#' See [approxfun][stats::approxfun] for possible arguments.
#' @return Object of class [prdfn]. If the function is called with parameters that
#' result from a parameter transformation (see [P]), the Jacobian of the parameter transformation
#' and the sensitivities of the ODE are multiplied according to the chain rule for
#' differentiation. The result is saved in the attributed "deriv", 
#' i.e. in this case the attibutes "deriv" and "sensitivities" do not coincide. 
#' @export
Xs <- function(odemodel, ...) {
  UseMethod("Xs", odemodel)
}

#' @export
Xs.deSolve <- function(odemodel, forcings = NULL, events = NULL, names = NULL, condition = NULL, 
                       optionsOde = list(method = "lsoda"), optionsSens = list(method = "lsodes"), 
                       fcontrol = NULL) {
  
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
  
  P2X <- function(times, pars, fixed = NULL, deriv = TRUE) {
    
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
      
      # Apply parameter transformation to the derivatives
      sensNames <- as.vector(outer(names, senspars, paste, sep = "."))
      mysensitivities <- aperm(array(outSens[, sensNames], dim = c(nrow(outSens), length(names), length(senspars))), c(2, 3, 1))
      
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        dPsub <- dP[senspars, , drop = FALSE]
        if(any(rownames(dP) %in% senspars)) {
          myderivs <- mysensitivities %bmm% dPsub
          dimnames(myderivs) <- list(names, colnames(dPsub), NULL)
        } else {
          myderivs <- NULL
        }
      } else {
        myderivs <- mysensitivities
        dimnames(myderivs) <- list(names, senspars, NULL)
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
Xs.CppODE <- function(odemodel, forcings = NULL, events = NULL, names = NULL, condition = NULL,
                      optionsOde = list(), optionsSens = list()) {
  
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
  
  optionsDefault <- list(atol = 1e-6, rtol = 1e-6, maxWithoutProgress = 50L, maxsteps = 1e6L,
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
  if (is.null(extended)) warning("Element 'extended' empty. ODE model does not contain sensitivities.")
  
  # Extract metadata
  paramNames <- c(attr(func, "variables"), attr(func, "parameters"))
  dim_names <- attr(func, "dim_names")
  dim_names_sens <- attr(extended, "dim_names")
  
  # Only a subset of all variables is returned
  if (is.null(names)) names <- dim_names$variable else names <- intersect(dim_names$variable, names)
  if (length(names) == 0) stop(paste("Valid names are:", paste(dim_names$variable, collapse = ", ")))
  
  # Controls to be modified from outside
  controls <- list(
    forcings = forcs,
    names = names,
    optionsOde = optionsOde,
    optionsSens = optionsSens,
    sensnames = dim_names_sens$sens
  )
  
  P2X <- function(times, pars, fixed = NULL, deriv = TRUE) {
    
    fixedNames <- intersect(names(fixed),controls$sensnames)
    params <- c(unclass(pars), unclass(fixed))
    forcings <- controls$forcings
    names <- controls$names
    sensnames <- setdiff(controls$sensnames, fixedNames)
    nvars <- length(names)
    nsens <- length(sensnames)
    optionsOde <- controls$optionsOde 
    optionsSens <- controls$optionsSens 
    
    dX <- NULL
    if (!deriv) {
      
      # Evaluate model without sensitivities
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

      out <- cbind(out$time, submatrix(t(out$variable), cols = names))
      colnames(out)[1] <- "time"

    } else {

      # Evaluate model with sensitivities
      outSens <- CppODE::solveODE(extended, times, params,
                                  sens1ini = NULL, sens2ini = NULL,
                                  fixed = fixedNames,
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
      
      out <- cbind(outSens$time, submatrix(t(outSens$variable), cols = names))
      colnames(out)[1] <- "time"
      
      ntimes <- nrow(out)
      
      # Apply parameter transformation to the derivatives (chain rule)
      mysensitivities <- outSens$sens1[names,,, drop = FALSE]
      
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        dPsub <- dP[sensnames, , drop = FALSE]
        dX <- mysensitivities %bmm% dPsub
        dimnames(dX) <- list(names, colnames(dPsub), NULL)
      } else {
        dX <- mysensitivities
      }
      
    }
    
    prdframe(out, deriv = dX, parameters = c(pars,fixed))
    
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
#' @description Interface to get an ODE 
#' into a model function `x(times, pars, forcings, events)` returning ODE output.
#' It is a reduced version of [Xs], missing the sensitivities. 
#' @param odemodel Object of class [odemodel].
#' @param forcings, see [Xs]
#' @param events, see [Xs]
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @param fcontrol list with additional fine-tuning arguments for the forcing interpolation. 
#' See [approxfun][stats::approxfun] for possible arguments.
#' @details Can be used to integrate additional quantities, e.g. fluxes, by adding them to `f`. 
#' All quantities that are not initialised by pars 
#' in `x(..., forcings, events)` are initialized with 0. For more details and
#' the return value see [Xs].
#' @export
Xf <- function(odemodel, forcings = NULL, events = NULL, condition = NULL, optionsOde=list(method = "lsoda"), fcontrol = NULL) {
  
  func <- odemodel$func
  
  myforcings <- forcings
  myevents <- events
  myfcontrol <- fcontrol
  
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  yini <- rep(0,length(variables))
  names(yini) <- variables
  
  # Controls to be modified from outside
  controls <- list(
    forcings = myforcings,
    events = myevents,
    optionsOde = optionsOde,
    fctonrol = myfcontrol
  )
  
  P2X <- function(times, pars, deriv = TRUE){
    
    events <- controls$events
    forcings <- controls$forcings
    optionsOde <- controls$optionsOde
    
    # Add event time points (required by integrator) 
    event.times <- unique(events$time)
    times <- sort(union(event.times, times))
    
    
    yini[names(pars[names(pars) %in% variables])] <- pars[names(pars) %in% variables]
    mypars <- pars[parameters]
    #alltimes <- unique(sort(c(times, forctimes)))
    
    # loadDLL(func)
    if(!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
    out <- suppressWarnings(do.call(odeC, c(list(y=yini, times=times, func=func, parms=mypars, forcings=forc,events = list(data = events), fcontrol = fcontrol), optionsOde)))
    #out <- cbind(out, out.inputs)      
    
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
  
  P2X <- function(times, pars, deriv=TRUE){
    
    
    predictions <- lapply(states, function(s) predL[[s]](times, pars)); names(predictions) <- states
    
    out <- cbind(times, do.call(cbind, predictions))
    colnames(out) <- c("time", states)
    
    mysensitivities <- NULL
    myderivs <- NULL
    if(deriv) {
      
      # Fill in sensitivities
      outSens <- matrix(0, nrow = length(times), ncol = length(sensNames), dimnames = list(NULL, c(sensNames)))
      for(s in states) {
        mysens <- attr(predictions[[s]], "sensitivities")
        mynames <- attr(predictions[[s]], "sensnames")
        outSens[, mynames] <- mysens
      }
      
      mysensitivities <- cbind(time = times, outSens)
      
      # Apply parameter transformation to the derivatives
      sensLong <- matrix(outSens, nrow = nrow(outSens)*length(states))
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        sensLong <- sensLong %*% submatrix(dP, rows = parameters)
        sensGrid <- expand.grid.alt(states, colnames(dP))
        sensNames <- paste(sensGrid[,1], sensGrid[,2], sep = ".")
      }
      outSens <- cbind(times, matrix(sensLong, nrow = dim(outSens)[1]))
      colnames(outSens) <- c("time", sensNames)
      
      myderivs <- outSens
      #attr(out, "deriv") <- outSens
    }
    
    #attr(out, "parameters") <- unique(sensGrid[,2])
    
    prdframe(out, deriv = myderivs, parameters = pars)
    
  }
  
  attr(P2X, "parameters") <- structure(parameters, names = NULL)
  attr(P2X, "pouter") <- pouter
  
  prdfn(P2X, attr(P2X, "parameters"), condition)
  
}


#' Observation functions. 
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
#' @param attach.input Logical, indicating whether the original model input should be 
#' included in the output.
#' @param deriv Logical, if `TRUE`, the function evaluates first-order derivatives
#' of observables with respect to parameters.
#' @param deriv2 Logical, if `TRUE`, the function also evaluates second derivatives 
#' of observables with respect to parameters.
#' @param compile Logical, if `TRUE`, the function is compiled (see [CppODE::funCpp]).
#' @param modelname Character, used if `compile = TRUE`, specifies a fixed filename
#' for the generated C file.
#' @param verbose Logical, print compiler output to the R console.
#' @param derivMode Character. Selects the derivative backend used by [funCpp]
#'   to evaluate observation Jacobians. One of `"ad"` (default, forward-mode
#'   automatic differentiation via `jac_chain` — typically faster when the
#'   number of fitted parameters is large; requires `compile = TRUE`),
#'   `"symbolic"` (classical SymPy Jacobian followed by an explicit chain
#'   rule against the upstream `dX`/`dP`), or `"none"` (no derivatives).
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
              derivMode = c("ad", "symbolic", "none")) {

  derivMode <- match.arg(derivMode)

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

  # Compile evaluator for g (value, Jacobian, AD chain)
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
      derivMode  = derivMode
    )
  )

  gfun       <- gEval$func
  gjac       <- gEval$jac
  gjac_chain <- gEval$jac_chain
  use_ad     <- derivMode == "ad"

  controls <- list(attach.input = attach.input)

  # Core observation mapping function
  X2Y <- function(out, pars, fixed = NULL, deriv = TRUE) {

    attach.input <- controls$attach.input
    fixedObsParams <- intersect(union(attr(pars, "fixed"), names(fixed)), obsParams)
    params <- c(unclass(pars), unclass(fixed))

    if (use_ad && deriv && !is.null(gjac_chain)) {
      # AD path: jac_chain returns y and dy already chain-ruled.
      dX <- attr(out, "deriv")
      # Gate dX the same way the symbolic path gates activeS: if no obsStates
      # appear in dX's row dim, it carries no upstream state sensitivity for
      # this observation function. Suppress to avoid spurious theta mismatches
      # against dP (e.g. Xt() returns a deriv array with unrelated layout).
      if (!is.null(dX) && length(intersect(obsStates, dimnames(dX)[[1]])) == 0)
        dX <- NULL
      ad_out <- gjac_chain(out[, obsStates, drop = FALSE], params[obsParams],
                           dX = dX, dP = attr(pars, "deriv"),
                           attach.input = attach.input, fixed = fixedObsParams)
      # Values: gjac_chain returns observables (and pass-through extras when
      # attach.input = TRUE) under attach.input semantics matching gfun.
      gAll <- ad_out$y
      gVal <- t(gAll[observables, , drop = FALSE])
      if (any(is.nan(gVal)))
        stop("Observable(s) evaluate to NaN: ",
             paste(observables[colSums(is.nan(gVal)) > 0], collapse = ", "),
             "\nLikely cause: division by zero or missing inputs.")
      values <- cbind(time = out[, "time"], gVal)
      if (attach.input) values <- cbind(values, submatrix(out, cols = -1))
      myderivs <- ad_out$dy
    } else {
      # Symbolic path (also serves "both" and "none").
      gVal <- t(gfun(out[, obsStates, drop = FALSE], params[obsParams], attach.input, fixedObsParams)[observables,, drop = FALSE])

      if (any(is.nan(gVal)))
        stop("Observable(s) evaluate to NaN: ",
             paste(observables[colSums(is.nan(gVal)) > 0], collapse = ", "),
             "\nLikely cause: division by zero or missing inputs.")

      values <- cbind(time = out[, "time"], gVal)
      if (attach.input) values <- cbind(values, submatrix(out, cols = -1))

      myderivs <- NULL
      if (deriv && !is.null(gjac)) {

        dX <- attr(out, "deriv")  # [states, theta, time] state sensitivities
        dP <- attr(pars, "deriv") # [p, theta] parameter transformation Jacobian
        dG <- gjac(out[, obsStates, drop = FALSE], params[obsParams]) # [obs, states+params, time]

        activeP <- setdiff(obsParams, fixedObsParams)
        activeS <- if (!is.null(dX)) intersect(obsStates, dimnames(dX)[[1]]) else character()
        theta <- if (!is.null(dP)) colnames(dP) else if (!is.null(dX)) dimnames(dX)[[2]] else NULL

        # Chain rule: dY/dtheta = dG/dX * dX/dtheta + dG/dP * dP/dtheta
        t1 <- if (length(activeS)) dG[,activeS,,drop=F] %bmm% dX[activeS,,,drop=F] else NULL
        t2 <- if (!is.null(dP) && length(activeP)) dG[,activeP,,drop=F] %bmm% dP[activeP,,drop=F] else NULL

        # Align by theta names before addition
        if (!is.null(t1)) dimnames(t1)[[2]] <- dimnames(dX)[[2]]
        if (!is.null(t2)) dimnames(t2)[[2]] <- colnames(dP)
        myderivs <- if (!is.null(t1) && !is.null(t2)) t1[,theta,,drop=F] + t2[,theta,,drop=F] else t1 %||% t2

        # Fallback: no upstream derivs, return dG/dp directly

        if (is.null(myderivs) && length(activeP)) myderivs <- dG[,activeP,,drop=F]
        if (!is.null(myderivs)) dimnames(myderivs) <- list(observables, theta, NULL)

        # Append original state sensitivities if attach.input
        if (attach.input && !is.null(myderivs) && !is.null(dX)) {
          outer_theta <- theta %||% dimnames(dX)[[2]]
          missing <- setdiff(outer_theta, dimnames(dX)[[2]])
          if (length(missing)) dX <- abind::abind(dX, array(0, c(dim(dX)[1], length(missing), dim(dX)[3]), dimnames = list(NULL, missing, NULL)), along = 2)
          myderivs <- abind::abind(myderivs, dX[, outer_theta, , drop = FALSE], along = 1)
        }
      }
    }

    prdframe(prediction = values, deriv = myderivs, parameters = c(pars, fixed))
  }

  attr(X2Y, "equations")  <- as.eqnvec(g)
  attr(X2Y, "parameters") <- parameters
  attr(X2Y, "states")     <- states
  attr(X2Y, "modelname")  <- modelname
  attr(X2Y, "compileInfo") <- collectCompileInfo(gfun, gjac, gjac_chain)

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
#' x <- Xt()
#' g <- Y(c(y = "a*time^2+b"), f = NULL, parameters = c("a", "b"))
#'
#' times <- seq(-1, 1, by = .05)
#' pars <- c(a = .1, b = 1)
#'
#' plot((g*x)(times, pars))
#' @export
Xt <- function(condition = NULL) {
  # Controls to be modified from outside
  controls <- list()
  P2X <- function(times, pars, deriv = TRUE, deriv2 = FALSE, ...) {
    n_times <- length(times)
    n_pars <- length(pars)
    par_names <- names(pars)
    
    # Output: matrix with time column
    out <- matrix(times, ncol = 1, dimnames = list(NULL, "time"))
    
    # Sensitivities (deriv): 3D array [n_times, n_states, n_pars]
    # time has no dependence on parameters, so all zeros
    sens <- array(0, 
                  dim = c(n_times, 1, n_pars),
                  dimnames = list(NULL, "time", par_names))
    
    prdframe(out, deriv = sens, parameters = pars)
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
