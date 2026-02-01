

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
#' @importFrom data.table CJ
Xs.deSolve <- function(odemodel, forcings=NULL, events=NULL, names = NULL, condition = NULL, optionsOde=list(method = "lsoda"), optionsSens=list(method = "lsodes"), fcontrol = NULL) {
  
  func <- odemodel$func
  extended <- odemodel$extended
  if (is.null(extended)) warning("Element 'extended' empty. ODE model does not contain sensitivities.")
  
  myforcings <- forcings
  myevents <- events
  myfcontrol <- fcontrol
  
  if (!is.null(attr(func, "events")) & !is.null(myevents))
    warning("Events already defined in odemodel. Additional events in Xs() will be ignored. Events need to be defined in either odemodel() or Xs().")
  if (is.null(attr(func, "events")) & !is.null(myevents))
    message("Events should be definend in odemodel(). If defined in Xs(), events will be applied, but sensitivities will not be reset accordingly.")
  
  
  
  # Variable and parameter names
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  forcnames <- attr(func, "forcings")
  
  # Variable and parameter names of sensitivities
  sensvar <- attr(extended, "variables")[!attr(extended, "variables")%in%variables]
  senssplit <- strsplit(sensvar, ".", fixed=TRUE)
  senssplit.1 <- unlist(lapply(senssplit, function(v) v[1]))
  senssplit.2 <- unlist(lapply(senssplit, function(v) paste(v[-1], collapse = ".")))
  svariables <- intersect(senssplit.2, variables)
  sparameters <- setdiff(senssplit.2, variables)
  senspars <- c(svariables, sparameters)
  
  # Initial values for sensitivities
  yiniSens <- as.numeric(senssplit.1 == senssplit.2)
  names(yiniSens) <- sensvar
  
  # Names for deriv output
  sensGrid <- data.table::CJ(variables, senspars, sort = FALSE)
  
  # Only a subset of all variables/forcings is returned
  if (is.null(names)) names <- c(variables, forcnames)
  
  # Update sensNames when names are set
  sensGrid <- sensGrid[sensGrid[[1]] %in% names]
  
  # Controls to be modified from outside
  controls <- list(
    forcings = myforcings,
    events = myevents,
    names = names,
    optionsOde = optionsOde,
    optionsSens = optionsSens,
    fcontrol = myfcontrol,
    sensGrid = sensGrid
  )
  
  P2X <- function(times, pars, fixed = NULL, deriv=TRUE) {
    
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
    sensGrid <- controls$sensGrid
    
    # Add event time points (required by integrator) 
    times <- sort(union(unique(events$time), times))
    
    # Sort event time points
    if (!is.null(events)) events <- events[order(events$time),]
    
    if (!is.null(fixedNames))
      sensGrid <- sensGrid[!sensGrid[[2]] %in% fixedNames]
      
    senspars <- unique(sensGrid[[2]])
    
    myderivs <- NULL
    if (!deriv) {
      
      # Evaluate model without sensitivities
      if (!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
      out <- suppressWarnings(do.call(odeC, c(list(y = unclass(yini), times = times, func = func, parms = mypars, forcings = forc, events = list(data = events), fcontrol = fcontrol), optionsOde)))
      out <- submatrix(out, cols = c("time", names))
      
    } else {
      
      # Evaluate extended model
      if (!is.null(forcings)) forc <- setForcings(extended, forcings) else forc <- NULL
      outSens <- suppressWarnings(do.call(odeC, c(list(y = c(unclass(yini), yiniSens), times = times, func = extended, parms = mypars, 
                                                       forcings = forc, fcontrol = fcontrol,
                                                       events = list(data = events)), optionsSens)))
      
      out <- submatrix(outSens, cols = c("time", names))
      
      # Apply parameter transformation to the derivatives
      sensNames <- paste(sensGrid[[1]], sensGrid[[2]], sep=".")  
      sensLong <- matrix(outSens[,sensNames], ncol = length(senspars))
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        myderivs <- array(
          data = sensLong %*% dP[senspars,],
          dim = c(nrow(outSens), length(svariables), ncol(dP)),
          dimnames = list(NULL, svariables, colnames(dP))
        )
      } else {
        myderivs <- array(
          data = sensLong,
          dim = c(nrow(outSens), length(svariables), length(senspars)),
          dimnames = list(NULL, svariables, senspars)
        )
      }
      
    }
    
    prdframe(out, deriv = myderivs, parameters = pars)
    
  }
  
  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  
  
  prdfn(P2X, c(variables, parameters), condition) 
  
}

#' @export
Xs.boost <- function(odemodel, forcings = NULL, events = NULL, names = NULL, condition = NULL, 
                     optionsOde = list(), optionsSens = list()) {
  
  if(!is.null(forcings)) {
    if (!inherits(forcings, "data.frame")) {
      stop("'forcings' must be a data.frame, data.table, or tibble")
    }
    
    if (!all(c("name", "time", "value") %in% names(forcings))) {
      stop(
        "'forcings' must contain columns: ",
        paste(c("name", "time", "value"), collapse = ", ")
      )
    }
    
    if (!is.character(forcings$name)) {
      stop("'name' must be a character")
    }
    
    if (!is.numeric(forcings$time)) {
      stop("'time' must be numeric")
    }
    
    if (!is.numeric(forcings$value)) {
      stop("'value' must be numeric")
    }
    
    if (anyNA(forcings[, required_cols])) {
      stop("'forcings' contains NA values")
    }
    forcs <- split(
      forcings[, c("time", "value")],
      forcings$name
    )
  } else {
    forcs <- NULL
  }
  
  if (!is.null(events)) {
    stop("Events should be passed to odemodel() when using solver = 'boost'")
  }
  
  optionsDefault  <- list(atol = 1e-6, rtol = 1e-6, maxattemps = 100, maxsteps = 1e6, hini = 0, roottol = 1e-6, maxroot = 1)
  
  ## --- Warn about unknown options
  warn_unknown <- function(user, defaults, label) {
    bad <- setdiff(names(user), names(defaults))
    if (length(bad) > 0)
      warning(sprintf("%s: Ignoring unknown option(s): %s", label, paste(bad, collapse=", ")))
  }
  warn_unknown(optionsOde,  optionsDefault, "optionsOde")
  warn_unknown(optionsSens, optionsDefault, "optionsSens")
  
  ## --- Merge user-supplied options with defaults
  optionsOde   <- modifyList(optionsDefault, optionsOde)
  optionsSens  <- modifyList(optionsDefault, optionsSens)
  
  func <- odemodel$func
  extended <- odemodel$extended
  if (is.null(extended)) warning("Element 'extended' empty. ODE model does not contain sensitivities.")
  
  # Extract metadata
  paramNames <- c(attr(func, "variables"), attr(func, "parameters"))
  dim_names <- attr(func, "dim_names")
  dim_names_sens <- attr(extended, "dim_names")
  
  
  # Only a subset of all variables is returned
  if (is.null(names)) names <- dim_names$variable else names <- intersect(dim_names$variable, names)
  if (is.null(names)) stop(paste("Valid names are:", dim_names$variable))
  
  # Controls to be modified from outside
  controls <- list(
    forcings = forcs,
    names = names,
    optionsOde = optionsOde,
    optionsSens = optionsSens,
    sensnames = dim_names_sens$sens
  )
  
  P2X <- function(times, pars, fixed = NULL, deriv=TRUE) {
    
    fixedNames <- names(fixed)
    params <- c(unclass(pars), unclass(fixed))
    forcings <- controls$forcs
    names <- controls$names
    sensnames <- setdiff(controls$sensnames, fixedNames)
    nvars <- length(names)
    nsens <- length(sensnames)
    optionsOde <- controls$optionsOde 
    optionsSens <- controls$optionsSens 
    
    dX <- NULL
    if (!deriv) {
      
      # Evaluate model without sensitivities
      out <- suppressWarnings(
        CppODE::solveODE(func, times, params, NULL, NULL, NULL, forcings, optionsOde$atol, optionsOde$rtol,
                         optionsOde$maxattemps, optionsOde$maxsteps, optionsOde$hini,
                         optionsOde$roottol, optionsOde$maxroot)
      )
      
      out <- cbind(out$time, submatrix(out$variable, cols = names))
      colnames(out)[1] <- "time"
      
    } else {
      # Evaluate model with sensitivities
      outSens <- suppressWarnings(
        CppODE::solveODE(extended, times, params, NULL, NULL, fixedNames, forcings, optionsSens$atol, optionsSens$rtol,
                         optionsSens$maxattemps, optionsSens$maxsteps, optionsSens$hini,
                         optionsSens$roottol, optionsSens$maxroot)
      )
      
      out <- cbind(outSens$time, submatrix(outSens$variable, cols = names))
      colnames(out)[1] <- "time"
      
      ntimes <- nrow(out)
      
      # Apply parameter transformation to the derivatives (chain rule)
      mysensitivities <- outSens$sens1[, names, , drop = FALSE]
      M <- matrix(mysensitivities, nrow = ntimes * nvars, ncol = nsens)
      
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        dPsub <- dP [sensnames, , drop = FALSE]
        dX <- array(M %*% dPsub, dim = c(ntimes, nvars, ncol(dPsub)))
        dimnames(dX) <- list(NULL, names, colnames(dPsub))
      } else {
        dX <- mysensitivities
      }
      
    }
    prdframe(out, deriv = dX, parameters = params)
    
  }
  
  attr(P2X, "parameters") <- paramNames
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  
  
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
#' 
#' @return 
#' An object of class [obsfn], i.e. a function 
#' `g(..., deriv = TRUE, deriv2 = FALSE, condition = NULL, verbose = F)` representing the evaluation of the 
#' observation function. The function returns observable values and, if requested, 
#' their first- and second-order derivatives with respect to the parameters.
#' 
#' @example inst/examples/prediction.R
#' 
#' @importFrom CppODE funCpp
#' @importFrom abind abind
#' @export
Y <- function(g, f = NULL, states = NULL, parameters = NULL, condition = NULL,
              attach.input = TRUE, deriv = TRUE, deriv2 = FALSE,
              compile = FALSE, modelname = NULL, verbose = FALSE) {
  
  if (deriv2 && !deriv) {
    warning("deriv2 = TRUE requires deriv = TRUE. Setting deriv = TRUE automatically.")
    deriv <- TRUE
  }
  
  if (is.null(f) && is.null(states) && is.null(parameters)) 
    stop("Not all three arguments f, states and parameters can be NULL")
  
  # --- Define model name with condition suffix (to avoid name collisions) ---
  if (is.null(modelname)) modelname <- "obsfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  # --- Identify symbols in g ---
  symbols <- getSymbols(unclass(g))
  
  # --- Infer states and parameters ---
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
  
  # --- Compile evaluator for g (value, Jacobian, Hessian) ---
  gEval <- CppODE::funCpp(
    g,
    variables  = obsStates,
    parameters = obsParams,
    compile    = compile,
    modelname  = modelname,
    outdir = getwd(),
    verbose    = verbose,
    convenient = FALSE,
    warnings   = FALSE,
    deriv      = deriv,
    deriv2     = deriv2
  )
  
  controls <- list(attach.input = attach.input)
  
  # =========================================================================
  # Helper: batched matrix multiply - einsum("ija,iak->ijk", A, B)
  # A[i,j,a], B[i,a,k] -> C[i,j,k] where C[i,,] = A[i,,] %*% B[i,,]
  # =========================================================================
  bmm <- function(A, B) {
    # A: [n_i, n_j, n_a], B: [n_i, n_a, n_k]
    # Strategy: Stack all i-slices into block-diagonal structure
    # A_block: (n_i * n_j) x (n_i * n_a) block-diagonal
    # B_stack: (n_i * n_a) x n_k
    # Result:  (n_i * n_j) x n_k, reshape to [n_i, n_j, n_k]
    n_i <- dim(A)[1]; n_j <- dim(A)[2]; n_a <- dim(A)[3]; n_k <- dim(B)[3]
    
    # Transpose A to [i, a, j] then reshape to (n_i * n_a) x n_j
    A_t <- aperm(A, c(1, 3, 2))  # [i, a, j]
    A_mat <- matrix(A_t, nrow = n_i * n_a, ncol = n_j)
    
    # Reshape B to (n_i * n_a) x n_k  
    B_mat <- matrix(B, nrow = n_i * n_a, ncol = n_k)
    
    # crossprod: t(A_mat) %*% B_mat won't work because we need per-i products
    # Instead: use row-wise Kronecker structure
    # 
    # Alternative: reshape to do n_j separate (n_i x n_a) @ (n_i x n_a -> n_i x n_k) ops
    # But that's still n_j ops.
    #
    # Best approach for pure R: use vapply over i
    # But we can be smarter: reshape A to (n_i, n_j * n_a) and use tcrossprod tricks
    
    # Actually, the cleanest vectorization:
    # C[i,j,k] = sum_a A[i,j,a] * B[i,a,k]
    # Reshape A to (n_i * n_j, n_a), B to (n_i, n_a * n_k)
    # But indices don't align for simple matmul...
    #
    # Use: for each k, C[,,k] = rowSums(A * B[,,k expanded])
    # That's n_k operations but each is vectorized
    
    C <- array(0, dim = c(n_i, n_j, n_k))
    for (k in seq_len(n_k)) {
      # B[i,a,k] expanded to [i,j,a] by replicating along j
      B_k <- B[, , k, drop = FALSE]  # [n_i, n_a, 1]
      dim(B_k) <- c(n_i, n_a)
      # A[i,j,a] * B_k[i,a] summed over a -> need to align dims
      # Expand B_k to [i, 1, a] then broadcast
      B_exp <- array(B_k, dim = c(n_i, 1, n_a))
      # Now A * B_exp[expanded to j] and sum over a
      # Use: C[,, k] = A %*% diag-like structure... 
      # Simpler: direct rowSums
      for (i in seq_len(n_i)) {
        C[i, , k] <- A[i, , ] %*% B_k[i, ]
      }
    }
    C
  }
  
  # =========================================================================
  # Better helper: truly vectorized batched matmul using single BLAS call
  # einsum("ija,iak->ijk", A, B) where A[i,j,a], B[i,a,k]
  # =========================================================================
  bmm_vec <- function(A, B) {
    n_i <- dim(A)[1]; n_j <- dim(A)[2]; n_a <- dim(A)[3]; n_k <- dim(B)[3]
    
    # Key insight: C[i,j,k] = sum_a A[i,j,a] * B[i,a,k]
    # 
    # Reshape A to (n_i * n_j) x n_a  (rows indexed by (i,j), cols by a)
    # Reshape B to n_a x (n_i * n_k)  (rows by a, cols indexed by (i,k))
    # 
    # But we need the i's to match! Standard matmul won't do this.
    #
    # Trick: use sparse block structure or loop over smallest dimension
    
    # If n_j is small (typical: few observables), loop over j:
    if (n_j <= n_k && n_j <= n_i) {
      C <- array(0, dim = c(n_i, n_j, n_k))
      for (j in seq_len(n_j)) {
        # A[i,j,a] for fixed j: matrix n_i x n_a
        A_j <- A[, j, , drop = TRUE]
        if (n_i == 1) dim(A_j) <- c(1, n_a)
        # B[i,a,k]: reshape to (n_i * n_a) x n_k, but we need per-i products
        # C[i,j,k] = sum_a A_j[i,a] * B[i,a,k] = row-wise dot products
        # = (A_j * B_slice) summed over a for each k
        # Use: C[,j,] = rowSums over a of (A_j[i,a] * B[i,a,k])
        # Vectorized: element-wise multiply and sum
        for (k in seq_len(n_k)) {
          B_k <- B[, , k, drop = TRUE]
          if (n_i == 1) dim(B_k) <- c(1, n_a)
          C[, j, k] <- rowSums(A_j * B_k)
        }
      }
      return(C)
    }
    
    # If n_k is smallest, loop over k:
    if (n_k <= n_j && n_k <= n_i) {
      C <- array(0, dim = c(n_i, n_j, n_k))
      for (k in seq_len(n_k)) {
        B_k <- B[, , k, drop = TRUE]
        if (n_i == 1) dim(B_k) <- c(1, n_a)
        # C[i,j,k] = sum_a A[i,j,a] * B_k[i,a]
        # For each i: C[i,,k] = A[i,,] %*% B_k[i,]
        for (i in seq_len(n_i)) {
          C[i, , k] <- A[i, , , drop = TRUE] %*% B_k[i, ]
        }
      }
      return(C)
    }
    
    # Default: loop over i (original approach)
    C <- array(0, dim = c(n_i, n_j, n_k))
    for (i in seq_len(n_i)) {
      C[i, , ] <- A[i, , , drop = TRUE] %*% B[i, , , drop = TRUE]
    }
    C
  }
  
  # =========================================================================
  # Even better: use rowSums trick for full vectorization
  # einsum("ija,iak->ijk", A, B)
  # =========================================================================
  bmm_full <- function(A, B) {
    n_i <- dim(A)[1]; n_j <- dim(A)[2]; n_a <- dim(A)[3]; n_k <- dim(B)[3]
    
    # C[i,j,k] = sum_a A[i,j,a] * B[i,a,k]
    #
    # Expand A to [i,j,a,1] and B to [i,1,a,k], multiply, sum over a
    # But R doesn't broadcast 4D arrays efficiently...
    #
    # Better: loop over the smallest dimension
    if (n_j == 1) {
      # A is [i,1,a], treat as [i,a]
      A_mat <- matrix(A, nrow = n_i, ncol = n_a)
      # B is [i,a,k]
      # C[i,1,k] = sum_a A_mat[i,a] * B[i,a,k]
      C <- array(0, dim = c(n_i, 1, n_k))
      for (k in seq_len(n_k)) {
        C[, 1, k] <- rowSums(A_mat * B[, , k])
      }
      return(C)
    }
    
    # General case: loop over j (usually small = n_observables)
    C <- array(0, dim = c(n_i, n_j, n_k))
    B_mat <- matrix(B, nrow = n_i * n_a, ncol = n_k)  # [i*a, k]
    
    for (j in seq_len(n_j)) {
      # A[,j,] is [i, a]
      A_j <- matrix(A[, j, ], nrow = n_i, ncol = n_a)
      # C[i,j,k] = sum_a A_j[i,a] * B[i,a,k]
      # Vectorize over k: for each k, rowSums(A_j * B[,,k])
      for (k in seq_len(n_k)) {
        C[, j, k] <- rowSums(A_j * matrix(B[, , k], nrow = n_i, ncol = n_a))
      }
    }
    C
  }
  
  # --- Core observation mapping function ---
  X2Y <- function(out, pars, deriv = TRUE, deriv2 = FALSE, env = parent.frame()) {
    
    if (deriv2 && !deriv) {
      warning("deriv2 = TRUE requires deriv = TRUE. Setting deriv = TRUE automatically.")
      deriv <- TRUE
    }
    
    attach.input <- controls$attach.input
    
    outEval <- gEval(out[, obsStates, drop = FALSE], 
                     pars[obsParams], 
                     deriv = deriv, 
                     deriv2 = deriv2)
    
    # --- Observable values ---
    values <- cbind(time = out[,"time"], outEval$out)
    if (attach.input) values <- cbind(values, submatrix(out, cols = -1))
    
    # --- Compute first and second derivatives ---
    myderivs <- myderivs2 <- NULL
    if (deriv && !deriv2) {
      
      dGdX <- outEval$jacobian[, , obsStates, drop = FALSE]
      dGdP <- outEval$jacobian[, , obsParams, drop = FALSE]
      dX   <- attr(out,  "deriv")
      dP   <- attr(pars, "deriv")
      
      if (!is.null(dX)) dXsub <- dX[, obsStates, , drop = FALSE]
      if (!is.null(dP)) dPsub <- dP[obsParams, , drop = FALSE]
      
      outer_pars <- character(0)
      
      # ---------------------------------------------------------------------
      # CASE 1: dX ≠ NULL, dP ≠ NULL → full chain rule
      # ---------------------------------------------------------------------
      if (!is.null(dX) && !is.null(dP)) {
        n_i <- dim(dGdX)[1]; n_j <- dim(dGdX)[2]; n_a <- dim(dGdX)[3]
        n_b <- dim(dGdP)[3]; n_k <- dim(dXsub)[3]
        
        # term11: einsum("ija,iak->ijk", dGdX, dXsub) - batched matmul
        term11 <- bmm_full(dGdX, dXsub)
        
        # term12: einsum("ijb,bk->ijk", dGdP, dPsub) - simple matmul
        M_dGdP <- matrix(dGdP, nrow = n_i * n_j, ncol = n_b)
        term12 <- array(M_dGdP %*% dPsub, dim = c(n_i, n_j, n_k))
        
        myderivs <- term11 + term12
        outer_pars <- colnames(dP)
      }
      
      # ---------------------------------------------------------------------
      # CASE 2: dX ≠ NULL, dP = NULL
      # ---------------------------------------------------------------------
      if (!is.null(dX) && is.null(dP)) {
        term11 <- bmm_full(dGdX, dXsub)
        
        dyn_params   <- intersect(obsParams, dimnames(dX)[[3]])
        local_params <- setdiff(obsParams, dyn_params)
        
        if (length(local_params) > 0) {
          term12 <- dGdP[, , local_params, drop = FALSE]
          myderivs <- abind::abind(term11, term12, along = 3)
          outer_pars <- c(dimnames(dX)[[3]], local_params)
        } else {
          myderivs <- term11
          outer_pars <- dimnames(dX)[[3]]
        }
      }
      
      # ---------------------------------------------------------------------
      # CASE 3: dX = NULL, dP ≠ NULL
      # ---------------------------------------------------------------------
      if (is.null(dX) && !is.null(dP)) {
        n_i <- dim(dGdP)[1]; n_j <- dim(dGdP)[2]; n_b <- dim(dGdP)[3]
        n_k <- ncol(dPsub)
        M_dGdP <- matrix(dGdP, nrow = n_i * n_j, ncol = n_b)
        myderivs <- array(M_dGdP %*% dPsub, dim = c(n_i, n_j, n_k))
        outer_pars <- colnames(dP)
      }
      
      # ---------------------------------------------------------------------
      # CASE 4: dX = NULL, dP = NULL
      # ---------------------------------------------------------------------
      if (is.null(dX) && is.null(dP)) {
        myderivs <- dGdP
        outer_pars <- obsParams
      }
      
      if (!is.null(myderivs))
        dimnames(myderivs) <- list(NULL, observables, outer_pars)
      
      if (attach.input && !is.null(myderivs) && !is.null(dX)) {
        dyn_params <- dimnames(dX)[[3]]
        all_params <- outer_pars
        
        if (length(extra <- setdiff(all_params, dyn_params)) > 0) {
          pad <- array(0, dim = c(dim(dX)[1:2], length(extra)),
                       dimnames = list(NULL, NULL, extra))
          dX <- abind::abind(dX, pad, along = 3)
        }
        dX <- dX[, , all_params, drop = FALSE]
        myderivs <- abind::abind(myderivs, dX, along = 2)
      }
      
    } else if (deriv && deriv2) {
      
      dGdX <- outEval$jacobian[, , obsStates, drop = FALSE]
      dGdP <- outEval$jacobian[, , obsParams, drop = FALSE]
      dG2dX2  <- outEval$hessian[, , obsStates,  obsStates,  drop = FALSE]
      dG2dXdP <- outEval$hessian[, , obsStates,  obsParams, drop = FALSE]
      dG2dPdX <- outEval$hessian[, , obsParams, obsStates,  drop = FALSE]
      dG2dP2  <- outEval$hessian[, , obsParams, obsParams, drop = FALSE]
      
      dX  <- attr(out,  "deriv")
      dP  <- attr(pars, "deriv")
      dX2 <- attr(out,  "deriv2")
      dP2 <- attr(pars, "deriv2")
      
      if (!is.null(dX))  dXsub  <- dX[,  obsStates, , drop = FALSE]
      if (!is.null(dX2)) dX2sub <- dX2[, obsStates, , , drop = FALSE]
      if (!is.null(dP))  dPsub  <- dP[obsParams, , drop = FALSE]
      if (!is.null(dP2)) dP2sub <- dP2[obsParams, , , drop = FALSE]
      
      outer_pars <- character(0)
      
      # ---------------------------------------------------------------------
      # CASE 1: dX ≠ NULL, dP ≠ NULL → full chain rule
      # ---------------------------------------------------------------------
      if (!is.null(dX) && !is.null(dP)) {
        n_i <- dim(dGdX)[1]; n_j <- dim(dGdX)[2]
        n_a <- dim(dGdX)[3]; n_b <- dim(dGdP)[3]
        n_k <- ncol(dPsub)
        
        # --- First derivatives ---
        term11 <- bmm_full(dGdX, dXsub)
        M_dGdP <- matrix(dGdP, nrow = n_i * n_j, ncol = n_b)
        term12 <- array(M_dGdP %*% dPsub, dim = c(n_i, n_j, n_k))
        myderivs <- term11 + term12
        outer_pars <- colnames(dP)
        
        # --- Second derivatives ---
        # Vectorize over j (observables), loop over i (time points)
        # Most terms are bilinear forms: t(L) %*% H %*% R
        
        # Precompute Kronecker products for fully vectorized terms
        dPkron <- dPsub %x% dPsub  # (n_b^2) x (n_k^2)
        
        # Reshape dP2sub and dX2sub for vectorized contraction
        dP2sub_flat <- matrix(dP2sub, nrow = n_b, ncol = n_k * n_k)
        
        myderivs2 <- array(0, dim = c(n_i, n_j, n_k, n_k))
        
        # term24: einsum("ijbc,bk,cl->ijkl", dG2dP2, dPsub, dPsub)
        # = vec(result[i,j,,]) = (dPsub %x% dPsub)^T %*% vec(dG2dP2[i,j,,])
        # Fully vectorized:
        M_H_P2 <- matrix(dG2dP2, nrow = n_i * n_j, ncol = n_b * n_b)
        term24_flat <- M_H_P2 %*% dPkron
        term24 <- array(term24_flat, dim = c(n_i, n_j, n_k, n_k))
        
        # term26: einsum("ijb,bkl->ijkl", dGdP, dP2sub)
        # M_dGdP is (n_i * n_j) x n_b, dP2sub_flat is n_b x (n_k^2)
        term26_flat <- M_dGdP %*% dP2sub_flat
        term26 <- array(term26_flat, dim = c(n_i, n_j, n_k, n_k))
        
        # Remaining terms need the i-dependent dXsub - loop over i
        for (i in seq_len(n_i)) {
          dX_i <- matrix(dXsub[i, , ], nrow = n_a, ncol = n_k)
          dX2_i_flat <- matrix(dX2sub[i, , , ], nrow = n_a, ncol = n_k * n_k)
          
          for (j in seq_len(n_j)) {
            # term21: t(dX_i) %*% dG2dX2[i,j,,] %*% dX_i
            H_X2 <- matrix(dG2dX2[i, j, , ], nrow = n_a, ncol = n_a)
            t21 <- crossprod(dX_i, H_X2 %*% dX_i)
            
            # term22: t(dX_i) %*% dG2dXdP[i,j,,] %*% dPsub
            H_XP <- matrix(dG2dXdP[i, j, , ], nrow = n_a, ncol = n_b)
            t22 <- crossprod(dX_i, H_XP %*% dPsub)
            
            # term23: t(dPsub) %*% dG2dPdX[i,j,,] %*% dX_i
            H_PX <- matrix(dG2dPdX[i, j, , ], nrow = n_b, ncol = n_a)
            t23 <- crossprod(dPsub, H_PX %*% dX_i)
            
            # term25: dGdX[i,j,] %*% dX2_i_flat -> [n_k^2], reshape to [n_k, n_k]
            t25 <- matrix(dGdX[i, j, ] %*% dX2_i_flat, nrow = n_k, ncol = n_k)
            
            myderivs2[i, j, , ] <- t21 + t22 + t23 + term24[i, j, , ] + t25 + term26[i, j, , ]
          }
        }
        outer_pars2 <- outer_pars
      }
      
      # ---------------------------------------------------------------------
      # CASE 2: dX ≠ NULL, dP = NULL
      # ---------------------------------------------------------------------
      if (!is.null(dX) && is.null(dP)) {
        n_i <- dim(dGdX)[1]; n_j <- dim(dGdX)[2]
        n_a <- dim(dGdX)[3]
        Kx <- dim(dXsub)[3]
        
        term11 <- bmm_full(dGdX, dXsub)
        
        dyn_params   <- dimnames(dX)[[3]]
        local_params <- setdiff(obsParams, dyn_params)
        Kl <- length(local_params)
        Ktot <- Kx + Kl
        
        if (Kl > 0) {
          term12 <- dGdP[, , local_params, drop = FALSE]
          myderivs <- abind::abind(term11, term12, along = 3)
          outer_pars <- c(dyn_params, local_params)
        } else {
          myderivs <- term11
          outer_pars <- dyn_params
        }
        
        myderivs2 <- array(0, dim = c(n_i, n_j, Ktot, Ktot))
        
        for (i in seq_len(n_i)) {
          dX_i <- matrix(dXsub[i, , ], nrow = n_a, ncol = Kx)
          dX2_i_flat <- matrix(dX2sub[i, , , ], nrow = n_a, ncol = Kx * Kx)
          
          for (j in seq_len(n_j)) {
            # xx-block: term21 + term25
            H_X2 <- matrix(dG2dX2[i, j, , ], nrow = n_a, ncol = n_a)
            t21 <- crossprod(dX_i, H_X2 %*% dX_i)
            t25 <- matrix(dGdX[i, j, ] %*% dX2_i_flat, nrow = Kx, ncol = Kx)
            myderivs2[i, j, 1:Kx, 1:Kx] <- t21 + t25
            
            if (Kl > 0) {
              # x-local block
              H_XP <- matrix(dG2dXdP[i, j, , local_params], nrow = n_a, ncol = Kl)
              myderivs2[i, j, 1:Kx, Kx + seq_len(Kl)] <- crossprod(dX_i, H_XP)
              
              # local-x block
              H_PX <- matrix(dG2dPdX[i, j, local_params, ], nrow = Kl, ncol = n_a)
              myderivs2[i, j, Kx + seq_len(Kl), 1:Kx] <- H_PX %*% dX_i
              
              # local-local block
              myderivs2[i, j, Kx + seq_len(Kl), Kx + seq_len(Kl)] <- 
                dG2dP2[i, j, local_params, local_params]
            }
          }
        }
        outer_pars2 <- outer_pars
      }
      
      # ---------------------------------------------------------------------
      # CASE 3: dX = NULL, dP ≠ NULL
      # ---------------------------------------------------------------------
      if (is.null(dX) && !is.null(dP)) {
        n_i <- dim(dGdP)[1]; n_j <- dim(dGdP)[2]; n_b <- dim(dGdP)[3]
        n_k <- ncol(dPsub)
        
        M_dGdP <- matrix(dGdP, nrow = n_i * n_j, ncol = n_b)
        myderivs <- array(M_dGdP %*% dPsub, dim = c(n_i, n_j, n_k))
        outer_pars <- colnames(dP)
        
        # term24 + term26 - fully vectorized
        dPkron <- dPsub %x% dPsub
        dP2sub_flat <- matrix(dP2sub, nrow = n_b, ncol = n_k * n_k)
        
        M_H_P2 <- matrix(dG2dP2, nrow = n_i * n_j, ncol = n_b * n_b)
        term24_flat <- M_H_P2 %*% dPkron
        term26_flat <- M_dGdP %*% dP2sub_flat
        
        myderivs2 <- array(term24_flat + term26_flat, dim = c(n_i, n_j, n_k, n_k))
        outer_pars2 <- outer_pars
      }
      
      # ---------------------------------------------------------------------
      # CASE 4: dX = NULL, dP = NULL
      # ---------------------------------------------------------------------
      if (is.null(dX) && is.null(dP)) {
        myderivs  <- dGdP
        myderivs2 <- dG2dP2
        outer_pars <- obsParams
        outer_pars2 <- obsParams
      }
      
      if (!is.null(myderivs))
        dimnames(myderivs) <- list(NULL, observables, outer_pars)
      if (!is.null(myderivs2))
        dimnames(myderivs2) <- list(NULL, observables, outer_pars2, outer_pars2)
      
      if (attach.input && !is.null(dX) && !is.null(myderivs)) {
        dyn_params <- dimnames(dX)[[3]]
        all_params <- outer_pars
        
        if (length(extra <- setdiff(all_params, dyn_params)) > 0L) {
          pad <- array(0, dim = c(dim(dX)[1:2], length(extra)),
                       dimnames = list(NULL, NULL, extra))
          dX <- abind::abind(dX, pad, along = 3)
        }
        dX <- dX[, , all_params, drop = FALSE]
        myderivs <- abind::abind(myderivs, dX, along = 2)
        
        if (!is.null(myderivs2) && !is.null(dX2)) {
          dyn_params2 <- dimnames(dX2)[[3]]
          all_params2 <- outer_pars2
          
          n_i <- dim(dX2)[1]
          n_a <- dim(dX2)[2]
          Ktot <- length(all_params2)
          
          dX2_full <- array(0, dim = c(n_i, n_a, Ktot, Ktot),
                            dimnames = list(NULL, dimnames(dX2)[[2]], all_params2, all_params2))
          
          idx <- match(dyn_params2, all_params2)
          dX2_full[, , idx, idx] <- dX2
          
          myderivs2 <- abind::abind(myderivs2, dX2_full, along = 2)
        }
      }
    }
    
    prdframe(prediction = values, deriv = myderivs, deriv2 = myderivs2, parameters = pars)
  }
  
  attr(X2Y, "equations")  <- as.eqnvec(g)
  attr(X2Y, "parameters") <- parameters
  attr(X2Y, "states")     <- states
  attr(X2Y, "modelname")  <- modelname
  
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
    
    # Second derivatives (deriv2): 4D array [n_times, n_states, n_pars, n_pars]
    # All zeros since time is independent of parameters
    deriv2_arr <- array(0,
                        dim = c(n_times, 1, n_pars, n_pars),
                        dimnames = list(NULL, "time", par_names, par_names))
    
    prdframe(out, deriv = sens, deriv2 = deriv2_arr, parameters = pars)
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
