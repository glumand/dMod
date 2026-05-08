
## Methods for class prdlist ------------------------------------------------



#' @export
#' @rdname prdlist
as.prdlist <- function(x, ...) {
  UseMethod("as.prdlist", x)
}

#' @export
#' @param x list of prediction frames
#' @param names character vector, the list names, e.g. the names of the experimental
#' @rdname prdlist
as.prdlist.list <- function(x = NULL, names = NULL, ...) {

  if (is.null(x)) x <- list()
  if (is.null(names)) mynames <- names(x) else mynames <- names 

  # if (length(mynames) != length(x)) stop("names argument has wrong length")

  ## Prepare output
  names(x) <- mynames
  class(x) <- c("prdlist", "list")

  return(x)

}


#' @export
c.prdlist <- function(...) {
  
  mylist <- list(...)
  mylist <- lapply(mylist, unclass)
  newlist <- do.call(c, mylist)
  
  as.prdlist(newlist)
  
}

#' @export
"[.prdlist" <- function(x, ...) {
  out <- unclass(x)[...]
  class(out) <- c("prdlist", "list")
  return(out)
}







#' @export
print.prdlist <- function(x, ...) {
  
  mynames <- names(x)
  if (is.null(mynames)) mynames <- rep("NULL", length(x))
  
  for (i in 1:length(x)) {
    cat(mynames[i], ":\n", sep = "")
    print(x[[i]])
  }
  
}


#' @export
#' @param data data list oject
#' @param errfn obsfn object, the error model function to predict sigma
#' @param ... not used right now
#' @rdname as.data.frame.dMod
as.data.frame.prdlist <- function(x, ..., data = NULL, errfn = NULL) {
  
  prediction <- x
  sigma <- NULL
  condition.grid <- attr(data, "condition.grid")
  
  if (!is.null(errfn)) {
    sigma <- as.prdlist(
      lapply(1:length(prediction), 
             function(i) errfn(prediction[[i]], 
                               getParameters(prediction[[i]]), 
                               conditions = names(prediction)[i])[[1]]),
      names = names(prediction)
    )
    sigma <- wide2long(sigma)
  }
  
  prediction <- wide2long(prediction)
  prediction$sigma <- NaN
  if (!is.null(sigma)) {
    common <- intersect(unique(prediction$name), unique(sigma$name))
    prediction$sigma[prediction$name %in% common] <- sigma$value[sigma$name %in% common]
  }
  
  if (!is.null(condition.grid)) {
    for (C in colnames(condition.grid)) {
      rows <- ifelse(is.na(prediction$condition), 1, as.character(prediction$condition))
      prediction[, C] <- condition.grid[rows, C]
    }
    n1 <- nrow(prediction)
  }
  
  
  return(prediction)
  
  
} 

#' @export
#' @param x prediction
#' @rdname plotCombined
plot.prdlist <- function(x, data = NULL, ..., scales = "free", facet = c("wrap", "grid", "wrap_plain"), transform = NULL) {
  
  prediction <- x
  
  if (is.null(names(prediction))) names(prediction) <- paste0("C", 1:length(prediction))
  if (!is.null(data) && is.null(names(data))) names(data) <- paste0("C", 1:length(data))
  
  plotCombined.prdlist(prediction = prediction, data = data, ..., scales = scales, facet = facet, transform = transform)
  
}


#' Plot combined prediction and data
#'
#' @description
#' Creates a combined plot of model predictions and observed data with flexible
#' faceting options. Supports error bars, below limit of quantification (BLoQ) 
#' indicators, and coordinate transformations.
#'
#' @param prediction A \code{prdlist} object containing model predictions.
#' @param data Optional data object with observed values. If provided, will be 
#'   merged with covariate information and displayed as points with error bars.
#' @param ... Filter expressions passed to \code{dplyr::filter} for subsetting
#'   both prediction and data.
#' @param scales Scale specification for facets, passed to facet functions.
#'   Default is \code{"free"}.
#' @param facet Faceting style. One of:
#'   \itemize{
#'     \item \code{"wrap"}: Facet by name, color by condition (default)
#'     \item \code{"grid"}: Facet grid with name as rows, condition as columns
#'     \item \code{"wrap_plain"}: Facet wrap by name and condition combined
#'   }
#' @param transform Optional transformation function applied to coordinates
#'   via \code{coordTransform}.
#' @param aesthetics Optional named list of aesthetic mappings (as strings) to
#'   override defaults. Default aesthetics include x, y, ymin, ymax, and 
#'   conditionally group and color.
#'
#' @return A \code{ggplot} object with an additional \code{"data"} attribute
#'   containing a list with the processed \code{data} and \code{prediction} 
#'   data frames.
#'
#' @examples
#' \dontrun{
#' plotCombined(pred, mydata, time < 100, facet = "grid")
#' plotCombined(pred, mydata, aesthetics = list(color = "treatment"))
#' }
#'
#' @export
#' @rdname plotCombined
#' @importFrom dplyr filter
#' @importFrom rlang parse_expr
plotCombined.prdlist <- function(prediction, data = NULL, ..., 
                                 scales = "free", 
                                 facet = c("wrap", "grid", "wrap_plain"), 
                                 transform = NULL, aesthetics = NULL) {
  
  facet <- match.arg(facet)
  
  make_aes <- function(mapping) {
    mapping <- lapply(mapping, function(col) {
      if (!is.null(col)) rlang::parse_expr(col)
    })
    do.call(aes, Filter(Negate(is.null), mapping))
  }
  
  mynames <- c("time", "name", "value", "sigma", "condition")
  covtable <- NULL
  
  # --- Prepare data ---
  if (!is.null(data)) {
    covtable <- covariates(data)
    covtable <- cbind(condition = rownames(covtable), covtable)
    covtable <- covtable[!duplicated(names(covtable))]
    
    data <- lbind(data)
    data <- base::merge(data, covtable, by = "condition", all.x = TRUE)
    data <- as.data.frame(dplyr::filter(data, ...))
    data$bloq <- ifelse(data$value <= data$lloq, "yes", "no")
    
    if (!is.null(transform)) data <- coordTransform(data, transform)
  }
  
  # --- Prepare prediction ---
  if (!is.null(prediction)) {
    prediction <- cbind(wide2long(prediction), sigma = NA)
    if (!is.null(covtable)) {
      prediction <- base::merge(prediction, covtable, by = "condition", all.x = TRUE)
    }
    prediction <- as.data.frame(dplyr::filter(prediction, ...))
    
    if (!is.null(transform)) prediction <- coordTransform(prediction, transform)
  }
  
  # --- Combine into single data frame ---
  keep_cols <- unique(c(mynames, names(covtable)))
  total <- rbind(
    if (!is.null(prediction)) prediction[, keep_cols] else NULL,
    if (!is.null(data)) data[, keep_cols] else NULL
  )
  
  # --- Build aesthetics ---
  aes_base <- list(x = "time", y = "value", 
                   ymin = "value - sigma", ymax = "value + sigma")
  if (facet == "wrap") {
    aes_base$group <- "condition"
    aes_base$color <- "condition"
  }
  aesthetics <- c(aes_base[setdiff(names(aes_base), names(aesthetics))], aesthetics)
  
  # --- Construct plot ---
  p <- ggplot(total, make_aes(aesthetics))
  
  p <- p + switch(facet,
                  wrap       = facet_wrap(~name, scales = scales),
                  grid       = facet_grid(name ~ condition, scales = scales),
                  wrap_plain = facet_wrap(~name * condition, scales = scales)
  )
  
  if (!is.null(prediction)) {
    p <- p + geom_line(data = prediction)
  }
  
  if (!is.null(data)) {
    p <- p + 
      geom_point(data = data, aes(pch = bloq)) + 
      geom_errorbar(data = data, width = 0) +
      scale_shape_manual(name = "BLoQ", values = c(yes = 4, no = 19))
    
    if (all(data$bloq == "no")) {
      p <- p + guides(shape = "none")
    }
  }
  
  attr(p, "data") <- list(data = data, prediction = prediction)
  p
}


#' Plot model predictions
#'
#' @description
#' Creates a plot of model predictions with optional error bands from an error
#' model function. Supports flexible faceting and coordinate transformations.
#'
#' @param prediction A \code{prdlist} object containing model predictions.
#' @param ... Filter expressions passed to \code{dplyr::filter} for subsetting
#'   the predictions.
#' @param errfn Optional error model function. If provided, predictions are
#'   augmented with sigma values and displayed with ribbon error bands.
#' @param scales Scale specification for facets, passed to facet functions.
#'   Default is \code{"free"}.
#' @param facet Faceting style. One of:
#'   \itemize{
#'     \item \code{"wrap"}: Facet by name, color by condition (default)
#'     \item \code{"grid"}: Facet grid with name as rows, condition as columns
#'   }
#' @param transform Optional transformation function applied to coordinates
#'   via \code{coordTransform}.
#'
#' @return A \code{ggplot} object with an additional \code{"data"} attribute
#'   containing the processed prediction data frame.
#'
#' @export
#' @rdname plotPrediction
#' @importFrom dplyr filter
plotPrediction.prdlist <- function(prediction, ..., errfn = NULL, 
                                   scales = "free", 
                                   facet = c("wrap", "grid"), 
                                   transform = NULL) {
  
  facet <- match.arg(facet)
  
  prediction <- as.data.frame(prediction, errfn = errfn)
  prediction <- dplyr::filter(prediction, ...)
  
  if (!is.null(transform)) prediction <- coordTransform(prediction, transform)
  
  # --- Construct plot ---
  p <- ggplot(prediction, aes(x = time, y = value))
  
  if (facet == "wrap") {
    p <- p + 
      aes(group = condition, color = condition) +
      facet_wrap(~name, scales = scales)
  } else {
    p <- p + facet_grid(name ~ condition, scales = scales)
  }
  
  if (!is.null(errfn)) {
    p <- p + geom_ribbon(
      aes(ymin = value - sigma, ymax = value + sigma, fill = condition), 
      lty = 0, alpha = 0.3
    )
  }
  
  p <- p + geom_line()
  
  attr(p, "data") <- prediction
  p
}



## Methods for class prdframe ----------------------------
#' @export
#' @rdname plotCombined
plot.prdframe <- function(x, data = NULL, ..., scales = "free", facet = c("wrap", "grid", "wrap_plain"), transform = NULL) {
  
  prediction <- x
  
  prediction <- list("C1" = prediction)
  if (!is.null(data) && is.data.frame(data))
    data <- list("C1" = data)
  
  
  plotCombined.prdlist(prediction = prediction, data = data, ..., scales = scales, facet = facet, transform = transform)
  
}

#' @export
print.prdframe <- function(x, ...) {

  d1 <- attr(x, "deriv")
  d2 <- attr(x, "deriv2")

  derivs <- if (!is.null(d1)) {
    sprintf("yes [%s]", paste(dim(d1), collapse = " x "))
  } else "no"
  derivs2 <- if (!is.null(d2)) {
    sprintf("yes [%s]", paste(dim(d2), collapse = " x "))
  } else "no"

  attr(x, "deriv")      <- NULL
  attr(x, "deriv2")     <- NULL
  attr(x, "parameters") <- NULL

  print(unclass(x))
  cat("\n")
  cat("The prediction contains 1st-order derivatives: ", derivs,  "\n", sep = "")
  cat("The prediction contains 2nd-order derivatives: ", derivs2, "\n", sep = "")

}


## Methods for class prdfn ----------------------------------

#' @export
print.prdfn <- function(x, ...) {
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Prediction function:\n")
  str(args(x))
  cat("\n")
  cat("... conditions:", paste0(conditions, collapse = ", "), "\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
 
}

#' @export
summary.prdfn <- function(object,...) {
  
  x <- object
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Details:\n")
  if (!inherits(x, "composed")) {
    
    output <- lapply(1:length(mappings), function(C) {
      
      list(
        equations = attr(mappings[[C]], "equations"),
        events = attr(mappings[[C]], "events"),
        forcings = attr(mappings[[C]], "forcings"),
        parameters = attr(mappings[[C]], "parameters")
      )
      
    })
    names(output) <- conditions
    
    #print(output, ...)
    output
    
  } else {
    
    cat("\nObject is composed. See original objects for more details.\n")
    
  }
}

#' @export
print.obsfn <- function(x, ...) {
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Observation function:\n")
  str(args(x))
  cat("\n")
  cat("... conditions:", paste0(conditions, collapse = ", "), "\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
 
}

#' @export
summary.obsfn <- function(object, ...) {
  
  x <- object
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Details:\n")
  if (!inherits(x, "composed")) {
    
    output <- lapply(1:length(mappings), function(C) {
      
      list(
        equations = attr(mappings[[C]], "equations"),
        states = attr(mappings[[C]], "states"),
        parameters = attr(mappings[[C]], "parameters")
      )
      
    })
    names(output) <- conditions
    
    #print(output, ...)
    output
    
  } else {
    
    cat("\nObject is composed. See original objects for more details.\n")
    
  }
}


#' Model Predictions
#' 
#' Make a model prediction for times and a parameter frame. The
#' function is a generalization of the standard prediction by a
#' prediction function object in that it allows to pass a parameter
#' frame instead of a single parameter vector.
#' 
#' @param object prediction function
#' @param ... Further arguments goint to the prediction function
#' @param times numeric vector of time points
#' @param pars parameter frame, e.g. output from [mstrust] or 
#' [profile]
#' @param data data list object. If data is passed, its condition.grid
#' attribute is used to augment the output dataframe by additional 
#' columns. `"data"` itself is returned as an attribute.
#' @return A data frame
#' @export
predict.prdfn <- function(object, ..., times, pars, data = NULL) {
  
  
  x <- object
  arglist <- list(...)
  if (any(names(arglist) == "conditions")) {
    C <- arglist[["conditions"]]
    if (!is.null(data)) {
      data <- data[C]
    }
  }
  if (is.null(data)) data <- data.frame()
  condition.grid.data <- attr(data, "condition.grid")
  
  prediction <- do.call(combine, lapply(1:nrow(pars), function(i) {
    
    mypar <- as.parvec(pars, i)
    prediction <- x(times, mypar, deriv = FALSE, ...)
    
    if (is.null(names(prediction))) {
      conditions <- 1
    } else {
      conditions <- names(prediction)
    }
    
    condition.grid <- data.frame(row.names = conditions)
    
    # Augment by parframe metanames and obj.attributes
    mygrid <- pars[i, !colnames(pars) %in% attr(pars, "parameters")]
    mynames <- colnames(mygrid)
    if (length(mynames) > 0) {
      mynames <- paste0(".", mynames)
      colnames(mygrid) <- mynames
      condition.grid <- cbind(condition.grid, mygrid)
    }
    
    # Augment by condition.grid of data
    if (!is.null(condition.grid.data) && ncol(condition.grid.data) > 1) 
      condition.grid <- cbind(condition.grid.data[conditions,], condition.grid)
    
    # Write condition.grid into data
    attr(data, "condition.grid") <- condition.grid

    # Return
    as.data.frame(prediction, data = data)
    
  }))
  
  n <- nrow(prediction)
  
  if (length(data) > 0) {
    attr(data, "condition.grid") <- condition.grid.data
    data <- as.data.frame(data)
    tmp <- combine(prediction, data)
    data <- tmp[-(1:n),]
  }
  
  attr(prediction, "data") <- data
  return(prediction)  
  
  
  
}

