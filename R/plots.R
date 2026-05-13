
# Pharmacometric abbreviations used in the NLME plot helpers below:
#   DV    observed value
#   IPRED individual prediction (per subject, with eta_i*)
#   PRED  population prediction (eta = 0)
#   IRES  individual residual = DV - IPRED
#   IWRES individual weighted residual = (DV - IPRED) / sigma

# ggplot2 / dplyr NSE column references; declared so R CMD check does not
# flag them as undefined globals.
utils::globalVariables(c("IPRED", "PRED", "predicted", "observed",
                         "sd_est", "iter", "level"))


# Custom interface to ggplot2 ---

#' Open last plot in external pdf viewer
#' 
#' @description Convenience function to show last plot in an external viewer.
#' @param plot `ggplot2` plot object.
#' @param command character, indicatig which pdf viewer is started.
#' @param ... arguments going to `ggsave`.
#' @export
ggopen <- function(plot = last_plot(), command = "xdg-open", ...) {
  filename <- tempfile(pattern = "Rplot", fileext = ".pdf")
  ggsave(filename = filename, plot = plot, ...)
  system(command = paste(command, filename))
}




#' Standard plotting theme of dMod
#' 
#' @param base_size numeric, font-size
#' @param base_family character, font-name
#' @export
theme_dMod <- function(base_size = 12, base_family = "") {
  colors <- list(
    medium = c(gray = '#737373', red = '#F15A60', green = '#7AC36A', blue = '#5A9BD4', orange = '#FAA75B', purple = '#9E67AB', maroon = '#CE7058', magenta = '#D77FB4'),
    dark = c(black = '#010202', red = '#EE2E2F', green = '#008C48', blue = '#185AA9', orange = '#F47D23', purple = '#662C91', maroon = '#A21D21', magenta = '#B43894'),
    light = c(gray = '#CCCCCC', red = '#F2AFAD', green = '#D9E4AA', blue = '#B8D2EC', orange = '#F3D1B0', purple = '#D5B2D4', maroon = '#DDB9A9', magenta = '#EBC0DA')
  )
  gray <- colors$medium["gray"]
  black <- colors$dark["black"]

  theme_bw(base_size = base_size, base_family = base_family) +
    theme(line = element_line(colour = "black"),
          rect = element_rect(fill = "white", colour = NA),
          text = element_text(colour = "black"),
          axis.text = element_text(size = rel(1.0), colour = "black"),
          axis.text.x = element_text(margin = margin(t = 4, r = 4, b = 0, l = 4, unit = "mm")),
          axis.text.y = element_text(margin = margin(t = 4, r = 4, b = 4, l = 0, unit = "mm")),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(-2, "mm"),
          legend.key = element_rect(colour = NA),
          panel.border = element_rect(colour = "black"),
          # panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(size = rel(1.0)))
  
}

dMod_colors <- c("#000000", "#C5000B", "#0084D1", "#579D1C", "#FF950E", "#4B1F6F", "#CC79A7","#006400", "#F0E442", "#8B4513", rep("gray", 100))

#' Standard dMod color palette
#'
#' @param ... arguments going to \code{scale_color_manual()}
#' @export
#' @examples
#' library(ggplot2)
#' times <- seq(0, 2*pi, 0.1)
#' values <- sin(times)
#' data <- data.frame(
#'    time = times, 
#'    value = c(values, 1.2*values, 1.4*values, 1.6*values), 
#'    group = rep(c("C1", "C2", "C3", "C4"), each = length(times))
#' )
#' qplot(time, value, data = data, color = group, geom = "line") + 
#'    theme_dMod() + scale_color_dMod()
#' @export
scale_color_dMod <- function(...) {
  scale_color_manual(..., values = dMod_colors)
}


#' Standard dMod color scheme
#'
#' @export
#' @param ... arguments going to \code{scale_color_manual()}
scale_fill_dMod <- function(...) {
  scale_fill_manual(..., values = dMod_colors)
}

ggplot <- function(...) ggplot2::ggplot(...) + scale_color_dMod() + theme_dMod()




# Other ---------------------------------------------

#' Coordinate transformation for data frames
#' 
#' Applies a symbolically defined transformation to the `value`
#' column of a data frame. Additionally, if a `sigma` column is
#' present, those values are transformed according to Gaussian error
#' propagation.
#' @param data data frame with at least columns "name" (character) and
#' "value" (numeric). Can optionally contain a column "sigma" (numeric).
#' @param transformations character (the transformation) or named list of
#' characters. In this case, the list names must be a subset of those 
#' contained in the "name" column.
#' @return The data frame with the transformed values and sigma uncertainties.
#' @export
#' 
#' @examples
#' mydata1 <- data.frame(name = c("A", "B"), time = 0:5, value = 0:5, sigma = .1)
#' coordTransform(mydata1, "log(value)")
#' coordTransform(mydata1, list(A = "exp(value)", B = "sqrt(value)"))
coordTransform <- function(data, transformations) {
  
  mynames <- unique(as.character(data$name))
  
  # Replicate transformation if not a list
  if (!is.list(transformations))
    transformations <- as.list(structure(rep(transformations, length(mynames)), names = mynames))
  
  out <- do.call(rbind, lapply(mynames, function(n) {
    
    subdata <- subset(data, name == n)
    
    if (n %in% names(transformations)) {
      
      mysymbol <- getSymbols(transformations[[n]])[1]
      mytrafo <- replaceSymbols(mysymbol, "value", transformations[[n]])
      mytrafo <- parse(text = mytrafo)
      
      if ("sigma" %in% colnames(subdata))
        subdata$sigma <- abs(with(subdata, eval(D(mytrafo, "value")))) * subdata$sigma
      subdata$value <- with(subdata, eval(mytrafo))
      
    }
    
    return(subdata)
    
  }))
  
  
  return(out)
  
  
}


# Method dispatch for plotX functions -------------



#' Plot a list of model predictions
#' 
#' @param prediction Named list of matrices or data.frames, usually the output of a prediction function
#' as generated by [Xs].
#' @param ... Further arguments going to `dplyr::filter`. 
#' @param scales The scales argument of `facet_wrap` or `facet_grid`, i.e. `"free"`, `"fixed"`, 
#' `"free_x"` or `"free_y"`
#' @param facet Either `"wrap"` or `"grid"`
#' @param transform list of transformation for the states, see [coordTransform].
#' @details The data.frame being plotted has columns `time`, `value`, `name` and `condition`.
#'  
#' 
#' @return A plot object of class `ggplot`.
#' @import ggplot2
#' @example inst/examples/plotting.R
#' @export
plotPrediction <- function(prediction,...) {
  UseMethod("plotPrediction", prediction)
}




#' Plot a list of model predictions and a list of data points in a combined plot
#' 
#' @param prediction Named list of matrices or data.frames, usually the output of a prediction function
#' as generated by [Xs].
#' @param data Named list of data.frames as being used in [res], i.e. with columns `name`, `time`, 
#' `value` and `sigma`.
#' @param ... Further arguments going to `dplyr::filter`. 
#' @param scales The scales argument of `facet_wrap` or `facet_grid`, i.e. `"free"`, `"fixed"`, 
#' `"free_x"` or `"free_y"`
#' @param facet `"wrap"` or `"grid"`. Try `"wrap_plain"` for high amounts of conditions and low amounts of observables.
#' @param transform list of transformation for the states, see [coordTransform].
#' @param aesthetics Named list of aesthetic mappings, specified as character, e.g. `list(linetype = "name")`. 
#' Can refer to variables in the condition.grid
#' @details The data.frame being plotted has columns `time`, `value`, `sigma`,
#' `name` and `condition`.
#'  
#' 
#' @return A plot object of class `ggplot`.
#' @example inst/examples/plotting.R
#' @importFrom graphics par
#' @export
plotCombined <- function(prediction,...) {
  UseMethod("plotCombined", prediction)
}


#' Plot a list data points
#' 
#' @param data Named list of data.frames as being used in [res], i.e. with columns `name`, `time`, 
#' `value` and `sigma`.
#' @param ... Further arguments going to `subset`. 
#' @param scales The scales argument of `facet_wrap` or `facet_grid`, i.e. `"free"`, `"fixed"`, 
#' `"free_x"` or `"free_y"`
#' @param facet Either `"wrap"` or `"grid"`
#' @param transform list of transformation for the states, see [coordTransform].
#' @details The data.frame being plotted has columns `time`, `value`, `sigma`,
#' `name` and `condition`.
#'  
#' 
#' @return A plot object of class `ggplot`.
#' @example inst/examples/plotting.R
#' @export
plotData  <- function(data,...) {
  UseMethod("plotData", data)
}

#' @export
#' @rdname plotData
plotData.data.frame <- function(data, ...) {
  plotData.datalist(as.datalist(data), ...)
}

#' Profile likelihood plot
#' 
#' @param profs Lists of profiles as being returned by [profile].
#' @param ... logical going to subset before plotting.
#' @param maxvalue Numeric, the value where profiles are cut off.
#' @param parlist Matrix or data.frame with columns for the parameters to be added to the plot as points.
#' If a "value" column is contained, deltas are calculated with respect to lowest chisquare of profiles.
#' @param ncol Number of columns in the resulting plot grid.
#' @return A plot object of class `ggplot`.
#' @details See [profile] for examples.
#' @export
plotProfile <- function(profs,...) {
  UseMethod("plotProfile", profs)
}


#' Profile likelihood: plot of the parameter paths.
#' 
#' @param profs profile or list of profiles as being returned by [profile]
#' @param ... arguments going to subset
#' @param whichPar Character or index vector, indicating the parameters that are taken as possible reference (x-axis)
#' @param sort Logical. If paths from different parameter profiles are plotted together, possible
#' combinations are either sorted or all combinations are taken as they are.
#' @param relative logical indicating whether the origin should be shifted.
#' @param scales character, either `"free"` or `"fixed"`.
#' @return A plot object of class `ggplot`.
#' @details See [profile] for examples.
#' @export
plotPaths <- function(profs, ..., whichPar = NULL, sort = FALSE, relative = TRUE, scales = "fixed") {
  
  if ("parframe" %in% class(profs)) 
    arglist <- list(profs)
  else
    arglist <- as.list(profs)
  
  
  if (is.null(names(arglist))) {
    profnames <- 1:length(arglist)
  } else {
    profnames <- names(arglist)
  }
  
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    # choose a proflist
    proflist <- as.data.frame(arglist[[i]])
    parameters <- attr(arglist[[i]], "parameters")
    
    if (is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }
    
    if (is.null(whichPar)) whichPar <- names(proflist)
    if (is.numeric(whichPar)) whichPar <- names(proflist)[whichPar]
    
    subdata <- do.call(rbind, lapply(whichPar, function(n) {
      # matirx
      paths <- as.matrix(proflist[[n]][, parameters])
      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      if (relative) 
        for(j in 1:ncol(paths)) paths[, j] <- as.numeric(paths[, j]) - as.numeric(paths[origin, j])
      
      combinations <- expand.grid.alt(whichPar, colnames(paths))
      if (sort) combinations <- apply(combinations, 1, sort) else combinations <- apply(combinations, 1, identity)
      combinations <- submatrix(combinations, cols = -which(combinations[1,] == combinations[2,]))
      combinations <- submatrix(combinations, cols = !duplicated(paste(combinations[1,], combinations[2,])))
      
      
      
      
      path.data <- do.call(rbind, lapply(1:dim(combinations)[2], function(j) {
        data.frame(chisquare = values, 
                   name = n,
                   proflist = profnames[i],
                   combination = paste(combinations[,j], collapse = " - \n"),
                   x = paths[, combinations[1,j]],
                   y = paths[, combinations[2,j]])
      }))
      
      return(path.data)
      
    }))
    
    return(subdata)
    
  }))
  
  data$proflist <- as.factor(data$proflist)
  
  
  if (relative)
    axis.labels <- c(expression(paste(Delta, "parameter 1")), expression(paste(Delta, "parameter 2")))  
  else
    axis.labels <- c("parameter 1", "parameter 2")
  
  
  data <- droplevels(subset(data, ...))
  data$y <- as.numeric(data$y)
  data$x <- as.numeric(data$x)
  
  suppressMessages(
    p <- ggplot(data, aes(x = x, y = y, group = interaction(name, proflist), color = name, lty = proflist)) + 
      facet_wrap(~combination, scales = scales) + 
      geom_path() + #geom_point(aes=aes(size=1), alpha=1/3) +
      xlab(axis.labels[1]) + ylab(axis.labels[2]) +
      scale_linetype_discrete(name = "profile\nlist") +
      scale_color_manual(name = "profiled\nparameter", values = dMod_colors)
  )
  
  attr(p, "data") <- data
  return(p)
  
}



#' Plot Fluxes given a list of flux Equations
#'
#' @param pouter parameters
#' @param x The model prediction function `x(times, pouter, fixed, ...)`
#' @param fluxEquations list of chars containing expressions for the fluxes,
#' if names are given, they are shown in the legend. Easy to obtain via [subset.eqnlist], see Examples.
#' @param nameFlux character, name of the legend.
#' @param times Numeric vector of time points for the model prediction
#' @param ... Further arguments going to x, such as `fixed` or `conditions`
#'
#'
#' @return A plot object of class `ggplot`.
#' @examples
#' \dontrun{
#'
#' plotFluxes(bestfit, x, times, subset(f, "B"%in%Product)$rates, nameFlux = "B production")
#' }
#' @export
plotFluxes <- function(pouter, x, times, fluxEquations, nameFlux = "Fluxes:", ...){

  if (is.null(names(fluxEquations))) names(fluxEquations) <- fluxEquations

  flux <- funCpp(fluxEquations, convenient = FALSE)$func
  prediction.all <- x(times, pouter, deriv = FALSE, ...)
  names.prediction.all <- names(prediction.all)
  if (is.null(names.prediction.all)) names.prediction.all <- paste0("C", 1:length(prediction.all))

  out <- lapply(1:length(prediction.all), function(cond) {
    prediction <- prediction.all[[cond]]
    pinner <- attr(prediction, "parameters")
    pinner.matrix <- matrix(pinner, nrow = length(pinner), ncol = nrow(prediction),
                            dimnames = list(names(pinner), NULL))
    fluxes <- cbind(time = prediction[, "time"], flux(cbind(prediction, t(pinner.matrix))))
    return(fluxes)
  }); names(out) <- names.prediction.all
  out <- wide2long(out)

  cbPalette <- c("#999999", "#E69F00", "#F0E442", "#56B4E9", "#009E73", "#0072B2",
                 "#D55E00", "#CC79A7","#CC6666", "#9999CC", "#66CC99","red", "blue", "green","black")

  P <- ggplot(out, aes(x = time, y = value, group = name, fill = name, log = "y")) +
    facet_wrap(~condition) + scale_fill_manual(values = cbPalette, name = nameFlux) +
    geom_density(stat = "identity", position = "stack", alpha = 0.3, color = "darkgrey", linewidth = 0.4) +
    xlab("time") + ylab("flux contribution")

  attr(P, "out") <- out

  return(P)

}


stepDetect <- function(x, tol) {
  
  jumps <- 1
  while (TRUE) {
    i <- which(x - x[1] > tol)[1]
    if (is.na(i)) break
    jumps <- c(jumps, tail(jumps, 1) - 1 + i)
    x <- x[-seq(1, i - 1, 1)]
  }
  
  return(jumps)
  
  
}

#' Plotting objective values of a collection of fits
#' 
#' @param x data.frame with columns "value", "converged" and "iterations", e.g. 
#' a [parframe].
#' @param ... arguments for subsetting of x
#' @param tol maximal allowed difference between neighboring objective values
#' to be recognized as one.
#' @export
plotValues <- function(x,...) {
  UseMethod("plotValues", x)
}



#' Plot parameter values for a fitlist
#' 
#' @param x parameter frame as obtained by as.parframe(mstrust)
#' @param tol maximal allowed difference between neighboring objective values
#' to be recognized as one.
#' @param ... arguments for subsetting of x
#' @export
plotPars <- function(x,...) {
  UseMethod("plotPars", x)
}



#' Plot residuals for a fitlist
#'
#' @description
#' Creates a plot of residuals from model fits, with flexible options for
#' grouping and faceting. Residuals can be summarized across different
#' dimensions (time, condition, observable, fit index).
#'
#' @param parframe Object of class \code{parframe}, e.g. returned by \link{mstrust}.
#' @param x Prediction function returning named list of data.frames with names 
#'   matching \code{data}.
#' @param data A \code{datalist} object, i.e. named list of data.frames with 
#'   columns \code{name}, \code{time}, \code{value}, and \code{sigma}.
#' @param split Character vector specifying how to summarize and display residuals.
#'   \itemize{
#'     \item \code{split[1]}: Variable for x-axis
#'     \item \code{split[2]}: Variable for grouping (color/line), defaults to \code{split[1]}
#'     \item \code{split[3+]}: Additional variables for \code{facet_wrap()}
#'   }
#' @param errmodel Optional error model function of type \code{prdfn}. If provided,
#'   residuals include the log-likelihood contribution from sigma.
#' @param ... Additional arguments passed to the prediction function \code{x}.
#'
#' @return A \code{ggplot} object with the summarized residual data frame
#'   attached as attribute \code{"out"}.
#'
#' @examples
#' \dontrun{
#' # Time on x-axis, faceted by condition and name
#' plotResiduals(myfitlist, g * x * p, data, 
#'               c("time", "index", "condition", "name"), 
#'               conditions = myconditions[1:4])
#'
#' # Condition on x-axis, residuals summed over time
#' plotResiduals(myfitlist, g * x * p, data, c("condition", "name", "index"))
#' }
#'
#' @export
#' @importFrom dplyr group_by summarise across
#' @importFrom rlang data_sym syms
plotResiduals <- function(parframe, x, data, split = "condition", errmodel = NULL, ...) {

  # Internal dispatch: NLME nlmeFit objects get their own residual diagnostic
  # (IWRES vs IPRED + vs TIME, see plotResiduals.nlmeFit below). The classical
  # parframe/x/data path below is preserved unchanged for back-compat.
  if (inherits(parframe, "nlmeFit")) return(plotResiduals.nlmeFit(parframe))

  timesD <- sort(unique(c(0, unlist(lapply(data, function(d) d$time)))))
  
  if (!("index" %in% colnames(parframe))) {
    parframe$index <- seq_len(nrow(parframe))
  }
  
  # --- Compute residuals for all fits and conditions ---
  out <- do.call(rbind, lapply(seq_len(nrow(parframe)), function(j) {
    pred <- x(timesD, as.parvec(parframe, j), deriv = FALSE, ...)
    
    out_con <- do.call(rbind, lapply(names(pred), function(con) {
      err <- NULL
      if (!is.null(errmodel)) {
        err <- errmodel(out = pred[[con]], pars = getParameters(pred[[con]]), conditions = con)
      }
      out <- res(data[[con]], pred[[con]], err[[con]])
      cbind(out, condition = con)
    }))
    
    cbind(index = as.character(parframe[j, "index"]), out_con)
  }))
  
  # --- Summarize residuals ---
  out <- dplyr::group_by(out, across(all_of(split)))
  
  if (!is.null(errmodel)) {
    out <- dplyr::summarise(out, res = sum(weighted.residual^2 + log(sigma^2)), .groups = "drop")
  } else {
    out <- dplyr::summarise(out, res = sum(weighted.residual^2), .groups = "drop")
  }
  
  out <- as.data.frame(out)
  
  # --- Build aesthetics ---
  groupvar <- if (length(split) > 1) split[2] else split[1]
  
  p <- ggplot(out, aes(x = !!rlang::data_sym(split[1]), 
                       y = res, 
                       color = !!rlang::data_sym(groupvar), 
                       group = !!rlang::data_sym(groupvar))) + 
    theme_dMod() + 
    geom_point() + 
    geom_line()
  
  if (length(split) > 2) {
    facet_vars <- rlang::syms(split[3:length(split)])
    p <- p + facet_wrap(vars(!!!facet_vars))
  }
  
  attr(p, "out") <- out
  p
}


# NLME plot helpers ---------------------------------------------------------

#' Predictions from an nlmeFit object
#'
#' @description
#' Returns a long-format `data.frame` of observed values vs population (`PRED`,
#' eta = 0) and individual (`IPRED`, eta at posterior modes) predictions per
#' condition, plus residuals. Used directly by the diagnostic plot helpers.
#'
#' @param object An [nlmeFit] (from [nlmeFit]).
#' @param times Optional numeric vector of additional times for the smooth
#'   IPRED/PRED curves. Defaults to the union of observed times.
#' @param ... Ignored.
#'
#' @return A long-format `data.frame` with columns
#'   `condition, time, name, observed, sigma, IPRED, PRED, IRES, IWRES, PRES,
#'   PWRES, source` where `source = "obs"` for rows aligned with data
#'   observations and `source = "grid"` for the dense IPRED/PRED smooth grid
#'   (observed/sigma NA for grid rows).
#' @export
predict.nlmeFit <- function(object, times = NULL, ...) {
  if (is.null(object$prdfn) || is.null(object$data) || is.null(object$omega))
    stop("predict.nlmeFit: fit is missing `prdfn`, `data`, or `omega` ",
         "(was the fit built by nlmeFit()?).")

  prdfn <- object$prdfn
  data  <- object$data
  om    <- object$omega
  etaModes <- object$etaModes
  pars      <- object$argument

  subjects <- rownames(om$subjectEtas)
  N <- length(subjects)
  obs_times <- sort(unique(unlist(lapply(data, `[[`, "time"))))
  grid_times <- if (is.null(times)) sort(unique(c(0, obs_times)))
                else sort(unique(c(0, obs_times, times)))

  eta_full <- as.vector(etaModes)
  names(eta_full) <- as.vector(om$subjectEtas)
  pars_ipred <- c(pars, eta_full)
  pars_pred  <- c(pars, setNames(rep(0, length(eta_full)), names(eta_full)))

  pred_ipred <- prdfn(times = grid_times, pars = pars_ipred, conditions = subjects)
  pred_pred  <- prdfn(times = grid_times, pars = pars_pred,  conditions = subjects)

  # Restrict to observables that actually appear in the data. Without this
  # filter, internal model states emitted by the prdfn (e.g. the depot Ag in a
  # PK model, present because the observation Y() is built with
  # attach.input = TRUE) would be plotted as ghost curves on the IPRED/PRED
  # axis next to the real observable, swamping the y-axis with the dose.
  data_names <- unique(unlist(lapply(data, `[[`, "name")))
  obs_names <- intersect(setdiff(colnames(pred_ipred[[1]]), "time"),
                         data_names)

  # If an errfn is attached and the data table has NA sigmas, evaluate
  # the errfn per-condition on the prediction so IWRES diagnostics get a
  # meaningful weight. Mirrors evalConditionResidual's contract.
  err <- object$errfn
  sigma_ipred <- NULL
  if (!is.null(err)) {
    sigma_ipred <- setNames(lapply(subjects, function(s) {
      pinner <- getParameters(pred_ipred[[s]])
      err(out = pred_ipred[[s]], pars = pinner, conditions = s)[[s]]
    }), subjects)
  }

  rows <- list()
  for (s in subjects) {
    d_s <- data[[s]]
    pi  <- pred_ipred[[s]]
    pp  <- pred_pred[[s]]
    si  <- if (!is.null(sigma_ipred)) sigma_ipred[[s]] else NULL
    # Coerce prdframe class away so plain matrix indexing works.
    if (!is.null(si)) si <- unclass(si)
    for (nm in obs_names) {
      d_nm <- d_s[d_s$name == nm, , drop = FALSE]
      get_sigma <- function(row_idx) {
        if (!is.null(si) && nm %in% colnames(si) && all(row_idx <= nrow(si)))
          as.numeric(si[row_idx, nm])
        else rep(NA_real_, length(row_idx))
      }
      if (nrow(d_nm)) {
        idx <- match(d_nm$time, pi[, "time"])
        sig_obs <- if (!is.null(si)) get_sigma(idx)
                   else if (all(is.na(d_nm$sigma))) rep(NA_real_, nrow(d_nm))
                   else d_nm$sigma
        rows[[length(rows) + 1L]] <- data.frame(
          condition = s, time = d_nm$time, name = nm,
          observed  = d_nm$value, sigma = sig_obs,
          IPRED     = pi[idx, nm], PRED  = pp[idx, nm],
          source    = "obs", stringsAsFactors = FALSE
        )
      }
      grid_idx <- setdiff(seq_len(nrow(pi)), match(d_nm$time, pi[, "time"]))
      if (length(grid_idx)) {
        rows[[length(rows) + 1L]] <- data.frame(
          condition = s, time = pi[grid_idx, "time"], name = nm,
          observed  = NA_real_, sigma = get_sigma(grid_idx),
          IPRED     = pi[grid_idx, nm], PRED = pp[grid_idx, nm],
          source    = "grid", stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$IRES  <- out$observed - out$IPRED
  out$PRES  <- out$observed - out$PRED
  out$IWRES <- out$IRES / out$sigma
  out$PWRES <- out$PRES / out$sigma
  rownames(out) <- NULL
  out
}



#' Per-subject individual fits (spaghetti plot)
#'
#' @description Faceted plot with one panel per subject: observed dots, IPRED
#'   curve, and (optionally) the population PRED curve overlaid dashed.
#' @param x Object to plot.
#' @param ... Method-specific arguments.
#' @return A ggplot.
#' @export
plotIndivs <- function(x, ...) UseMethod("plotIndivs", x)


#' Per-subject individual fits for an nlmeFit
#'
#' @description Per-subject IPRED curve with an IPRED ± sigma ribbon derived
#'   from the fit's attached errfn, optional population PRED overlay dashed,
#'   and observed values as points. When the fit has a single observable the
#'   layout is `facet_wrap(~ condition)`; with multiple observables it switches
#'   to `facet_grid(name ~ condition)` (observables in rows, subjects in
#'   columns). Use `subjectsPerPage` to split very large cohorts into several
#'   plots.
#' @param x An [nlmeFit].
#' @param times Optional grid of additional times for the smooth IPRED/PRED.
#' @param ncol Facet column count for the single-observable
#'   `facet_wrap(~ condition)` layout (default 4). Ignored in
#'   `facet_grid(name ~ condition)` mode.
#' @param showPred Logical; overlay population PRED dashed (default TRUE).
#' @param showBand Logical; draw the IPRED ± sigma ribbon if the fit carries
#'   an errfn (default TRUE).
#' @param subjectsPerPage Optional integer. If set, subjects are split into
#'   pages of at most `subjectsPerPage` each and the function returns a list
#'   of ggplots (one per page, page index appended to the title). `NULL`
#'   (default) keeps the single-plot behaviour.
#' @param ... Ignored.
#' @return A ggplot, or a list of ggplots when `subjectsPerPage` is set.
#' @export
plotIndivs.nlmeFit <- function(x, times = NULL, ncol = 4L,
                              showPred = TRUE, showBand = TRUE,
                              subjectsPerPage = NULL, ...) {
  fit <- x
  obs_times <- sort(unique(unlist(lapply(fit$data, `[[`, "time"))))
  if (is.null(times))
    times <- seq(min(obs_times), max(obs_times), length.out = 200L)

  pf <- predict(fit, times = times)
  # Trim grid rows to the observed time range. predict.nlmeFit prepends t=0 to
  # the dense grid for the ODE solver, but if no subject has an observation at
  # that time the grid IPRED/PRED there can be far outside the data range
  # (e.g. log(Cc + eps) at Cc(0)=0 sits at log(eps)), which would dominate the
  # y-axis. Observation rows are kept verbatim.
  pf <- pf[pf$source == "obs" |
             (pf$time >= min(obs_times) & pf$time <= max(obs_times)),
           , drop = FALSE]

  n_obs_names <- length(unique(pf$name))
  cond_levels <- if (is.factor(pf$condition))
    levels(droplevels(pf$condition)) else unique(as.character(pf$condition))

  build_page <- function(page_pf, page_label = NULL) {
    obs <- page_pf[page_pf$source == "obs", , drop = FALSE]
    grd <- page_pf[order(page_pf$condition, page_pf$name, page_pf$time), ,
                   drop = FALSE]
    has_band <- showBand && any(is.finite(grd$sigma))

    p <- ggplot2::ggplot()
    if (has_band) {
      p <- p + ggplot2::geom_ribbon(
        data = grd,
        ggplot2::aes(x = time,
                     ymin = IPRED - sigma, ymax = IPRED + sigma),
        fill = dMod_colors[3], alpha = 0.2, linetype = 0)
    }
    p <- p +
      ggplot2::geom_line(data = grd,
                         ggplot2::aes(x = time, y = IPRED, color = "IPRED"),
                         linewidth = 0.7)
    if (showPred) {
      p <- p + ggplot2::geom_line(
        data = grd,
        ggplot2::aes(x = time, y = PRED, color = "PRED"),
        linewidth = 0.5, linetype = "dashed")
    }
    p <- p +
      ggplot2::geom_point(data = obs,
                          ggplot2::aes(x = time, y = observed),
                          size = 1.6, alpha = 0.85)

    if (n_obs_names > 1L) {
      p <- p + ggplot2::facet_grid(name ~ condition, scales = "free_y")
    } else {
      p <- p + ggplot2::facet_wrap(~ condition, ncol = ncol, scales = "free_y")
    }

    title <- sprintf("Individual fits (method = %s)", fit$method)
    if (!is.null(page_label))
      title <- paste0(title, " ", page_label)

    p +
      ggplot2::scale_color_manual(values = c(IPRED = dMod_colors[3],
                                             PRED  = dMod_colors[2])) +
      ggplot2::labs(x = "Time", y = "Value", color = NULL, title = title) +
      theme_dMod(base_size = 11)
  }

  if (is.null(subjectsPerPage))
    return(build_page(pf))

  subjectsPerPage <- as.integer(subjectsPerPage)
  if (length(subjectsPerPage) != 1L || is.na(subjectsPerPage) ||
      subjectsPerPage < 1L)
    stop("`subjectsPerPage` must be a positive integer or NULL.")

  pages <- split(cond_levels,
                 ceiling(seq_along(cond_levels) / subjectsPerPage))
  n_pages <- length(pages)
  lapply(seq_along(pages), function(i) {
    page_conds <- pages[[i]]
    sub <- pf[as.character(pf$condition) %in% page_conds, , drop = FALSE]
    if (is.factor(sub$condition))
      sub$condition <- factor(sub$condition, levels = page_conds)
    label <- sprintf("(page %d/%d)", i, n_pages)
    build_page(sub, page_label = label)
  })
}



#' Observed vs predicted scatter (DV vs IPRED and DV vs PRED)
#'
#' @description S3 plot method for [nlmeFit]. Two-panel scatter: DV vs IPRED
#'   on the left, DV vs PRED on the right, identity line shown dashed. With a
#'   multi-observable fit the observable becomes a faceting row
#'   (`name ~ panel`); otherwise one row with `facet_wrap(~ panel)`.
#' @param x An [nlmeFit].
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plot.nlmeFit <- function(x, ...) {
  fit <- x
  pf <- predict(fit)
  pf <- pf[pf$source == "obs", , drop = FALSE]
  long <- rbind(
    data.frame(condition = pf$condition, name = pf$name,
               observed = pf$observed, predicted = pf$IPRED,
               panel = "IPRED", stringsAsFactors = FALSE),
    data.frame(condition = pf$condition, name = pf$name,
               observed = pf$observed, predicted = pf$PRED,
               panel = "PRED", stringsAsFactors = FALSE))
  multi_obs <- length(unique(long$name)) > 1L
  p <- ggplot2::ggplot(long, ggplot2::aes(x = predicted, y = observed,
                                          color = condition)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey50") +
    ggplot2::geom_point(alpha = 0.75)
  if (multi_obs) {
    # Per-observable scales differ by orders of magnitude (e.g. log-cp vs
    # linear PCA); free scales are required and coord_equal is incompatible
    # with that.
    p <- p + ggplot2::facet_grid(name ~ panel, scales = "free")
  } else {
    rng <- range(c(long$observed, long$predicted), na.rm = TRUE)
    p <- p +
      ggplot2::coord_equal(xlim = rng, ylim = rng) +
      ggplot2::facet_wrap(~ panel)
  }
  p +
    scale_color_dMod() +
    ggplot2::labs(x = "Prediction", y = "Observed",
                  title = sprintf("Observed vs predicted (method = %s)", fit$method)) +
    theme_dMod(base_size = 11) +
    ggplot2::theme(legend.position = "none")
}



#' Weighted-residual diagnostics for an nlmeFit
#'
#' @description Two-panel scatter: IWRES vs IPRED and IWRES vs TIME with a
#'   loess smoother. Sourced from `plotResiduals(fit, ...)` when the first
#'   argument inherits from `nlmeFit` (internal type dispatch, see
#'   [plotResiduals]).
#' @param fit An [nlmeFit].
#' @param ... Ignored.
#' @return A ggplot.
#' @keywords internal
plotResiduals.nlmeFit <- function(fit, ...) {
  pf <- predict(fit)
  pf <- pf[pf$source == "obs", , drop = FALSE]
  long <- rbind(
    data.frame(x = pf$IPRED, y = pf$IWRES, name = pf$name,
               panel = "IWRES vs IPRED", stringsAsFactors = FALSE),
    data.frame(x = pf$time,  y = pf$IWRES, name = pf$name,
               panel = "IWRES vs TIME",  stringsAsFactors = FALSE))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = x, y = y, color = name)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_smooth(se = FALSE, method = "loess", formula = y ~ x,
                         color = dMod_colors[2], linewidth = 0.6)
  if (length(unique(long$name)) > 1L) {
    p <- p + ggplot2::facet_grid(name ~ panel, scales = "free_x")
  } else {
    p <- p + ggplot2::facet_wrap(~ panel, scales = "free_x")
  }
  p +
    scale_color_dMod() +
    ggplot2::labs(x = NULL, y = "IWRES",
                  title = "Weighted residual diagnostics") +
    theme_dMod(base_size = 11)
}



#' Random-effect distribution diagnostics
#'
#' @description Per-eta histogram against the estimated `N(0, Omega_kk)`
#'   density plus a QQ-plot against the estimated normal. Detects systematic
#'   shrinkage, bimodality, or distributional misfit. Generic to leave room
#'   for non-nlmeFit methods in the future.
#' @param x Object to plot.
#' @param ... Method-specific arguments.
#' @return A ggplot (or a list of ggplots if `cowplot` is unavailable).
#' @export
plotHistIndivs <- function(x, ...) UseMethod("plotHistIndivs", x)


#' @rdname plotHistIndivs
#' @param x An [nlmeFit].
#' @param ... Ignored.
#' @return A ggplot (or a list with `hist` and `qq` if cowplot is unavailable).
#' @export
plotHistIndivs.nlmeFit <- function(x, ...) {
  fit <- x
  etaModes <- fit$etaModes
  om <- fit$omega
  if (is.null(etaModes) || is.null(om))
    stop("plotHistIndivs.nlmeFit: fit has no etaModes or omega.")
  Omega <- if (!is.null(fit$Omega)) fit$Omega else {
    L <- om$buildL(fit$argument[om$cholPars])
    tcrossprod(L)
  }
  K <- ncol(etaModes)
  eta_long <- do.call(rbind, lapply(seq_len(K), function(k) {
    sd_k <- sqrt(Omega[k, k])
    data.frame(eta_name = om$eta[k],
               value    = etaModes[, k],
               sd_est   = sd_k,
               stringsAsFactors = FALSE)
  }))
  # Per-eta Gaussian curve on a wide grid so the full N(0, Omega_kk) shape is
  # visible; using stat_function with a single mean(sd_est) painted all panels
  # with the same density (wrong when omega varies across etas) and was
  # clipped to the data extent by facet_wrap(scales = "free"). geom_line on
  # the explicit grid both picks up the panel-specific SD and widens the
  # x-range to ~ +/- 3.5 sigma.
  dens_long <- do.call(rbind, lapply(seq_len(K), function(k) {
    sd_k  <- sqrt(Omega[k, k])
    xlim  <- max(abs(etaModes[, k]), 3.5 * sd_k)
    xseq  <- seq(-xlim, xlim, length.out = 200L)
    data.frame(eta_name = om$eta[k],
               x        = xseq,
               density  = dnorm(xseq, 0, sd_k),
               stringsAsFactors = FALSE)
  }))
  p_hist <- ggplot2::ggplot(eta_long, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 12, fill = "grey80", color = "grey40") +
    ggplot2::geom_line(data = dens_long,
                       ggplot2::aes(x = x, y = density),
                       inherit.aes = FALSE,
                       color = dMod_colors[2], linewidth = 0.8) +
    ggplot2::facet_wrap(~ eta_name, scales = "free") +
    ggplot2::labs(x = "eta value", y = "density",
                  title = "Eta distribution vs N(0, Omega_kk)") +
    theme_dMod(base_size = 11)
  p_qq <- ggplot2::ggplot(eta_long,
                          ggplot2::aes(sample = value / sd_est)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line(color = dMod_colors[2]) +
    ggplot2::facet_wrap(~ eta_name) +
    ggplot2::labs(x = "Theoretical quantile (N(0,1))",
                  y = "Standardised eta",
                  title = "Eta QQ vs N(0,1)") +
    theme_dMod(base_size = 11)
  if (requireNamespace("cowplot", quietly = TRUE))
    cowplot::plot_grid(p_hist, p_qq, ncol = 1)
  else list(hist = p_hist, qq = p_qq)
}



#' ECM convergence trace (OFV, |delta-psi|) per stage
#'
#' @description Four-panel trace of OFV, structural-parameter step
#'   `|delta psi|`, max softmax weight, and minimum effective node count
#'   across ECM iterations. Quadrature-method nlmeFit only.
#' @param x Object to plot.
#' @param ... Method-specific arguments.
#' @return A ggplot.
#' @export
plotTrace <- function(x, ...) UseMethod("plotTrace", x)


#' @rdname plotTrace
#' @param x An [nlmeFit] with `method = "quadrature"`. Errors otherwise.
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plotTrace.nlmeFit <- function(x, ...) {
  fit <- x
  if (!fit$method %in% c("quadrature", "foceiQuadrature") ||
      is.null(fit$stageTrace))
    stop("plotTrace requires an nlmeFit fit with method 'quadrature' or 'foceiQuadrature'.")
  tr <- fit$stageTrace
  tr$iter <- seq_len(nrow(tr))
  long <- rbind(
    data.frame(iter = tr$iter, value = tr$OFV,        panel = "OFV",
               level = factor(tr$level)),
    data.frame(iter = tr$iter, value = tr$deltaPsi,  panel = "|delta psi|",
               level = factor(tr$level)),
    data.frame(iter = tr$iter, value = tr$maxSoftmax, panel = "max softmax",
               level = factor(tr$level)),
    data.frame(iter = tr$iter, value = tr$nEffMin,   panel = "n_eff min",
               level = factor(tr$level)))
  ggplot2::ggplot(long, ggplot2::aes(x = iter, y = value, color = level)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    scale_color_dMod() +
    ggplot2::labs(x = "ECM iteration",
                  title = "ECM convergence trace") +
    theme_dMod(base_size = 11)
}
