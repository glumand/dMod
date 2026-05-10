#' Adapter from a dMod prediction function to an `nlme`-compatible model
#'
#' Wraps a [prdfn] (typically `g * x * p`) so it can be passed to `nlme::nlme`.
#' The returned closure expects the per-observation arguments
#' `time, name, <parameters>, <covariates>`, evaluates the prdfn for each unique
#' subject (defined by the covariate columns), and returns the residual prediction
#' with the gradient against parameters attached as `"gradient"`.
#'
#' @param prdfn A `prdfn` (e.g. `g * x * p`) describing the deterministic part of
#'   the mixed-effects model.
#' @param covtable Optional `data.frame` whose column names label the subject
#'   covariates. Used only to read the names — values are passed through `...` of
#'   the returned closure.
#' @param cores Number of cores for `parallel::mclapply`. Forced to 1 on Windows.
#'
#' @return A function `model(time, name, ...)` that returns the prediction with a
#'   `"gradient"` attribute.
#' @export
modelNLME <- function(prdfn, covtable = NULL, cores = 1) {


  covnames <- names(covtable)
  parnames <- getParameters(prdfn)


  model <- function(time, name, ...) {

    pars <- as.data.frame(c(list(time, name), list(...)))

    names(pars) <- c("time", "name", parnames, covnames)

    id <- cumsum(Reduce("|", lapply(pars[-(1:2)], function(x) !duplicated(x))))
    pars <- split(pars, id)

    output <- parallel::mclapply(pars, function(sub) {
      timesD <- unique(sub$time)
      parsD <- unlist(sub[1, parnames])
      condition <- paste(unlist(sub[1, covnames]), collapse = "_")

      prediction <- prdfn(timesD, parsD, conditions = condition)[[1]]
      template <- data.frame(name = sub$name, time = sub$time, value = 0, sigma = 1, lloq = -Inf)

      myres <- res(template, prediction)

      output <- myres$prediction
      deriv <- as.matrix(attr(myres, "deriv")[, -(1:2)])

      list(output, deriv)
    }, mc.cores = cores)

    gradient <- do.call(rbind, lapply(output, function(x) x[[2]]))
    output <- unlist(lapply(output, function(x) x[[1]]))

    attr(output, "gradient") <- gradient

    return(output)

  }

  return(model)

}


#' Adapter from a dMod prediction function to a `saemix`-compatible model
#'
#' Wraps a [prdfn] for use as the `model` argument of `saemix::saemixModel`.
#' The returned closure has the saemix signature `function(psi, id, xidep)`,
#' splitting subjects by `id` and evaluating the prdfn per subject.
#'
#' @param prdfn A `prdfn` (e.g. `g * x * p`).
#' @param cores Number of cores for `parallel::mclapply`. Forced to 1 on Windows.
#'
#' @return A function `model(psi, id, xidep)` returning a numeric vector of
#'   predictions, in `id`-major order.
#' @export
modelSAEMIX <- function(prdfn, cores = 1) {

  parnames <- getParameters(prdfn)

  model <- function(psi, id, xidep) {

    pars <- split(as.data.frame(psi[id, ]), id)

    output <- do.call(c, parallel::mclapply(1:nrow(psi), function(i) {

      parsD <- unlist(psi[i,])
      names(parsD) <- parnames
      timesD <- as.numeric(xidep[id == i, 1])
      namesD <- as.character(xidep[id == i, 2])

      prediction <- prdfn(timesD, parsD, deriv = FALSE)[[1]]
      template <- data.frame(name = namesD, time = timesD, value = 0, sigma = 1)

      myres <- res(template, prediction)

      return(myres$prediction)


    }, mc.cores = cores))

    return(output)


  }

  return(model)


}
