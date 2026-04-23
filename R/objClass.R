## Methods for class "objfn" -----------------------------------------------



## Class "objlist" and its constructors ------------------------------------



#' Generate objective list from numeric vector
#' 
#' @param p Named numeric vector
#' @return list with entries value (\code{0}), 
#' gradient (\code{rep(0, length(p))}) and 
#' hessian (\code{matrix(0, length(p), length(p))}) of class \code{obj}.
#' @examples
#' p <- c(A = 1, B = 2)
#' as.objlist(p)
#' @export
as.objlist <- function(p) {
  
  objlist(value = 0,
          gradient = structure(rep(0, length(p)), names = names(p)),
          hessian = matrix(0, length(p), length(p), dimnames = list(names(p), names(p))))
  
}


#' Compute a differentiable box prior
#' 
#' @param p Named numeric, the parameter value
#' @param mu Named numeric, the prior values, means of boxes
#' @param sigma Named numeric, half box width
#' @param k Named numeric, shape of box; if 0 a quadratic prior is obtained, the higher k the more box shape, gradient at border of the box (-sigma, sigma) is equal to sigma*k
#' @param fixed Named numeric with fixed parameter values (contribute to the prior value but not to gradient and Hessian)
#' @return list with entries: value (numeric, the weighted residual sum of squares), 
#' gradient (numeric, gradient) and 
#' hessian (matrix of type numeric). Object of class \code{objlist}.
constraintExp2 <- function(p, mu, sigma = 1, k = 0.05, fixed=NULL) {
  
  
  ##
  ## This function need to be extended according to constraintL2()
  ## The parameters sigma and k need to be replaced by more
  ## meaningful parameters.
  ##
  
  kmin <- 1e-5
  
  ## Augment sigma if length = 1
  if(length(sigma) == 1) 
    sigma <- structure(rep(sigma, length(mu)), names = names(mu)) 
  ## Augment k if length = 1
  if(length(k) == 1) 
    k <- structure(rep(k, length(mu)), names = names(mu))
  
  k <- sapply(k, function(ki){
    if(ki < kmin){
      kmin
    } else ki
  })
  
  
  ## Extract contribution of fixed pars and delete names for calculation of gr and hs  
  par.fixed <- intersect(names(mu), names(fixed))
  sumOfFixed <- 0
  if(!is.null(par.fixed)) sumOfFixed <- sum(0.5*(exp(k[par.fixed]*((fixed[par.fixed] - mu[par.fixed])/sigma[par.fixed])^2)-1)/(exp(k[par.fixed])-1))
  
  
  par <- intersect(names(mu), names(p))
  t <- p[par]
  mu <- mu[par]
  s <- sigma[par]
  k <- k[par]
  
  # Compute prior value and derivatives 
  
  gr <- rep(0, length(t)); names(gr) <- names(t)
  hs <- matrix(0, length(t), length(t), dimnames = list(names(t), names(t)))
  
  val <- sum(0.5*(exp(k*((t-mu)/s)^2)-1)/(exp(k)-1)) + sumOfFixed
  gr <- (k*(t-mu)/(s^2)*exp(k*((t-mu)/s)^2)/(exp(k)-1))
  diag(hs)[par] <- k/(s*s)*exp(k*((t-mu)/s)^2)/(exp(k)-1)*(1+2*k*(t-mu)/(s^2))
  
  dP <- attr(p, "deriv")
  if(!is.null(dP)) {
    gr <- as.vector(gr%*%dP); names(gr) <- colnames(dP)
    hs <- t(dP)%*%hs%*%dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
  }
  
  objlist(value=val,gradient=gr,hessian=hs)
  
}


#' L2 norm between data and model prediction
#'
#' @description
#' Creates an objective function for parameter estimation based on the
#' (negative log-likelihood) L2 norm between observed data and model predictions.
#' The returned objective function can be used with optimizers such as
#' [mstrust] and supports aggregation over multiple experimental conditions.
#'
#' @param data Object of class [datalist].
#' @param x Object of class [prdfn].
#' @param errmodel Optional object of class [obsfn]. The error model may be
#'   defined only for a subset of conditions.
#' @param times Optional numeric vector of additional time points at which the
#'   prediction function is evaluated. If NULL, time points are taken from the
#'   data. Event times should be included here if the prediction model uses events.
#' @param attr.name Character string. The objective value is additionally returned
#'   as an attribute with this name.
#' @param use.bessel Logical. If TRUE and an error model is provided, applies a
#'   global Bessel correction to variance estimates to account for finite-sample
#'   bias. Defaults to TRUE if an error model is supplied, FALSE otherwise.
#' @param cores Integer. Number of CPU cores used for parallel evaluation over
#'   conditions. Must be >= 1. Parallelization is configured once when the
#'   objective function is created.
#'
#' @return
#' An object of class `objfn`, i.e. a function
#' \code{obj(pars, fixed, deriv, env)} returning an [objlist].
#'
#' @details
#' Objective functions can be combined using the \code{+} operator, see
#' [sumobjfn].
#'
#' The Bessel correction is applied globally across all conditions and is given by
#' \deqn{\sqrt{n / (n - p)}}
#' where \eqn{n} is the total number of data points and \eqn{p} is the number of
#' structural (non-error-model) parameters.
#'
#' Parallelization is performed over experimental conditions if
#' \code{cores > 1}. The number of cores is fixed when calling \code{normL2()}
#' and cannot be changed at evaluation time.
#'
#' @example inst/examples/normL2.R
#' @export
normL2 <- function(data, x, errmodel = NULL, times = NULL,
                   attr.name = "data", use.bessel = !is.null(errmodel), cores = 1L) {
  
  stopifnot(cores >= 1L)
  
  timesD <- sort(unique(c(0, unlist(lapply(data, `[[`, "time")), times)))
  
  x.cond <- names(attr(x, "mappings"))
  d.cond <- names(data)
  stopifnot(all(d.cond %in% x.cond))
  
  e.cond <- if (!is.null(errmodel)) names(attr(errmodel, "mappings")) else NULL
  conditions <- intersect(x.cond, d.cond)
  
  # Precompute Bessel correction
  bessel <- 1
  if (use.bessel && !is.null(errmodel)) {
    n <- sum(vapply(data, nrow, 0L))
    p.all <- union(getParameters(x), getParameters(errmodel))
    p.err <- setdiff(getSymbols(unlist(getEquations(errmodel))),
                     names(unlist(getEquations(errmodel))))
    bessel <- sqrt(n / (n - length(p.all) + length(p.err)))
  }
  
  # Force early binding
  force(errmodel); force(bessel); force(conditions); force(timesD)
  
  # Core evaluation function
  eval_condition <- function(cn, prediction, pars, deriv) {
    err_cn <- NULL
    if (!is.null(errmodel) && (is.null(e.cond) || cn %in% e.cond)) {
      
      pinner <- getParameters(prediction[[cn]])
      fixedinner <- pinner[attr(pinner, "fixed")]
      pinner  <- as.parvec(pinner[setdiff(names(pinner), names(fixed))])
      fixedinner <- as.parvec(fixedinner, deriv = FALSE)
      
      err_cn <- errmodel(out = prediction[[cn]], pars = pinner,
                         fixed = fixedinner, conditions = cn)[[cn]]
    }
    nll(res(data[[cn]], prediction[[cn]], err_cn), pars = pars, deriv = deriv,
        bessel.correction = bessel)
  }
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, env = NULL) {
    pars <- ..1
    if (is.null(env)) env <- new.env()
    
    prediction <- x(times = timesD, pars = pars, fixed = fixed,
                    deriv = deriv, conditions = conditions)
    
    out <- if (cores == 1L) {
      Reduce(`+`, lapply(conditions, eval_condition, prediction, pars, deriv))
    } else if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      Reduce(`+`, parallel::parLapply(cl, conditions, eval_condition, prediction, pars, deriv))
    } else {
      Reduce(`+`, parallel::mclapply(conditions, eval_condition, prediction, pars, deriv, mc.cores = cores))
    }
    
    attr(out, attr.name) <- out$value
    env$prediction <- prediction
    attr(out, "env") <- env
    out
  }
  
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- d.cond
  attr(myfn, "parameters") <- attr(x, "parameters")
  attr(myfn, "modelname") <- modelname(x, errmodel)
  myfn
}



#' Soft L2 constraint on parameters
#'
#' @param mu Named numeric vector of prior means.
#' @param sigma Named numeric or character vector. Character entries indicate
#'   log-scale sigma parameters to be estimated.
#' @param attr.name Character. Name of the attribute storing the constraint value.
#' @param condition Optional character vector of conditions.
#'
#' @details
#' Computes
#' \deqn{(p-\mu)^2 / \sigma^2}
#' or, if sigma is estimated,
#' \deqn{(p-\mu)^2 / \sigma^2 + 2\log(\sigma)},
#' with sigma internally transformed via \code{exp()}.
#'
#' @return Object of class \code{objfn}.
#' @export
constraintL2 <- function(mu, sigma = 1, attr.name = "prior", condition = NULL) {
  
  est <- is.character(sigma)
  if (length(sigma) == 1) sigma <- setNames(rep(sigma, length(mu)), names(mu))
  if (is.null(names(sigma))) names(sigma) <- names(mu)
  sigma <- sigma[names(mu)]
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = condition, env = NULL) {
    
    p <- list(...)[[match.fnargs(list(...), "pars")]]
    dP <- attr(p, "deriv", exact = TRUE)
    
    allp <- c(p, fixed)
    avail <- intersect(names(mu), names(allp))
    if (!length(avail))
      return(objlist(value = 0))
    
    pa <- allp[avail]
    sg <- if (est) exp(allp[sigma[avail]]) else sigma[avail]
    r <- pa - mu[avail]
    
    val <- sum(r^2 / sg^2) + est * sum(2 * log(sg))
    
    if (!deriv)
      return(objlist(value = val))
    
    gr <- setNames(numeric(length(p)), names(p))
    hs <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
    
    p1 <- intersect(avail, names(p))
    gr[p1] <- 2 * r[p1] / sg[p1]^2
    diag(hs)[p1] <- 2 / sg[p1]^2
    
    if (est) for (sp in intersect(unique(sigma[avail]), names(p))) {
      idx <- sigma[avail] == sp
      gr[sp] <- sum(-2 * r[idx]^2 / sg[idx]^2 + 2)
      hs[sp, sp] <- sum(4 * r[idx]^2 / sg[idx]^2)
      cm <- intersect(names(idx)[idx], p1)
      hs[cm, sp] <- hs[sp, cm] <- -4 * r[cm] / sg[cm]^2
    }
    
    if (!is.null(dP)) {
      gi <- gr
      gr <- drop(gi %*% dP); names(gr) <- colnames(dP)
      hs <- t(dP) %*% hs %*% dP
      dimnames(hs) <- list(colnames(dP), colnames(dP))
    }
    
    out <- objlist(value = val, gradient = gr, hessian = hs)
    attr(out, attr.name) <- out$value
    attr(out, "env") <- env
    out
  }
  
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- names(mu)
  myfn
}




#' L2 objective function for validation data point
#' 
#' @param name character, the name of the prediction, e.g. a state name.
#' @param time numeric, the time-point associated to the prediction
#' @param value character, the name of the parameter which contains the
#' prediction value.
#' @param sigma numeric, the uncertainty of the introduced test data point
#' @param attr.name character. The constraint value is additionally returned in an 
#' attributed with this name
#' @param condition character, the condition for which the prediction is made.
#' @return List of class \code{objlist}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}, \link{constraintL2}
#' @details Computes the constraint value 
#' \deqn{\left(\frac{x(t)-\mu}{\sigma}\right)^2}{(pred-p[names(mu)])^2/sigma^2}
#' and its derivatives with respect to p.
#' @examples
#' prediction <- list(a = matrix(c(0, 1), nrow = 1, dimnames = list(NULL, c("time", "A"))))
#' derivs <- matrix(c(0, 1, 0.1), nrow = 1, dimnames = list(NULL, c("time", "A.A", "A.k1")))
#' attr(prediction$a, "deriv") <- derivs
#' p0 <- c(A = 1, k1 = 2)
#' 
#' vali <- datapointL2(name = "A", time = 0, value = "newpoint", sigma = 1, condition = "a")
#' vali(pars = c(p0, newpoint = 1), env = .GlobalEnv)
#' @export
datapointL2 <- function(name, time, value, sigma = 1, attr.name = "validation", condition) {
  
  controls <- list(
    mu        = structure(name, names = value)[1], # only one data point is allowed
    time      = time[1],
    sigma     = sigma[1],
    attr.name = attr.name
  )
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL, env = NULL) {
    mu        <- controls$mu
    t         <- controls$time
    sigma     <- controls$sigma
    attr.name <- controls$attr.name
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter  <- arglist[[1]]
    if (is.null(env)) {
      stop("No prediction available. Use the argument env to pass an environment that contains the prediction.")
    }
    prediction <- as.list(env)$prediction
    
    if (!is.null(conditions) && !condition %in% conditions)
      return()
    if (is.null(conditions) && !condition %in% names(prediction))
      stop("datapointL2 requests unavailable condition. Call the objective function explicitly stating the conditions argument.")
    
    datapar <- setdiff(names(mu), names(fixed))
    parapar <- setdiff(names(pouter), c(datapar, names(fixed)))
    
    time.index <- which(prediction[[condition]][, "time"] == t)
    if (!length(time.index))
      stop("datapointL2() requests time point for which no prediction is available. Please add missing time point by the times argument in normL2()")
    withDeriv <- !is.null(attr(prediction[[condition]], "deriv"))
    
    pred  <- prediction[[condition]][time.index, ][mu]
    deriv <- NULL
    if (withDeriv) {
      dfull <- attr(prediction[[condition]], "deriv")
      if (length(dim(dfull)) == 3L) {
        # new format: [time x variable x parameter]
        avail_pars <- dimnames(dfull)[[3]]
        use_pars   <- intersect(parapar, avail_pars)
        if (length(use_pars)) {
          dtmp  <- dfull[time.index, mu, use_pars, drop = TRUE]
          deriv <- setNames(as.numeric(dtmp), use_pars)
        }
      } else {
        # fallback to old matrix format with "var.par" column names
        mu.para <- intersect(paste(mu, parapar, sep = "."), names(dfull))
        deriv   <- dfull[mu.para]
      }
    }
    
    res <- as.numeric(pred - c(fixed, pouter)[names(mu)])
    val <- (res / sigma)^2
    
    gr <- hs <- NULL
    if (withDeriv) {
      dres.dp <- setNames(numeric(length(pouter)), names(pouter))
      if (length(deriv))    dres.dp[names(deriv)] <- deriv
      if (length(datapar))  dres.dp[datapar] <- -1
      gr <- 2 * res * dres.dp / sigma^2
      hs <- 2 * outer(dres.dp, dres.dp, "*") / sigma^2
      colnames(hs) <- rownames(hs) <- names(pouter)
    }
    
    out <- objlist(value = val, gradient = gr, hessian = hs)
    attr(out, attr.name)   <- out$value
    attr(out, "prediction") <- pred
    attr(out, "env")       <- env
    class(out)             <- NULL
    out
  }
  class(myfn)             <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- value[1]
  myfn
}


#' L2 objective function for prior value
#' 
#' @description As a prior function, it returns derivatives with respect to
#' the penalty parameter in addition to parameter derivatives.
#' 
#' @param mu Named numeric, the prior values
#' @param lambda Character of length one. The name of the penalty paramter in \code{p}.
#' @param attr.name character. The constraint value is additionally returned in an 
#' attributed with this name
#' @param condition character, the condition for which the constraint should apply. If
#' \code{NULL}, applies to any condition.
#' @return List of class \code{objlist}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}
#' @details Computes the constraint value 
#' \deqn{e^{\lambda} \| p-\mu \|^2}{exp(lambda)*sum((p-mu)^2)}
#' and its derivatives with respect to p and lambda.
#' @examples
#' p <- c(A = 1, B = 2, C = 3, lambda = 0)
#' mu <- c(A = 0, B = 0)
#' obj <- priorL2(mu = mu, lambda = "lambda")
#' obj(pars = p + rnorm(length(p), 0, .1))
#' @export
priorL2 <- function(mu, lambda = "lambda", attr.name = "prior", condition = NULL) {
  
  
  controls <- list(mu = mu, lambda = lambda, attr.name = attr.name)
  
  myfn <- function(..., fixed = NULL, deriv=TRUE, conditions = condition, env = NULL) {
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    # Import from controls 
    mu <- controls$mu
    lambda <- controls$lambda
    attr.name <- controls$attr.name
    
    # pouter can be a list (if result from a parameter transformation)
    # In this case match with conditions and evaluate only those
    # If there is no overlap, return NULL
    # If pouter is not a list, evaluate the constraint function 
    # for this pouter.
    
    if (is.list(pouter) && !is.null(conditions)) {
      available <- intersect(names(pouter), conditions)
      defined <- ifelse(is.null(condition), TRUE, condition %in% conditions)
      
      if (length(available) == 0 | !defined) return()
      pouter <- pouter[intersect(available, condition)]
    }
    if (!is.list(pouter)) pouter <- list(pouter)
    
    outlist <- lapply(pouter, function(p) {
      
      
      ## Extract contribution of fixed pars and delete names for calculation of gr and hs  
      par.fixed <- intersect(names(mu), names(fixed))
      sumOfFixed <- 0
      if (!is.null(par.fixed)) sumOfFixed <- sum(exp(c(fixed, p)[lambda])*(fixed[par.fixed] - mu[par.fixed]) ^ 2)
      
      # Compute prior value and derivatives
      par <- intersect(names(mu), names(p))
      par0 <- setdiff(par, lambda)
      
      val <- sum(exp(c(fixed, p)[lambda]) * (p[par] - mu[par]) ^ 2) + sumOfFixed
      
      gr <- hs <- NULL
      if (deriv) {
        gr <- rep(0, length(p)); names(gr) <- names(p)
        gr[par] <- 2*exp(c(fixed, p)[lambda])*(p[par] - mu[par])
        if (lambda %in% names(p)) {
          gr[lambda] <- sum(exp(c(fixed, p)[lambda]) * (p[par0] - mu[par0]) ^ 2) + 
            sum(exp(c(fixed, p)[lambda]) * (fixed[par.fixed] - mu[par.fixed]) ^ 2)
        }
        
        hs <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
        diag(hs)[par] <- 2*exp(c(fixed, p)[lambda])
        if (lambda %in% names(p)) {
          hs[lambda, lambda] <- gr[lambda] 
          hs[lambda, par0] <- hs[par0, lambda] <- gr[par0]
        }
        
        dP <- attr(p, "deriv")
        if (!is.null(dP)) {
          gr <- as.vector(gr %*% dP); names(gr) <- colnames(dP)
          hs <- t(dP) %*% hs %*% dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
        }
      }
      
      objlist(value = val, gradient = gr, hessian = hs)
      
    })
    
    out <- Reduce("+", outlist)
    attr(out, controls$attr.name) <- out$value
    attr(out, "env") <- env
    
    return(out)
    
    
  }
  
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- names(mu)
  return(myfn)
  
}


#' Compute the negative log-likelihood
#' 
#' @description Gaussian Log-likelihood. Supports NONMEM-like BLOQ handling methods M1, M3 and M4 
#' and estimation of error models with optional Bessel correction for variance parameter bias.
#' The Hessian is approximated via the Jacobian of the residuals (Gauss-Newton approximation).
#' Supports different parameter sets per condition - gradients and Hessians are merged by parameter name.
#' 
#' @param nout data.frame (result of [res]) or object of class [res].
#' @param pars Named vector of ALL outer parameters (union across conditions)
#' @param deriv Logical. If TRUE, compute gradient and hessian
#' @param opt.BLOQ Character denoting the method to deal with BLOQ data. 
#' One of "M1", "M3", "M4NM", or "M4BEAL".
#' @param opt.hessian Named logical vector to include or exclude various 
#' summands of the hessian matrix.
#' @param bessel.correction Numeric. Bessel correction factor for variance estimation.
#' 
#' @md
#' @return list with entries value, gradient, and hessian (Gauss-Newton approximation).
#' @export
nll <- function(nout, pars, deriv, opt.BLOQ = "M3", opt.hessian = c(
  ALOQ_part1 = TRUE, ALOQ_part2 = TRUE, ALOQ_part3 = TRUE,
  BLOQ_part1 = TRUE, BLOQ_part2 = TRUE, BLOQ_part3 = TRUE,
  PD = TRUE), bessel.correction = 1) {
  
  is.bloq <- nout$bloq
  nout.bloq <- nout[is.bloq, , drop = FALSE]
  nout.aloq <- nout[!is.bloq, , drop = FALSE]
  
  derivs <- attr(nout, "deriv")
  derivs.bloq <- if (!is.null(derivs)) derivs[is.bloq, , drop = FALSE] else NULL
  derivs.aloq <- if (!is.null(derivs)) derivs[!is.bloq, , drop = FALSE] else NULL
  
  derivs.err <- attr(nout, "deriv.err")
  derivs.err.bloq <- if (!is.null(derivs.err)) derivs.err[is.bloq, , drop = FALSE] else NULL
  derivs.err.aloq <- if (!is.null(derivs.err)) derivs.err[!is.bloq, , drop = FALSE] else NULL
  
  n_pars <- length(pars)
  par_names <- names(pars)
  
  mywrss <- {
    gr <- if (deriv) setNames(numeric(n_pars), par_names) else NULL
    he <- if (deriv) matrix(0, n_pars, n_pars, dimnames = list(par_names, par_names)) else NULL
    objlist(value = 0, gradient = gr, hessian = he)
  }
  
  nll_ALOQ_result <- NULL
  if (!all(is.bloq)) {
    nll_ALOQ_result <- nll_ALOQ(nout.aloq, derivs.aloq, derivs.err.aloq,
                                par_names = par_names,
                                opt.BLOQ = opt.BLOQ, opt.hessian = opt.hessian,
                                bessel.correction = bessel.correction)
  }
  mywrss <- mywrss + nll_ALOQ_result
  
  if (any(is.bloq) && opt.BLOQ != "M1") {
    mywrss <- mywrss + nll_BLOQ(nout.bloq, derivs.bloq, derivs.err.bloq,
                                par_names = par_names,
                                opt.BLOQ = opt.BLOQ, opt.hessian = opt.hessian)
  }
  
  chisquare <- attr(nll_ALOQ_result, "chisquare")
  nll_val <- attr(nll_ALOQ_result, "nll")
  attr(mywrss, "chisquare") <- if (length(chisquare)) chisquare else 0
  attr(mywrss, "nll") <- if (length(nll_val)) nll_val else 0
  
  mywrss
}


#' Non-linear log likelihood for the ALOQ part of the data
#' 
#' @param nout output of [res()]
#' @param derivs,derivs.err matrix of first derivatives (may have subset of parameters)
#' @param par_names Character vector of ALL parameter names (full set)
#' @param opt.BLOQ Character denoting the method to deal with BLOQ data
#' @param opt.hessian Named logical vector for hessian components
#' @param bessel.correction Numeric. Bessel correction factor.
#' @md
#' @importFrom stats pnorm dnorm
nll_ALOQ <- function(nout, derivs, derivs.err,
                     par_names,
                     opt.BLOQ = c("M3", "M4NM", "M4BEAL", "M1"),
                     opt.hessian = c(ALOQ_part1 = TRUE, ALOQ_part2 = TRUE, ALOQ_part3 = TRUE),
                     bessel.correction = 1) {
  
  wr <- nout$weighted.residual
  w0 <- nout$weighted.0
  s  <- nout$sigma
  
  chisquare_ml <- sum(wr^2)
  neg2ll_ml <- chisquare_ml + sum(log(2 * pi * s^2))
  
  use_bessel <- bessel.correction != 1
  if (use_bessel) {
    wr <- wr * bessel.correction
    w0 <- w0 * bessel.correction
  }
  
  chisquare <- sum(wr^2)
  obj <- chisquare + sum(log(2 * pi * s^2))
  
  if (opt.BLOQ[1] == "M4BEAL") {
    bloq_term <- 2 * sum(stats::pnorm(w0, log.p = TRUE))
    obj <- obj + bloq_term
    neg2ll_ml <- neg2ll_ml + bloq_term
  }
  
  n_pars_full <- length(par_names)
  grad <- NULL
  hessian <- NULL
  
  if (!is.null(derivs) && nrow(derivs) > 0) {
    local_pars <- colnames(derivs)
    local_pars_err <- if (!is.null(derivs.err)) colnames(derivs.err) else character(0)
    n_local <- length(local_pars)
    n_data <- nrow(derivs)
    
    idx_map <- match(local_pars, par_names)
    
    dxdp <- derivs
    inv_s <- 1 / s
    
    # Build aligned dsdp matrix (same columns as dxdp)
    dsdp <- matrix(0, n_data, n_local)
    colnames(dsdp) <- local_pars
    if (length(local_pars_err) > 0) {
      common <- intersect(local_pars, local_pars_err)
      if (length(common) > 0) {
        dsdp[, common] <- derivs.err[, common, drop = FALSE]
      }
    }
    
    # Compute derivatives
    dwrdp <- inv_s * dxdp - (wr * inv_s) * dsdp
    dw0dp <- inv_s * dxdp - (w0 * inv_s) * dsdp
    dlogsdp <- inv_s * dsdp
    
    if (use_bessel) {
      dwrdp <- dwrdp * bessel.correction
      dw0dp <- dw0dp * bessel.correction
    }
    
    # Local gradient
    grad_local <- 2 * (colSums(wr * dwrdp) + colSums(dlogsdp))
    
    if (opt.BLOQ[1] == "M4BEAL") {
      G_by_Phi <- exp(stats::dnorm(w0, log = TRUE) - stats::pnorm(w0, log.p = TRUE))
      grad_local <- grad_local + 2 * colSums(G_by_Phi * dw0dp)
    }
    
    # Map to full gradient
    grad <- setNames(numeric(n_pars_full), par_names)
    grad[idx_map] <- grad_local
    
    # Local Hessian (Gauss-Newton)
    hessian_local <- 2 * crossprod(dwrdp)
    
    if (opt.hessian["ALOQ_part1"]) {
      tmp <- (-wr * inv_s^2) * dxdp
      hessian_local <- hessian_local + 2 * (crossprod(tmp, dsdp) + crossprod((-wr * inv_s^2) * dsdp, dxdp))
    }
    
    if (opt.hessian["ALOQ_part2"]) {
      hessian_local <- hessian_local + 4 * crossprod((wr * inv_s) * dsdp)
    }
    
    if (opt.hessian["ALOQ_part3"]) {
      hessian_local <- hessian_local - 2 * crossprod(dlogsdp)
    }
    
    if (opt.BLOQ[1] == "M4BEAL") {
      G_w0 <- exp(stats::dnorm(w0, log = TRUE) - stats::pnorm(w0, log.p = TRUE))
      coef <- pmax(0, -w0 * G_w0 - G_w0^2)
      hessian_local <- hessian_local + 2 * crossprod(sqrt(coef) * dw0dp)
      
      tmp_G <- G_w0 * (-inv_s^2)
      hessian_local <- hessian_local + 2 * (crossprod(tmp_G * dxdp, dsdp) + crossprod(tmp_G * dsdp, dxdp))
      
      if (opt.hessian["ALOQ_part1"]) {
        hessian_local <- hessian_local + 4 * crossprod(sqrt(pmax(0, G_w0 * w0) * inv_s) * dsdp)
      }
    }
    
    # Map to full Hessian
    hessian <- matrix(0, n_pars_full, n_pars_full, dimnames = list(par_names, par_names))
    hessian[idx_map, idx_map] <- hessian_local
  }
  
  out <- objlist(value = obj, gradient = grad, hessian = hessian)
  attr(out, "chisquare") <- chisquare_ml
  attr(out, "nll") <- neg2ll_ml
  attr(out, "besselcorrected") <- use_bessel
  out
}


#' Non-linear log likelihood for the BLOQ part of the data
#' @md
#' @param nout.bloq The bloq output of [res()]
#' @param derivs.bloq,derivs.err.bloq matrix of first derivatives
#' @param par_names Character vector of ALL parameter names (full set)
#' @param opt.BLOQ Character denoting the method to deal with BLOQ data
#' @param opt.hessian Named logical vector for hessian components
#' @importFrom stats pnorm dnorm
nll_BLOQ <- function(nout.bloq, derivs.bloq, derivs.err.bloq,
                     par_names,
                     opt.BLOQ = c("M3", "M4NM", "M4BEAL", "M1"),
                     opt.hessian = c(BLOQ_part1 = TRUE, BLOQ_part2 = TRUE, BLOQ_part3 = TRUE)) {
  
  if (opt.BLOQ[1] %in% c("M4NM", "M4BEAL") && any(nout.bloq$value < 0)) {
    stop("M4-Method cannot handle LLOQ < 0. Possible solutions:\n",
         "  * Use M3 which allows negative LLOQ (recommended)\n",
         "  * If you are working with log-transformed DV, exponentiate DV and LLOQ\n")
  }
  
  wr <- nout.bloq$weighted.residual
  w0 <- nout.bloq$weighted.0
  s  <- nout.bloq$sigma
  inv_s <- 1 / s
  
  n_pars_full <- length(par_names)
  
  if (opt.BLOQ[1] == "M3") {
    obj.bloq <- -2 * sum(stats::pnorm(-wr, log.p = TRUE))
  } else {
    objvals <- -2 * log(1 - stats::pnorm(wr) / stats::pnorm(w0))
    bad <- !is.finite(objvals)
    if (any(bad)) {
      diff_w <- w0[bad] - wr[bad]
      intercept <- ifelse(log(diff_w) > 0, 1.8, -1.9 * log(diff_w) + 0.9)
      lin <- ifelse(log(diff_w) > 0, 0.9, 0.5)
      objvals[bad] <- intercept + lin * w0[bad] + 0.95 * w0[bad]^2
    }
    obj.bloq <- sum(objvals)
  }
  
  grad.bloq <- NULL
  hessian.bloq <- NULL
  
  if (!is.null(derivs.bloq) && nrow(derivs.bloq) > 0) {
    local_pars <- colnames(derivs.bloq)
    local_pars_err <- if (!is.null(derivs.err.bloq)) colnames(derivs.err.bloq) else character(0)
    n_local <- length(local_pars)
    n_data <- nrow(derivs.bloq)
    
    idx_map <- match(local_pars, par_names)
    
    dxdp <- derivs.bloq
    
    # Build aligned dsdp
    dsdp <- matrix(0, n_data, n_local)
    colnames(dsdp) <- local_pars
    if (length(local_pars_err) > 0) {
      common <- intersect(local_pars, local_pars_err)
      if (length(common) > 0) {
        dsdp[, common] <- derivs.err.bloq[, common, drop = FALSE]
      }
    }
    
    dwrdp <- inv_s * dxdp - (wr * inv_s) * dsdp
    dw0dp <- inv_s * dxdp - (w0 * inv_s) * dsdp
    
    G_by_Phi <- function(w1, w2 = w1) {
      exp(stats::dnorm(w1, log = TRUE) - stats::pnorm(w2, log.p = TRUE))
    }
    
    if (opt.BLOQ[1] == "M3") {
      G_neg_wr <- G_by_Phi(-wr)
      grad_local <- 2 * colSums(G_neg_wr * dwrdp)
    } else {
      c1 <- 1 / (1/G_by_Phi(wr, w0) - 1/G_by_Phi(wr, wr))
      c2 <- 1 / (1/G_by_Phi(w0, w0) - 1/G_by_Phi(w0, wr))
      c3 <- G_by_Phi(w0)
      grad_local <- 2 * colSums(c1 * dwrdp - c2 * dw0dp + c3 * dw0dp)
    }
    
    grad.bloq <- setNames(numeric(n_pars_full), par_names)
    grad.bloq[idx_map] <- grad_local
    
    hessian_local <- matrix(0, n_local, n_local, dimnames = list(local_pars, local_pars))
    
    if (opt.BLOQ[1] == "M3") {
      G_neg_wr <- G_by_Phi(-wr)
      
      if (opt.hessian["BLOQ_part1"]) {
        coef <- -wr * G_neg_wr + G_neg_wr^2
        hessian_local <- hessian_local + 2 * crossprod(dwrdp, coef * dwrdp)
      }
      
      if (opt.hessian["BLOQ_part2"]) {
        tmp <- G_neg_wr * inv_s^2
        hessian_local <- hessian_local - 2 * (crossprod(tmp * dxdp, dsdp) + crossprod(tmp * dsdp, dxdp))
      }
      
      if (opt.hessian["BLOQ_part3"]) {
        hessian_local <- hessian_local - 2 * crossprod(dsdp, (G_neg_wr * 2 * (-wr) * inv_s^2) * dsdp)
      }
      
    } else {
      stable <- function(wn, w0, wr) {
        out <- stats::dnorm(wn) / (stats::pnorm(w0) - stats::pnorm(wr))
        if (identical(wn, w0)) { out[is.infinite(out)] <- 0; return(out) }
        if (identical(wn, wr)) { out[is.infinite(out)] <- 1/(w0 - wr) + wr; return(out) }
        out
      }
      
      A1 <- -wr * stable(wr, w0, wr)
      A2 <- stable(wr, w0, wr)
      A3 <- -w0 * stable(w0, w0, wr)
      A4 <- stable(w0, w0, wr)
      G_w0 <- G_by_Phi(w0)
      A5 <- -w0 * G_w0 - G_w0^2
      A6 <- G_w0
      
      if (opt.hessian["BLOQ_part1"]) {
        hessian_local <- hessian_local + 2 * (
          crossprod(dwrdp, A1 * dwrdp) +
            crossprod(dw0dp, A3 * dw0dp) +
            crossprod(dw0dp, A5 * dw0dp)
        )
      }
      
      if (opt.hessian["BLOQ_part2"]) {
        part2_vec <- A2 * dwrdp - A4 * dw0dp
        hessian_local <- hessian_local - 2 * crossprod(part2_vec)
        
        hessian_local <- hessian_local + 2 * (
          crossprod(A2 * (-inv_s^2) * dxdp, dsdp) + crossprod(A2 * (-inv_s^2) * dsdp, dxdp) +
            crossprod(A4 * (-inv_s^2) * dxdp, dsdp) + crossprod(A4 * (-inv_s^2) * dsdp, dxdp) +
            crossprod(A6 * (-inv_s^2) * dxdp, dsdp) + crossprod(A6 * (-inv_s^2) * dsdp, dxdp)
        )
      }
      
      if (opt.hessian["BLOQ_part3"]) {
        hessian_local <- hessian_local + 2 * (
          crossprod(dsdp, (A2 * 2 * wr * inv_s^2) * dsdp) +
            crossprod(dsdp, (A4 * 2 * w0 * inv_s^2) * dsdp) +
            crossprod(dsdp, (A6 * 2 * w0 * inv_s^2) * dsdp)
        )
      }
    }
    
    hessian.bloq <- matrix(0, n_pars_full, n_pars_full, dimnames = list(par_names, par_names))
    hessian.bloq[idx_map, idx_map] <- hessian_local
  }
  
  objlist(value = obj.bloq, gradient = grad.bloq, hessian = hessian.bloq)
}


## Methods for class objlist ------------------------------------------------

#' Add two lists element by element
#' 
#' @param out1 List of numerics or matrices
#' @param out2 List with the same structure as out1 (there will be no warning when mismatching)
#' @details If out1 has names, out2 is assumed to share these names. Each element of the list out1
#' is inspected. If it has a \code{names} attributed, it is used to do a matching between out1 and out2.
#' The same holds for the attributed \code{dimnames}. In all other cases, the "+" operator is applied
#' the corresponding elements of out1 and out2 as they are.
#' @return List of length of out1. 
#' @aliases sumobjlist
#' @export
#' 
"+.objlist" <- function(out1, out2) {
  
  if (is.null(out1)) return(out2)
  if (is.null(out2)) return(out1)
  
  what <- intersect(c("value", "gradient", "hessian"), c(names(out1), names(out2)))
  
  add_vector <- function(a,b) {
    # add vector b to a by names
    i <- intersect(names(a), names(b))
    a[i] <- a[i] + b[i]
    a}
  add_matrix <- function(a,b) {
    i <- intersect(rownames(a), rownames(b))
    a[i,i] <- a[i,i] + b[i,i]
    a}
  
  gn1 <- names(out1$gradient)
  gn2 <- names(out2$gradient)
  
  one_includes_two <- all(gn2 %in% gn1) 
  two_includes_one <- all(gn1 %in% gn2)
  neither_included <- !(one_includes_two | two_includes_one)
  
  out12 <- lapply(what, function(w) {
    v1 <- out1[[w]]
    v2 <- out2[[w]]
    if (w == "value") 
      return(v1 + v2)
    if (w == "gradient"){
      if (neither_included) return(add_vector(add_vector(setNames(rep(0, length(union(gn1, gn2))), union(gn1, gn2)),v1),v2))
      if (one_includes_two) return(add_vector(v1,v2))
      if (two_includes_one) return(add_vector(v2,v1))
    }
    if (w == "hessian") {
      if (neither_included) return(add_matrix(add_matrix(matrix(0, length(union(gn1,gn2)),length(union(gn1,gn2)),
                                                                dimnames = list(union(gn1,gn2), union(gn1,gn2))
      ),v1),v2))
      if (one_includes_two) return(add_matrix(v1,v2))
      if (two_includes_one) return(add_matrix(v2,v1))
    }
  })
  names(out12) <- what
  
  # Summation of numeric attributes 
  out1.attributes <- attributes(out1)[sapply(attributes(out1), is.numeric)]
  out2.attributes <- attributes(out2)[sapply(attributes(out2), is.numeric)]
  attr.names <- union(names(out1.attributes), names(out2.attributes))
  out12.attributes <- lapply(attr.names, function(n) {
    x1 <- ifelse(is.null(out1.attributes[[n]]), 0, out1.attributes[[n]])
    x2 <- ifelse(is.null(out2.attributes[[n]]), 0, out2.attributes[[n]])
    x1 + x2
  })
  attributes(out12)[attr.names] <- out12.attributes
  
  class(out12) <- "objlist"
  return(out12)
}


#' @export
print.objlist <- function(x, n1 = 20, n2 = 6, ...) {
  n1 <- min(n1,length(x$gradient))
  n2 <- min(n2,length(x$gradient))
  cat("value\n", "==================\n",x$value, "\n")
  cat("gradient[1:",n1,"] (full length = ",length(x$gradient),")\n", "==================\n", sep = "")
  print(x$gradient[1:n1])
  cat("\n")
  cat("hessian[1:",n2,",1:",n2,"]","\n", "==================\n", sep = "")
  print(x$hessian[1:n2,1:n2])
  cat("\n\n")
  cat("attributes\n", "==================\n")
  cat(capture.output(str(attributes(x), max.level = 1)), sep = "\n")
  
}



#' @export
print.objfn <- function(x, ...) {
  
  parameters <- attr(x, "parameters")
  
  cat("Objective function:\n")
  str(args(x))
  cat("\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
  
}
