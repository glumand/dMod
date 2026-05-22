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
#' @keywords internal
#' @noRd
constraintExp2 <- function(p, mu, sigma = 1, k = 0.05, fixed=NULL) {
  
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


#' Per-condition residual contribution to an L2 objective
#'
#' @description
#' Computes the negative-log-likelihood residual contribution of a single
#' condition (with optional error model). Exposed so quadrature node-loops
#' can evaluate one condition without paying the per-call cost of
#' [normL2]'s multi-condition setup.
#'
#' @param dataI datalist entry for one condition (data.frame with
#'   `name`, `time`, `value`, `sigma` columns).
#' @param predictionI prdframe for that condition (typically `prediction[[cn]]`
#'   from a prdfn call).
#' @param pars Named numeric parameter vector at which to evaluate.
#' @param errfn Optional obsfn defining a parameter-dependent error model.
#' @param fixed Optional fixed-parameter vector (passed through to `errfn`).
#' @param cn Character condition name. Required when `errfn` is set, used
#'   for errmodel condition routing.
#' @param eCondNames Optional character vector of condition names that have
#'   an errmodel mapping. NULL means `errfn` applies to all.
#' @param bessel Bessel correction factor (default 1, matching the
#'   `use.bessel = FALSE` branch of `normL2`).
#' @param deriv,deriv2 Logical. Whether to return gradient/Hessian.
#' @param opt.BLOQ Character. BLOQ likelihood treatment forwarded to [nll()].
#'   One of `"M1"`, `"M3"` (default), `"M4NM"`, `"M4BEAL"`.
#'
#' @return An [objlist] for the single condition's contribution.
#' @export
evalConditionResidual <- function(dataI, predictionI, pars,
                                  errfn      = NULL,
                                  fixed      = NULL,
                                  cn         = NULL,
                                  eCondNames = NULL,
                                  bessel     = 1,
                                  deriv      = TRUE,
                                  deriv2     = FALSE,
                                  opt.BLOQ   = c("M3", "M1", "M4NM", "M4BEAL")) {
  opt.BLOQ <- match.arg(opt.BLOQ)
  err_cn <- NULL
  if (!is.null(errfn) && (is.null(eCondNames) || cn %in% eCondNames)) {
    if (is.null(cn))
      stop("evalConditionResidual: `cn` must be supplied when `errfn` is set.")
    pinner     <- getParameters(predictionI)
    fixedinner <- pinner[attr(pinner, "fixed")]
    pinner     <- as.parvec(pinner[setdiff(names(pinner), names(fixed))])
    fixedinner <- as.parvec(fixedinner, deriv = FALSE, deriv2 = FALSE)
    err_cn <- errfn(out = predictionI, pars = pinner,
                    fixed = fixedinner, conditions = cn)[[cn]]
  }
  nll(res(dataI, predictionI, err_cn),
      pars = pars, deriv = deriv, deriv2 = deriv2,
      opt.BLOQ = opt.BLOQ,
      bessel.correction = bessel)
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
#' @param threads Integer. Per-call OpenMP threads passed to the C++ residual
#'   kernel (used when the `dMod.objfn.cpp` path is active). Capped by
#'   `getOption("dMod.objfn.threads")`. Default 1.
#' @param opt.BLOQ Character. NONMEM-style treatment of below-LOQ rows
#'   (those with `value <= lloq` in the data). One of `"M1"` (drop BLOQ rows
#'   from the objective), `"M3"` (censored log-likelihood, default), `"M4NM"`
#'   or `"M4BEAL"` (truncated variants; require non-negative LOQ). Forwarded
#'   to [nll()] for the R path and to the C++ kernel.
#'
#' @return
#' An object of class `objfn`, i.e. a function
#' \code{obj(pars, fixed, deriv, env)} returning an [objlist].
#'
#' @details
#' Combine objectives with `+` (see [sumobjfn]). The Bessel correction
#' \eqn{\sqrt{n/(n-p)}} is applied globally (\eqn{n} = total data points,
#' \eqn{p} = structural parameters). When `cores > 1`, conditions are
#' evaluated in parallel; the core count is fixed at construction.
#'
#' @example inst/examples/normL2.R
#' @export
normL2 <- function(data, x, errmodel = NULL, times = NULL,
                   attr.name = "data", use.bessel = !is.null(errmodel),
                   cores = 1L, threads = 1L,
                   opt.BLOQ = c("M3", "M1", "M4NM", "M4BEAL")) {

  stopifnot(cores >= 1L)
  stopifnot(threads >= 1L)
  opt.BLOQ <- match.arg(opt.BLOQ)

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
  force(threads)

  # Lazy meta cache for the C++ kernel path. Built on first call; rebuilt
  # if the deriv column set changes (e.g. when `fixed` toggles between
  # calls — uncommon, but cheap to detect via length+name compare).
  .meta_cache <- new.env(parent = emptyenv())
  .meta_cache$meta_list        <- NULL
  .meta_cache$par_names_global <- NULL
  .meta_cache$signature        <- NULL  # used to invalidate on shape change

  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, env = NULL) {
    pars <- ..1
    if (is.null(env)) env <- new.env()

    prediction <- x(times = timesD, pars = pars, fixed = fixed,
                    deriv = deriv, deriv2 = deriv2, conditions = conditions)

    use_cpp <- isTRUE(getOption("dMod.objfn.cpp", FALSE)) && deriv

    if (use_cpp) {
      # Build errmodel output per condition (if any).
      err_list <- NULL
      if (!is.null(errmodel)) {
        err_list <- lapply(conditions, function(cn) {
          if (!is.null(e.cond) && !(cn %in% e.cond)) return(NULL)
          pinner     <- getParameters(prediction[[cn]])
          fixedinner <- pinner[attr(pinner, "fixed")]
          pinner     <- as.parvec(pinner[setdiff(names(pinner), names(fixed))])
          fixedinner <- as.parvec(fixedinner, deriv = FALSE, deriv2 = FALSE)
          errmodel(out = prediction[[cn]], pars = pinner,
                   fixed = fixedinner, conditions = cn)[[cn]]
        })
      }

      # Determine current deriv signature (per-condition local par names).
      cur_sig <- lapply(prediction, function(pr) dimnames(attr(pr, "deriv"))[[3]])
      if (is.null(.meta_cache$meta_list) ||
          !identical(.meta_cache$signature, cur_sig)) {
        .meta_cache$par_names_global <- unique(unlist(cur_sig))
        .meta_cache$meta_list <- .build_normL2_meta(
          data, prediction, err_list, conditions, e.cond)
        .meta_cache$signature <- cur_sig
      }

      eff_threads <- max(
        1L, min(as.integer(threads),
                as.integer(getOption("dMod.objfn.threads", threads))))

      kr <- normL2_kernel(
        prediction       = prediction,
        err_list_opt     = err_list,
        meta_list        = .meta_cache$meta_list,
        par_names_global = .meta_cache$par_names_global,
        bessel           = bessel,
        deriv2_requested = isTRUE(deriv2),
        threads          = eff_threads,
        bloq_mode        = opt.BLOQ
      )
      out <- objlist(value    = kr$value,
                     gradient = kr$gradient,
                     hessian  = kr$hessian)
      attr(out, attr.name) <- out$value
      env$prediction <- prediction
      attr(out, "env") <- env
      return(out)
    }

    # ---- R fallback path (existing behaviour) ----
    one <- function(cn) {
      evalConditionResidual(dataI = data[[cn]], predictionI = prediction[[cn]],
                            pars = pars, errfn = errmodel, fixed = fixed,
                            cn = cn, eCondNames = e.cond, bessel = bessel,
                            deriv = deriv, deriv2 = deriv2,
                            opt.BLOQ = opt.BLOQ)
    }
    out <- if (cores == 1L) {
      Reduce(`+`, lapply(conditions, one))
    } else if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      Reduce(`+`, parallel::parLapply(cl, conditions, one))
    } else {
      Reduce(`+`, parallel::mclapply(conditions, one, mc.cores = cores))
    }

    attr(out, attr.name) <- out$value
    env$prediction <- prediction
    attr(out, "env") <- env
    out
  }

  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- d.cond
  # Union of prediction-fn and errmodel parameters so the errmodel's sigma
  # parameters survive when the inner solver reads `full_pars`.
  err_pars <- if (!is.null(errmodel)) attr(errmodel, "parameters") else character(0)
  attr(myfn, "parameters") <- union(attr(x, "parameters"), err_pars)
  attr(myfn, "modelname") <- modelname(x, errmodel)
  myfn
}


# Build per-condition metadata for the C++ normL2_kernel. Indexes data rows
# into prediction/errmodel matrices, encodes ALOQ/BLOQ partition, and stores
# the LOQ-substituted y values (matching res()'s `pmax(value, lloq)`).
.build_normL2_meta <- function(data, prediction, err_list, conditions, e_cond) {
  err_list_named <- if (!is.null(err_list)) {
    setNames(err_list, conditions)
  } else {
    NULL
  }
  lapply(conditions, function(cn) {
    dataI <- data[[cn]]
    dataI$name <- as.character(dataI$name)
    prdfI <- prediction[[cn]]
    pcols <- colnames(prdfI)
    d_dn  <- dimnames(attr(prdfI, "deriv"))

    t_idx_in_pred  <- match(dataI$time, prdfI[, "time"])
    o_idx_in_pred  <- match(dataI$name, pcols)
    o_idx_in_deriv <- match(dataI$name, d_dn[[2]])

    if (anyNA(t_idx_in_pred) || anyNA(o_idx_in_pred) || anyNA(o_idx_in_deriv)) {
      stop(".build_normL2_meta: data point not found in prediction for condition '",
           cn, "'.", call. = FALSE)
    }

    sig <- if (!is.null(dataI$sigma)) dataI$sigma else rep(NA_real_, nrow(dataI))
    sigma_is_na <- is.na(sig)
    sigma_fixed <- ifelse(sigma_is_na, 0, sig)

    t_idx_in_err <- rep(0L, nrow(dataI))
    o_idx_in_err <- rep(0L, nrow(dataI))
    o_idx_in_err_deriv <- rep(0L, nrow(dataI))
    if (any(sigma_is_na) && !is.null(err_list_named)) {
      erm <- err_list_named[[cn]]
      if (!is.null(erm)) {
        t_idx_in_err <- match(dataI$time, erm[, "time"])
        o_idx_in_err <- match(dataI$name, colnames(erm))
        e_dn <- dimnames(attr(erm, "deriv"))
        if (!is.null(e_dn)) {
          o_idx_in_err_deriv <- match(dataI$name, e_dn[[2]])
          o_idx_in_err_deriv[is.na(o_idx_in_err_deriv)] <- 0L
        }
      }
    }

    lloq <- if (!is.null(dataI$lloq)) dataI$lloq else rep(-Inf, nrow(dataI))
    val  <- pmax(dataI$value, lloq)
    bloq_mask <- as.integer(val <= lloq)

    list(
      t_idx_in_pred       = as.integer(t_idx_in_pred),
      o_idx_in_pred       = as.integer(o_idx_in_pred),
      o_idx_in_deriv      = as.integer(o_idx_in_deriv),
      t_idx_in_err        = as.integer(t_idx_in_err),
      o_idx_in_err        = as.integer(o_idx_in_err),
      o_idx_in_err_deriv  = as.integer(o_idx_in_err_deriv),
      sigma_is_na         = as.integer(sigma_is_na),
      sigma_fixed         = as.numeric(sigma_fixed),
      y_data              = as.numeric(val),
      lloq                = as.numeric(lloq),
      bloq_mask           = bloq_mask
    )
  })
}



#' Soft L2 constraint on parameters
#'
#' @param mu Named numeric vector of prior means. For the MVN path
#'   (`Omega` set), `mu` may be a scalar (broadcast across all etas) or a
#'   length-K named vector matching `Omega$eta`.
#' @param sigma Named numeric or character vector. Character entries indicate
#'   log-scale sigma parameters to be estimated. Used only when `Omega` is NULL.
#' @param Omega Optional `omegaSpec` object (see [omega]) describing a
#'   multivariate Gaussian prior over subject-level random effects with full
#'   Cholesky-parametrised covariance. When supplied, the function switches to
#'   the MVN path and ignores `sigma`.
#' @param attr.name Character. Name of the attribute storing the constraint value.
#' @param condition Optional character vector of conditions.
#' @param threads Integer. Per-call OpenMP threads passed to the C++ constraint
#'   kernel. Default 1.
#'
#' @details
#' Computes, depending on which path is selected,
#' \deqn{(p-\mu)^2 / \sigma^2}
#' or, if sigma is estimated,
#' \deqn{(p-\mu)^2 / \sigma^2 + 2\log(\sigma)},
#' with sigma internally transformed via \code{exp()}.
#'
#' When `Omega` is set, computes the multivariate-normal prior
#' \deqn{\sum_i (\eta_i - \mu)^T \Omega^{-1} (\eta_i - \mu) + N \log|\Omega|}
#' over subject-level random effects \eqn{\eta_i}, where \eqn{\Omega = L L^T}
#' with \eqn{L} lower-triangular and log-parametrised on the diagonal. The
#' parameter vector at evaluation time must contain all subject-level eta
#' parameters listed in `Omega$subjectEtas` and all Cholesky parameters in
#' `Omega$cholPars`.
#'
#' @return Object of class \code{objfn}.
#' @export
constraintL2 <- function(mu, sigma = 1, Omega = NULL,
                         attr.name = "prior", condition = NULL,
                         threads = 1L) {

  if (!is.null(Omega)) {
    if (missing(mu)) mu <- 0
    return(constraintL2_mvn(mu = mu, Omega = Omega,
                            attr.name = attr.name, condition = condition,
                            threads = threads))
  }

  est <- is.character(sigma)
  if (length(sigma) == 1) sigma <- setNames(rep(sigma, length(mu)), names(mu))
  if (is.null(names(sigma))) names(sigma) <- names(mu)
  sigma <- sigma[names(mu)]
  force(threads)

  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = condition, env = NULL) {

    p <- list(...)[[match.fnargs(list(...), "pars")]]
    dP <- attr(p, "deriv", exact = TRUE)
    dP2 <- if (deriv2) attr(p, "deriv2", exact = TRUE) else NULL

    use_cpp <- isTRUE(getOption("dMod.objfn.cpp", FALSE)) && deriv

    if (use_cpp) {
      inner_par_names <- names(p)
      # Build sigma_par names (only meaningful if est==TRUE)
      sigma_pars <- if (est) sigma[names(mu)] else rep("", length(mu))
      sigma_vec  <- if (est) rep(0.0, length(mu)) else as.numeric(sigma[names(mu)])
      kr <- constraintL2_scalar_kernel(
        pars = p,
        dP_opt = if (!is.null(dP)) dP else NULL,
        dP2_opt = if (!is.null(dP2)) dP2 else NULL,
        inner_par_names = inner_par_names,
        fixed_opt = fixed,
        mu_names = names(mu),
        mu = as.numeric(mu),
        sigma = sigma_vec,
        sigma_pars = as.character(sigma_pars),
        est = est
      )
      out <- objlist(value = kr$value, gradient = kr$gradient,
                     hessian = kr$hessian)
      attr(out, attr.name) <- out$value
      attr(out, "env") <- env
      return(out)
    }

    # ---- R fallback (existing path) ----
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

      # Exact Hessian addition: gi . dP2 contributes the (dL/dp) * (d^2 p/dtheta^2)
      # term that the sandwich (dP^T H dP) drops.
      if (!is.null(dP2)) {
        common <- intersect(names(gi), dimnames(dP2)[[1]])
        if (length(common) > 0L) {
          theta_names <- colnames(dP)
          dP2_sub <- dP2[common, theta_names, theta_names, drop = FALSE]
          gi_sub <- gi[common]
          # H_add[k1, k2] = sum_p gi_p * dP2[p, k1, k2]
          flat <- matrix(dP2_sub, nrow = length(common), ncol = length(theta_names)^2)
          h_add_flat <- crossprod(flat, gi_sub)
          h_add <- matrix(h_add_flat, length(theta_names), length(theta_names))
          hs <- hs + h_add
        }
      }
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



# Multivariate-normal Gaussian prior over subject-level random effects.
# Internal: dispatched from constraintL2() when its `Omega` argument is set.
# Returns sum_i (eta_i - mu)^T Omega^-1 (eta_i - mu) + N log|Omega| with exact
# value/gradient and Gauss-Newton Hessian (block-diagonal in eta, exact crosses
# to the Cholesky parameters; sandwich via dP for chain rule).
constraintL2_mvn <- function(mu, Omega, attr.name = "prior", condition = NULL,
                              threads = 1L) {

  if (!inherits(Omega, "omegaSpec"))
    stop("`Omega` must be an omegaSpec object built by omega().")
  if (is.null(Omega$subjectEtas))
    stop("`Omega` must have subject expansion. Call omega(..., subjects = ...).")

  K            <- Omega$K
  subject_etas <- Omega$subjectEtas
  N            <- nrow(subject_etas)
  chol_pars    <- Omega$cholPars
  is_diag      <- Omega$isDiag
  chol_loc     <- Omega$cholLoc
  build_L      <- Omega$buildL

  if (length(mu) == 1L) mu <- rep(mu, K)
  if (length(mu) != K)
    stop("`mu` must have length 1 or K = ", K, " for the MVN constraintL2 path.")
  if (is.null(names(mu))) names(mu) <- Omega$eta

  all_eta_names <- as.vector(subject_etas)
  parnames      <- c(all_eta_names, chol_pars)

  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = condition, env = NULL) {

    p  <- list(...)[[match.fnargs(list(...), "pars")]]
    dP <- attr(p, "deriv", exact = TRUE)
    dP2 <- if (deriv2) attr(p, "deriv2", exact = TRUE) else NULL

    allp <- c(p, fixed)

    if (!all(parnames %in% names(allp)))
      return(objlist(value = 0,
                     gradient = setNames(numeric(length(p)), names(p)),
                     hessian  = matrix(0, length(p), length(p),
                                       dimnames = list(names(p), names(p)))))

    # --- value -------------------------------------------------------------
    chol_vec <- allp[chol_pars]
    L        <- build_L(chol_vec)

    use_cpp <- isTRUE(getOption("dMod.objfn.cpp", FALSE)) && deriv
    # The C++ path is correct only when the chol params are HELD FIXED (i.e.
    # none of `chol_pars` is in `names(p)` as a free parameter). When the
    # caller is estimating chol params (ECM workflow), fall back to R.
    if (use_cpp && length(intersect(chol_pars, names(p))) > 0L) use_cpp <- FALSE

    if (use_cpp) {
      kr <- constraintL2_mvn_kernel(
        pars                = p,
        fixed_opt           = fixed,
        dP_opt              = if (!is.null(dP)) dP else NULL,
        dP2_opt             = if (!is.null(dP2)) dP2 else NULL,
        inner_par_names     = names(p),
        K                   = K,
        N                   = N,
        all_eta_names       = all_eta_names,
        mu                  = as.numeric(mu),
        L_lower             = L,
        include_chol_block  = FALSE
      )
      out <- objlist(value = kr$value, gradient = kr$gradient,
                     hessian = kr$hessian)
      attr(out, attr.name) <- kr$value
      attr(out, "env") <- env
      return(out)
    }

    eta_mat  <- matrix(allp[all_eta_names], nrow = N, ncol = K,
                       dimnames = dimnames(subject_etas))
    R        <- t(sweep(eta_mat, 2, mu, "-"))   # K x N
    Z        <- forwardsolve(L, R)              # K x N
    W        <- backsolve(t(L), Z)              # K x N

    quad     <- sum(Z * Z)
    logdetO  <- 2 * sum(log(diag(L)))
    val      <- quad + N * logdetO

    if (!deriv)
      return(objlist(value    = val,
                     gradient = setNames(numeric(length(p)), names(p)),
                     hessian  = matrix(0, length(p), length(p),
                                       dimnames = list(names(p), names(p)))))

    np <- length(p)
    gr <- setNames(numeric(np), names(p))
    hs <- matrix(0, np, np, dimnames = list(names(p), names(p)))

    free_etas  <- intersect(all_eta_names, names(p))
    free_chols <- intersect(chol_pars,    names(p))

    # --- gradient: eta block -----------------------------------------------
    if (length(free_etas) > 0L) {
      idx_mat <- match(free_etas, all_eta_names)
      sub_idx <- ((idx_mat - 1L) %% N) + 1L
      eta_idx <- ((idx_mat - 1L) %/% N) + 1L
      gr[free_etas] <- 2 * W[cbind(eta_idx, sub_idx)]
    }

    # --- gradient: chol block ----------------------------------------------
    if (length(free_chols) > 0L) {
      WZt <- W %*% t(Z)
      for (nm in free_chols) {
        m <- match(nm, chol_pars)
        k <- chol_loc[m, 1L]; l <- chol_loc[m, 2L]
        if (is_diag[m]) {
          gr[nm] <- -2 * WZt[k, k] * L[k, k] + 2 * N
        } else {
          gr[nm] <- -2 * WZt[k, l]
        }
      }
    }

    # --- Hessian: Gauss-Newton on z_i --------------------------------------
    Linv      <- forwardsolve(L, diag(K))   # K x K, lower-triangular inverse
    Omega_inv <- crossprod(Linv)            # = L^{-T} L^{-1}

    if (length(free_chols) > 0L) {
      J_chol_template <- matrix(0, K, length(free_chols),
                                dimnames = list(NULL, free_chols))
      chol_meta <- vapply(free_chols, function(nm) {
        m <- match(nm, chol_pars); c(m, chol_loc[m, 1L], chol_loc[m, 2L],
                                     as.integer(is_diag[m]))
      }, integer(4))
      colnames(chol_meta) <- free_chols
    } else {
      J_chol_template <- NULL
    }

    for (i in seq_len(N)) {
      eta_i_names  <- subject_etas[i, ]
      eta_i_active <- eta_i_names %in% names(p)

      # eta-eta block: 2 * Omega^{-1}[active, active]
      if (any(eta_i_active)) {
        idx <- eta_i_names[eta_i_active]
        hs[idx, idx] <- hs[idx, idx] + 2 * Omega_inv[eta_i_active, eta_i_active]
      }

      if (!is.null(J_chol_template)) {
        J_chol <- J_chol_template
        for (j_idx in seq_along(free_chols)) {
          k <- chol_meta[2L, j_idx]; l <- chol_meta[3L, j_idx]
          if (chol_meta[4L, j_idx] == 1L) {
            J_chol[, j_idx] <- -L[k, k] * Z[k, i] * Linv[, k]
          } else {
            J_chol[, j_idx] <- -Z[l, i] * Linv[, k]
          }
        }
        hs[free_chols, free_chols] <- hs[free_chols, free_chols] + 2 * crossprod(J_chol)

        if (any(eta_i_active)) {
          idx   <- eta_i_names[eta_i_active]
          J_eta <- Linv[, eta_i_active, drop = FALSE]
          colnames(J_eta) <- idx
          cross <- 2 * crossprod(J_eta, J_chol)
          hs[idx, free_chols] <- hs[idx, free_chols] + cross
          hs[free_chols, idx] <- hs[free_chols, idx] + t(cross)
        }
      }
    }

    # --- chain rule via dP -------------------------------------------------
    if (!is.null(dP)) {
      gi <- gr
      gr <- drop(gi %*% dP); names(gr) <- colnames(dP)
      hs <- t(dP) %*% hs %*% dP
      dimnames(hs) <- list(colnames(dP), colnames(dP))

      # Exact Hessian addition: gi . dP2.
      if (!is.null(dP2)) {
        common <- intersect(names(gi), dimnames(dP2)[[1]])
        if (length(common) > 0L) {
          theta_names <- colnames(dP)
          dP2_sub <- dP2[common, theta_names, theta_names, drop = FALSE]
          gi_sub <- gi[common]
          flat <- matrix(dP2_sub, nrow = length(common), ncol = length(theta_names)^2)
          h_add_flat <- crossprod(flat, gi_sub)
          h_add <- matrix(h_add_flat, length(theta_names), length(theta_names))
          hs <- hs + h_add
        }
      }
    }

    out <- objlist(value = val, gradient = gr, hessian = hs)
    attr(out, attr.name) <- val
    attr(out, "env") <- env
    out
  }

  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- parnames
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
#' @param threads Integer. Per-call OpenMP threads passed to the C++ kernel. Default 1.
#' @return List of class \code{objlist}, i.e. objective value, gradient and Hessian as list.
#' @seealso [normL2], [constraintL2]
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
datapointL2 <- function(name, time, value, sigma = 1, attr.name = "validation", condition,
                         threads = 1L) {

  controls <- list(
    mu        = structure(name, names = value)[1], # only one data point is allowed
    time      = time[1],
    sigma     = sigma[1],
    attr.name = attr.name
  )
  force(threads)
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = NULL, env = NULL) {
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

    use_cpp <- isTRUE(getOption("dMod.objfn.cpp", FALSE)) && withDeriv && deriv
    if (use_cpp) {
      prdf <- prediction[[condition]]
      dpred_attr  <- attr(prdf, "deriv")
      d2pred_attr <- if (deriv2) attr(prdf, "deriv2") else NULL
      kr <- datapointL2_kernel(
        pouter           = pouter,
        fixed_opt        = fixed,
        prdf             = prdf,
        dpred_attr_opt   = dpred_attr,
        d2pred_attr_opt  = d2pred_attr,
        obs_name         = as.character(mu),
        t                = as.numeric(t),
        sigma            = as.numeric(sigma),
        value_par        = names(mu)[1]
      )
      out <- objlist(value = kr$value, gradient = kr$gradient,
                     hessian = kr$hessian)
      attr(out, attr.name)    <- out$value
      attr(out, "prediction") <- kr$prediction
      attr(out, "env")        <- env
      class(out) <- NULL
      return(out)
    }

    pred  <- prediction[[condition]][time.index, ][mu]
    deriv <- NULL
    deriv2_pred <- NULL
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
      if (deriv2) {
        d2full <- attr(prediction[[condition]], "deriv2")
        if (!is.null(d2full) && length(dim(d2full)) == 4L) {
          avail_pars <- dimnames(d2full)[[3]]
          use_pars   <- intersect(parapar, avail_pars)
          if (length(use_pars)) {
            d2tmp <- d2full[time.index, mu, use_pars, use_pars, drop = TRUE]
            deriv2_pred <- matrix(d2tmp, length(use_pars), length(use_pars),
                                  dimnames = list(use_pars, use_pars))
          }
        }
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

      if (!is.null(deriv2_pred)) {
        # Exact second-order term: 2 * res * d^2 pred / sigma^2.
        rn <- rownames(deriv2_pred)
        hs[rn, rn] <- hs[rn, rn] + 2 * res * deriv2_pred / sigma^2
      }
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
#' @seealso [normL2]
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
#' @param deriv Logical. If TRUE, compute gradient and hessian.
#' @param deriv2 Logical. If TRUE, also propagate second-order derivatives.
#' @param opt.BLOQ Character denoting the method to deal with BLOQ data.
#' One of "M1", "M3", "M4NM", or "M4BEAL".
#' @param opt.hessian Named logical vector to include or exclude various 
#' summands of the hessian matrix.
#' @param bessel.correction Numeric. Bessel correction factor for variance estimation.
#' 
#' @md
#' @return list with entries value, gradient, and hessian (Gauss-Newton approximation).
#' @export
nll <- function(nout, pars, deriv, deriv2 = FALSE,
                opt.BLOQ = "M3", opt.hessian = c(
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

  derivs2 <- if (deriv2) attr(nout, "deriv2") else NULL
  derivs2.aloq <- if (!is.null(derivs2)) derivs2[!is.bloq, , , drop = FALSE] else NULL
  derivs2.bloq <- if (!is.null(derivs2)) derivs2[is.bloq, , , drop = FALSE] else NULL

  derivs2.err <- if (deriv2) attr(nout, "deriv2.err") else NULL
  derivs2.err.aloq <- if (!is.null(derivs2.err))
    derivs2.err[!is.bloq, , , drop = FALSE] else NULL
  derivs2.err.bloq <- if (!is.null(derivs2.err))
    derivs2.err[is.bloq, , , drop = FALSE] else NULL

  n_pars <- length(pars)
  par_names <- names(pars)

  mywrss <- {
    gr <- if (deriv) setNames(numeric(n_pars), par_names) else NULL
    he <- if (deriv) matrix(0, n_pars, n_pars, dimnames = list(par_names, par_names)) else NULL
    objlist(value = 0, gradient = gr, hessian = he)
  }

  # When the caller asked only for the value, drop the deriv matrices so
  # nll_ALOQ / nll_BLOQ skip the gradient branch entirely. Important for
  # error models where derivs.err can have NULL colnames (e.g. PEtab case
  # 0015 with a single symbolic noise parameter), which otherwise trigger
  # a non-conformable-arrays error in the dwrdp/dlogsdp arithmetic.
  if (!isTRUE(deriv)) {
    derivs.aloq <- NULL; derivs.bloq <- NULL
    derivs.err.aloq <- NULL; derivs.err.bloq <- NULL
    derivs2.aloq <- NULL; derivs2.bloq <- NULL
    derivs2.err.aloq <- NULL; derivs2.err.bloq <- NULL
  }

  nll_ALOQ_result <- NULL
  if (!all(is.bloq)) {
    nll_ALOQ_result <- nll_ALOQ(nout.aloq, derivs.aloq, derivs.err.aloq,
                                derivs2 = derivs2.aloq,
                                derivs2.err = derivs2.err.aloq,
                                par_names = par_names,
                                opt.BLOQ = opt.BLOQ, opt.hessian = opt.hessian,
                                bessel.correction = bessel.correction)
  }
  mywrss <- mywrss + nll_ALOQ_result

  if (any(is.bloq) && opt.BLOQ != "M1") {
    mywrss <- mywrss + nll_BLOQ(nout.bloq, derivs.bloq, derivs.err.bloq,
                                derivs2 = derivs2.bloq,
                                derivs2.err = derivs2.err.bloq,
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
#' @param derivs2 Optional 3D array of second prediction derivatives. If `NULL`,
#'   only first-order contributions through the structural model are propagated.
#' @param derivs2.err Optional 3D array of second error-model derivatives
#'   (`d^2 sigma / d theta^2`). If `NULL`, the exact d2sigma Hessian term is
#'   skipped (mathematically zero when sigma does not depend on parameters).
#' @param par_names Character vector of ALL parameter names (full set)
#' @param opt.BLOQ Character denoting the method to deal with BLOQ data
#' @param opt.hessian Named logical vector for hessian components
#' @param bessel.correction Numeric. Bessel correction factor.
#' @md
#' @importFrom stats pnorm dnorm
nll_ALOQ <- function(nout, derivs, derivs.err,
                     derivs2 = NULL,
                     derivs2.err = NULL,
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
    
    # Compute derivatives of the (possibly bessel-scaled) weighted residual.
    # We want d wr_scaled / d theta where wr_scaled = bessel * wr_orig:
    #   d wr_scaled / d theta = bessel * (inv_s * dxdp - wr_orig * inv_s * dsdp)
    #                         = bessel * inv_s * dxdp - (wr_scaled * inv_s) * dsdp
    # `wr` carries the SCALED value at this point (see use_bessel block above),
    # so the dsdp term is correct as is; only the dxdp term needs the bessel
    # factor. The earlier implementation multiplied the ENTIRE dwrdp by bessel
    # AFTER the wr-was-scaled formula, which double-applied bessel to the dsdp
    # term and broke the gradient for parameters that affect sigma (e.g.
    # `sigma_add` in an additive errmodel). Same fix for dw0dp / BLOQ.
    if (use_bessel) {
      dwrdp <- bessel.correction * inv_s * dxdp - (wr * inv_s) * dsdp
      dw0dp <- bessel.correction * inv_s * dxdp - (w0 * inv_s) * dsdp
    } else {
      dwrdp <- inv_s * dxdp - (wr * inv_s) * dsdp
      dw0dp <- inv_s * dxdp - (w0 * inv_s) * dsdp
    }
    dlogsdp <- inv_s * dsdp
    
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

    # Exact Hessian: residual second-derivative contribution.
    # For wr = (pred - val) / s and constant s in theta:
    # d^2(wr_i)/dtheta^2 = inv_s_i * d^2(pred_i)/dtheta^2
    # H_exact_addition = 2 * sum_i wr_i * d^2(wr_i)/dtheta^2
    # (We keep the GN crossprod above; this term is the contribution that
    # GN drops. Sigma-theta cross terms remain GN-approximated when the
    # error model shares parameters with the structural model.)
    if (!is.null(derivs2)) {
      d2_local_pars <- dimnames(derivs2)[[2]]
      common_d2 <- intersect(local_pars, d2_local_pars)
      if (length(common_d2) > 0L) {
        idx_d2_local <- match(common_d2, local_pars)
        # weight per residual: 2 * wr_i * inv_s_i
        wts <- 2 * wr * inv_s
        # contract first axis (residual index) of derivs2 with weights:
        # H_add[j,k] = sum_i wts_i * derivs2[i, j, k]
        d2_sub <- derivs2[, common_d2, common_d2, drop = FALSE]
        n_common <- length(common_d2)
        # Use einsum-style contraction via apply or matrix product on flattened
        # second/third axes. Reshape d2_sub to [n_data, n_common^2] then matmul
        # against wts and reshape back.
        flat <- matrix(d2_sub, nrow = nrow(d2_sub), ncol = n_common * n_common)
        h_add_flat <- crossprod(flat, wts)
        h_add <- matrix(h_add_flat, n_common, n_common)
        hessian_local[idx_d2_local, idx_d2_local] <-
          hessian_local[idx_d2_local, idx_d2_local] + h_add
      }
    }

    # Exact Hessian: error-model second-derivative contribution.
    # From 2*wr*d2wr (which contributes -2 wr^2 / sigma * d2sigma) plus
    # 2*d2(log sigma) (which contributes +2/sigma * d2sigma), the net weight
    # on d2sigma is (2/sigma)(1 - wr^2). Skipped when sigma does not depend
    # on theta (derivs2.err is then either NULL or all zero).
    if (!is.null(derivs2.err) && length(local_pars_err) > 0L) {
      d2e_local_pars <- dimnames(derivs2.err)[[2]]
      common_d2e <- intersect(local_pars, d2e_local_pars)
      if (length(common_d2e) > 0L) {
        idx_d2e_local <- match(common_d2e, local_pars)
        wts_sig <- 2 * inv_s * (1 - wr^2)
        d2e_sub <- derivs2.err[, common_d2e, common_d2e, drop = FALSE]
        n_common_e <- length(common_d2e)
        flat_e <- matrix(d2e_sub, nrow = nrow(d2e_sub),
                         ncol = n_common_e * n_common_e)
        h_add_flat_e <- crossprod(flat_e, wts_sig)
        h_add_e <- matrix(h_add_flat_e, n_common_e, n_common_e)
        hessian_local[idx_d2e_local, idx_d2e_local] <-
          hessian_local[idx_d2e_local, idx_d2e_local] + h_add_e
      }
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
#' @param derivs2,derivs2.err Optional 3D arrays of second prediction and
#'   error-model derivatives over the BLOQ rows. When non-`NULL`, the exact
#'   d^2 pred / d^2 sigma contributions to the Hessian are added on top of
#'   the Gauss-Newton-plus-Parts form.
#' @param par_names Character vector of ALL parameter names (full set)
#' @param opt.BLOQ Character denoting the method to deal with BLOQ data
#' @param opt.hessian Named logical vector for hessian components
#' @importFrom stats pnorm dnorm
nll_BLOQ <- function(nout.bloq, derivs.bloq, derivs.err.bloq,
                     derivs2 = NULL,
                     derivs2.err = NULL,
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

    # Exact Hessian: prediction second-derivative contribution.
    # See math reference in src/residual_kernel.h.
    if (!is.null(derivs2)) {
      d2_local_pars <- dimnames(derivs2)[[2]]
      common_d2 <- intersect(local_pars, d2_local_pars)
      if (length(common_d2) > 0L) {
        idx_d2_local <- match(common_d2, local_pars)
        if (opt.BLOQ[1] == "M3") {
          G_neg_wr <- exp(stats::dnorm(-wr, log = TRUE) -
                          stats::pnorm(-wr, log.p = TRUE))
          wts <- 2 * G_neg_wr * inv_s
        } else {
          # M4*: c1 + c3 - c2 weight on d2pred (both d2wr and d2w0 contribute
          # inv_s * d2pred).
          dP    <- stats::pnorm(w0) - stats::pnorm(wr)
          c1    <- stats::dnorm(wr) / dP
          c2    <- stats::dnorm(w0) / dP
          c3    <- exp(stats::dnorm(w0, log = TRUE) -
                       stats::pnorm(w0, log.p = TRUE))
          wts   <- 2 * inv_s * (c1 + c3 - c2)
        }
        d2_sub <- derivs2[, common_d2, common_d2, drop = FALSE]
        n_common <- length(common_d2)
        flat <- matrix(d2_sub, nrow = nrow(d2_sub),
                       ncol = n_common * n_common)
        h_add_flat <- crossprod(flat, wts)
        h_add <- matrix(h_add_flat, n_common, n_common)
        hessian_local[idx_d2_local, idx_d2_local] <-
          hessian_local[idx_d2_local, idx_d2_local] + h_add
      }
    }

    # Exact Hessian: error-model second-derivative contribution.
    if (!is.null(derivs2.err) && length(local_pars_err) > 0L) {
      d2e_local_pars <- dimnames(derivs2.err)[[2]]
      common_d2e <- intersect(local_pars, d2e_local_pars)
      if (length(common_d2e) > 0L) {
        idx_d2e_local <- match(common_d2e, local_pars)
        if (opt.BLOQ[1] == "M3") {
          G_neg_wr <- exp(stats::dnorm(-wr, log = TRUE) -
                          stats::pnorm(-wr, log.p = TRUE))
          wts_sig <- -2 * wr * G_neg_wr * inv_s
        } else {
          dP    <- stats::pnorm(w0) - stats::pnorm(wr)
          c1    <- stats::dnorm(wr) / dP
          c2    <- stats::dnorm(w0) / dP
          c3    <- exp(stats::dnorm(w0, log = TRUE) -
                       stats::pnorm(w0, log.p = TRUE))
          # 2 (c1 d2wr + (c3 - c2) d2w0) projects -wr/sigma and -w0/sigma
          # respectively onto d2sigma: net -2/sigma * (c1*wr + (c3 - c2)*w0).
          wts_sig <- -2 * inv_s * (c1 * wr + (c3 - c2) * w0)
        }
        d2e_sub <- derivs2.err[, common_d2e, common_d2e, drop = FALSE]
        n_common_e <- length(common_d2e)
        flat_e <- matrix(d2e_sub, nrow = nrow(d2e_sub),
                         ncol = n_common_e * n_common_e)
        h_add_flat_e <- crossprod(flat_e, wts_sig)
        h_add_e <- matrix(h_add_flat_e, n_common_e, n_common_e)
        hessian_local[idx_d2e_local, idx_d2e_local] <-
          hessian_local[idx_d2e_local, idx_d2e_local] + h_add_e
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
