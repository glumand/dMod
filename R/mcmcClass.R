#' Construct a flat MCMC target
#'
#' Bundles a likelihood objfn, an optional prior objfn, and an optional
#' prior-sample generator into a single target object consumable by
#' [mcmc()]. Both objfns must follow dMod's `-2 log p` convention; the
#' kernel converts to log-posterior space internally.
#'
#' @param likObj An `objfn` returning \eqn{-2 \log p(y \mid \theta)},
#'   typically from [normL2()].
#' @param priorObj Optional `objfn` returning \eqn{-2 \log p(\theta)}.
#'   May be `NULL` (improper flat prior).
#' @param priorSample Optional function `function(n)` returning an
#'   `n x K` matrix of prior draws with column names matching the
#'   likelihood's parameter names. Required when [mcmc()] is called
#'   with `sequenceType = "sequential"`.
#'
#' @return An object of class `c("mcmcTarget", "flatTarget", "list")`.
#'
#' @seealso [mcmc()], [bayesNLMEMarginal()], [bayesNLMEJoint()]
#' @export
flatTarget <- function(likObj, priorObj = NULL, priorSample = NULL) {
  if (!inherits(likObj, "objfn"))
    stop("flatTarget: likObj must be an objfn.")
  if (!is.null(priorObj) && !inherits(priorObj, "objfn"))
    stop("flatTarget: priorObj must be an objfn or NULL.")
  if (!is.null(priorSample) && !is.function(priorSample))
    stop("flatTarget: priorSample must be a function(n) or NULL.")
  parNames <- attr(likObj, "parameters")
  structure(list(likObj      = likObj,
                 priorObj    = priorObj,
                 priorSample = priorSample,
                 parNames    = parNames,
                 kind        = "flat"),
            class = c("mcmcTarget", "flatTarget", "list"))
}


.buildBayesSubjectMeta <- function(omegaSpec, initFull, prdfn, data,
                                   errfn        = NULL,
                                   innerControl = list(),
                                   trustControl = list()) {
  if (is.null(errfn)) errfn <- .makeStaticErr(data)

  K        <- omegaSpec$K
  subjects <- rownames(omegaSpec$subjectEtas)
  if (is.null(subjects))
    stop(".buildBayesSubjectMeta: omegaSpec needs subject expansion ",
         "(omega(..., subjects = ...)).")
  N <- length(subjects)
  outer_names   <- names(initFull)
  eta_names_all <- as.vector(omegaSpec$subjectEtas)
  pars_full_names <- c(outer_names, eta_names_all)
  eta_idx_global  <- matrix(0L, N, K)
  eta_names_list  <- vector("list", N)
  for (i in seq_len(N)) {
    nms <- omegaSpec$subjectEtas[i, ]
    eta_names_list[[i]] <- nms
    eta_idx_global[i, ] <- match(nms, pars_full_names)
  }
  outer_idx_full  <- match(outer_names, pars_full_names)
  other_etas_init <- setNames(rep(0, length(eta_names_all)), eta_names_all)

  subject_meta <- list(
    subjects          = subjects,
    eta_idx_global    = eta_idx_global,
    eta_names         = eta_names_list,
    K                 = K,
    outer_names       = outer_names,
    outer_idx_in_full = outer_idx_full,
    other_etas_init   = other_etas_init,
    pars_full_names   = pars_full_names)

  ic <- modifyList(list(rinit = 1, rmax = 10, iterlim = 30,
                        fterm = 1e-7, mterm = 1e-7,
                        eigen_floor_relative = 1e-10),
                   innerControl)
  oc <- modifyList(list(rinit = 1, rmax = 10, iterlim = 200,
                        fterm = 1e-7, mterm = 1e-7),
                   trustControl)

  pars_probe <- setNames(numeric(length(pars_full_names)), pars_full_names)
  pars_probe[outer_names] <- initFull
  fast_meta <- .buildFastMeta(prdfn, errfn, data, subjects,
                              eta_names_list, pars_full_names, pars_probe)
  om_meta <- list(
    chol_pars = omegaSpec$cholPars,
    chol_loc  = matrix(as.integer(omegaSpec$cholLoc), ncol = 2, dimnames = NULL),
    is_diag   = as.logical(omegaSpec$isDiag))
  subject_meta$fast_meta  <- fast_meta
  subject_meta$omega_meta <- om_meta

  list(subjectMeta  = subject_meta,
       innerControl = ic,
       outerControl = oc,
       errfn        = errfn,
       N            = N,
       K            = K,
       outerNames   = outer_names)
}


#' Bayesian NLME target with Laplace-marginalised etas
#'
#' MCMC target for the posterior
#' \eqn{p(\theta_{\text{struct}}, \omega_{\text{chol}} \mid y)} with the
#' per-subject \eqn{\eta_i} integrated out via FOCEI's Laplace
#' approximation. The marginal likelihood matches the frequentist
#' [nlmeFit()] path with `method = "focei"`.
#'
#' @param obj A joint `objfn`, conventionally
#'   `normL2(data, prdfn) + constraintL2(mu = 0, Omega = omegaSpec)`.
#' @param omegaSpec An `omegaSpec` from [omega()] with subject expansion.
#' @param prdfn The prediction function used inside `obj`.
#' @param data The [datalist] used inside `obj`.
#' @param init Optional named numeric naming the outer parameters.
#'   Must contain structural names plus all of `omegaSpec$cholPars`.
#'   When `NULL`, structural names are inferred from
#'   `attr(obj, "parameters")` minus subject etas and cholPars.
#' @param priorTheta Optional `objfn` returning
#'   \eqn{-2 \log p(\theta_{\text{struct}})}.
#' @param priorOmegaObj Optional `objfn` from [priorOmega()] over
#'   `omegaSpec$cholPars`.
#' @param priorSample Function `function(n)` returning an `n x P` matrix
#'   of prior draws with column names matching the structural and
#'   Cholesky parameter names. Required for sequential SMC sampling.
#' @param errfn Optional error model (`obsfn`).
#' @param control Named list with optional `innerControl` and
#'   `trustControl` lists forwarded to focei_inner_trust.
#'
#' @return An object of class
#'   `c("mcmcTarget", "bayesNLMEMarginal", "flatTarget", "list")`,
#'   consumed by [mcmc()].
#'
#' @seealso [mcmc()], [bayesNLMEJoint()], [priorOmega()], [omega()]
#' @export
bayesNLMEMarginal <- function(obj, omegaSpec, prdfn, data,
                              init          = NULL,
                              priorTheta    = NULL,
                              priorOmegaObj = NULL,
                              priorSample   = NULL,
                              errfn         = NULL,
                              control       = list()) {
  likObj <- .buildBayesNLMEMarginalLik(obj, omegaSpec, prdfn, data,
                                       init = init, errfn = errfn,
                                       control = control)
  priorObj <- .combinePriors(priorTheta, priorOmegaObj)
  parNames <- attr(likObj, "parameters")

  structure(list(likObj          = likObj,
                 priorObj        = priorObj,
                 priorSample     = priorSample,
                 parNames        = parNames,
                 omegaSpec       = omegaSpec,
                 structuralNames = attr(likObj, "structuralNames"),
                 kind            = "flat"),
            class = c("mcmcTarget", "bayesNLMEMarginal", "flatTarget", "list"))
}


.buildBayesNLMEMarginalLik <- function(obj, omegaSpec, prdfn, data,
                                       init    = NULL,
                                       errfn   = NULL,
                                       control = list()) {
  if (!inherits(obj, "objfn"))
    stop("bayesNLMEMarginal: obj must be an objfn (normL2 + constraintL2(Omega = ...)).")
  if (!inherits(omegaSpec, "omegaSpec"))
    stop("bayesNLMEMarginal: omegaSpec must come from omega().")
  if (is.null(omegaSpec$subjectEtas))
    stop("bayesNLMEMarginal: omegaSpec must have subject expansion ",
         "(omega(..., subjects = ...)).")

  cholPars <- omegaSpec$cholPars
  if (!is.null(init)) {
    if (!all(cholPars %in% names(init)))
      stop("bayesNLMEMarginal: init must contain all of omegaSpec$cholPars: ",
           paste(setdiff(cholPars, names(init)), collapse = ", "))
    outer_names_full <- names(init)
    structural_names <- setdiff(outer_names_full, cholPars)
    init_full        <- init
  } else {
    joint_pars  <- attr(obj, "parameters")
    eta_names   <- as.vector(omegaSpec$subjectEtas)
    structural_names <- setdiff(joint_pars, c(eta_names, cholPars))
    outer_names_full <- c(structural_names, cholPars)
    init_full        <- setNames(numeric(length(outer_names_full)),
                                 outer_names_full)
  }

  meta <- .buildBayesSubjectMeta(
    omegaSpec    = omegaSpec,
    initFull     = init_full,
    prdfn        = prdfn,
    data         = data,
    errfn        = errfn,
    innerControl = control$innerControl %||% list(),
    trustControl = control$trustControl %||% list())

  N <- meta$N; K <- meta$K
  subjects   <- meta$subjectMeta$subjects
  outer_names_meta <- meta$outerNames
  data_per_subject <- lapply(subjects, function(s) data[[s]])
  names(data_per_subject) <- subjects
  times_union <- sort(unique(c(0, unlist(lapply(data_per_subject, `[[`,
                                                "time")))))
  correction_cb <- function(full_pars, joint_hessian, H_inv_list) {
    .computeFoceiCorrection(
      full_pars = full_pars, joint_hessian = joint_hessian,
      fixed = NULL, outer_names = outer_names_meta,
      H_inv_list = H_inv_list,
      prdfn = prdfn, errfn = meta$errfn, omega = omegaSpec,
      subjects = subjects, subject_etas = omegaSpec$subjectEtas,
      K = K, N = N,
      data_per_subject = data_per_subject, times_union = times_union)
  }

  eta_warm_cache <- new.env(parent = emptyenv())
  eta_warm_cache$mat <- matrix(0, N, K)

  bayesFn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                      conditions = NULL, env = NULL) {
    pars <- ..1
    if (!all(outer_names_full %in% names(pars)))
      stop("bayesNLMEMarginal: missing parameters: ",
           paste(setdiff(outer_names_full, names(pars)), collapse = ", "))

    omegaVec  <- pars[cholPars]
    L         <- omegaSpec$buildL(omegaVec)
    Omega_inv <- chol2inv(t(L))
    omega_log_det <- 2 * sum(log(diag(L)))

    outer_vec <- pars[outer_names_full]
    inner <- focei_outer_objfn(
      model_cb        = prdfn,
      err_cb          = meta$errfn,
      joint_cb        = obj,
      outer_pars      = outer_vec,
      eta_warmstart   = eta_warm_cache$mat,
      subject_meta    = meta$subjectMeta,
      Omega_inv_mat   = Omega_inv,
      Omega_log_det   = omega_log_det,
      fixed           = fixed,
      inner_ctrl      = meta$innerControl,
      correction_mode = "eager",
      correction_cb_opt = correction_cb)

    eta_warm_cache$mat <- as.matrix(inner$eta_modes)

    gradient <- setNames(numeric(length(pars)), names(pars))
    g_names <- intersect(names(inner$gradient), names(pars))
    gradient[g_names] <- inner$gradient[g_names]

    H <- matrix(0, length(pars), length(pars),
                dimnames = list(names(pars), names(pars)))
    if (!is.null(inner$hessian) && nrow(inner$hessian) > 0L) {
      h_names <- intersect(rownames(inner$hessian), names(pars))
      H[h_names, h_names] <- inner$hessian[h_names, h_names, drop = FALSE]
    }

    out <- structure(list(value    = as.numeric(inner$value),
                          gradient = gradient,
                          hessian  = H),
                     class = c("objlist", "list"))
    if (is.null(env)) env <- new.env()
    env$etaModes <- inner$eta_modes
    env$HInvList <- inner$H_inv
    attr(out, "env") <- env
    out
  }

  class(bayesFn) <- c("bayesNLMEMarginalLik", "objfn", "fn")
  attr(bayesFn, "parameters")      <- outer_names_full
  attr(bayesFn, "conditions")      <- subjects
  attr(bayesFn, "modelname")       <- attr(obj, "modelname")
  attr(bayesFn, "omegaSpec")       <- omegaSpec
  attr(bayesFn, "subjectMeta")     <- meta
  attr(bayesFn, "structuralNames") <- structural_names
  attr(bayesFn, "prdfn")           <- prdfn
  attr(bayesFn, "data")            <- data
  attr(bayesFn, "errfn")           <- meta$errfn
  attr(bayesFn, "joint_obj")       <- obj
  bayesFn
}


.makeSubjectEtaObj <- function(i, subjectMeta, modelCb, errCb,
                                parsFullTemplate, etaIndices) {
  meta_i  <- subjectMeta$fast_meta[[i]]
  eta_idx <- as.integer(etaIndices)
  fn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                 conditions = NULL, env = NULL, .pars_full = NULL,
                 .Omega_inv = NULL, .Omega_log_det = NULL) {
    eta_i <- ..1
    p_full <- if (is.null(.pars_full)) parsFullTemplate else .pars_full
    p_full[eta_idx] <- as.numeric(eta_i)
    kr <- focei_eval_one_subject(
      model_cb       = modelCb,
      err_cb         = errCb,
      pars_full      = p_full,
      fixed          = NULL,
      meta_i         = meta_i,
      eta_block      = as.numeric(eta_i),
      Omega_inv      = .Omega_inv,
      Omega_log_det  = .Omega_log_det)
    out <- structure(list(value = kr$value, gradient = kr$gradient,
                          hessian = kr$hessian),
                     class = c("objlist", "list"))
    attr(out, "env") <- env
    out
  }
  class(fn) <- c("objfn", "fn")
  attr(fn, "parameters") <- meta_i$eta_names
  fn
}


.makeJointThetaOmegaObj <- function(jointObj, priorTheta, priorOmegaObj,
                                     parsFullTemplate, outerNames) {
  outer_idx <- match(outerNames, names(parsFullTemplate))
  if (anyNA(outer_idx))
    stop(".makeJointThetaOmegaObj: outerNames not in parsFullTemplate.")
  fn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                 conditions = NULL, env = NULL, .etas = NULL) {
    p_outer <- ..1
    p_full <- parsFullTemplate
    if (!is.null(.etas)) p_full[names(.etas)] <- as.numeric(.etas)
    p_full[names(p_outer)] <- as.numeric(p_outer)

    out_lik <- jointObj(pars = p_full, fixed = fixed,
                        deriv = deriv, conditions = conditions, env = env)
    g <- setNames(numeric(length(p_outer)), names(p_outer))
    if (length(out_lik$gradient)) {
      g_names <- intersect(names(p_outer), names(out_lik$gradient))
      g[g_names] <- out_lik$gradient[g_names]
    }
    H <- matrix(0, length(p_outer), length(p_outer),
                dimnames = list(names(p_outer), names(p_outer)))
    if (!is.null(out_lik$hessian) && nrow(out_lik$hessian) > 0L) {
      h_names <- intersect(names(p_outer), rownames(out_lik$hessian))
      H[h_names, h_names] <- out_lik$hessian[h_names, h_names, drop = FALSE]
    }
    out <- structure(list(value = out_lik$value, gradient = g, hessian = H),
                     class = c("objlist", "list"))
    if (!is.null(priorTheta)) {
      pt <- priorTheta(pars = p_outer, deriv = deriv)
      out <- out + pt
    }
    if (!is.null(priorOmegaObj)) {
      po <- priorOmegaObj(pars = p_outer, deriv = deriv)
      out <- out + po
    }
    out
  }
  class(fn) <- c("objfn", "fn")
  attr(fn, "parameters") <- outerNames
  fn
}


#' Bayesian NLME target with joint sampling over etas
#'
#' MCMC target for the full hierarchical posterior
#' \eqn{p(\theta, \omega, \{\eta_i\} \mid y)} without marginalisation.
#' [mcmc()] routes this target through a Particle-Gibbs orchestrator that
#' alternates per-subject \eqn{\eta_i} blocks with a joint
#' \eqn{(\theta, \omega)} block.
#'
#' @param obj A joint `objfn`.
#' @param omegaSpec An `omegaSpec` from [omega()] with subject expansion.
#' @param prdfn The prediction function used inside `obj`.
#' @param data The [datalist] used inside `obj`.
#' @param init Optional named numeric naming the outer parameters.
#' @param priorTheta Optional structural prior.
#' @param priorOmegaObj Optional Omega-Cholesky prior from [priorOmega()].
#' @param priorSample Required function `function(n)` returning an
#'   `n x P` matrix of prior draws over the outer parameter set.
#' @param errfn Optional error model (`obsfn`).
#'
#' @return An object of class
#'   `c("mcmcTarget", "bayesNLMEJoint", "list")`, consumed by [mcmc()].
#'
#' @seealso [mcmc()], [bayesNLMEMarginal()], [priorOmega()]
#' @export
bayesNLMEJoint <- function(obj, omegaSpec, prdfn, data,
                           init          = NULL,
                           priorTheta    = NULL,
                           priorOmegaObj = NULL,
                           priorSample   = NULL,
                           errfn         = NULL) {
  if (!inherits(obj, "objfn"))
    stop("bayesNLMEJoint: obj must be an objfn.")
  if (!inherits(omegaSpec, "omegaSpec"))
    stop("bayesNLMEJoint: omegaSpec must come from omega().")
  if (is.null(omegaSpec$subjectEtas))
    stop("bayesNLMEJoint: omegaSpec must have subject expansion.")
  if (is.null(priorSample))
    stop("bayesNLMEJoint: priorSample is required.")
  if (is.null(errfn)) errfn <- .makeStaticErr(data)

  eta_names_all <- as.vector(omegaSpec$subjectEtas)
  cholPars      <- omegaSpec$cholPars
  if (!is.null(init)) {
    if (!all(cholPars %in% names(init)))
      stop("bayesNLMEJoint: init must contain all of omegaSpec$cholPars.")
    outer_names_full <- names(init)
    structural_names <- setdiff(outer_names_full, cholPars)
    init_full        <- init
  } else {
    joint_pars       <- attr(obj, "parameters")
    structural_names <- setdiff(joint_pars, c(eta_names_all, cholPars))
    outer_names_full <- c(structural_names, cholPars)
    init_full        <- setNames(numeric(length(outer_names_full)),
                                 outer_names_full)
  }

  meta_pkg <- .buildBayesSubjectMeta(omegaSpec, init_full, prdfn, data, errfn)
  meta <- meta_pkg$subjectMeta
  N    <- meta_pkg$N
  K    <- meta_pkg$K
  subjects <- meta$subjects

  structure(list(jointObj         = obj,
                 omegaSpec        = omegaSpec,
                 prdfn            = prdfn,
                 data             = data,
                 errfn            = errfn,
                 init             = init_full,
                 priorTheta       = priorTheta,
                 priorOmegaObj    = priorOmegaObj,
                 priorSample      = priorSample,
                 subjectMeta      = meta,
                 N                = N,
                 K                = K,
                 subjects         = subjects,
                 etaNamesAll      = eta_names_all,
                 cholPars         = cholPars,
                 outerNames       = outer_names_full,
                 structuralNames  = structural_names,
                 kind             = "blocked"),
            class = c("mcmcTarget", "bayesNLMEJoint", "list"))
}


# Compute partial G[a,b] / partial theta[c] for a single condition.
.metric_derivative_one_condition <- function(prediction, data, err,
                                             par_names) {
  rf <- res(data, prediction, err = err)
  J  <- attr(rf, "deriv",  exact = TRUE)
  D2 <- attr(rf, "deriv2", exact = TRUE)
  if (is.null(J) || is.null(D2))
    stop("metric derivative needs both 'deriv' and 'deriv2' on the prediction.")

  K  <- length(par_names)
  n  <- nrow(rf)
  s  <- rf$sigma
  w  <- 1 / (s * s)

  Jp_full  <- matrix(0, n, K, dimnames = list(NULL, par_names))
  jp_loc   <- colnames(J)
  Jp_full[, intersect(jp_loc, par_names)] <-
    J[, intersect(jp_loc, par_names), drop = FALSE]

  d2_loc   <- dimnames(D2)[[2]]
  D2p_full <- array(0, c(n, K, K),
                    dimnames = list(NULL, par_names, par_names))
  keep <- intersect(d2_loc, par_names)
  D2p_full[, keep, keep] <- D2[, keep, keep, drop = FALSE]

  dG <- array(0, c(K, K, K),
              dimnames = list(par_names, par_names, par_names))
  for (c_idx in seq_len(K)) {
    D2_c <- D2p_full[, c_idx, , drop = TRUE]
    if (n == 1L) D2_c <- matrix(D2_c, nrow = 1L)
    term1 <- crossprod(D2_c * w, Jp_full)
    term2 <- crossprod(Jp_full * w, D2_c)
    dG[, , c_idx] <- term1 + term2
  }
  dG
}


# Compute the full dG/dtheta tensor by summing per-condition contributions.
.computeMetricDerivative <- function(metric_context, theta) {

  data     <- metric_context$data
  prdfn    <- metric_context$prdfn
  errmodel <- metric_context$errmodel
  par_names <- names(theta)

  conditions <- intersect(names(data), names(attr(prdfn, "mappings")))
  if (!length(conditions))
    stop("metricContext: no conditions shared between data and prdfn.")

  timesD <- sort(unique(c(0,
    unlist(lapply(data[conditions], function(d) d$time)))))

  prediction <- prdfn(times = timesD, pars = theta,
                      deriv = TRUE, deriv2 = TRUE,
                      conditions = conditions)

  K <- length(par_names)
  dG <- array(0, c(K, K, K),
              dimnames = list(par_names, par_names, par_names))

  for (cn in conditions) {
    pi   <- prediction[[cn]]
    erri <- NULL
    if (!is.null(errmodel) && cn %in% names(attr(errmodel, "mappings"))) {
      pinner <- getParameters(pi)
      pinner <- as.parvec(pinner)
      erri <- errmodel(out = pi, pars = pinner,
                       conditions = cn)[[cn]]
    }
    dG <- dG + .metric_derivative_one_condition(pi, data[[cn]], erri,
                                                par_names)
  }
  dG
}


#' Context for the full-RMALA path of [mcmc()]
#'
#' Bundles `data`, `prdfn` and an optional `errmodel` so the Christoffel
#' tensor \eqn{\partial G / \partial \theta} can be evaluated from the
#' deriv2 chain. Required by `langevinControl(correction = "full")`. The
#' prediction chain must have been built with `deriv2 = TRUE`.
#'
#' @param data A [datalist].
#' @param prdfn A [prdfn], typically `g * x * p`.
#' @param errmodel Optional [obsfn] error model.
#'
#' @return A list of class `metricContext` consumed by `mcmc`.
#' @seealso [mcmc()], [langevinControl()], [metricControl()]
#' @export
metricContext <- function(data, prdfn, errmodel = NULL) {
  if (!inherits(prdfn, "prdfn"))
    stop("`prdfn` must be of class 'prdfn'.")
  if (!inherits(data, "datalist"))
    data <- as.datalist(data)
  out <- list(data = data, prdfn = prdfn, errmodel = errmodel)
  class(out) <- c("metricContext", "list")
  out
}


#' Hyperprior on a random-effects covariance Omega
#'
#' Returns an `objfn` over `omegaSpec$cholPars` representing
#' \eqn{-2 \log p(\Omega)}: LKJ(\eqn{\eta}) on the correlation matrix
#' combined with a half-Normal or half-Cauchy prior on the marginal
#' standard deviations \eqn{\sigma_k = \sqrt{\Omega_{kk}}}. Gradient and
#' Hessian are evaluated analytically (the Hessian is block-diagonal by row
#' of \eqn{L}).
#'
#' @param omegaSpec An `omegaSpec` from [omega()].
#' @param kind Prior family. `"LKJHalfNormal"` and `"LKJHalfCauchy"` use the
#'   indicated marginal-scale prior together with LKJ on the correlation.
#'   `"inverseWishart"` is not yet implemented.
#' @param lkjEta LKJ shape \eqn{\eta > 0}. Default 2 (mildly favours the
#'   identity correlation). Ignored when `omegaSpec` is purely diagonal.
#' @param scaleSD Scale of the marginal-SD prior. Default 1.
#' @param df Degrees of freedom for the unimplemented `inverseWishart`
#'   variant.
#'
#' @return An `objfn` returning an [objlist] with `value = -2 log p`,
#'   `gradient` and `hessian` over `omegaSpec$cholPars`.
#'
#' @references
#'   Lewandowski, D., Kurowicka, D. and Joe, H. (2009). Generating random
#'   correlation matrices based on vines and extended onion method. J. Mult.
#'   Anal. 100(9), 1989-2001.
#'
#'   Gelman, A. (2006). Prior distributions for variance parameters in
#'   hierarchical models. Bayesian Analysis 1(3), 515-534.
#'
#' @seealso [omega()], [constraintL2()], [mcmc()]
#' @export
priorOmega <- function(omegaSpec,
                       kind    = c("LKJHalfNormal", "LKJHalfCauchy",
                                   "inverseWishart"),
                       lkjEta  = 2.0,
                       scaleSD = 1.0,
                       df      = NULL) {

  if (!inherits(omegaSpec, "omegaSpec"))
    stop("`omegaSpec` must be of class 'omegaSpec'.")
  kind <- match.arg(kind)
  stopifnot(lkjEta > 0, scaleSD > 0)
  if (kind == "inverseWishart")
    stop("kind = \"inverseWishart\" is not implemented yet; ",
         "use LKJHalfNormal or LKJHalfCauchy.")

  K        <- omegaSpec$K
  cholPars <- omegaSpec$cholPars
  cholLoc  <- omegaSpec$cholLoc
  isDiag   <- omegaSpec$isDiag
  P        <- length(cholPars)
  kindFlag <- if (kind == "LKJHalfNormal") 0L else 1L

  cholLocInt <- matrix(as.integer(cholLoc), nrow = P, ncol = 2L,
                       dimnames = list(NULL, c("row", "col")))
  isDiagLog  <- as.logical(isDiag)

  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                   conditions = NULL, env = NULL) {
    pars <- ..1
    if (!all(cholPars %in% names(pars)))
      stop("priorOmega: missing Cholesky parameters: ",
           paste(setdiff(cholPars, names(pars)), collapse = ", "))
    omegaVec <- pars[cholPars]

    kr <- priorOmegaKernel(omegaVec   = as.numeric(omegaVec),
                           cholLoc    = cholLocInt,
                           isDiag     = isDiagLog,
                           K          = as.integer(K),
                           lkjEta     = as.numeric(lkjEta),
                           scaleSD    = as.numeric(scaleSD),
                           kindFlag   = kindFlag)

    if (!is.finite(kr$value)) {
      out <- structure(list(value    = Inf,
                            gradient = setNames(numeric(length(pars)),
                                                 names(pars)),
                            hessian  = matrix(0, length(pars), length(pars),
                                              dimnames = list(names(pars),
                                                              names(pars)))),
                       class = c("objlist", "list"))
      attr(out, "env") <- env
      return(out)
    }

    if (!isTRUE(deriv))
      return(structure(list(value = kr$value),
                       class = c("objlist", "list")))

    gradient <- setNames(numeric(length(pars)), names(pars))
    gradient[cholPars] <- kr$gradient

    H <- matrix(0, length(pars), length(pars),
                dimnames = list(names(pars), names(pars)))
    H[cholPars, cholPars] <- kr$hessian

    out <- structure(list(value = kr$value, gradient = gradient,
                          hessian = H),
                     class = c("objlist", "list"))
    attr(out, "env") <- env
    out
  }

  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- NULL
  attr(myfn, "parameters") <- cholPars
  attr(myfn, "modelname")  <- NULL
  myfn
}


#' Non-centered reparameterisation of subject etas
#'
#' Deterministic mapping \eqn{\eta_i = L_\Omega z_i} from standard-normal
#' latent variables to centered etas, suitable as a parfn-style
#' preprocessing step for the joint Particle-Gibbs path of [mcmc()].
#' Removes Neal's-funnel geometry. Returns mapping closures plus the
#' generated z parameter names.
#'
#' @param omegaSpec An `omegaSpec` from [omega()] with subject expansion.
#' @param zPrefix Character prefix for the latent z parameter names.
#'   Defaults to `"z"`. Each subject's K z names are
#'   `paste0(zPrefix, "_", eta_label, "_", subject)`.
#'
#' @return A list with
#'   * `zNames`: an `N x K` character matrix of z parameter names;
#'   * `mapZtoEta(zVec, omegaChol)`: function returning a named eta
#'     numeric vector;
#'   * `mapEtaToZ(etaVec, omegaChol)`: inverse mapping
#'     \eqn{z_i = L_\Omega^{-1} \eta_i}.
#'
#' @seealso [omega()], [mcmc()], [bayesNLMEJoint()]
#' @export
nonCenteredOmega <- function(omegaSpec, zPrefix = "z") {

  if (!inherits(omegaSpec, "omegaSpec"))
    stop("`omegaSpec` must come from omega().")
  if (is.null(omegaSpec$subjectEtas))
    stop("`omegaSpec` must have subject expansion (omega(..., subjects = ...)).")

  K        <- omegaSpec$K
  subjects <- rownames(omegaSpec$subjectEtas)
  N        <- length(subjects)
  eta_lab  <- sub("^eta_", "", omegaSpec$eta)
  z_mat    <- outer(subjects, eta_lab,
                    function(s, e) paste0(zPrefix, "_", e, "_", s))
  dimnames(z_mat) <- list(subjects, omegaSpec$eta)

  buildL <- omegaSpec$buildL
  subj_eta_names <- omegaSpec$subjectEtas
  z_names_flat   <- as.vector(z_mat)
  eta_names_flat <- as.vector(subj_eta_names)

  mapZtoEta <- function(zVec, omegaChol) {
    if (!all(z_names_flat %in% names(zVec)))
      stop("nonCenteredOmega$mapZtoEta: missing z names.")
    L <- buildL(omegaChol)
    Z <- matrix(zVec[z_names_flat], nrow = N, ncol = K,
                dimnames = dimnames(z_mat))
    E <- Z %*% t(L)
    out <- as.numeric(E)
    names(out) <- eta_names_flat
    out
  }

  mapEtaToZ <- function(etaVec, omegaChol) {
    if (!all(eta_names_flat %in% names(etaVec)))
      stop("nonCenteredOmega$mapEtaToZ: missing eta names.")
    L <- buildL(omegaChol)
    E <- matrix(etaVec[eta_names_flat], nrow = N, ncol = K,
                dimnames = dimnames(subj_eta_names))
    Z <- t(forwardsolve(L, t(E)))
    out <- as.numeric(Z)
    names(out) <- z_names_flat
    out
  }

  structure(list(zNames    = z_mat,
                 mapZtoEta = mapZtoEta,
                 mapEtaToZ = mapEtaToZ,
                 omegaSpec = omegaSpec),
            class = c("nonCenteredOmega", "list"))
}


#' @export
print.nonCenteredOmega <- function(x, ...) {
  cat("nonCenteredOmega (",
      nrow(x$zNames), " subjects, ",
      ncol(x$zNames), " etas)\n", sep = "")
  cat(" z names per subject:\n")
  print(x$zNames)
  invisible(x)
}


#' Analytical d OFV / d(omegaChol) for the FOCEI marginal likelihood
#'
#' Closed-form gradient of the Laplace-approximated marginal
#' \eqn{-2 \log p(y \mid \theta, \Omega)} with respect to
#' `omegaSpec$cholPars`, given the per-subject modes \eqn{\hat\eta_i} and
#' inverse Hessians \eqn{H_i^{-1}} at the FOCEI inner optimum.
#' Fully vectorised over subjects and Cholesky parameters.
#'
#' @param omegaChol Named numeric vector of Cholesky parameter values
#'   (matching `omegaSpec$cholPars`).
#' @param omegaSpec An `omegaSpec` from [omega()].
#' @param etaModes Numeric matrix `[N x K]` with row \eqn{i} equal to
#'   \eqn{\hat\eta_i}.
#' @param HInvList List of `N` symmetric `[K x K]` matrices, one per
#'   subject, representing \eqn{H_i^{-1}} at the inner optimum.
#'
#' @return A named numeric vector of length `length(omegaSpec$cholPars)`
#'   carrying \eqn{\partial \mathrm{OFV} / \partial \omega_{cholPar}}.
#'
#' @seealso [priorOmega()] (which supplies the Hessian metric on the Omega
#'   axis in Bayesian use)
#' @export
foceiOmegaGradient <- function(omegaChol, omegaSpec, etaModes, HInvList) {

  if (!inherits(omegaSpec, "omegaSpec"))
    stop("`omegaSpec` must be of class 'omegaSpec'.")
  stopifnot(all(omegaSpec$cholPars %in% names(omegaChol)))
  stopifnot(is.matrix(etaModes))
  K_eta <- omegaSpec$K
  if (ncol(etaModes) != K_eta)
    stop(sprintf("etaModes must have %d columns (K_eta), got %d.",
                 K_eta, ncol(etaModes)))
  N <- nrow(etaModes)
  if (length(HInvList) != N)
    stop("HInvList must have one entry per subject (row of etaModes).")

  L         <- omegaSpec$buildL(omegaChol[omegaSpec$cholPars])
  Omega_inv <- chol2inv(t(L))
  Linv_T    <- backsolve(t(L), diag(K_eta))

  U <- t(forwardsolve(L, t(etaModes)))
  V <- etaModes %*% Omega_inv

  sum_H_inv <- Reduce(`+`, HInvList)
  M_sum_L   <- Omega_inv %*% sum_H_inv %*% Linv_T

  rows   <- omegaSpec$cholLoc[, "row"]
  cols   <- omegaSpec$cholLoc[, "col"]
  isDiag <- as.logical(omegaSpec$isDiag)

  quadform <- -2 * colSums(U[, cols, drop = FALSE] *
                           V[, rows, drop = FALSE])
  logdetH  <- -4 * M_sum_L[cbind(rows, cols)]
  logdetOmega <- ifelse(isDiag, 2 * N, 0)
  chain <- ifelse(isDiag, L[cbind(rows, rows)], 1)

  grad <- chain * (quadform + logdetH) + logdetOmega
  setNames(grad, omegaSpec$cholPars)
}


# Sum two objfns, returning NULL when both inputs are NULL.
.combinePriors <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  a + b
}


#' @export
print.mcmcResult <- function(x, ...) {
  type <- if (inherits(x, "mcmcResultMulti"))      "multi-chain"
          else if (inherits(x, "mcmcResultSequential")) "SMC"
          else if (inherits(x, "mcmcResultBlocked"))    "blocked Particle-Gibbs"
          else "single-chain"
  cat("mcmcResult [", type, "]: ", nrow(x$samples), " samples, ",
      ncol(x$samples), " parameters\n", sep = "")
  if (!is.null(x$moveType))
    cat(" moveType: ", x$moveType, "    metric: ",
        x$metric %||% "n/a", "\n", sep = "")
  if (!is.null(x$acceptRate))
    cat(" acceptance rate: ", format(x$acceptRate, digits = 3), "\n", sep = "")
  if (!is.null(x$logEvidence))
    cat(" log-evidence:    ", format(x$logEvidence, digits = 6), "\n", sep = "")
  if (!is.null(x$rHat) && any(is.finite(x$rHat))) {
    rh <- x$rHat[is.finite(x$rHat)]
    cat(" R-hat range:     ", format(min(rh), digits = 3), " - ",
        format(max(rh), digits = 3), "\n", sep = "")
  }
  invisible(x)
}
