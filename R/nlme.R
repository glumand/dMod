## Stage-2 d log|H_GN| / d theta correction. Called as an Rcpp::Function from
## the C++ FOCEI kernel once per outer iter (after the inner trust has
## converged at the modes), then added into the outer gradient.
##
## Math identities: envelope theorem for d/d theta at eta = eta*(theta),
## implicit chain via Newton hessian for d eta*/d theta, and the sigma-driven
## contribution when the errfn depends on theta or eta.
.computeFoceiCorrection <- function(full_pars, joint_hessian, fixed,
                                    outer_names, H_inv_list,
                                    prdfn, errfn, omega,
                                    subjects, subject_etas, K, N,
                                    data_per_subject, times_union,
                                    conv2 = 2) {
  Q <- length(outer_names)
  pred <- prdfn(times = times_union, pars = full_pars, fixed = fixed,
                deriv = TRUE, deriv2 = TRUE, conditions = subjects)
  correction <- setNames(numeric(Q), outer_names)

  err_pred <- NULL
  err_par_names <- character()
  if (!is.null(errfn)) {
    err_pred <- vector("list", N); names(err_pred) <- subjects
    for (i_e in seq_len(N)) {
      cn_e   <- subjects[i_e]
      pred_e <- pred[[cn_e]]
      pinner     <- getParameters(pred_e)
      fixedinner <- pinner[attr(pinner, "fixed")]
      pinner     <- as.parvec(pinner[setdiff(names(pinner), names(fixed))])
      fixedinner <- as.parvec(fixedinner, deriv = FALSE, deriv2 = FALSE)
      err_pred[[cn_e]] <- errfn(out = pred_e, pars = pinner,
                                   fixed = fixedinner,
                                   conditions = cn_e)[[cn_e]]
    }
    err_par_names <- dimnames(attr(err_pred[[1]], "deriv"))[[3]]
  }
  err_in_outer <- intersect(outer_names, err_par_names)

  chol_in_outer <- intersect(outer_names, omega$cholPars)
  dOmega_inv_dchol <- list()
  if (length(chol_in_outer) > 0L) {
    chol_vec <- full_pars[omega$cholPars]
    L0       <- omega$buildL(chol_vec)
    h_chol   <- sqrt(.Machine$double.eps) * pmax(abs(chol_vec), 1)
    for (cp in chol_in_outer) {
      cv_p <- chol_vec; cv_p[cp] <- cv_p[cp] + h_chol[cp]
      cv_m <- chol_vec; cv_m[cp] <- cv_m[cp] - h_chol[cp]
      Lp <- omega$buildL(cv_p); Lm <- omega$buildL(cv_m)
      Op_inv <- chol2inv(chol(Lp %*% t(Lp)))
      Om_inv <- chol2inv(chol(Lm %*% t(Lm)))
      dOmega_inv_dchol[[cp]] <- (Op_inv - Om_inv) / (2 * h_chol[cp])
      dimnames(dOmega_inv_dchol[[cp]]) <-
        list(omega$eta, omega$eta)
    }
  }

  for (i in seq_len(N)) {
    cn         <- subjects[i]
    eta_i_nm   <- subject_etas[i, ]
    d_i        <- data_per_subject[[cn]]
    pred_i     <- pred[[cn]]
    d_full     <- attr(pred_i, "deriv")
    d2_full    <- attr(pred_i, "deriv2")
    if (is.null(d_full) || is.null(d2_full))
      stop("`prdfn` did not produce 'deriv'/'deriv2' attributes for ",
           "condition '", cn, "'. Did you build it with deriv2-capable ",
           "Xs/Y/P?")

    avail_pars <- dimnames(d_full)[[3]]
    th_avail   <- intersect(outer_names, avail_pars)
    eta_avail  <- intersect(eta_i_nm,    avail_pars)
    if (!length(eta_avail)) next

    times_i  <- d_i$time
    names_i  <- as.character(d_i$name)
    sigma_i  <- d_i$sigma
    ti_e <- oi_e <- NULL
    if (!is.null(err_pred)) {
      err_i  <- err_pred[[cn]]
      pt     <- as.numeric(err_i[, "time"])
      ti_e   <- match.num(times_i, pt)
      ni_e   <- match(names_i, colnames(err_i))
      if (anyNA(ti_e) || anyNA(ni_e))
        stop("compute_correction: cannot align data point to errfn grid.")
      sigma_from_err <- err_i[cbind(ti_e, ni_e)]
      sigma_i <- ifelse(is.na(sigma_i), sigma_from_err, sigma_i)
      d_err  <- attr(err_i, "deriv")
      if (!is.null(d_err)) {
        oi_e <- match(names_i, dimnames(d_err)[[2]])
        if (anyNA(oi_e))
          stop("compute_correction: data name not in errfn deriv axis.")
      }
    }
    time_idx <- match(times_i, times_union)
    if (anyNA(time_idx))
      stop("compute_correction: data times missing from prediction grid.")

    Tn <- length(times_i)
    Kn <- length(eta_avail)
    Qn <- length(th_avail)

    G    <- matrix(0,    Tn, Kn, dimnames = list(NULL, eta_avail))
    Wmix <- array (0, c(Tn, Kn, Qn))
    Veta <- array (0, c(Tn, Kn, Kn))
    for (jr in seq_len(Tn)) {
      ti <- time_idx[jr]; nm <- names_i[jr]
      G[jr, ]      <- d_full [ti, nm, eta_avail]
      Wmix[jr, , ] <- d2_full[ti, nm, eta_avail, th_avail]
      Veta[jr, , ] <- d2_full[ti, nm, eta_avail, eta_avail]
    }
    Sinv <- 1 / sigma_i^2

    H_inv_full <- H_inv_list[[i]]
    H_inv      <- H_inv_full[eta_avail, eta_avail, drop = FALSE]

    HG   <- G %*% H_inv
    qf_t <- rowSums(HG * G)

    if (Qn > 0L) {
      explicit <- numeric(Qn)
      for (k in seq_len(Qn)) {
        Mk <- crossprod(G * Sinv, Wmix[, , k, drop = TRUE])
        Mk <- Mk + t(Mk)
        explicit[k] <- conv2 * sum(H_inv * Mk)
      }
      correction[th_avail] <- correction[th_avail] + explicit
    }

    if (length(chol_in_outer) > 0L) {
      eta_pos  <- match(eta_avail, eta_i_nm)
      eta_base <- omega$eta[eta_pos]
      for (cp in chol_in_outer) {
        dOmega_inv_sub <-
          dOmega_inv_dchol[[cp]][eta_base, eta_base, drop = FALSE]
        correction[cp] <- correction[cp] +
          conv2 * sum(H_inv * dOmega_inv_sub)
      }
    }

    if (!is.null(err_pred)) {
      err_i <- err_pred[[cn]]
      d_err <- attr(err_i, "deriv")
      if (!is.null(d_err)) {
        err_par_set  <- dimnames(d_err)[[3]]
        theta_in_err <- intersect(outer_names, err_par_set)
        for (q in theta_in_err) {
          dsigma_q_t <- d_err[cbind(ti_e, oi_e, match(q, err_par_set))]
          correction[q] <- correction[q] +
            (-4) * sum(dsigma_q_t * qf_t / sigma_i^3)
        }
      }
    }

    th_for_implicit <- intersect(outer_names, colnames(joint_hessian))
    if (length(th_for_implicit) > 0L) {
      cross_block <- joint_hessian[eta_avail, th_for_implicit, drop = FALSE]
      H_newton    <- joint_hessian[eta_avail, eta_avail, drop = FALSE]
      eigN          <- eigen(H_newton, symmetric = TRUE)
      trN           <- sum(eigN$values)
      epsN          <- 1e-10 * abs(trN) / Kn
      lamN          <- pmax(eigN$values, epsN)
      H_inv_Newton  <- eigN$vectors %*% (t(eigN$vectors) / lamN)
      dimnames(H_inv_Newton) <- dimnames(H_newton)

      deta_dth <- -H_inv_Newton %*% cross_block

      dlogH_deta <- numeric(Kn)
      for (l in seq_len(Kn)) {
        Ml <- crossprod(G * Sinv, Veta[, , l, drop = TRUE])
        Ml <- Ml + t(Ml)
        dlogH_deta[l] <- conv2 * sum(H_inv * Ml)
      }
      if (!is.null(err_pred)) {
        err_i <- err_pred[[cn]]
        d_err <- attr(err_i, "deriv")
        if (!is.null(d_err)) {
          err_par_set <- dimnames(d_err)[[3]]
          for (l in seq_len(Kn)) {
            et <- eta_avail[l]
            if (et %in% err_par_set) {
              dsigma_et_t <- d_err[cbind(ti_e, oi_e,
                                          match(et, err_par_set))]
              dlogH_deta[l] <- dlogH_deta[l] +
                (-4) * sum(dsigma_et_t * qf_t / sigma_i^3)
            }
          }
        }
      }
      implicit <- as.numeric(dlogH_deta %*% deta_dth)
      correction[th_for_implicit] <-
        correction[th_for_implicit] + implicit
    }
  }
  correction
}



#' Build a quadrature-based NLME marginal-likelihood objective function
#'
#' Returns a callable `em(pars)` that integrates out the per-subject random
#' effects using a sparse-grid Gauss-Hermite rule. The ECM E-step refreshes
#' the integration nodes between outer iterations via
#' `attr(em, "rebuildQuadrature")`. Most users want [nlmeFit], which builds
#' the right `em` internally and runs a solver; use `emObjfn` directly only
#' if you need to evaluate the objective by hand.
#'
#' The Laplace approximation that was historically available here as
#' `method = "laplace"` has been folded into the C++ FOCEI kernel exposed by
#' [nlmeFit] with `method = "focei"`. Call `nlmeFit` directly.
#'
#' @param obj An \code{objfn}, typically
#'   `normL2(data, g*x*p, errmodel = err) + constraintL2(mu = 0, Omega = om)`.
#' @param omega An [omega] spec with subject expansion.
#' @param prdfn A `prdfn` (`g*x*p`). Required.
#' @param data A [datalist]. Required.
#' @param errfn Optional obsfn defining a parameter-dependent error model.
#' @param method Character. Only `"quadrature"` is supported. Retained for
#'   forward compatibility; the historical `"laplace"` path now lives inside
#'   the C++ FOCEI kernel reached via [nlmeFit] with `method = "focei"`.
#' @param control Named list with `level` (Smolyak depth, default 4) and
#'   `cores` (default 1).
#'
#' @return A callable of class `c("emObjfn", "objfn", "fn")`. Calling it on
#'   `pars` returns an [objlist] with an `emDiag` attribute carrying
#'   quadrature diagnostics.
#'
#' @seealso [nlmeFit], [omega]
#' @export
emObjfn <- function(obj, omega,
                    prdfn    = NULL,
                    data     = NULL,
                    errfn = NULL,
                    method   = "quadrature",
                    control  = list()) {
  if (!identical(method, "quadrature"))
    stop("emObjfn(): only `method = \"quadrature\"` is supported. The Laplace ",
         "approximation was folded into the C++ FOCEI kernel; call ",
         "`nlmeFit(method = \"focei\")` directly.")
  .emObjfn_quadrature(obj, omega,
                      prdfn    = prdfn,
                      data     = data,
                      errfn = errfn,
                      level    = control$level %||% 4L,
                      cores    = control$cores %||% 1L)
}



## Internal quadrature-method emObjfn constructor (Phase 4b).
##
## Closure state separates frozen E-step (nodes_per_subject, etaModes,
## chol_value, current_level) from the trust-varying structural pars. The
## orchestrator (nlmeFit) calls attr(em, "rebuildQuadrature")(psiFull, level)
## to update the frozen state, then runs trust(em, init = psi_structural) to
## step over structural params while the integration grid stays fixed.
.emObjfn_quadrature <- function(obj, omega, prdfn, data, errfn,
                                level, cores) {
  if (!inherits(obj, "objfn"))
    stop("`obj` must be an objfn.")
  if (!inherits(omega, "omegaSpec"))
    stop("`omega` must be built by omega().")
  if (is.null(omega$subjectEtas))
    stop("`omega` must have subject expansion (call omega(..., subjects = ...)).")
  if (is.null(prdfn))
    stop("`prdfn` (the prdfn used to build `obj`) is required for ",
         "method = \"quadrature\".")
  if (is.null(data))
    stop("`data` (the datalist used for obj) is required for ",
         "method = \"quadrature\".")

  K            <- omega$K
  subject_etas <- omega$subjectEtas
  subjects     <- rownames(subject_etas)
  N            <- length(subjects)
  all_eta_nm   <- as.vector(subject_etas)
  chol_pars    <- omega$cholPars
  joint_pars   <- attr(obj, "parameters")

  # Frozen E-step state (mutable via <<- inside rebuildQuadrature only).
  nodes_per_subject <- NULL
  eta_modes_state   <- NULL
  H_i_list_state    <- NULL
  chol_value_state  <- NULL
  current_level     <- as.integer(level)
  n_floored_state   <- integer(N)
  converged_state   <- logical(N)
  iter_state        <- integer(N)

  # Mode-finder for one subject. Per-subject inner trust over eta, with the
  # other subjects' etas held at their cached values.
  find_mode_one <- function(subjIdx, outer_input_with_chol, fixed,
                            eta_init_full) {
    eta_i_names <- subject_etas[subjIdx, ]
    eta_i_init  <- setNames(eta_init_full[subjIdx, ], eta_i_names)
    other_nm    <- as.vector(subject_etas[-subjIdx, , drop = FALSE])
    other_vals  <- as.vector(eta_init_full[-subjIdx, , drop = FALSE])
    names(other_vals) <- other_nm

    inner_objfn <- function(eta_i_in, ...) {
      full_pars <- c(outer_input_with_chol, other_vals, eta_i_in)
      out <- obj(pars = full_pars, fixed = fixed, deriv = TRUE,
                   conditions = subjects[subjIdx])
      gr <- out$gradient[eta_i_names]
      hs <- out$hessian[eta_i_names, eta_i_names, drop = FALSE]
      objlist(value = out$value, gradient = gr, hessian = hs)
    }

    fit <- try(suppressMessages(trust(inner_objfn, parinit = eta_i_init,
                                      rinit = 1, rmax = 10, iterlim = 50L,
                                      fterm = 1e-7, mterm = 1e-7)), silent = TRUE)
    if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
      H_fallback <- diag(K)
      dimnames(H_fallback) <- list(eta_i_names, eta_i_names)
      list(etaStar  = eta_i_init,
           H_i       = H_fallback,
           iter      = NA_integer_,
           converged = FALSE)
    } else {
      list(etaStar  = setNames(as.numeric(fit$argument), eta_i_names),
           H_i       = fit$hessian,
           iter      = fit$iterations,
           converged = isTRUE(fit$converged))
    }
  }

  # E-step: rebuild per-subject modes, eigen-floored H_i, and Smolyak nodes.
  rebuildQuadrature <- function(psiFull, level_new = NULL,
                                  fixed = NULL,
                                  eta_init = NULL) {
    if (!is.null(level_new)) current_level <<- as.integer(level_new)
    if (!all(chol_pars %in% names(psiFull)))
      stop("rebuildQuadrature: `psiFull` is missing omega$cholPars.")
    outer_with_chol_names <- intersect(names(psiFull),
                                       setdiff(joint_pars, all_eta_nm))
    outer_with_chol <- psiFull[outer_with_chol_names]

    if (is.null(eta_init)) {
      eta_init <- if (!is.null(eta_modes_state)) eta_modes_state
                  else matrix(0, N, K, dimnames = dimnames(subject_etas))
    }

    new_modes <- matrix(0, N, K, dimnames = dimnames(subject_etas))
    new_H     <- vector("list", N)
    new_nodes <- vector("list", N)
    new_nfloor <- integer(N)
    new_conv   <- logical(N)
    new_iter   <- integer(N)

    for (i in seq_len(N)) {
      m <- find_mode_one(i, outer_with_chol, fixed, eta_init)
      new_modes[i, ] <- m$etaStar
      new_conv[i]    <- m$converged
      new_iter[i]    <- if (is.null(m$iter)) NA_integer_ else m$iter
      H_i <- m$H_i
      eig <- eigen(H_i, symmetric = TRUE)
      tr_H <- sum(eig$values)
      eps <- 1e-10 * abs(tr_H) / K
      lambda_floored <- pmax(eig$values, eps)
      new_nfloor[i] <- sum(eig$values < eps)
      H_safe <- eig$vectors %*% (lambda_floored * t(eig$vectors))
      dimnames(H_safe) <- dimnames(H_i)
      new_H[[i]] <- H_safe
      new_nodes[[i]] <- makeSubjectNodes(new_modes[i, ], H_safe, current_level)
    }

    nodes_per_subject <<- new_nodes
    eta_modes_state   <<- new_modes
    H_i_list_state    <<- new_H
    chol_value_state  <<- psiFull[chol_pars]
    n_floored_state   <<- new_nfloor
    converged_state   <<- new_conv
    iter_state        <<- new_iter

    invisible(list(etaModes  = new_modes,
                   HiList   = new_H,
                   level      = current_level,
                   nFloored  = new_nfloor,
                   converged  = new_conv,
                   iterations = new_iter))
  }

  myfn <- function(..., fixed = NULL, deriv = TRUE, env = NULL) {
    p <- list(...)[[match.fnargs(list(...), "pars")]]
    if (is.null(env)) env <- new.env()

    if (is.null(nodes_per_subject))
      stop("emObjfn: no quadrature grid built. Call ",
           "`attr(em, 'rebuildQuadrature')(psiFull, level)` first.")

    # Build psiFull from p + frozen chol_value_state (CM-1 holds chol fixed).
    structural_names <- setdiff(names(p), c(chol_pars, all_eta_nm))
    chol_in_p <- intersect(chol_pars, names(p))
    psiFull <- c(p[structural_names],
                  if (length(chol_in_p)) p[chol_in_p] else chol_value_state)

    # Outer-active params = those varied by the trust caller = names(p).
    outer_active <- structural_names

    if (cores == 1L) {
      per_subj <- lapply(seq_len(N), function(i) ecmEvaluateSubject(
        i, psiFull, eta_modes_state, omega, nodes_per_subject[[i]],
        xPred = prdfn, datalist = data, errfn = errfn,
        fixed = fixed, outerActiveNames = outer_active,
        mode = if (deriv) "with_grad" else "moments_only"))
    } else if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      per_subj <- parallel::parLapply(cl, seq_len(N), function(i) ecmEvaluateSubject(
        i, psiFull, eta_modes_state, omega, nodes_per_subject[[i]],
        xPred = prdfn, datalist = data, errfn = errfn,
        fixed = fixed, outerActiveNames = outer_active,
        mode = if (deriv) "with_grad" else "moments_only"))
    } else {
      per_subj <- parallel::mclapply(seq_len(N), function(i) ecmEvaluateSubject(
        i, psiFull, eta_modes_state, omega, nodes_per_subject[[i]],
        xPred = prdfn, datalist = data, errfn = errfn,
        fixed = fixed, outerActiveNames = outer_active,
        mode = if (deriv) "with_grad" else "moments_only"),
        mc.cores = cores)
    }

    tot_logL <- 0
    tot_gr   <- setNames(numeric(length(p)), names(p))
    tot_he   <- matrix(0, length(p), length(p),
                       dimnames = list(names(p), names(p)))
    MHatList      <- vector("list", N)
    mHatList      <- vector("list", N)
    max_softmax_per <- numeric(N)
    n_eff_per       <- numeric(N)

    for (i in seq_len(N)) {
      o <- per_subj[[i]]
      tot_logL          <- tot_logL + o$logLhat
      MHatList[[i]]   <- o$M_hat
      mHatList[[i]]   <- o$m_hat
      max_softmax_per[i] <- o$maxSoftmax
      n_eff_per[i]      <- o$n_eff
      if (deriv && !is.null(o$gradient)) {
        nm <- intersect(outer_active, names(tot_gr))
        tot_gr[nm] <- tot_gr[nm] + o$gradient[nm]
        tot_he[nm, nm] <- tot_he[nm, nm] + o$hessian[nm, nm]
      }
    }

    OFV <- -2 * tot_logL
    out <- objlist(value = OFV, gradient = tot_gr, hessian = tot_he)
    emDiag <- list(method          = "quadrature",
                    level           = current_level,
                    etaModes       = eta_modes_state,
                    HiList        = H_i_list_state,
                    mHatList      = mHatList,
                    MHatList      = MHatList,
                    maxSoftmax     = max_softmax_per,
                    n_eff           = n_eff_per,
                    nFloored       = n_floored_state,
                    innerConverged = converged_state,
                    innerIter      = iter_state)
    attr(out, "emDiag") <- emDiag
    attr(out, "env")     <- env
    out
  }

  class(myfn) <- c("emObjfn", "objfn", "fn")
  attr(myfn, "method")     <- "quadrature"
  attr(myfn, "conditions") <- attr(obj, "conditions")
  attr(myfn, "parameters") <- setdiff(joint_pars, all_eta_nm)
  attr(myfn, "omega")  <- omega
  attr(myfn, "obj")      <- obj
  attr(myfn, "prdfn")      <- prdfn
  attr(myfn, "data")       <- data
  attr(myfn, "errfn")   <- errfn
  attr(myfn, "rebuildQuadrature") <- rebuildQuadrature
  attr(myfn, "control")    <- list(level = level, cores = cores)
  myfn
}



# nlmeFit S3 constructor. Bundles solver output with the prdfn / data /
# omega references that predict.nlmeFit and the diagnostic plot
# helpers (in plots.R) consume.
nlmeFit_make <- function(argument, value, gradient, hessian, Omega, etaModes,
                         converged, iterations, emDiag, method,
                         foceiStart = NULL, stageTrace = NULL,
                         prdfn = NULL, data = NULL, omega = NULL,
                         errfn = NULL) {
  etaInfo  <- .computeEtaInfo(emDiag, Omega, etaModes)
  out <- list(argument   = argument,
              value      = value,
              gradient   = gradient,
              hessian    = hessian,
              Omega      = Omega,
              etaModes   = etaModes,
              etaSE      = etaInfo$etaSE,
              shrinkage  = etaInfo$shrinkage,
              converged  = converged,
              iterations = iterations,
              emDiag     = emDiag,
              method     = method,
              foceiStart = foceiStart,
              stageTrace = stageTrace,
              prdfn      = prdfn,
              data       = data,
              omega      = omega,
              errfn      = errfn)
  class(out) <- c("nlmeFit", "list")
  out
}

# Posterior-mode standard errors and shrinkage diagnostics for the per-subject
# random effects. Caller passes the full emDiag (which carries HInvList from
# either the R or C++ Laplace path); returns NULLs when the inverse Hessian
# list is unavailable (e.g. quadrature method).
.computeEtaInfo <- function(emDiag, Omega, etaModes) {
  out <- list(etaSE = NULL, shrinkage = NULL)
  if (is.null(emDiag) || is.null(emDiag$HInvList) || is.null(etaModes))
    return(out)
  H_inv_list <- emDiag$HInvList
  N <- nrow(etaModes); K <- ncol(etaModes)
  if (length(H_inv_list) != N) return(out)

  diag_template <- rep(NA_real_, K)
  diags <- vapply(H_inv_list, function(Hi) {
    if (is.null(Hi) || !is.matrix(Hi) || any(dim(Hi) != c(K, K)))
      return(diag_template)
    di <- diag(Hi)
    di[di < 0] <- NA_real_
    di
  }, numeric(K))
  dim(diags) <- c(K, N)
  etaSE <- sqrt(t(diags))
  dimnames(etaSE) <- dimnames(etaModes)
  out$etaSE <- etaSE

  if (!is.null(Omega) && all(dim(Omega) == c(K, K))) {
    Omega_sd <- sqrt(pmax(diag(Omega), 0))
    sd_mat <- matrix(Omega_sd, N, K, byrow = TRUE)
    shrink <- 1 - etaSE / sd_mat
    shrink[, Omega_sd <= 0] <- NA_real_
    dimnames(shrink) <- dimnames(etaModes)
    out$shrinkage <- shrink
  }
  out
}


# Synthetic err callable for the case where the user supplied no errfn
# and instead set `data$sigma` directly. Returns a prdlist whose matrix has
# the per-observation sigmas (from the data) padded onto the prediction's
# time grid. No `deriv` attribute, so the C++ kernel treats `dsigma/deta = 0`.
.makeStaticErr <- function(data) {
  data_per_subject <- lapply(seq_along(data), function(i) {
    d <- data[[i]]
    d$name <- as.character(d$name)
    d
  })
  names(data_per_subject) <- names(data)

  function(out, pars = NULL, conditions = NULL, ...) {
    s <- if (!is.null(conditions)) conditions[1] else names(data_per_subject)[1]
    out_mat <- if (inherits(out, "prdlist")) out[[s]] else out
    pred_times <- as.numeric(out_mat[, "time"])
    obs_names  <- setdiff(colnames(out_mat), "time")
    d <- data_per_subject[[s]]
    if (is.null(d))
      stop(".makeStaticErr: condition '", s,
           "' missing from data; supply an errfn.")
    if (anyNA(d$sigma))
      stop(".makeStaticErr: condition '", s,
           "' has NA in data$sigma; supply an errfn.")
    err_mat <- matrix(NA_real_, length(pred_times), length(obs_names) + 1L)
    colnames(err_mat) <- c("time", obs_names)
    err_mat[, "time"] <- pred_times
    for (o in obs_names) {
      drow <- d[d$name == o, , drop = FALSE]
      if (!nrow(drow)) next
      idx <- match(pred_times, drow$time)
      err_mat[, o] <- ifelse(is.na(idx), median(drow$sigma), drow$sigma[idx])
    }
    setNames(list(structure(err_mat, class = c("prdframe", "matrix"))), s)
  }
}


# Builds the long-format subject_meta consumed by focei_inner_trust / the
# C++ kernel. One row in t_idx_in_pred / o_idx_in_pred / y_data / ... per
# observed data point of the subject, covering all observables. eta_idx_in_*
# arrays are length K (one entry per random effect). Values 0 mark "no
# contribution" (e.g. an eta does not appear in the err prdfn's deriv).
.buildFastMeta <- function(prdfn, errfn, data, subjects,
                           eta_names_list, pars_full_names, pars_probe) {
  N <- length(subjects)
  fast_meta <- vector("list", N)
  for (i in seq_len(N)) {
    s      <- subjects[i]
    data_i <- data[[s]]
    data_i$name <- as.character(data_i$name)
    data_i <- data_i[order(data_i$time, data_i$name), , drop = FALSE]
    times_union <- sort(unique(data_i$time))

    # Probe prdfn + err once to learn deriv array structure.
    pred_probe <- prdfn(times = times_union, pars = pars_probe,
                        deriv = TRUE, conditions = s)
    prdf  <- pred_probe[[1]]
    pcols <- colnames(prdf)
    d_dn  <- dimnames(attr(prdf, "deriv"))
    eta_names_i      <- eta_names_list[[i]]
    eta_idx_in_deriv <- match(eta_names_i, d_dn[[3]])
    if (anyNA(eta_idx_in_deriv))
      stop(".buildFastMeta: eta names not present in prdfn deriv array for ",
           "subject ", s, ".")

    pinner    <- attr(prdf, "parameters")
    err_probe <- errfn(out = prdf, pars = pinner, conditions = s)
    erm       <- err_probe[[1]]
    ecols     <- colnames(erm)
    e_attr    <- attr(erm, "deriv")
    if (!is.null(e_attr)) {
      e_dn               <- dimnames(e_attr)
      eta_idx_in_err_deriv <- match(eta_names_i, e_dn[[3]])
      eta_idx_in_err_deriv[is.na(eta_idx_in_err_deriv)] <- 0L
    } else {
      e_dn               <- NULL
      eta_idx_in_err_deriv <- rep(0L, length(eta_names_i))
    }

    # Long-format per-row indices.
    t_idx_in_pred  <- match(data_i$time, prdf[, "time"])
    o_idx_in_pred  <- match(data_i$name, pcols)
    if (anyNA(o_idx_in_pred))
      stop(".buildFastMeta: observable(s) not found in prdfn output for ",
           "subject ", s, ": ",
           paste(setdiff(unique(data_i$name), pcols), collapse = ", "))
    o_idx_in_deriv <- match(data_i$name, d_dn[[2]])
    t_idx_in_err   <- match(data_i$time, erm[, "time"])
    o_idx_in_err   <- match(data_i$name, ecols)
    if (!is.null(e_dn)) {
      o_idx_in_err_deriv <- match(data_i$name, e_dn[[2]])
      o_idx_in_err_deriv[is.na(o_idx_in_err_deriv)] <- 0L
    } else {
      o_idx_in_err_deriv <- rep(0L, nrow(data_i))
    }

    fast_meta[[i]] <- list(
      times                = as.numeric(times_union),
      eta_idx_in_pars      = as.integer(match(eta_names_i, pars_full_names)),
      t_idx_in_pred        = as.integer(t_idx_in_pred),
      y_data               = as.numeric(data_i$value),
      o_idx_in_pred        = as.integer(o_idx_in_pred),
      eta_idx_in_deriv     = as.integer(eta_idx_in_deriv),
      o_idx_in_deriv       = as.integer(o_idx_in_deriv),
      t_idx_in_err         = as.integer(t_idx_in_err),
      o_idx_in_err         = as.integer(o_idx_in_err),
      eta_idx_in_err_deriv = as.integer(eta_idx_in_err_deriv),
      o_idx_in_err_deriv   = as.integer(o_idx_in_err_deriv),
      condition            = s,
      eta_names            = eta_names_i)
  }
  fast_meta
}


# Internal: run the C++ FOCEI kernel and package its output as nlmeFit.
# Used by nlmeFit(method = "focei"). Always uses the fast-inner C++ path
# with eager Stage-2 correction, calling .computeFoceiCorrection as an
# Rcpp::Function once per outer iter.
.runFoceiCpp <- function(obj, omega, init, prdfn, data, errfn,
                         fixed = NULL,
                         innerControl = list(), trustControl = list(),
                         methodLabel = "focei") {
  if (is.null(prdfn))
    stop(".runFoceiCpp: `prdfn` (the prdfn) is required.")
  if (is.null(data))
    stop(".runFoceiCpp: `data` (the datalist) is required.")
  if (is.null(omega$subjectEtas))
    stop(".runFoceiCpp: omega has no subject expansion. Call ",
         "omega(..., subjects = ...).")

  # When the user did not supply an errfn but recorded sigma in the data
  # itself, wrap the per-row sigmas into a synthetic obsfn-like callable.
  # The fast-inner kernel treats this as "sigma constant in eta" (no err
  # deriv attribute -> Js = 0).
  if (is.null(errfn)) errfn <- .makeStaticErr(data)

  K        <- omega$K
  subjects <- rownames(omega$subjectEtas)
  N        <- length(subjects)
  outer_names <- names(init)
  eta_names_all <- as.vector(omega$subjectEtas)
  pars_full_names <- c(outer_names, eta_names_all)
  eta_idx_global <- matrix(0L, N, K)
  eta_names_list <- vector("list", N)
  for (i in seq_len(N)) {
    nms <- omega$subjectEtas[i, ]
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
    pars_full_names   = pars_full_names
  )
  ic <- modifyList(list(rinit = 1, rmax = 10, iterlim = 30,
                        fterm = 1e-7, mterm = 1e-7,
                        eigen_floor_relative = 1e-10),
                   innerControl)
  oc <- modifyList(list(rinit = 1, rmax = 10, iterlim = 200,
                        fterm = 1e-7, mterm = 1e-7),
                   trustControl)
  control_cpp <- list(inner = ic, outer = oc)

  pars_probe <- setNames(numeric(length(pars_full_names)), pars_full_names)
  pars_probe[outer_names] <- init
  fast_meta <- .buildFastMeta(prdfn, errfn, data, subjects,
                              eta_names_list, pars_full_names, pars_probe)
  om_meta <- list(
    chol_pars = omega$cholPars,
    chol_loc  = matrix(as.integer(omega$cholLoc), ncol = 2,
                       dimnames = NULL),
    is_diag   = as.logical(omega$isDiag))
  subject_meta$fast_meta  <- fast_meta
  subject_meta$omega_meta <- om_meta

  data_per_subject <- lapply(subjects, function(s) data[[s]])
  names(data_per_subject) <- subjects
  times_union <- sort(unique(c(0, unlist(lapply(data_per_subject, `[[`,
                                                "time")))))
  correction_cb <- function(full_pars, joint_hessian, H_inv_list) {
    .computeFoceiCorrection(
      full_pars = full_pars, joint_hessian = joint_hessian,
      fixed = fixed, outer_names = outer_names, H_inv_list = H_inv_list,
      prdfn = prdfn, errfn = errfn, omega = omega,
      subjects = subjects, subject_etas = omega$subjectEtas,
      K = K, N = N,
      data_per_subject = data_per_subject, times_union = times_union)
  }

  fit <- focei_run(model_cb = prdfn, err_cb = errfn,
                        joint_cb = obj, init = init,
                        subject_meta = subject_meta,
                        fixed = fixed, control = control_cpp,
                        correction_mode = "eager",
                        correction_cb = correction_cb)

  L_omega <- if (all(omega$cholPars %in% names(fit$argument)))
    omega$buildL(fit$argument[omega$cholPars]) else NULL
  Omega <- if (!is.null(L_omega)) tcrossprod(L_omega) else NULL
  etaModes <- fit$etaModes
  rownames(etaModes) <- subjects
  colnames(etaModes) <- omega$eta
  emDiag <- list(etaStar = etaModes, logdet = fit$log_det_H,
                 sum_logdetH = fit$sum_logdetH, trace = fit$trace,
                 backend = "cpp", HInvList = fit$H_inv)
  nlmeFit_make(argument = fit$argument, value = fit$value,
               gradient = fit$gradient, hessian = fit$hessian,
               Omega = Omega, etaModes = etaModes,
               converged = isTRUE(fit$converged),
               iterations = fit$iterations, emDiag = emDiag,
               method = methodLabel,
               prdfn = prdfn, data = data,
               omega = omega, errfn = errfn)
}


# Internal: run ECM (E-step / CM-1 / CM-2) polish on a quadrature-method
# emObjfn, package as nlmeFit.
.runQuadratureEcm <- function(em, init, fixed = NULL, foceiStart = NULL,
                              epsQuadLevels = NULL, epsEcm = 1e-4,
                              epsOfvRel = 1e-5, maxEcmPerStage = 5L,
                              maxCm1Iter = 30L, cm1Control = list(),
                              methodLabel = "quadrature", verbose = TRUE) {
  om <- attr(em, "omega")
  K  <- om$K
  chol_pars <- om$cholPars
  if (is.null(epsQuadLevels)) epsQuadLevels <- K + 1:3
  cm1 <- modifyList(list(rinit = 1, rmax = 10,
                         fterm = 1e-6, mterm = 1e-6), cm1Control)
  psi <- init
  if (!all(chol_pars %in% names(psi)))
    stop("nlmeFit: `init` is missing omega$cholPars (",
         paste(setdiff(chol_pars, names(psi)), collapse = ", "), ").")
  structural_names <- setdiff(names(psi), chol_pars)
  rebuild <- attr(em, "rebuildQuadrature")
  stage_rows <- list(); prev_psi <- psi; prev_ofv <- NA_real_

  for (stage in seq_along(epsQuadLevels)) {
    level <- epsQuadLevels[stage]
    if (verbose) message(sprintf("nlmeFit(%s): stage %d / %d (level = %d)",
                                 methodLabel, stage,
                                 length(epsQuadLevels), level))
    e_info <- rebuild(psi, level_new = level, fixed = fixed,
                      eta_init = if (!is.null(foceiStart))
                                   foceiStart$emDiag$etaStar else NULL)
    for (ecmIter in seq_len(maxEcmPerStage)) {
      if (ecmIter > 1L)
        e_info <- rebuild(psi, level_new = level, fixed = fixed,
                          eta_init = e_info$etaModes)
      cm1_fit <- suppressMessages(trust(
        em, parinit = psi[structural_names], fixed = fixed,
        rinit = cm1$rinit, rmax = cm1$rmax,
        iterlim = maxCm1Iter,
        fterm = cm1$fterm, mterm = cm1$mterm, on_step = NULL))
      psi[structural_names] <- cm1_fit$argument
      out_after_cm1 <- em(psi[structural_names], fixed = fixed, deriv = FALSE)
      diag_after    <- attr(out_after_cm1, "emDiag")
      psi[chol_pars] <- updateOmegaChol(diag_after$MHatList, om)
      ofv         <- out_after_cm1$value
      deltaPsi    <- max(abs(psi - prev_psi))
      deltaOfvRel <- if (is.na(prev_ofv) || abs(prev_ofv) < .Machine$double.eps)
                       Inf else abs(ofv - prev_ofv) / abs(prev_ofv)
      stage_rows[[length(stage_rows) + 1L]] <- data.frame(
        stage = stage, ecmIter = ecmIter, level = level,
        OFV = ofv, deltaPsi = deltaPsi, deltaOfvRel = deltaOfvRel,
        maxSoftmax = max(diag_after$maxSoftmax),
        nEffMin = min(diag_after$n_eff),
        cm1TrustIter = cm1_fit$iterations)
      if (verbose) message(sprintf(
        "  ecm %d : OFV=%.6f  |dpsi|=%.2e  |dOFV/OFV|=%.2e  max_smax=%.3f  nEffMin=%.1f",
        ecmIter, ofv, deltaPsi, deltaOfvRel,
        max(diag_after$maxSoftmax), min(diag_after$n_eff)))
      prev_psi <- psi; prev_ofv <- ofv
      if (deltaPsi < epsEcm || deltaOfvRel < epsOfvRel) break
    }
  }

  rebuild(psi, fixed = fixed,
          eta_init = if (length(stage_rows))
                       attr(em(psi[structural_names], deriv = FALSE),
                            "emDiag")$etaModes else NULL)
  final_out  <- em(psi[structural_names], fixed = fixed)
  final_diag <- attr(final_out, "emDiag")
  L_omega    <- om$buildL(psi[chol_pars])
  Omega      <- tcrossprod(L_omega)
  conv       <- length(stage_rows) > 0L && {
    last <- tail(stage_rows, 1)[[1]]
    last$deltaPsi < epsEcm || last$deltaOfvRel < epsOfvRel
  }
  nlmeFit_make(argument   = psi,
               value      = final_out$value,
               gradient   = final_out$gradient,
               hessian    = final_out$hessian,
               Omega      = Omega,
               etaModes   = final_diag$etaModes,
               converged  = conv,
               iterations = length(stage_rows),
               emDiag     = final_diag,
               method     = methodLabel,
               foceiStart = foceiStart,
               stageTrace = do.call(rbind, stage_rows),
               prdfn      = attr(em, "prdfn"),
               data       = attr(em, "data"),
               omega      = attr(em, "omega"),
               errfn      = attr(em, "errfn"))
}


#' Fit a nonlinear mixed-effects prdfn
#'
#' Builds the marginal-likelihood objective via [emObjfn] and runs the
#' selected estimator. Returns an `nlmeFit` S3 object consumable by the
#' diagnostic helpers ([predict.nlmeFit], [plot.nlmeFit], [plotIndivs] etc.).
#'
#' @param obj An \code{objfn} of the form
#'   `normL2(data, prdfn, errmodel = err) + constraintL2(mu = 0, Omega = om)`.
#' @param omega An [omega] spec with subject expansion.
#' @param init Named numeric starting parameter vector. Must contain all
#'   structural parameters and all `omega$cholPars`.
#' @param prdfn The prediction function `g * x * p` used to build `obj`.
#'   Required.
#' @param data The [datalist] used for `obj`. Required.
#' @param errfn Optional obsfn defining a parameter-dependent error model.
#' @param fixed Optional named-numeric of fixed parameters.
#' @param method Estimator. \code{"focei"} runs the C++ FOCEI kernel
#'   (Laplace + trust + eager Stage-2 correction); \code{"quadrature"} runs
#'   adaptive sparse-grid Gauss-Hermite + ECM with a cold start;
#'   \code{"foceiQuadrature"} runs FOCEI first and uses the converged
#'   structural pars and modes as warmstart for the quadrature polish.
#' @param control Nested list of method-specific knobs. Entries:
#'   \describe{
#'     \item{`$focei`}{Recognised keys: `innerControl`, `trustControl`.}
#'     \item{`$quadrature`}{Passed to the quadrature [emObjfn] and ECM solver.
#'       Recognised keys: `level`, `cores`, `epsQuadLevels`, `epsEcm`,
#'       `epsOfvRel`, `maxEcmPerStage`, `maxCm1Iter`, `cm1Control`.}
#'   }
#' @param verbose Logical. If TRUE prints solver progress.
#'
#' @return An `nlmeFit` S3 list with fields `argument`, `value`, `gradient`,
#'   `hessian`, `omega`, `etaModes`, `converged`, `iterations`, `emDiag`,
#'   `method`, `foceiStart`, `stageTrace`, `prdfn`, `data`, `omega`,
#'   `errfn`.
#'
#' @seealso [emObjfn], [omega], [predict.nlmeFit]
#' @export
nlmeFit <- function(obj, omega, init,
                    prdfn    = NULL,
                    data     = NULL,
                    errfn = NULL,
                    fixed    = NULL,
                    method   = c("focei", "quadrature", "foceiQuadrature"),
                    control  = list(),
                    verbose  = TRUE) {
  method <- match.arg(method)
  fc <- control$focei      %||% list()
  qc <- control$quadrature %||% list()

  if (method == "focei") {
    return(.runFoceiCpp(obj, omega, init,
                        prdfn = prdfn, data = data, errfn = errfn,
                        fixed = fixed,
                        innerControl = fc$innerControl %||% list(),
                        trustControl = fc$trustControl %||% list(),
                        methodLabel  = "focei"))
  }

  if (method == "quadrature") {
    em <- emObjfn(obj, omega, prdfn = prdfn, data = data,
                  errfn = errfn, control = qc)
    return(.runQuadratureEcm(em, init, fixed = fixed,
                             foceiStart     = NULL,
                             epsQuadLevels  = qc$epsQuadLevels,
                             epsEcm         = qc$epsEcm         %||% 1e-4,
                             epsOfvRel      = qc$epsOfvRel      %||% 1e-5,
                             maxEcmPerStage = qc$maxEcmPerStage %||% 5L,
                             maxCm1Iter     = qc$maxCm1Iter     %||% 30L,
                             cm1Control     = qc$cm1Control     %||% list(),
                             methodLabel    = "quadrature",
                             verbose        = verbose))
  }

  # foceiQuadrature: FOCEI warmstart + quadrature polish.
  if (verbose) message("nlmeFit: running FOCEI warmstart ...")
  foceiStart <- .runFoceiCpp(obj, omega, init,
                             prdfn = prdfn, data = data, errfn = errfn,
                             fixed = fixed,
                             innerControl = fc$innerControl %||% list(),
                             trustControl = fc$trustControl %||% list(),
                             methodLabel  = "focei")
  if (verbose) message(sprintf("  warmstart OFV = %.6f", foceiStart$value))
  em_qd <- emObjfn(obj, omega, prdfn = prdfn, data = data,
                   errfn = errfn, control = qc)
  .runQuadratureEcm(em_qd, foceiStart$argument, fixed = fixed,
                    foceiStart     = foceiStart,
                    epsQuadLevels  = qc$epsQuadLevels,
                    epsEcm         = qc$epsEcm         %||% 1e-4,
                    epsOfvRel      = qc$epsOfvRel      %||% 1e-5,
                    maxEcmPerStage = qc$maxEcmPerStage %||% 5L,
                    maxCm1Iter     = qc$maxCm1Iter     %||% 30L,
                    cm1Control     = qc$cm1Control     %||% list(),
                    methodLabel    = "foceiQuadrature",
                    verbose        = verbose)
}


#' Print an nlmeFit object
#'
#' @param x An `nlmeFit` object (see [nlmeFit]).
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.nlmeFit <- function(x, ...) {
  cat("nlmeFit (method = ", x$method, ")\n", sep = "")
  cat(sprintf("  OFV (-2 log L): %.6f\n", x$value))
  cat(sprintf("  converged    : %s   iterations: %s\n",
              x$converged, format(x$iterations %||% NA_integer_)))
  cat("  argument     :\n")
  print(x$argument)
  if (!is.null(x$Omega)) {
    cat("  Omega        :\n")
    print(round(x$Omega, 4))
  }
  if (!is.null(x$etaModes) && !is.null(x$etaSE)) {
    cat("  eta (mode +/- SE, shrinkage):\n")
    print(.formatEtaTable(x$etaModes, x$etaSE, x$shrinkage))
  }
  invisible(x)
}

# Build a numeric data.frame with one column per quantity (mode, SE,
# optional shrinkage) per eta. Relies on print.data.frame for alignment.
.formatEtaTable <- function(etaModes, etaSE, shrinkage) {
  K <- ncol(etaModes)
  eta_names <- colnames(etaModes)
  cols <- c(rbind(eta_names, paste0("SE.", eta_names)))
  vals <- cbind(etaModes, etaSE)[, c(rbind(seq_len(K), seq_len(K) + K)),
                                 drop = FALSE]
  colnames(vals) <- cols

  if (!is.null(shrinkage)) {
    shr <- shrinkage
    colnames(shr) <- paste0("shrink.", eta_names)
    vals <- cbind(vals, shr)
  }
  round(as.data.frame(vals), 3)
}



#' Per-subject adaptive-quadrature marginal-likelihood evaluator
#'
#' @description
#' Evaluates the per-subject marginal likelihood
#' \eqn{\hat L_i = \int p(y_i \mid \eta) \, N(\eta \mid 0, \Omega)\, d\eta}
#' via sparse-grid Gauss-Hermite quadrature at a precomputed set of nodes
#' (output of [makeSubjectNodes]). Returns the log-likelihood, posterior
#' first and second moments of \eqn{\eta_i}, plus (optionally) the gradient
#' and Hessian of \eqn{-2 \log \hat L_i} with respect to a chosen subset of
#' outer (population) parameters.
#'
#' This helper bypasses [normL2] / [constraintL2] entirely: the data
#' likelihood contribution comes from the lifted [evalConditionResidual]
#' (one prediction call per node, single-condition), the MVN prior on the
#' subject's \eqn{\eta_b} is added in closed form via `omega$buildL`.
#' Avoids the multi-condition parameter rebinding and full-population MVN
#' contribution that calling `obj(..., conditions = subjects[i])` per node
#' would incur.
#'
#' @param subjIdx Integer in `1..N` selecting the subject.
#' @param psiFull Named numeric, the full outer parameter vector at which to
#'   evaluate (structural + chol_pars; chol_pars are frozen during CM-1).
#' @param etaModes N x K matrix of all subjects' eta values; the row at
#'   `subjIdx` is ignored. Other rows are passed through to the prediction
#'   call (required for parameter completeness in joint models).
#' @param omega An [omega] spec object with subject expansion.
#' @param nodesSubj Output of `makeSubjectNodes(eta_hat_i, H_i, level)` for
#'   the active subject.
#' @param xPred The prediction function (e.g. `g * x * p`).
#' @param datalist The [datalist] used for the joint objective.
#' @param errfn Optional obsfn (passed through to
#'   [evalConditionResidual]).
#' @param fixed Optional fixed-parameter vector.
#' @param outerActiveNames Character vector of parameter names to track
#'   gradient/Hessian for. Defaults to `names(psiFull)`. For CM-1 with frozen
#'   chol_pars, pass just the structural names.
#' @param mode `"moments_only"` (no gradient/Hessian, cheap) or
#'   `"with_grad"` (gradient + Hessian of `-2 log L_i`).
#'
#' @return A list with components:
#' \describe{
#'   \item{`logLhat`}{`log L_i` (scalar).}
#'   \item{`m_hat`}{Length-K named numeric, posterior mean of \eqn{\eta_i}.}
#'   \item{`M_hat`}{K x K named matrix, posterior 2nd moment
#'     \eqn{\hat E[\eta_i \eta_i^T \mid y_i]}.}
#'   \item{`maxSoftmax`}{Max |softmax weight| across nodes; diagnostic of
#'     grid concentration.}
#'   \item{`n_eff`}{Effective node count, 1/sum(softmax^2); diagnostic.}
#'   \item{`gradient`, `hessian`}{When `mode = "with_grad"`: gradient and
#'     Hessian of `-2 log L_i` w.r.t. `outerActiveNames`. NULL otherwise.}
#' }
#'
#' @seealso [makeSubjectNodes], [evalConditionResidual]
#' @keywords internal
ecmEvaluateSubject <- function(subjIdx, psiFull, etaModes,
                                 omega, nodesSubj,
                                 xPred, datalist,
                                 errfn           = NULL,
                                 fixed              = NULL,
                                 outerActiveNames = NULL,
                                 mode               = c("moments_only", "with_grad")) {
  mode <- match.arg(mode)
  with_grad <- (mode == "with_grad")

  K            <- omega$K
  subject_etas <- omega$subjectEtas
  subjects     <- rownames(subject_etas)
  cn           <- subjects[subjIdx]
  eta_i_names  <- subject_etas[subjIdx, ]
  chol_pars    <- omega$cholPars

  if (is.null(outerActiveNames)) outerActiveNames <- names(psiFull)
  if (!all(chol_pars %in% names(psiFull)))
    stop("ecmEvaluateSubject: `psiFull` must contain all omega$cholPars.")

  # Build closed-form prior log-normalisation: log N(eta_b|0,Omega)
  # = -K/2 log(2 pi) - log|L_omega| - 0.5 * z^T z, z = L_omega^{-1} eta_b.
  L_omega         <- omega$buildL(psiFull[chol_pars])
  log_det_L_omega <- sum(log(diag(L_omega)))
  log_norm_prior  <- -K / 2 * log(2 * pi) - log_det_L_omega

  # Other subjects' eta values, named with their per-subject eta names.
  if (nrow(etaModes) != length(subjects))
    stop("ecmEvaluateSubject: `etaModes` must have one row per subject.")
  other_eta_nm   <- as.vector(subject_etas[-subjIdx, , drop = FALSE])
  other_eta_vals <- as.vector(etaModes[-subjIdx, , drop = FALSE])
  names(other_eta_vals) <- other_eta_nm

  dataI  <- datalist[[cn]]
  times_i <- dataI$time

  B           <- nrow(nodesSubj$etaNodes)
  Q           <- length(outerActiveNames)
  log_int     <- numeric(B)
  sign_int    <- nodesSubj$weightSigns
  per_node_gr <- if (with_grad) matrix(0, B, Q,
                                       dimnames = list(NULL, outerActiveNames))
                 else NULL
  per_node_he <- if (with_grad) array(0, c(B, Q, Q),
                                       dimnames = list(NULL, outerActiveNames,
                                                       outerActiveNames))
                 else NULL

  for (b in seq_len(B)) {
    eta_b     <- setNames(nodesSubj$etaNodes[b, ], eta_i_names)
    full_pars <- c(psiFull, eta_b, other_eta_vals)

    pred_b <- xPred(times = times_i, pars = full_pars, fixed = fixed,
                     deriv = with_grad, conditions = cn)
    res_b  <- evalConditionResidual(
      dataI        = dataI,
      predictionI  = pred_b[[cn]],
      pars          = full_pars,
      errfn      = errfn,
      fixed         = fixed,
      cn            = cn,
      eCondNames  = NULL,
      bessel        = 1,
      deriv         = with_grad,
      deriv2        = FALSE)

    # Closed-form log-prior on eta_b alone.
    z_prior     <- forwardsolve(L_omega, eta_b)
    log_prior_b <- log_norm_prior - 0.5 * sum(z_prior^2)

    # log integrand = log|W_b| + log p(y|eta_b) + log p(eta_b|Omega).
    # log p(y|eta_b) = -0.5 * res_b$value (the value carries -2 log p form).
    log_int[b] <- nodesSubj$logAbsWeights[b] - 0.5 * res_b$value + log_prior_b

    if (with_grad) {
      gr <- res_b$gradient[outerActiveNames]
      gr[is.na(gr)] <- 0
      per_node_gr[b, ] <- -0.5 * gr
      he <- res_b$hessian[outerActiveNames, outerActiveNames, drop = FALSE]
      he[is.na(he)] <- 0
      per_node_he[b, , ] <- -0.5 * he
    }
  }

  # Signed log-sum-exp: log(sum_b sign_b * exp(log_int[b])).
  M       <- max(log_int)
  shifted <- sign_int * exp(log_int - M)
  s       <- sum(shifted)
  if (!is.finite(s) || s <= 0) {
    warning(sprintf(paste0("ecmEvaluateSubject(subjIdx=%d): signed-LSE ",
                           "sum = %.3e; grid is under-resolved (raise `level`)."),
                    subjIdx, s))
    logLhat <- -Inf
    softmax <- numeric(B)
  } else {
    logLhat <- M + log(s)
    softmax <- sign_int * exp(log_int - logLhat)
  }

  m_hat <- as.numeric(softmax %*% nodesSubj$etaNodes)
  M_hat <- crossprod(nodesSubj$etaNodes,
                     softmax * nodesSubj$etaNodes)
  names(m_hat) <- omega$eta
  dimnames(M_hat) <- list(omega$eta, omega$eta)

  out <- list(logLhat     = logLhat,
              m_hat       = m_hat,
              M_hat       = M_hat,
              maxSoftmax = max(abs(softmax)),
              n_eff       = if (any(softmax != 0)) 1 / sum(softmax^2) else 0)

  if (with_grad) {
    # Gradient of log L_i: sum_b softmax_b * (d log integrand_b / d theta).
    gr_lse <- as.numeric(softmax %*% per_node_gr)
    names(gr_lse) <- outerActiveNames

    # Hessian of LSE = sum softmax_b H_b + Cov_softmax(grad_b).
    # Cov = E[gg^T] - E[g] E[g]^T (with signed softmax this is an algebraic identity).
    H_avg <- matrix(0, Q, Q, dimnames = list(outerActiveNames, outerActiveNames))
    for (b in seq_len(B)) H_avg <- H_avg + softmax[b] * per_node_he[b, , ]
    Eggt   <- crossprod(per_node_gr, softmax * per_node_gr)
    cov_g  <- Eggt - tcrossprod(gr_lse)
    H_lse  <- H_avg + cov_g

    out$gradient <- -2 * gr_lse
    out$hessian  <- -2 * H_lse
  }
  out
}


# Drop heavy state (emDiag, prdfn, data, omega, errfn, foceiStart,
# stageTrace) from an nlmeFit so a parlist of msnlmeFit results stays small.
# Keeps everything as.parframe.parlist + summary.parlist + downstream
# diagnostics consume: argument, value, gradient, hessian, omega, etaModes,
# etaSE, shrinkage, converged, iterations, method.
.stripNlmeFit <- function(fit) {
  keep <- c("argument", "value", "gradient", "hessian",
            "omega", "etaModes", "etaSE", "shrinkage",
            "converged", "iterations", "method")
  out <- fit[intersect(keep, names(fit))]
  class(out) <- class(fit)
  out
}


#' Multi-start nonlinear mixed-effects fit
#'
#' Runs [nlmeFit()] from many starting points in parallel, returning a
#' [parlist] of fits sorted-ready for [as.parframe()] / [summary()]. The
#' design mirrors [mstrust()]: `center` is either a named-numeric (perturbed
#' by `samplefun` per fit) or a [parframe] (each row used as a starting
#' point). Use this to characterise the multi-modality of the marginal
#' likelihood and pick the best optimum.
#'
#' @param obj An `objfn` passed straight to [nlmeFit()] (typically
#'   `normL2(data, g*x*p) + constraintL2(mu = 0, Omega = om)`).
#' @param omega An [omega] spec with subject expansion.
#' @param center Named numeric or [parframe]. If numeric, the population
#'   parameter vector around which random starts are sampled (structural pars
#'   plus `omega$cholPars`). If a parframe, each row is used as a fixed
#'   starting point and `fits` is overridden by `nrow(center)`.
#' @param prdfn The prediction function `g * x * p` used to build `obj`.
#' @param data The [datalist] used for `obj`.
#' @param errfn Optional obsfn defining a parameter-dependent error model.
#' @param fixed Optional named-numeric of fixed parameters.
#' @param method Estimator. Passed through to [nlmeFit()]: `"focei"`,
#'   `"quadrature"`, or `"foceiQuadrature"`.
#' @param control Nested control list, passed through to [nlmeFit()].
#' @param fits Integer, number of random starts. Ignored when `center` is a
#'   parframe.
#' @param cores Integer, number of parallel workers. On Unix uses
#'   `parallel::mclapply` (fork); on Windows a PSOCK cluster + `foreach`.
#'   Outer parallelism multiplies with the inner OpenMP threads of the C++
#'   objective kernels (`getOption("dMod.objfn.threads")`); keep
#'   `cores * dMod.objfn.threads` below your core count.
#' @param samplefun Name of a random-number generator (default `"rnorm"`)
#'   used to perturb `center`. Extra args in `...` whose names match
#'   `formals(samplefun)` are forwarded.
#' @param start1stfromCenter Logical. If `TRUE`, the first fit starts at
#'   `center` itself (no perturbation). Ignored when `center` is a parframe.
#' @param keepFull Logical. If `FALSE` (the default), each returned fit is
#'   stripped of heavy state (`emDiag`, `prdfn`, `data`, `omega`,
#'   `errfn`, `foceiStart`, `stageTrace`) so the result stays small.
#'   Set `TRUE` if you need to call [predict.nlmeFit()] / [plot.nlmeFit()]
#'   etc. on individual fits.
#' @param studyname Optional character. If `output = TRUE`, fits are written
#'   to `<resultPath>/<studyname>/trial-N-<timestamp>/interRes/`. Defaults to
#'   `"msnlmeFit"`.
#' @param resultPath Character, base directory for the on-disk dump.
#' @param output Logical. If `TRUE`, each fit is saved as it completes
#'   (crash-resilient) and the full parlist is written at the end.
#' @param verbose Logical. If `TRUE`, prints per-fit progress and forwards
#'   `verbose = TRUE` into [nlmeFit()].
#' @param ... Forwarded to `samplefun` (e.g. `sd = 0.3`).
#'
#' @return A [parlist] of length `fits`, each element an `nlmeFit` (stripped
#'   per `keepFull`) with `parinit` and `index` attached. Pass to
#'   [as.parframe()] for a sorted-by-`value` table. Failed fits are stored
#'   as `list(error = ..., value = NA, converged = FALSE, ...)`.
#'
#' @seealso [nlmeFit()], [mstrust()], [msParframe()], [parlist].
#' @export
msnlmeFit <- function(obj, omega, center,
                      prdfn    = NULL,
                      data     = NULL,
                      errfn = NULL,
                      fixed    = NULL,
                      method   = c("focei", "quadrature", "foceiQuadrature"),
                      control  = list(),
                      fits     = 20,
                      cores    = 1,
                      samplefun = "rnorm",
                      start1stfromCenter = FALSE,
                      keepFull   = FALSE,
                      studyname  = NULL,
                      resultPath = ".",
                      output     = FALSE,
                      verbose    = FALSE,
                      ...) {

  method <- match.arg(method)
  cores  <- sanitizeCores(cores)

  # Build the per-fit starting points. parframe input fixes the number of
  # fits and pulls starts directly from rows (sorted-by-value via
  # as.parvec.parframe). Numeric input perturbs by samplefun().
  varargslist <- list(...)
  if (is.parframe(center)) {
    fits <- nrow(center)
    parInitList <- lapply(seq_len(fits), function(i) as.parvec(center, i))
  } else {
    if (is.null(names(center)) || any(!nzchar(names(center))))
      stop("`center` must be a fully named numeric vector or a parframe.")
    namessample <- intersect(names(formals(samplefun)), names(varargslist))
    argssample  <- varargslist[namessample]
    argssample$n <- length(center)
    parInitList <- lapply(seq_len(fits), function(i) {
      if (i == 1L && start1stfromCenter) {
        center
      } else {
        perturb <- do.call(samplefun, argssample)
        out <- center + perturb
        names(out) <- names(center)
        out
      }
    })
  }
  cores <- min(fits, cores)

  # Optional on-disk dump (crash-resilient): one .Rda per fit + a final
  # parameterList.Rda. Mirrors mstrust()'s folder layout.
  interResultFolder <- NULL
  resultFolder      <- NULL
  if (output) {
    if (is.null(studyname)) studyname <- "msnlmeFit"
    m_timeStamp <- format(Sys.time(), "%d-%m-%Y-%H%M%S")
    resultFolderBase <- file.path(resultPath, studyname)
    n_existing <- length(dir(resultFolderBase, pattern = "trial*"))
    m_trial <- paste0("trial-", n_existing + 1L)
    resultFolder <- file.path(resultFolderBase,
                              paste0(m_trial, "-", m_timeStamp))
    interResultFolder <- file.path(resultFolder, "interRes")
    dir.create(interResultFolder, showWarnings = FALSE, recursive = TRUE)
  }

  digits <- if (fits >= 10L) floor(log10(fits)) + 1L else 1L

  doOne <- function(i) {
    init_i <- parInitList[[i]]
    # Invalidate Pimpl warm-start caches per fit so we do not inherit roots
    # from a neighbouring start.
    options(.dMod.fit_token = paste0("msnlme_", i, "_", as.numeric(Sys.time())))
    t0 <- Sys.time()
    fit <- try(suppressMessages(
      nlmeFit(obj, omega, init_i,
              prdfn = prdfn, data = data, errfn = errfn,
              fixed = fixed, method = method, control = control,
              verbose = verbose)),
      silent = !verbose)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (inherits(fit, "try-error")) {
      fit <- list(error      = as.character(fit),
                  value      = NA_real_,
                  converged  = FALSE,
                  iterations = NA_integer_,
                  method     = method)
      class(fit) <- c("nlmeFit", "list")
    } else if (!keepFull) {
      fit <- .stripNlmeFit(fit)
    }
    fit$parinit <- init_i
    fit$index   <- i
    fit$elapsed <- elapsed

    if (verbose) {
      msg <- if (!is.null(fit$error)) {
        sprintf("[msnlmeFit %s/%d] FAILED (%.1fs): %s",
                formatC(i, width = digits, flag = "0"), fits,
                elapsed, fit$error)
      } else {
        sprintf("[msnlmeFit %s/%d] OFV=%.6f  conv=%s  iter=%s  (%.1fs)",
                formatC(i, width = digits, flag = "0"), fits,
                fit$value, fit$converged,
                format(fit$iterations %||% NA_integer_), elapsed)
      }
      message(msg)
    }

    if (output) {
      saveRDS(fit, file = file.path(interResultFolder,
                                    sprintf("fit-%d.Rda", i)))
    }
    fit
  }

  # Parallel dispatch. Fork on Unix; PSOCK + foreach on Windows so the
  # workers can pick up obj / prdfn / data via clusterExport. Falls back to
  # serial when cores == 1 or fits == 1.
  if (cores > 1L) {
    if (Sys.info()[['sysname']] == "Windows") {
      cluster <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
      doParallel::registerDoParallel(cluster)
      parallel::clusterCall(cl = cluster,
                            function(x) .libPaths(x), .libPaths())
      parallel::clusterExport(
        cluster, envir = environment(),
        varlist = c("obj", "omega", "parInitList", "prdfn", "data",
                    "errfn", "fixed", "method", "control", "verbose",
                    "keepFull", "output", "interResultFolder", "fits",
                    "digits"))
      `%mydo%` <- foreach::`%dopar%`
      i <- NULL
      results <- foreach::foreach(
        i = seq_len(fits),
        .packages = .packages(),
        .inorder = TRUE,
        .options.multicore = list(preschedule = FALSE)) %mydo% doOne(i)
    } else {
      results <- parallel::mclapply(seq_len(fits), doOne,
                                    mc.cores = cores,
                                    mc.preschedule = FALSE)
    }
  } else {
    results <- lapply(seq_len(fits), doOne)
  }

  if (output) {
    saveRDS(results, file = file.path(resultFolder, "parameterList.Rda"))
  }

  as.parlist(results)
}

