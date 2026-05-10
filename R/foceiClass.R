## FOCEI: nonlinear mixed-effects estimation -----------------------------------
##
## Two user-facing entries: omega() builds an omegaSpec describing the random-
## effect covariance structure, and focei() wraps a joint objfn into the bilevel
## FOCEI estimator. The MVN path of constraintL2 (in objClass.R) consumes
## omegaSpec for the prior contribution sum_i eta_i^T Omega^-1 eta_i + N log|Omega|.
##
## The naming convention generated here is:
##   diagonal Cholesky entries:  paste0(prefix, "_", short, "_", short)
##                               with L_kk = exp(par)  (log-Cholesky on the diagonal)
##   off-diagonal entries:       paste0(prefix, "_", short_k, "_", short_l)  for k > l
##                               with L_kl = par       (free, signed)
## where short = sub("^eta_", "", eta) so eta_Cl -> Cl, eta_V -> V, etc.



#' Specify random-effect covariance structure for FOCEI
#'
#' @description
#' Builds an `omegaSpec` object describing a Cholesky-parametrised random-effect
#' covariance matrix `Omega = L L^T`, with `L` lower-triangular, log-diagonal,
#' and arbitrary sparsity below the diagonal. The result is consumed by
#' [constraintL2] (when its `Omega` argument is set) and by [focei].
#'
#' @param eta Character vector of base eta names, one per random effect, e.g.
#'   `c("eta_Cl", "eta_V", "eta_Ka")`. Determines the dimension `K`.
#' @param structure One of `"diag"` (default, only diagonal Cholesky entries
#'   estimated, K free parameters) or `"full"` (full lower-triangular Cholesky,
#'   K(K+1)/2 free parameters).
#' @param correlate Optional list of length-2 character vectors selecting
#'   off-diagonal Cholesky entries to be estimated. Alternative to
#'   `structure = "full"`. Each entry is a pair of eta names whose correlation
#'   is freed.
#' @param subjects Optional character vector of subject identifiers. When given,
#'   the per-subject random-effect parameter names are constructed as
#'   `paste0(eta[k], "_", subjects[i])` and stored in `subject_etas`. Required
#'   downstream by [focei] but optional here.
#' @param prefix Character string used as a name stem for Cholesky parameters.
#'   Default `"omega"`. With `eta = c("eta_Cl", "eta_V")` this yields
#'   `omega_Cl_Cl`, `omega_V_V`, `omega_V_Cl`.
#'
#' @return An `omegaSpec` (S3 list) with components:
#' \describe{
#'   \item{`K`}{Dimension.}
#'   \item{`eta`}{Base eta names.}
#'   \item{`chol_pars`}{Character vector of Cholesky parameter names in
#'     column-major order.}
#'   \item{`name_matrix`}{`K x K` character matrix; entry `[k, l]` is the
#'     parameter name of `L[k, l]` or `NA` if structurally zero.}
#'   \item{`is_diag`}{Named logical vector: TRUE for log-diagonal Cholesky
#'     parameters, FALSE for free off-diagonal entries.}
#'   \item{`chol_loc`}{`length(chol_pars) x 2` integer matrix giving the
#'     `(row, col)` location of each Cholesky parameter in `L`.}
#'   \item{`subject_etas`}{`N x K` character matrix of per-subject eta names,
#'     or `NULL` if `subjects` was not supplied.}
#'   \item{`build_L(chol_vec)`}{Function returning the lower-triangular `L`
#'     given a named numeric vector containing at least `chol_pars`.}
#' }
#'
#' @details
#' The Cholesky diagonal is on the log scale (enforces positive definiteness
#' of \eqn{\Omega}); off-diagonal entries are unrestricted reals. Matches
#' `lme4`/`TMB` conventions.
#'
#' @examples
#' om1 <- omega(eta = c("eta_Cl", "eta_V"))
#' om1$chol_pars
#' om1$build_L(c(omega_Cl_Cl = log(0.3), omega_V_V = log(0.2)))
#'
#' om2 <- omega(eta = c("eta_Cl", "eta_V", "eta_Ka"), structure = "full",
#'              subjects = paste0("subj", 1:3))
#' om2$subject_etas
#' om2$chol_pars
#'
#' @seealso [constraintL2], [focei]
#' @export
omega <- function(eta,
                  structure = c("diag", "full"),
                  correlate = NULL,
                  subjects  = NULL,
                  prefix    = "omega") {

  structure <- match.arg(structure)
  if (!is.character(eta) || length(eta) < 1L)
    stop("`eta` must be a non-empty character vector of random-effect names.")
  if (anyDuplicated(eta))
    stop("`eta` names must be unique.")

  K <- length(eta)
  short <- sub("^eta_", "", eta)

  # --- name_matrix: K x K character (NA = structurally zero) -----------------
  name_matrix <- matrix(NA_character_, K, K, dimnames = list(eta, eta))
  for (k in seq_len(K)) {
    name_matrix[k, k] <- paste(prefix, short[k], short[k], sep = "_")
  }

  if (structure == "full") {
    if (K >= 2L) for (k in 2:K) for (l in 1:(k - 1)) {
      name_matrix[k, l] <- paste(prefix, short[k], short[l], sep = "_")
    }
  } else if (!is.null(correlate)) {
    if (!is.list(correlate))
      stop("`correlate` must be a list of length-2 character vectors.")
    for (pair in correlate) {
      if (length(pair) != 2L || !is.character(pair))
        stop("Each entry of `correlate` must be a length-2 character vector.")
      i1 <- match(pair[1], eta)
      i2 <- match(pair[2], eta)
      if (is.na(i1) || is.na(i2))
        stop("`correlate` refers to unknown eta names: ",
             paste(pair[is.na(c(i1, i2))], collapse = ", "))
      kk <- max(i1, i2); ll <- min(i1, i2)
      if (kk == ll) next
      name_matrix[kk, ll] <- paste(prefix, short[kk], short[ll], sep = "_")
    }
  }

  # --- chol_pars in column-major order (matches forwardsolve traversal) ------
  chol_pars <- character(0)
  chol_loc  <- matrix(integer(0), nrow = 0, ncol = 2,
                      dimnames = list(NULL, c("row", "col")))
  for (l in seq_len(K)) for (k in l:K) {
    nm <- name_matrix[k, l]
    if (is.na(nm)) next
    chol_pars <- c(chol_pars, nm)
    chol_loc  <- rbind(chol_loc, c(k, l))
  }
  rownames(chol_loc) <- chol_pars

  if (anyDuplicated(chol_pars))
    stop("Generated duplicate Cholesky parameter names. Check your `prefix`.")

  is_diag <- chol_loc[, "row"] == chol_loc[, "col"]
  names(is_diag) <- chol_pars

  # --- subject expansion -----------------------------------------------------
  subject_etas <- NULL
  if (!is.null(subjects)) {
    if (!is.character(subjects) || anyDuplicated(subjects))
      stop("`subjects` must be a non-empty character vector with unique entries.")
    subject_etas <- outer(subjects, eta, function(s, e) paste(e, s, sep = "_"))
    dimnames(subject_etas) <- list(subjects, eta)
  }

  # --- build_L closure -------------------------------------------------------
  build_L <- function(chol_vec) {
    if (!all(chol_pars %in% names(chol_vec)))
      stop("`chol_vec` is missing Cholesky parameters: ",
           paste(setdiff(chol_pars, names(chol_vec)), collapse = ", "))
    L <- matrix(0, K, K, dimnames = list(eta, eta))
    for (m in seq_along(chol_pars)) {
      k <- chol_loc[m, 1L]; l <- chol_loc[m, 2L]
      v <- chol_vec[[chol_pars[m]]]
      L[k, l] <- if (is_diag[m]) exp(v) else v
    }
    L
  }

  out <- list(
    K            = K,
    eta          = eta,
    short        = short,
    structure    = structure,
    chol_pars    = chol_pars,
    name_matrix  = name_matrix,
    is_diag      = is_diag,
    chol_loc     = chol_loc,
    subject_etas = subject_etas,
    build_L      = build_L,
    prefix       = prefix
  )
  class(out) <- c("omegaSpec", "list")
  out
}



#' Print method for omegaSpec
#'
#' @param x An `omegaSpec` object.
#' @param ... Ignored.
#' @export
print.omegaSpec <- function(x, ...) {
  cat("omegaSpec\n")
  cat(sprintf("  K            : %d\n", x$K))
  cat(sprintf("  eta          : %s\n", paste(x$eta, collapse = ", ")))
  cat(sprintf("  structure    : %s\n", x$structure))
  cat(sprintf("  chol_pars    : %d (%s)\n",
              length(x$chol_pars), paste(x$chol_pars, collapse = ", ")))
  if (!is.null(x$subject_etas)) {
    cat(sprintf("  subjects     : %d (%s)\n",
                nrow(x$subject_etas),
                paste(rownames(x$subject_etas), collapse = ", ")))
    cat(sprintf("  subject_etas : %d unique random-effect parameters\n",
                length(unique(as.vector(x$subject_etas)))))
  } else {
    cat("  subjects     : none (call omega(..., subjects = ...) before focei())\n")
  }
  invisible(x)
}



#' All random-effect parameter names of an omegaSpec
#'
#' @description
#' Returns the union of `chol_pars` and all subject-level eta names in a single
#' character vector. Useful for constructing initial parameter vectors and for
#' partitioning the outer/inner parameter space in [focei].
#'
#' @param x An `omegaSpec` object.
#' @param what One of `"all"` (default, both eta and Cholesky parameters),
#'   `"eta"` (subject-level random effects only), or `"chol"` (Cholesky
#'   parameters only).
#' @return Character vector of parameter names.
#' @export
parnames.omegaSpec <- function(x, what = c("all", "eta", "chol")) {
  what <- match.arg(what)
  eta_names <- if (is.null(x$subject_etas)) character(0) else as.vector(x$subject_etas)
  switch(what,
         all  = c(x$chol_pars, eta_names),
         eta  = eta_names,
         chol = x$chol_pars)
}



#' FOCEI bilevel estimator
#'
#' @description
#' Wraps a user-built joint objective into a FOCEI (First-Order Conditional
#' Estimation with Interaction) marginal-likelihood approximation. The returned
#' objfn lives over the population-level parameters
#' \eqn{(\theta, \Omega, \Sigma)}; the subject-level random effects \eqn{\eta_i}
#' are integrated out by Laplace approximation at each evaluation.
#'
#' @param joint An `objfn` typically built as
#'   `normL2(data, g*x*p, errmodel = err) + constraintL2(mu = 0, Omega = om)`.
#'   It must accept the full parameter set including all subject-level eta
#'   parameters listed in `omegaSpec$subject_etas`. Subjects must be encoded as
#'   dMod conditions (one condition per subject), typically created via
#'   [branch] with one row per subject in the substitution table.
#' @param omegaSpec An `omegaSpec` object with subject expansion (see [omega]).
#'   Determines the subjects, the K random-effect names per subject, and the
#'   Cholesky parametrisation of \eqn{\Omega}.
#' @param model Optional `prdfn` (typically `g * x * p`) supporting
#'   `deriv2 = TRUE`. Required when `correction != "none"` so the analytical
#'   `\partial \log|H_i|/\partial \theta` can be assembled from
#'   `\partial^2 g / \partial \eta \partial \theta`. Pass the same prediction
#'   function that was used to build `joint`.
#' @param data Optional [datalist] supplying per-observation `sigma`.
#'   Required when `correction != "none"`. `sigma` may depend on `theta` but
#'   not on `eta`.
#' @param innerControl Named list of inner-trust options: `rtol`
#'   (default `1e-6`), `maxit` (default `50`), `rinit`/`rmax`
#'   (default `1`/`10`).
#' @param correction One of `"none"` (default, Stage 1 envelope gradient
#'   only), `"lagged"` (Stage 2 analytical correction with cached refresh),
#'   or `"eager"` (Stage 2 recomputed every outer iteration).
#' @param correctionControl Named list for `correction = "lagged"`: `tau`
#'   (anchor-distance threshold, default `0.1`) and `M` (periodic refresh
#'   interval, default `5L`).
#' @param cores Integer. Cores for per-subject inner solves (default `1`).
#'   `mclapply` on POSIX, on-the-fly `parLapply` cluster on Windows.
#'
#' @return An object of class `c("focei", "objfn", "fn")`, callable as
#'   `obj(pars, fixed = NULL, deriv = TRUE, env = NULL)`. The returned
#'   `objlist` carries per-subject diagnostics under `attr(., "focei_diag")`.
#'
#' @seealso [omega], [constraintL2], [trust]
#' @examples
#' \dontrun{
#' # ... build model, data, joint, omega spec ...
#' focei_obj <- focei(joint, om, model = g*x*p, data = dlist,
#'                    correction = "lagged")
#' fit <- trust(focei_obj, init,
#'              rinit = 1, rmax = 10, iterlim = 200,
#'              on_step = attr(focei_obj, "on_step"))
#' }
#' @export
focei <- function(joint, omegaSpec,
                  model             = NULL,
                  data              = NULL,
                  innerControl      = list(),
                  correction        = c("none", "lagged", "eager"),
                  correctionControl = list(tau = 0.1, M = 5L),
                  cores             = 1L) {

  if (!inherits(joint, "objfn"))
    stop("`joint` must be an objfn (typically normL2 + constraintL2 MVN).")
  if (!inherits(omegaSpec, "omegaSpec"))
    stop("`omegaSpec` must be an omegaSpec object built by omega().")
  if (is.null(omegaSpec$subject_etas))
    stop("`omegaSpec` must have subject expansion. Call omega(..., subjects = ...).")

  correction <- match.arg(correction)
  if (correction != "none") {
    if (is.null(model))
      stop("`model` (the prdfn used to build `joint`) is required when ",
           "correction != \"none\".")
    if (is.null(data))
      stop("`data` (the datalist used to build `joint`) is required when ",
           "correction != \"none\".")
  }

  K            <- omegaSpec$K
  subject_etas <- omegaSpec$subject_etas
  subjects     <- rownames(subject_etas)
  N            <- length(subjects)
  all_eta_nm   <- as.vector(subject_etas)

  joint_pars <- attr(joint, "parameters")
  if (!all(all_eta_nm %in% joint_pars))
    warning("`joint` does not expose all subject-level eta names as parameters; ",
            "make sure your parameter transformation references them.")

  rtol  <- innerControl$rtol  %||% 1e-6
  maxit <- innerControl$maxit %||% 50L
  rinit <- innerControl$rinit %||% 1
  rmax  <- innerControl$rmax  %||% 10

  ctl_tau <- correctionControl$tau %||% 0.1
  ctl_M   <- correctionControl$M   %||% 5L

  # Warmstart cache (mutable via <<-)
  cache_etas <- matrix(0, N, K, dimnames = dimnames(subject_etas))

  # Lagging cache for the log|H_i| correction (mutable via <<-)
  cache_correction   <- NULL
  cache_anchor       <- NULL
  iter_since_refresh <- 0L
  last_step_rejected <- FALSE
  last_rho           <- NA_real_
  refresh_count      <- 0L

  # Pre-extract per-subject data slices once.
  data_per_subject <- NULL
  times_union      <- NULL
  if (correction != "none") {
    data_per_subject <- lapply(subjects, function(s) {
      d <- data[[s]]
      if (is.null(d))
        stop("`data` is missing condition '", s,
             "' required by omegaSpec$subject_etas.")
      d
    })
    names(data_per_subject) <- subjects
    times_union <- sort(unique(c(0, unlist(lapply(data_per_subject, `[[`, "time")))))
  }

  # --- analytical log|H_i|-correction (Stage 2) ---------------------------
  # GN form throughout: H_i = 2 * (sum_t (1/sigma_t^2) G_t G_t^T + Omega^-1).
  # The leading factor 2 is dMod's normL2/constraintL2 convention -- both
  # carry value = (residual/sigma)^2 (no 1/2 prefactor), so their second
  # derivatives are 2x the conventional Gaussian negative-log-likelihood
  # Hessian. We mirror that here so dH/d theta and H^-1 (read from the inner
  # solver via H_inv_list) live in the same convention; the trace
  # tr(H^-1 dH/d theta) is independent of an overall scale of H, but we keep
  # both pieces consistent so the formula stays auditable.
  # The correction itself is
  #   c_k = sum_i [ tr(H_i^-1 dH_i/d theta_k|_eta)
  #               + sum_l tr(H_i^-1 dH_i/d eta_l) * (d eta_l*/d theta_k) ]
  # with d eta*/d theta = - H_i^-1 * (d^2 J / d eta d theta) read from the
  # joint's GN Hessian cross block. For Cholesky parameters the prediction
  # carries no sensitivity (they live in the prior); the explicit term
  # reduces to tr(H_i^-1 d Omega^-1/d chol_par) computed from a small FD
  # over omegaSpec$build_L.
  conv2 <- 2  # dMod (res/sigma)^2 convention factor -- see comment above.
  compute_correction <- function(full_pars, joint_hessian, fixed,
                                 outer_names, H_inv_list) {
    Q <- length(outer_names)
    pred <- model(times = times_union, pars = full_pars, fixed = fixed,
                  deriv = TRUE, deriv2 = TRUE, conditions = subjects)
    correction <- setNames(numeric(Q), outer_names)

    # d Omega^-1 / d chol_par via central FD on omegaSpec$build_L: K is small
    # (single-digit), the build_L closure is pure R, so the cost is negligible
    # vs the model evaluation above.
    chol_in_outer <- intersect(outer_names, omegaSpec$chol_pars)
    dOmega_inv_dchol <- list()
    if (length(chol_in_outer) > 0L) {
      chol_vec <- full_pars[omegaSpec$chol_pars]
      L0       <- omegaSpec$build_L(chol_vec)
      h_chol   <- sqrt(.Machine$double.eps) * pmax(abs(chol_vec), 1)
      for (cp in chol_in_outer) {
        cv_p <- chol_vec; cv_p[cp] <- cv_p[cp] + h_chol[cp]
        cv_m <- chol_vec; cv_m[cp] <- cv_m[cp] - h_chol[cp]
        Lp <- omegaSpec$build_L(cv_p); Lm <- omegaSpec$build_L(cv_m)
        Op_inv <- chol2inv(chol(Lp %*% t(Lp)))
        Om_inv <- chol2inv(chol(Lm %*% t(Lm)))
        dOmega_inv_dchol[[cp]] <- (Op_inv - Om_inv) / (2 * h_chol[cp])
        dimnames(dOmega_inv_dchol[[cp]]) <-
          list(omegaSpec$eta, omegaSpec$eta)
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
        stop("`model` did not produce 'deriv'/'deriv2' attributes for ",
             "condition '", cn, "'. Did you build it with deriv2-capable ",
             "Xs/Y/P?")

      avail_pars <- dimnames(d_full)[[3]]
      th_avail   <- intersect(outer_names, avail_pars)
      eta_avail  <- intersect(eta_i_nm,    avail_pars)
      if (!length(eta_avail)) next

      times_i  <- d_i$time
      names_i  <- as.character(d_i$name)
      sigma_i  <- d_i$sigma
      time_idx <- match(times_i, times_union)
      if (anyNA(time_idx))
        stop("compute_correction: data times missing from prediction grid.")

      Tn <- length(times_i)
      Kn <- length(eta_avail)
      Qn <- length(th_avail)

      # Per-row slices: G[t,] = d g(t,name_t)/d eta, Wmix[t,k,m] = d^2 g/d eta d theta.
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

      # Explicit term, structural-parameter contribution: dH/d theta_k =
      # conv2 * sum_t Sinv[t] * (Wmix[t,,k] outer G[t,] + G[t,] outer Wmix[t,,k])
      if (Qn > 0L) {
        explicit <- numeric(Qn)
        for (k in seq_len(Qn)) {
          Mk <- crossprod(G * Sinv, Wmix[, , k, drop = TRUE])  # K x K
          Mk <- Mk + t(Mk)
          explicit[k] <- conv2 * sum(H_inv * Mk)
        }
        correction[th_avail] <- correction[th_avail] + explicit
      }

      # Explicit term, Cholesky-parameter contribution: dH/d chol_par
      # = d Omega^-1 / d chol_par (the predictions are blind to chol pars).
      # dOmega_inv_dchol is indexed by BASE eta names (omegaSpec$eta);
      # eta_avail uses SUBJECT-specific eta names. Map via column index.
      if (length(chol_in_outer) > 0L) {
        eta_pos  <- match(eta_avail, eta_i_nm)        # 1..K positions
        eta_base <- omegaSpec$eta[eta_pos]
        for (cp in chol_in_outer) {
          dOmega_inv_sub <-
            dOmega_inv_dchol[[cp]][eta_base, eta_base, drop = FALSE]
          correction[cp] <- correction[cp] +
            conv2 * sum(H_inv * dOmega_inv_sub)
        }
      }

      # Implicit chain via eta*(theta) -- covers BOTH structural and chol pars
      # because joint$hessian carries the cross block for the whole outer set.
      # The deta*/dtheta product H_inv %*% cross_block is convention-neutral
      # (both factors carry the same conv2), so we only scale dlogH_deta.
      th_for_implicit <- intersect(outer_names, colnames(joint_hessian))
      if (length(th_for_implicit) > 0L) {
        cross_block <- joint_hessian[eta_avail, th_for_implicit, drop = FALSE]
        deta_dth <- -H_inv %*% cross_block                     # K x Q'

        dlogH_deta <- numeric(Kn)
        for (l in seq_len(Kn)) {
          Ml <- crossprod(G * Sinv, Veta[, , l, drop = TRUE])
          Ml <- Ml + t(Ml)
          dlogH_deta[l] <- conv2 * sum(H_inv * Ml)
        }
        implicit <- as.numeric(dlogH_deta %*% deta_dth)
        correction[th_for_implicit] <-
          correction[th_for_implicit] + implicit
      }
    }
    correction
  }

  # on_step closure for trust(): flips the rejection flag the focei wrapper
  # consumes on its next refresh decision. Stored on attr(myfn, "on_step").
  on_step_fn <- function(rho, accepted, iter, r) {
    last_rho           <<- rho
    last_step_rejected <<- !isTRUE(accepted)
    invisible(NULL)
  }

  inner_solve_one <- function(i, outer_input, fixed) {
    eta_i_names  <- subject_etas[i, ]
    eta_i_init   <- cache_etas[i, ]
    names(eta_i_init) <- eta_i_names
    other_etas <- as.vector(subject_etas[-i, , drop = FALSE])
    other_vals <- as.vector(cache_etas[-i, , drop = FALSE])
    names(other_vals) <- other_etas

    inner_objfn <- function(eta_i_in, ...) {
      full_pars <- c(outer_input, other_vals, eta_i_in)
      out <- joint(pars = full_pars, fixed = fixed, deriv = TRUE,
                   conditions = subjects[i])
      gr <- out$gradient[eta_i_names]
      hs <- out$hessian[eta_i_names, eta_i_names, drop = FALSE]
      objlist(value = out$value, gradient = gr, hessian = hs)
    }

    fit <- try(suppressMessages(trust(
      inner_objfn, parinit = eta_i_init,
      rinit = rinit, rmax = rmax,
      iterlim = maxit, fterm = rtol, mterm = rtol)),
      silent = TRUE)

    if (inherits(fit, "try-error") || isTRUE(fit$converged) == FALSE) {
      list(eta_star = eta_i_init, H_i = diag(K),
           value = NA_real_, iter = NA_integer_, converged = FALSE)
    } else {
      list(eta_star  = setNames(as.numeric(fit$argument), eta_i_names),
           H_i       = fit$hessian,
           value     = fit$value,
           iter      = fit$iterations,
           converged = isTRUE(fit$converged))
    }
  }

  myfn <- function(..., fixed = NULL, deriv = TRUE, env = NULL) {

    p <- list(...)[[match.fnargs(list(...), "pars")]]
    if (is.null(env)) env <- new.env()

    outer_input_names <- intersect(names(p), setdiff(joint_pars, all_eta_nm))
    outer_input       <- p[outer_input_names]

    # --- per-subject inner solves ------------------------------------------
    if (cores == 1L) {
      results <- lapply(seq_len(N), inner_solve_one, outer_input, fixed)
    } else if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      results <- parallel::parLapply(cl, seq_len(N), inner_solve_one,
                                     outer_input, fixed)
    } else {
      results <- parallel::mclapply(seq_len(N), inner_solve_one,
                                    outer_input, fixed, mc.cores = cores)
    }

    # --- update warmstart cache --------------------------------------------
    for (i in seq_len(N)) {
      cache_etas[i, ] <<- results[[i]]$eta_star
    }

    # --- assemble full pars at modes --------------------------------------
    eta_full           <- as.vector(cache_etas)
    names(eta_full)    <- as.vector(subject_etas)
    full_pars          <- c(outer_input, eta_full)

    # --- one full joint evaluation at modes -------------------------------
    # GN (deriv2 = FALSE) keeps the inner Hessian consistent with the GN H_i
    # used inside the log|H_i|-correction trace, so eigen-floor + Schur block
    # are coherent. The cross block joint$hessian[eta, theta] read out below
    # for `d eta*/d theta` is the GN chain rule -- matches H_i^GN.
    joint_at_modes <- joint(pars = full_pars, fixed = fixed,
                            deriv = deriv, conditions = subjects, env = env)

    # --- log|H_i| with eigen floor + collect diagnostics -------------------
    sum_logdetH <- 0
    diag_eta    <- matrix(NA_real_, N, K, dimnames = dimnames(subject_etas))
    diag_eigs   <- vector("list", N)
    diag_logdet <- numeric(N)
    diag_nfloor <- integer(N)
    diag_iter   <- integer(N)
    diag_conv   <- logical(N)
    H_i_inv_list <- vector("list", N)

    for (i in seq_len(N)) {
      diag_eta[i, ]     <- results[[i]]$eta_star
      diag_iter[i]      <- if (is.null(results[[i]]$iter)) NA_integer_ else results[[i]]$iter
      diag_conv[i]      <- isTRUE(results[[i]]$converged)
      H_i               <- results[[i]]$H_i
      eig               <- eigen(H_i, symmetric = TRUE)
      tr_H              <- sum(eig$values)
      epsilon           <- 1e-10 * abs(tr_H) / K
      lambda_floored    <- pmax(eig$values, epsilon)
      diag_eigs[[i]]    <- eig$values
      diag_nfloor[i]    <- sum(eig$values < epsilon)
      diag_logdet[i]    <- sum(log(lambda_floored))
      sum_logdetH       <- sum_logdetH + diag_logdet[i]
      H_inv_i           <- eig$vectors %*% (t(eig$vectors) / lambda_floored)
      dimnames(H_inv_i) <- dimnames(H_i)
      H_i_inv_list[[i]] <- H_inv_i
    }

    OFV <- joint_at_modes$value + sum_logdetH

    if (!deriv) {
      out <- objlist(value    = OFV,
                     gradient = setNames(numeric(length(p)), names(p)),
                     hessian  = matrix(0, length(p), length(p),
                                       dimnames = list(names(p), names(p))))
      attr(out, "focei_diag") <- list(eta_star         = diag_eta,
                                      eigs             = diag_eigs,
                                      logdet           = diag_logdet,
                                      n_floored        = diag_nfloor,
                                      iter             = diag_iter,
                                      converged        = diag_conv,
                                      correction_value = cache_correction,
                                      refresh_count    = refresh_count,
                                      iter_since_refresh = iter_since_refresh,
                                      last_rho         = last_rho)
      attr(out, "env") <- env
      return(out)
    }

    # --- outer gradient: envelope theorem from joint_at_modes -------------
    grad_full <- joint_at_modes$gradient
    outer_active <- intersect(outer_input_names, names(grad_full))
    grad_outer <- setNames(numeric(length(p)), names(p))
    grad_outer[outer_active] <- grad_full[outer_active]

    # --- log|H_i| correction with BDF-style lagging (Stage 2) -------------
    if (correction != "none" && length(outer_active) > 0L) {
      theta_now <- as.numeric(outer_input[outer_active])
      do_refresh <- correction == "eager" ||
                    is.null(cache_correction) ||
                    isTRUE(last_step_rejected) ||
                    (!is.na(last_rho) && last_rho < 0.25) ||
                    iter_since_refresh >= ctl_M ||
                    (!is.null(cache_anchor) &&
                     max(abs(theta_now - cache_anchor) /
                         pmax(abs(cache_anchor), 1)) > ctl_tau)

      if (do_refresh) {
        cache_correction   <<- compute_correction(
          full_pars     = full_pars,
          joint_hessian = joint_at_modes$hessian,
          fixed         = fixed,
          outer_names   = outer_active,
          H_inv_list    = H_i_inv_list)
        cache_anchor       <<- theta_now
        iter_since_refresh <<- 0L
        last_step_rejected <<- FALSE
        refresh_count      <<- refresh_count + 1L
      } else {
        iter_since_refresh <<- iter_since_refresh + 1L
      }

      take <- intersect(outer_active, names(cache_correction))
      grad_outer[take] <- grad_outer[take] + cache_correction[take]
    }

    # --- outer Hessian: Schur complement over eta blocks -------------------
    H_full <- joint_at_modes$hessian
    Hess <- matrix(0, length(p), length(p),
                   dimnames = list(names(p), names(p)))
    if (length(outer_active) > 0L) {
      Hess[outer_active, outer_active] <-
        H_full[outer_active, outer_active, drop = FALSE]
      for (i in seq_len(N)) {
        eta_i_names <- subject_etas[i, ]
        if (!all(eta_i_names %in% rownames(H_full))) next
        H_oi <- H_full[outer_active, eta_i_names, drop = FALSE]
        if (sum(abs(H_oi)) == 0) next
        Hess[outer_active, outer_active] <-
          Hess[outer_active, outer_active] - H_oi %*% H_i_inv_list[[i]] %*% t(H_oi)
      }
    }

    out <- objlist(value = OFV, gradient = grad_outer, hessian = Hess)
    attr(out, "focei_diag") <- list(eta_star          = diag_eta,
                                    eigs              = diag_eigs,
                                    logdet            = diag_logdet,
                                    n_floored         = diag_nfloor,
                                    iter              = diag_iter,
                                    converged         = diag_conv,
                                    correction_value  = cache_correction,
                                    refresh_count     = refresh_count,
                                    iter_since_refresh = iter_since_refresh,
                                    last_rho          = last_rho)
    attr(out, "env") <- env
    out
  }

  class(myfn) <- c("focei", "objfn", "fn")
  attr(myfn, "conditions") <- attr(joint, "conditions")
  attr(myfn, "parameters") <- setdiff(joint_pars, all_eta_nm)
  attr(myfn, "omegaSpec")  <- omegaSpec
  attr(myfn, "correction") <- correction
  attr(myfn, "on_step")    <- on_step_fn
  myfn
}
