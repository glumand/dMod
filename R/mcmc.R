## Sampler entry, control constructors, dispatcher and run drivers for
## mcmc(). Target / metric / Omega class definitions live in R/mcmcClass.R;
## plot methods on mcmcResult* are in R/plots.R.


#' Control constructors for [mcmc()]
#'
#' Typed parameter bundles for each axis of the [mcmc()] sampler.
#' `mhControl` / `langevinControl` / `hmcControl` / `nutsControl` /
#' `pathIntegralControl` configure the inner move kernel; `smcControl`
#' configures the SMC tempering loop; `metricControl` configures
#' geometry options shared by Langevin / HMC / NUTS.
#'
#' @param stepsize Initial step size. `NULL` triggers a per-move heuristic
#'   refined by dual averaging during warmup.
#' @param acceptTarget Target Metropolis acceptance probability for dual
#'   averaging adaption.
#' @param proposalCov For [mhControl()]: one of `"identity"`, `"fisher"`
#'   (use Fisher information from the first objfn call), or `"fixed"`
#'   (use `covFixed`).
#' @param covFixed Fixed K x K proposal covariance for `proposalCov = "fixed"`.
#' @param preconditioner For [langevinControl()]: one of `"local"`
#'   (Fisher recomputed each step), `"fixed"` (use `GFixed`), or `"identity"`.
#' @param correction For [langevinControl()]: `"none"` (simplified mMALA)
#'   or `"full"` (RMALA with Christoffel drift; requires `metricContext`).
#' @param GFixed K x K Fisher information used when
#'   `preconditioner = "fixed"`. Expected in `-2 log L` units (matches
#'   `trust()$hessian`); the kernel divides by 2 internally.
#' @param ridge Adaptive-ridge floor for the metric Cholesky.
#' @param leapfrogSteps Number of leapfrog steps per HMC iteration.
#' @param maxTreeDepth Maximum NUTS tree depth.
#' @param deltaMax NUTS divergence threshold on the Hamiltonian.
#' @param massMatrix For [hmcControl()] / [nutsControl()]: `"identity"`,
#'   `"fisher"` (frozen at parinit), or `"fixed"` (use `MFixed`).
#' @param MFixed Fixed K x K mass matrix for `massMatrix = "fixed"`.
#' @param P Path-integral replica count (slot only; not yet implemented).
#' @param hbar Path-integral coupling (slot only; not yet implemented).
#' @param nParticles Number of SMC particles.
#' @param essThreshold ESS-fraction trigger for the adaptive-beta bisection.
#' @param malaSteps Inner move iterations per SMC level.
#' @param schedule Optional fixed beta ladder for
#'   `sequenceSchedule = "fixed"`. Numeric vector starting at 0, ending at 1.
#' @param parallelStrategy Particle-parallelism: `"forks"` (mclapply on Unix,
#'   serial on Windows) or `"serial"`.
#' @param continuousAdaption If `TRUE`, dual-averaging state persists across
#'   SMC levels rather than resetting per level.
#' @param verbose Print one line per SMC level.
#' @param blather Attach full per-level particle / weight history.
#' @param daGamma,daT0,daKappa Dual-averaging hyperparameters
#'   (Hoffman-Gelman 2014).
#' @param metricContext A [metricContext()] object, required when
#'   `correction = "full"` in [langevinControl()].
#' @param fisherMode Reserved for future Fisher-vs-exact-Hessian selection.
#'
#' @return A tagged control list consumed by [mcmc()].
#'
#' @name mcmcControl
#' @rdname mcmcControl
NULL


#' @export
#' @rdname mcmcControl
mhControl <- function(stepsize     = NULL,
                      acceptTarget = 0.234,
                      proposalCov  = c("identity", "fisher", "fixed"),
                      covFixed     = NULL,
                      ridge        = 1e-8,
                      daGamma      = 0.05, daT0 = 10, daKappa = 0.75) {
  proposalCov <- match.arg(proposalCov)
  stopifnot(is.null(stepsize) || (is.numeric(stepsize) && stepsize > 0),
            acceptTarget > 0, acceptTarget < 1,
            ridge >= 0, daGamma > 0, daT0 >= 0,
            daKappa > 0.5, daKappa <= 1)
  if (proposalCov == "fixed" && is.null(covFixed))
    stop("mhControl: proposalCov = 'fixed' requires covFixed.")
  structure(list(stepsize     = stepsize,
                 acceptTarget = acceptTarget,
                 proposalCov  = proposalCov,
                 covFixed     = covFixed,
                 ridge        = ridge,
                 daGamma      = daGamma, daT0 = daT0, daKappa = daKappa),
            class = c("mhControl", "list"))
}


#' @export
#' @rdname mcmcControl
langevinControl <- function(stepsize       = NULL,
                            acceptTarget   = 0.574,
                            preconditioner = c("local", "fixed", "identity"),
                            correction     = c("none", "full"),
                            GFixed         = NULL,
                            ridge          = 1e-8,
                            daGamma = 0.05, daT0 = 10, daKappa = 0.75) {
  preconditioner <- match.arg(preconditioner)
  correction     <- match.arg(correction)
  stopifnot(is.null(stepsize) || (is.numeric(stepsize) && stepsize > 0),
            acceptTarget > 0, acceptTarget < 1, ridge >= 0,
            daGamma > 0, daT0 >= 0, daKappa > 0.5, daKappa <= 1)
  if (preconditioner == "fixed" && is.null(GFixed))
    stop("langevinControl: preconditioner = 'fixed' requires GFixed.")
  if (correction == "full" && preconditioner != "local")
    stop("langevinControl: correction = 'full' requires preconditioner = 'local'.")
  structure(list(stepsize       = stepsize,
                 acceptTarget   = acceptTarget,
                 preconditioner = preconditioner,
                 correction     = correction,
                 GFixed         = GFixed,
                 ridge          = ridge,
                 daGamma        = daGamma, daT0 = daT0, daKappa = daKappa),
            class = c("langevinControl", "list"))
}


#' @export
#' @rdname mcmcControl
hmcControl <- function(stepsize     = NULL,
                       leapfrogSteps = 10L,
                       massMatrix    = c("identity", "fisher", "fixed"),
                       MFixed        = NULL,
                       acceptTarget  = 0.8,
                       ridge         = 1e-8,
                       daGamma = 0.05, daT0 = 10, daKappa = 0.75) {
  massMatrix <- match.arg(massMatrix)
  stopifnot(is.null(stepsize) || (is.numeric(stepsize) && stepsize > 0),
            leapfrogSteps >= 1, acceptTarget > 0, acceptTarget < 1,
            ridge >= 0)
  if (massMatrix == "fixed" && is.null(MFixed))
    stop("hmcControl: massMatrix = 'fixed' requires MFixed.")
  structure(list(stepsize      = stepsize,
                 leapfrogSteps = as.integer(leapfrogSteps),
                 massMatrix    = massMatrix,
                 MFixed        = MFixed,
                 acceptTarget  = acceptTarget,
                 ridge         = ridge,
                 daGamma       = daGamma, daT0 = daT0, daKappa = daKappa),
            class = c("hmcControl", "list"))
}


#' @export
#' @rdname mcmcControl
nutsControl <- function(stepsize     = NULL,
                        maxTreeDepth = 10L,
                        deltaMax     = 1000.0,
                        massMatrix   = c("identity", "fisher", "fixed"),
                        MFixed       = NULL,
                        acceptTarget = 0.8,
                        ridge        = 1e-8,
                        daGamma = 0.05, daT0 = 10, daKappa = 0.75) {
  massMatrix <- match.arg(massMatrix)
  stopifnot(is.null(stepsize) || (is.numeric(stepsize) && stepsize > 0),
            maxTreeDepth >= 1, deltaMax > 0,
            acceptTarget > 0, acceptTarget < 1, ridge >= 0)
  if (massMatrix == "fixed" && is.null(MFixed))
    stop("nutsControl: massMatrix = 'fixed' requires MFixed.")
  structure(list(stepsize     = stepsize,
                 maxTreeDepth = as.integer(maxTreeDepth),
                 deltaMax     = deltaMax,
                 massMatrix   = massMatrix,
                 MFixed       = MFixed,
                 acceptTarget = acceptTarget,
                 ridge        = ridge,
                 daGamma      = daGamma, daT0 = daT0, daKappa = daKappa),
            class = c("nutsControl", "list"))
}


#' @export
#' @rdname mcmcControl
pathIntegralControl <- function(P = 8L, hbar = 1.0) {
  structure(list(P = as.integer(P), hbar = as.numeric(hbar)),
            class = c("pathIntegralControl", "list"))
}


#' @export
#' @rdname mcmcControl
smcControl <- function(nParticles       = 1000L,
                       essThreshold     = 0.5,
                       malaSteps        = 5L,
                       schedule         = list(),
                       parallelStrategy = c("forks", "serial"),
                       continuousAdaption = TRUE,
                       verbose = FALSE, blather = FALSE,
                       stepsize = NULL) {
  parallelStrategy <- match.arg(parallelStrategy)
  stopifnot(nParticles >= 4L, essThreshold > 0, essThreshold < 1,
            malaSteps >= 0L)
  structure(list(nParticles         = as.integer(nParticles),
                 essThreshold       = essThreshold,
                 malaSteps          = as.integer(malaSteps),
                 schedule           = schedule,
                 parallelStrategy   = parallelStrategy,
                 continuousAdaption = isTRUE(continuousAdaption),
                 verbose            = isTRUE(verbose),
                 blather            = isTRUE(blather),
                 stepsize           = stepsize),
            class = c("smcControl", "list"))
}


#' @export
#' @rdname mcmcControl
metricControl <- function(metricContext = NULL,
                          ridge         = 1e-8,
                          fisherMode    = c("gn", "exact")) {
  fisherMode <- match.arg(fisherMode)
  if (!is.null(metricContext) && !inherits(metricContext, "metricContext"))
    stop("metricControl: metricContext must come from metricContext().")
  structure(list(metricContext = metricContext,
                 ridge         = ridge,
                 fisherMode    = fisherMode),
            class = c("metricControl", "list"))
}


.resolve_bounds <- function(par_names, parupper, parlower) {
  K <- length(par_names)
  upper <- rep(Inf,  K); lower <- rep(-Inf, K)
  names(upper) <- names(lower) <- par_names
  if (!is.null(parupper)) {
    if (length(parupper) == 1L && is.null(names(parupper))) {
      upper[] <- parupper
    } else if (!is.null(names(parupper))) {
      idx <- intersect(names(parupper), par_names)
      upper[idx] <- parupper[idx]
    } else if (length(parupper) == K) {
      upper[] <- parupper
    }
  }
  if (!is.null(parlower)) {
    if (length(parlower) == 1L && is.null(names(parlower))) {
      lower[] <- parlower
    } else if (!is.null(names(parlower))) {
      idx <- intersect(names(parlower), par_names)
      lower[idx] <- parlower[idx]
    } else if (length(parlower) == K) {
      lower[] <- parlower
    }
  }
  list(upper = upper, lower = lower)
}


.resolve_parscale <- function(par_names, parscale) {
  K <- length(par_names)
  if (is.null(parscale)) return(rep(1, K))
  if (length(parscale) == 1L) return(rep(as.numeric(parscale), K))
  if (!is.null(names(parscale))) {
    idx <- match(par_names, names(parscale))
    if (anyNA(idx)) stop("parscale names do not cover parinit.")
    return(as.numeric(parscale[idx]))
  }
  if (length(parscale) == K) return(as.numeric(parscale))
  stop("parscale length does not match parinit.")
}


# Gelman-Rubin R-hat over the common tail of `n` samples per chain.
.computeRHat <- function(chains, par_names, n) {
  M <- length(chains)
  rHat <- setNames(numeric(length(par_names)), par_names)
  for (pj in par_names) {
    chain_means <- vapply(chains, function(c) {
      x <- c$samples[, pj]
      mean(utils::tail(x, n))
    }, 0.0)
    chain_vars <- vapply(chains, function(c) {
      x <- utils::tail(c$samples[, pj], n)
      stats::var(x)
    }, 0.0)
    grand_mean <- mean(chain_means)
    W <- mean(chain_vars)
    B <- (n / (M - 1)) * sum((chain_means - grand_mean)^2)
    if (!is.finite(W) || W <= 0) { rHat[pj] <- NA_real_; next }
    var_hat <- ((n - 1) / n) * W + B / n
    rHat[pj] <- sqrt(var_hat / W)
  }
  rHat
}


# Bake dots + per-call extras into a closure callable as
# objfun(pars, deriv, deriv2). `extra` threads per-block locals
# (.pars_full / .Omega_inv for the joint NLME path) that must be
# re-evaluated by the caller's R closure each iteration.
.bake_objfun <- function(raw_objfun, dots = list(), extra = list()) {
  force(raw_objfun)
  force(dots)
  function(pars, deriv = TRUE, deriv2 = FALSE) {
    args <- c(list(pars = pars), dots, extra,
              list(deriv = deriv, deriv2 = deriv2))
    do.call(raw_objfun, args)
  }
}


.build_mass <- function(M_user, K, par_names, ridge = 1e-8) {
  if (is.null(M_user)) M_user <- diag(K)
  if (!is.matrix(M_user) || nrow(M_user) != K || ncol(M_user) != K)
    stop(".build_mass: mass matrix must be K x K (K = ", K, ").")
  dimnames(M_user) <- list(par_names, par_names)
  cholM <- tryCatch(chol(M_user + ridge * diag(K)),
                    error = function(e) NULL)
  if (is.null(cholM))
    stop(".build_mass: mass matrix not positive-definite even after ridge.")
  Minv <- chol2inv(cholM)
  cholMi <- chol(Minv + ridge * diag(K))
  list(MCholUpper    = cholM,
       MinvCholUpper = cholMi)
}


.run_single_chain <- function(objfun_baked, parinit, n, warmup, moveType,
                              moveControl, metricControl, bounds, parscale,
                              dG_cb) {
  K <- length(parinit)
  par_names <- names(parinit)
  ridge <- moveControl$ridge %||% (metricControl$ridge %||% 1e-8)
  control <- list(ridge        = ridge,
                  daGamma      = moveControl$daGamma %||% 0.05,
                  daT0         = moveControl$daT0    %||% 10,
                  daKappa      = moveControl$daKappa %||% 0.75,
                  acceptTarget = moveControl$acceptTarget %||% 0.574,
                  stepsize     = moveControl$stepsize)

  if (moveType == "mh") {
    move_code <- 0L
    pc <- switch(moveControl$proposalCov, identity = 2L, fixed = 1L,
                                            fisher = 0L)
    control$preconditioner <- pc
    control$correctionFull <- FALSE
    control$GFixed <- if (pc == 1L)
      moveControl$covFixed * 1
    else NULL
  } else if (moveType == "langevin") {
    move_code <- 1L
    pc <- switch(moveControl$preconditioner, local = 0L, fixed = 1L,
                                              identity = 2L)
    control$preconditioner <- pc
    control$correctionFull <- (moveControl$correction == "full")
    control$GFixed <- if (pc == 1L) moveControl$GFixed else NULL
  } else if (moveType == "hmc") {
    move_code <- 2L
    control$preconditioner <- 2L
    control$correctionFull <- FALSE
    M <- switch(moveControl$massMatrix,
                identity = diag(K),
                fixed    = moveControl$MFixed,
                fisher   = .fisher_at(objfun_baked, parinit, K))
    mass <- .build_mass(M, K, par_names, ridge)
    control$MCholUpper    <- mass$MCholUpper
    control$MinvCholUpper <- mass$MinvCholUpper
    control$leapfrogSteps <- moveControl$leapfrogSteps
  } else if (moveType == "nuts") {
    move_code <- 3L
    control$preconditioner <- 2L
    control$correctionFull <- FALSE
    M <- switch(moveControl$massMatrix,
                identity = diag(K),
                fixed    = moveControl$MFixed,
                fisher   = .fisher_at(objfun_baked, parinit, K))
    mass <- .build_mass(M, K, par_names, ridge)
    control$MCholUpper    <- mass$MCholUpper
    control$MinvCholUpper <- mass$MinvCholUpper
    control$maxTreeDepth  <- moveControl$maxTreeDepth
    control$deltaMax      <- moveControl$deltaMax
  } else {
    stop(".run_single_chain: unknown moveType = '", moveType, "'.")
  }

  raw <- mcmcChainRun(objfun_baked, parinit, as.integer(n),
                       as.integer(warmup), move_code, control,
                       bounds, parscale, dG_cb)
  raw
}


# Single-call Fisher info evaluation for HMC/NUTS "fisher" mass matrix
# initialisation. Returns G = hessian/2 at parinit.
.fisher_at <- function(objfun, parinit, K) {
  out <- tryCatch(objfun(pars = parinit, deriv = TRUE, deriv2 = FALSE),
                  error = function(e) NULL)
  if (is.null(out$hessian))
    stop(".fisher_at: objfn returned no hessian. Use massMatrix = ",
         "'identity' or 'fixed' instead of 'fisher'.")
  H <- as.matrix(out$hessian) / 2
  dimnames(H) <- list(names(parinit), names(parinit))
  H
}


.run_smc <- function(target, n_arg, sequenceSchedule, populationSize,
                     reSampling, moveType, moveControl, metricControl,
                     sequenceControl, bounds, parscale, dG_cb, cores) {

  if (sequenceSchedule == "learned")
    stop("sequenceSchedule = 'learned' is not yet implemented.")
  if (moveType == "pathIntegral")
    stop("moveType = 'pathIntegral' is not yet implemented.")

  likObj      <- target$likObj
  priorObj    <- target$priorObj
  priorSample <- target$priorSample
  if (!is.function(priorSample))
    stop(".run_smc: target$priorSample must be a function(n) -> n x K matrix.")

  N <- populationSize

  particles <- as.matrix(priorSample(N))
  if (nrow(particles) != N)
    stop(".run_smc: priorSample(N) must return exactly N rows.")
  par_names <- colnames(particles)
  if (is.null(par_names))
    stop(".run_smc: priorSample output must carry parameter names as colnames.")
  K <- ncol(particles)

  call_lik <- function(theta) {
    tryCatch(likObj(pars = theta, deriv = TRUE, deriv2 = FALSE),
             error = function(e) NULL)
  }

  init_one <- function(i) {
    out <- call_lik(setNames(particles[i, ], par_names))
    if (is.null(out) || !is.finite(out$value)) -Inf
    else -as.numeric(out$value) / 2
  }
  logL <- unlist(.parallelLapply(seq_len(N), init_one, cores = cores,
                                  extraExports = c("particles", "par_names",
                                                    "call_lik")))
  if (sum(!is.finite(logL)) == N)
    stop(".run_smc: initial likelihood non-finite for all particles.")

  logw         <- rep(0.0, N)
  beta         <- 0.0
  betaPath     <- 0.0
  ESSPath      <- as.numeric(N)
  acceptRates  <- numeric(0)
  stepsizePath <- numeric(0)
  logEvidence  <- 0.0
  level        <- 0L
  eps          <- sequenceControl$stepsize
  if (is.null(eps)) eps <- .smc_initial_stepsize(particles)

  history_particles <- if (sequenceControl$blather) list(particles) else NULL
  history_weights   <- if (sequenceControl$blather)
                         list(rep(1 / N, N)) else NULL

  targetESS  <- sequenceControl$essThreshold * N
  malaSteps  <- sequenceControl$malaSteps

  # Fixed schedule (if requested).
  fixed_betas <- NULL
  if (sequenceSchedule == "fixed") {
    if (!length(sequenceControl$schedule))
      stop("sequenceControl$schedule must be a non-empty numeric vector ",
           "for sequenceSchedule = 'fixed'.")
    fixed_betas <- as.numeric(sequenceControl$schedule)
    if (fixed_betas[1] != 0 || tail(fixed_betas, 1) != 1)
      stop("sequenceControl$schedule must start at 0 and end at 1.")
  }
  beta_step_idx <- 1L

  while (beta < 1.0) {
    level <- level + 1L

    if (sequenceSchedule == "adaptiveEss") {
      delta <- smcBetaBisect(logL = logL, logwPrev = logw,
                             betaOld = beta, targetESS = targetESS)
      if (delta <= 0) delta <- min(1 - beta, 1e-3)
      betaNew <- min(1.0, beta + delta)
    } else {
      beta_step_idx <- beta_step_idx + 1L
      betaNew <- fixed_betas[beta_step_idx]
    }

    rw <- mcmcSmcReweight(logL, logw, beta, betaNew)
    logEvidence <- logEvidence + rw$logZinc
    logw <- rw$logw
    beta <- betaNew

    w <- exp(logw); w <- w / sum(w)
    ix <- switch(reSampling,
                 systematic   = smcSystematicResample(weights = w, u = stats::runif(1L)),
                 stratified   = mcmcStratifiedResample(weights = w),
                 residual     = mcmcResidualResample(weights = w),
                 multinomial  = mcmcMultinomialResample(weights = w),
                 none         = seq_len(N))
    particles <- particles[ix, , drop = FALSE]
    logL      <- logL[ix]
    logw      <- rep(0.0, N)

    # ---- MCMC move per particle ---------------------------------------
    target_temp <- .smc_make_tempered_obj(likObj, priorObj, beta)
    move_one <- function(i) {
      theta_i <- setNames(particles[i, ], par_names)
      bake <- .bake_objfun(target_temp)
      raw <- tryCatch(
        .run_single_chain(bake, theta_i, n = malaSteps, warmup = 0L,
                          moveType, moveControl, metricControl,
                          bounds, parscale, dG_cb),
        error = function(e) NULL)
      if (is.null(raw)) return(list(theta = theta_i, accepts = 0L))
      list(theta = raw$samples[nrow(raw$samples), ],
           accepts = sum(raw$accept))
    }
    moves <- .parallelLapply(seq_len(N), move_one, cores = cores,
                              extraExports = c("particles", "par_names",
                                                "target_temp", "moveType",
                                                "moveControl", "metricControl",
                                                "bounds", "parscale", "dG_cb",
                                                "malaSteps"))
    for (i in seq_len(N)) particles[i, ] <- moves[[i]]$theta

    refresh_one <- function(i) {
      out <- call_lik(setNames(particles[i, ], par_names))
      if (is.null(out) || !is.finite(out$value)) -Inf
      else -as.numeric(out$value) / 2
    }
    logL <- unlist(.parallelLapply(seq_len(N), refresh_one, cores = cores,
                                    extraExports = c("particles", "par_names",
                                                      "call_lik")))

    accCount <- sum(vapply(moves, `[[`, integer(1), "accepts"))
    accRate  <- if (malaSteps > 0L) accCount / (N * malaSteps) else NA_real_
    if (is.finite(accRate) && sequenceControl$continuousAdaption) {
      if (accRate < 0.30)      eps <- eps / 1.5
      else if (accRate > 0.75) eps <- eps * 1.5
      eps <- max(eps, 1e-8)
      moveControl$stepsize <- eps
    }

    betaPath     <- c(betaPath, beta)
    ESSPath      <- c(ESSPath, N)
    acceptRates  <- c(acceptRates, accRate)
    stepsizePath <- c(stepsizePath, eps)

    if (sequenceControl$blather) {
      history_particles <- c(history_particles, list(particles))
      history_weights   <- c(history_weights, list(rep(1 / N, N)))
    }
    if (sequenceControl$verbose)
      message(sprintf("[mcmc/smc] level %d: beta=%.4f, acc=%.3f, eps=%.3g, logZ=%.3f",
                      level, beta, accRate, eps, logEvidence))
  }

  out <- list(samples      = particles,
              weights      = rep(1 / N, N),
              logEvidence  = logEvidence,
              betaPath     = betaPath,
              ESSPath      = ESSPath,
              acceptRates  = acceptRates,
              stepsizePath = stepsizePath,
              nLevels      = level)
  if (sequenceControl$blather) {
    out$particlesHistory <- history_particles
    out$weightsHistory   <- history_weights
  }
  out
}


# Tempered objfn at temperature beta: value/grad/hess = beta * lik + prior.
.smc_make_tempered_obj <- function(likObj, priorObj, beta) {
  force(likObj); force(priorObj); force(beta)
  fn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                 conditions = NULL, env = NULL) {
    pars <- ..1
    args_lik <- list(pars = pars, fixed = fixed, deriv = deriv,
                     deriv2 = deriv2, conditions = conditions, env = env)
    out_lik <- do.call(likObj, args_lik)
    out <- list(value    = beta * as.numeric(out_lik$value),
                gradient = beta * as.numeric(out_lik$gradient),
                hessian  = beta * as.matrix(out_lik$hessian))
    names(out$gradient)   <- names(out_lik$gradient)
    dimnames(out$hessian) <- dimnames(out_lik$hessian)
    if (!is.null(priorObj)) {
      args_pr <- list(pars = pars, fixed = fixed, deriv = deriv,
                      deriv2 = deriv2, conditions = conditions, env = env)
      out_pr <- do.call(priorObj, args_pr)
      out <- .smc_add_obj(out, out_pr)
    }
    class(out) <- c("objlist", "list")
    attr(out, "env") <- attr(out_lik, "env")
    out
  }
  class(fn) <- c("objfn", "fn")
  attr(fn, "conditions") <- attr(likObj, "conditions")
  attr(fn, "parameters") <- union(attr(likObj, "parameters"),
                                   attr(priorObj, "parameters"))
  fn
}


.smc_add_obj <- function(a, b) {
  if (is.null(a$gradient) || is.null(b$gradient)) {
    a$value <- a$value + b$value
    return(a)
  }
  gn <- union(names(a$gradient), names(b$gradient))
  g  <- setNames(numeric(length(gn)), gn)
  g[names(a$gradient)] <- g[names(a$gradient)] + a$gradient
  g[names(b$gradient)] <- g[names(b$gradient)] + b$gradient
  H <- matrix(0, length(gn), length(gn), dimnames = list(gn, gn))
  if (!is.null(a$hessian))
    H[rownames(a$hessian), colnames(a$hessian)] <-
      H[rownames(a$hessian), colnames(a$hessian)] + a$hessian
  if (!is.null(b$hessian))
    H[rownames(b$hessian), colnames(b$hessian)] <-
      H[rownames(b$hessian), colnames(b$hessian)] + b$hessian
  list(value = a$value + b$value, gradient = g, hessian = H)
}


.smc_initial_stepsize <- function(particles) {
  K <- ncol(particles)
  v <- apply(particles, 2L, stats::var)
  v[!is.finite(v) | v <= 0] <- 1
  max(1e-6, 0.5 * mean(v) * K^(-1 / 3))
}


.run_pgibbs <- function(target, n_arg, moveType, moveControl, metricControl,
                        sequenceControl, bounds, parscale, dG_cb, cores,
                        warmup, n_keep) {

  omegaSpec     <- target$omegaSpec
  prdfn         <- target$prdfn
  data          <- target$data
  errfn         <- target$errfn
  obj           <- target$jointObj
  priorTheta    <- target$priorTheta
  priorOmegaObj <- target$priorOmegaObj
  priorSample   <- target$priorSample
  init          <- target$init
  N             <- target$N
  K             <- target$K
  subjects      <- target$subjects
  meta          <- target$subjectMeta
  eta_names_all <- target$etaNamesAll
  cholPars      <- target$cholPars
  outer_names   <- target$outerNames
  structural_names <- target$structuralNames

  parsFull <- setNames(numeric(length(meta$pars_full_names)),
                       meta$pars_full_names)
  init_draw <- as.numeric(priorSample(1L)[1, ])
  names(init_draw) <- colnames(priorSample(1L))
  if (!all(outer_names %in% names(init_draw)))
    stop(".run_pgibbs: priorSample must include all outer names.")
  parsFull[outer_names] <- init_draw[outer_names]
  # etas start at 0

  subjEtaObjList <- lapply(seq_len(N), function(i) {
    .makeSubjectEtaObj(i, meta, prdfn, errfn, parsFull,
                       meta$eta_idx_global[i, ])
  })
  jointThetaOmegaObj <- .makeJointThetaOmegaObj(
      obj, priorTheta, priorOmegaObj, parsFull, outer_names)

  total <- as.integer(n_keep) + as.integer(warmup)
  burnin <- as.integer(warmup)
  P_outer <- length(outer_names)
  malaSteps <- sequenceControl$malaSteps %||% 5L

  eta_step    <- rep(0.1, N)
  outer_step  <- 0.05

  samples_outer <- matrix(NA_real_, total, P_outer,
                          dimnames = list(NULL, outer_names))
  samples_eta   <- array(NA_real_, c(total, N, K),
                         dimnames = list(NULL, subjects, omegaSpec$eta))
  logp_path    <- numeric(total)
  accept_outer <- numeric(total)
  accept_eta   <- matrix(NA_real_, total, N,
                         dimnames = list(NULL, subjects))

  # Per-block move controls (clone the user's moveControl per block).
  mc_eta   <- moveControl
  mc_outer <- moveControl

  for (sweep in seq_len(total)) {
    omegaChol_vec <- parsFull[cholPars]
    L_om          <- omegaSpec$buildL(omegaChol_vec)
    Omega_inv     <- chol2inv(t(L_om))
    Omega_log_det <- 2 * sum(log(diag(L_om)))

    move_subj <- function(i) {
      eta_init <- parsFull[meta$eta_names[[i]]]
      objSubj  <- subjEtaObjList[[i]]
      mc_eta_i <- mc_eta
      mc_eta_i$stepsize <- eta_step[i]
      bake <- .bake_objfun(objSubj, dots = list(),
                           extra = list(.pars_full = parsFull,
                                         .Omega_inv = Omega_inv,
                                         .Omega_log_det = Omega_log_det))
      par_names_i <- meta$eta_names[[i]]
      bounds_i <- list(upper = rep(Inf, K),  lower = rep(-Inf, K))
      raw <- tryCatch(
        .run_single_chain(bake, eta_init, n = malaSteps, warmup = 0L,
                          moveType, mc_eta_i, metricControl,
                          bounds_i, rep(1, K), dG_cb),
        error = function(e) NULL)
      if (is.null(raw)) return(list(eta = eta_init, accept = 0))
      list(eta    = raw$samples[malaSteps, ],
           accept = mean(raw$accept))
    }
    eta_results <- .parallelLapply(seq_len(N), move_subj, cores = cores,
                                   extraExports = c("parsFull", "meta",
                                                     "subjEtaObjList",
                                                     "Omega_inv",
                                                     "Omega_log_det",
                                                     "eta_step", "malaSteps",
                                                     "mc_eta", "moveType",
                                                     "metricControl", "K"))
    for (i in seq_len(N)) {
      parsFull[meta$eta_names[[i]]] <- eta_results[[i]]$eta
      accept_eta[sweep, i] <- eta_results[[i]]$accept
      a <- eta_results[[i]]$accept
      if (is.finite(a)) {
        if (a < 0.3)      eta_step[i] <- eta_step[i] / 1.5
        else if (a > 0.75) eta_step[i] <- eta_step[i] * 1.5
        eta_step[i] <- max(eta_step[i], 1e-6)
      }
    }

    # Outer block
    outer_init   <- parsFull[outer_names]
    current_etas <- parsFull[eta_names_all]
    mc_outer_step <- mc_outer
    mc_outer_step$stepsize <- outer_step
    bake_outer <- .bake_objfun(jointThetaOmegaObj, dots = list(),
                                extra = list(.etas = current_etas))
    bounds_outer <- list(upper = rep(Inf, P_outer),
                         lower = rep(-Inf, P_outer))
    chainB <- tryCatch(
      .run_single_chain(bake_outer, outer_init, n = malaSteps, warmup = 0L,
                        moveType, mc_outer_step, metricControl,
                        bounds_outer, rep(1, P_outer), dG_cb),
      error = function(e) NULL)
    if (!is.null(chainB)) {
      parsFull[outer_names] <- chainB$samples[malaSteps, ]
      a_outer <- mean(chainB$accept)
    } else {
      a_outer <- 0
    }
    accept_outer[sweep] <- a_outer
    if (is.finite(a_outer)) {
      if (a_outer < 0.3)      outer_step <- outer_step / 1.5
      else if (a_outer > 0.75) outer_step <- outer_step * 1.5
      outer_step <- max(outer_step, 1e-6)
    }

    samples_outer[sweep, ] <- parsFull[outer_names]
    for (i in seq_len(N))
      samples_eta[sweep, i, ] <- parsFull[meta$eta_names[[i]]]
    logp_path[sweep] <- if (!is.null(chainB))
      -as.numeric(chainB$logp[malaSteps]) else NA_real_
  }

  keep <- if (burnin >= total) seq_len(total) else (burnin + 1L):total
  list(samples         = samples_outer[keep, , drop = FALSE],
       etaSamples      = samples_eta[keep, , , drop = FALSE],
       logp            = logp_path[keep],
       acceptOuter     = accept_outer[keep],
       acceptEta       = accept_eta[keep, , drop = FALSE],
       nSweeps         = total,
       burnin          = burnin,
       omegaSpec       = omegaSpec,
       structuralNames = structural_names)
}


.run_chains <- function(runOne, nChains, seeds, chainCores, par_names) {
  if (is.null(seeds)) seeds <- seq_len(nChains)
  if (length(seeds) != nChains)
    stop(".run_chains: seeds must have length nChains.")
  chains <- .parallelLapply(seq_len(nChains), function(idx) {
    set.seed(seeds[idx])
    runOne(idx)
  }, cores = chainCores,
     extraExports = c("seeds", "runOne"))

  per_chain_n <- vapply(chains, function(c) nrow(c$samples), 0L)
  if (length(unique(per_chain_n)) != 1L)
    warning(".run_chains: chains differ in sample count; padding with NA.")
  all_samples <- do.call(rbind, lapply(chains, `[[`, "samples"))
  chain_ids   <- rep(seq_len(nChains), times = per_chain_n)

  n_per <- min(per_chain_n)
  rHat <- if (n_per >= 4L) .computeRHat(chains, par_names, n_per)
          else setNames(rep(NA_real_, length(par_names)), par_names)
  logEvidence <- vapply(chains, function(c) {
    v <- c$logEvidence
    if (is.null(v)) NA_real_ else as.numeric(v)
  }, 0.0)

  list(samples     = all_samples,
       chainId     = chain_ids,
       rHat        = rHat,
       logEvidence = logEvidence,
       seeds       = seeds,
       nChains     = nChains,
       chains      = chains)
}


#' Posterior sampling for dMod objective functions
#'
#' MCMC sampler with single-chain (`sequenceType = "single"`) and SMC
#' (`sequenceType = "sequential"`) modes. A bare `objfn` passed as `target`
#' is auto-wrapped to a [flatTarget()]; hierarchical NLME posteriors are
#' built via [bayesNLMEMarginal()] or [bayesNLMEJoint()].
#'
#' @param target Either an `objfn`, a [flatTarget()], or one of the
#'   NLME target constructors. An `objfn` is silently wrapped to a flat
#'   target (no prior); pass an explicit [flatTarget()] when you want to
#'   supply `priorObj` / `priorSample` for sequential SMC.
#' @param sequenceType `"single"` (one or more independent chains, each
#'   moving with the chosen `moveType`) or `"sequential"` (SMC tempering
#'   along an annealed-likelihood path, requires `priorSample` in the
#'   target).
#' @param sequenceSchedule For sequential sampling: `"fixed"` (uses
#'   `smcControl()$schedule`), `"adaptiveEss"` (ESS-bisection;
#'   recommended default), or `"learned"` (slot only; not implemented).
#' @param populationSize Number of particles for SMC.
#' @param reSampling Per-level particle resampling scheme.
#' @param moveType Inner move kernel: `"mh"` / `"langevin"` / `"hmc"` /
#'   `"nuts"` are real implementations; `"pathIntegral"` is a slot.
#' @param metric Metric tensor for `langevin` / `hmc` / `nuts`:
#'   `"euclidean"` (G = I), `"riemannFisher"` (G = local Fisher info),
#'   `"riemannHessian"` (placeholder, falls back to Fisher), `"fixed"`
#'   (user-supplied via the move control).
#' @param parinit Optional named numeric. Required for `"single"` mode;
#'   unused in `"sequential"` (initial particles come from `priorSample`).
#' @param nIter Post-warmup samples per chain (single mode) or per
#'   particle (joint Particle-Gibbs).
#' @param warmup Warmup iterations (single mode) or burn-in sweeps
#'   (joint Particle-Gibbs). Ignored for `"sequential"` (tempering plays
#'   the role of warmup); a warning is emitted if non-zero.
#' @param chains Number of independent chains. Setting `chains > 1`
#'   wraps the sampler with multi-chain orchestration and computes
#'   Gelman-Rubin R-hat.
#' @param cores Parallel workers for particle / chain moves.
#' @param chainCores Outer-axis cores when `chains > 1`. Defaults to
#'   `min(chains, parallel::detectCores())`.
#' @param seeds Optional integer vector of length `chains` to seed each
#'   chain. Defaults to `seq_len(chains)`.
#' @param parlower,parupper,parscale Forwarded to the inner kernel.
#'   Bounds are enforced by immediate rejection of out-of-box proposals.
#' @param sequenceControl,moveControl,metricControl Typed control objects
#'   from the per-axis constructors. Sensible defaults are used when the
#'   user passes `NULL` (which is the default).
#' @param ... Additional arguments forwarded as dots into the objfn.
#'
#' @return An object of class `c("mcmcResult", ...)`, subclassed to
#'   `mcmcResultSingle`, `mcmcResultSequential` (carries `logEvidence`,
#'   `betaPath`, `ESSPath`, `acceptRates`, `stepsizePath`, `nLevels`),
#'   `mcmcResultBlocked` (adds `etaSamples`, `acceptOuter`, `acceptEta`),
#'   or `mcmcResultMulti` (adds `chainId`, `rHat`) depending on the run
#'   configuration.
#'
#' @seealso [flatTarget()], [bayesNLMEMarginal()], [bayesNLMEJoint()],
#'   [mcmcControl], [priorOmega()], [metricContext()]
#' @export
mcmc <- function(target,
                 sequenceType     = c("single", "sequential"),
                 sequenceSchedule = c("adaptiveEss", "fixed", "learned"),
                 populationSize   = 1L,
                 reSampling       = c("systematic", "stratified",
                                      "residual", "multinomial", "none"),
                 moveType         = c("langevin", "mh", "hmc", "nuts",
                                       "pathIntegral"),
                 metric           = c("euclidean", "riemannFisher",
                                      "riemannHessian", "fixed"),
                 parinit          = NULL,
                 nIter            = 1000L,
                 warmup           = NULL,
                 chains           = 1L,
                 cores            = 1L,
                 chainCores       = NULL,
                 seeds            = NULL,
                 parlower         = NULL,
                 parupper         = NULL,
                 parscale         = NULL,
                 sequenceControl  = NULL,
                 moveControl      = NULL,
                 metricControl    = NULL,
                 ...) {

  sequenceType     <- match.arg(sequenceType)
  sequenceSchedule <- match.arg(sequenceSchedule)
  reSampling       <- match.arg(reSampling)
  moveType         <- match.arg(moveType)
  metric           <- match.arg(metric)
  call_capture     <- sys.call()

  # Default warmup: 500 for single chains, 0 for SMC tempering.
  if (is.null(warmup))
    warmup <- if (sequenceType == "sequential") 0L else 500L

  if (moveType == "pathIntegral")
    stop("moveType = 'pathIntegral' is a slot for future work and is not ",
         "yet implemented. Use 'mh' / 'langevin' / 'hmc' / 'nuts'.")
  if (sequenceSchedule == "learned")
    stop("sequenceSchedule = 'learned' is a slot for future work and is ",
         "not yet implemented. Use 'fixed' or 'adaptiveEss'.")

  # Resolve default controls.
  if (is.null(moveControl)) {
    moveControl <- switch(moveType,
                          mh        = mhControl(),
                          langevin  = langevinControl(),
                          hmc       = hmcControl(),
                          nuts      = nutsControl())
  }
  expected_class <- paste0(moveType, "Control")
  if (!inherits(moveControl, expected_class))
    stop("mcmc: moveControl is not of class '", expected_class, "'.")
  if (is.null(metricControl)) metricControl <- metricControl()
  if (is.null(sequenceControl)) sequenceControl <- smcControl()

  # Wire metric into the move control where it changes the kernel's
  # preconditioner / mass-matrix choice.
  if (moveType == "langevin") {
    if (metric == "euclidean")        moveControl$preconditioner <- "identity"
    else if (metric == "fixed")       moveControl$preconditioner <- "fixed"
    else                              moveControl$preconditioner <- "local"
    if (metric != "fixed") moveControl$GFixed <- NULL
  } else if (moveType %in% c("hmc", "nuts")) {
    if (metric == "euclidean")       moveControl$massMatrix <- "identity"
    else if (metric == "fixed")      moveControl$massMatrix <- "fixed"
    else                             moveControl$massMatrix <- "fisher"
  }

  # Normalise the target.
  target_obj <- .normalise_target(target)

  # Branch: blocked (Particle-Gibbs) target.
  if (inherits(target_obj, "bayesNLMEJoint")) {
    if (sequenceType != "single")
      stop("mcmc: bayesNLMEJoint targets only support sequenceType = 'single'.")
    if (!is.null(chains) && chains > 1L) {
      out <- .run_chains(function(idx) {
        raw <- .run_pgibbs(target_obj, nIter, moveType, moveControl,
                            metricControl, sequenceControl,
                            list(upper = rep(Inf, length(target_obj$outerNames)),
                                 lower = rep(-Inf, length(target_obj$outerNames))),
                            rep(1, length(target_obj$outerNames)),
                            NULL, cores, warmup, nIter)
        .finish_pgibbs(raw, call_capture)
      }, chains, seeds, chainCores %||% min(chains, parallel::detectCores()),
         target_obj$outerNames)
      class(out) <- c("mcmcResultMulti", "mcmcResultBlocked",
                       "mcmcResult", "list")
      out$method <- "blocked"
      out$call <- call_capture
      return(out)
    }
    raw <- .run_pgibbs(target_obj, nIter, moveType, moveControl,
                        metricControl, sequenceControl,
                        list(upper = rep(Inf, length(target_obj$outerNames)),
                             lower = rep(-Inf, length(target_obj$outerNames))),
                        rep(1, length(target_obj$outerNames)),
                        NULL, cores, warmup, nIter)
    return(.finish_pgibbs(raw, call_capture))
  }

  # Flat target -> single chain or sequential.
  par_names <- target_obj$parNames
  if (is.null(par_names) && !is.null(parinit)) par_names <- names(parinit)

  if (sequenceType == "sequential") {
    if (is.null(target_obj$priorSample))
      stop("mcmc: sequenceType = 'sequential' requires target$priorSample.")
    bounds <- if (!is.null(par_names))
      .resolve_bounds(par_names, parupper, parlower)
    else list(upper = Inf, lower = -Inf)
    parsc <- if (!is.null(par_names))
      .resolve_parscale(par_names, parscale)
    else 1
    dG_cb <- .build_dG_cb(metricControl, moveType, moveControl)

    runOne <- function(idx) {
      raw <- .run_smc(target_obj, nIter, sequenceSchedule, populationSize,
                      reSampling, moveType, moveControl, metricControl,
                      sequenceControl, bounds, parsc, dG_cb, cores)
      .finish_smc(raw, call_capture, target_obj)
    }
    if (chains > 1L) {
      cc <- chainCores %||% min(chains, parallel::detectCores())
      out <- .run_chains(runOne, chains, seeds, cc,
                          par_names %||% colnames(target_obj$priorSample(1L)))
      class(out) <- c("mcmcResultMulti", "mcmcResultSequential",
                       "mcmcResult", "list")
      out$method <- "sequential"
      out$call <- call_capture
      return(out)
    }
    return(runOne(1L))
  }

  # Single-chain run with a flat target.
  if (is.null(parinit))
    stop("mcmc: parinit is required for sequenceType = 'single'.")
  par_names <- names(parinit)
  bounds <- .resolve_bounds(par_names, parupper, parlower)
  parsc  <- .resolve_parscale(par_names, parscale)
  dG_cb  <- .build_dG_cb(metricControl, moveType, moveControl)
  raw_obj <- if (!is.null(target_obj$priorObj))
    .smc_make_tempered_obj(target_obj$likObj, target_obj$priorObj, 1.0)
  else target_obj$likObj
  bake <- .bake_objfun(raw_obj, dots = list(...))

  runOne <- function(idx) {
    raw <- .run_single_chain(bake, parinit, nIter, warmup, moveType,
                              moveControl, metricControl, bounds, parsc,
                              dG_cb)
    .finish_single(raw, parinit, moveType, metric, sequenceType, call_capture)
  }
  if (chains > 1L) {
    cc <- chainCores %||% min(chains, parallel::detectCores())
    out <- .run_chains(runOne, chains, seeds, cc, par_names)
    class(out) <- c("mcmcResultMulti", "mcmcResultSingle",
                     "mcmcResult", "list")
    out$method <- "single"
    out$call <- call_capture
    return(out)
  }
  runOne(1L)
}


.normalise_target <- function(target) {
  if (inherits(target, "mcmcTarget")) return(target)
  if (inherits(target, "objfn"))
    return(flatTarget(likObj = target, priorObj = NULL, priorSample = NULL))
  stop("mcmc: target must be an objfn, a flatTarget(), or a ",
       "bayesNLMEMarginal()/bayesNLMEJoint() result.")
}


.build_dG_cb <- function(metricControl, moveType, moveControl) {
  if (moveType != "langevin") return(NULL)
  if (is.null(moveControl$correction) || moveControl$correction != "full")
    return(NULL)
  ctx <- metricControl$metricContext
  if (is.null(ctx))
    stop("mcmc: correction = 'full' requires metricControl(metricContext = ...).")
  force(ctx)
  function(theta) .computeMetricDerivative(ctx, theta)
}


.finish_single <- function(raw, parinit, moveType, metric, sequenceType,
                           call_capture) {
  out <- list(samples       = raw$samples,
              logp          = raw$logp,
              accept        = raw$accept,
              acceptRate    = raw$acceptRate,
              stepsize      = raw$stepsize,
              ess           = raw$ess,
              parinit       = parinit,
              moveType      = moveType,
              metric        = metric,
              sequenceType  = sequenceType,
              call          = call_capture)
  if (moveType == "nuts") out$treedepth <- raw$treedepth
  class(out) <- c("mcmcResultSingle", "mcmcResult", "list")
  out
}


.finish_smc <- function(raw, call_capture, target_obj) {
  out <- list(samples      = raw$samples,
              weights      = raw$weights,
              logEvidence  = raw$logEvidence,
              betaPath     = raw$betaPath,
              ESSPath      = raw$ESSPath,
              acceptRates  = raw$acceptRates,
              stepsizePath = raw$stepsizePath,
              nLevels      = raw$nLevels,
              call         = call_capture)
  if (!is.null(raw$particlesHistory)) out$particlesHistory <- raw$particlesHistory
  if (!is.null(raw$weightsHistory))   out$weightsHistory   <- raw$weightsHistory
  out_classes <- c("mcmcResultSequential", "mcmcResult", "list")
  if (inherits(target_obj, "bayesNLMEMarginal")) {
    out$omegaSpec       <- target_obj$omegaSpec
    out$structuralNames <- target_obj$structuralNames
    out_classes <- c("bayesNLMEMarginal", out_classes)
  }
  class(out) <- out_classes
  out
}


.finish_pgibbs <- function(raw, call_capture) {
  out <- list(samples         = raw$samples,
              etaSamples      = raw$etaSamples,
              logp            = raw$logp,
              acceptOuter     = raw$acceptOuter,
              acceptEta       = raw$acceptEta,
              nSweeps         = raw$nSweeps,
              burnin          = raw$burnin,
              omegaSpec       = raw$omegaSpec,
              structuralNames = raw$structuralNames,
              call            = call_capture)
  class(out) <- c("bayesNLMEJoint", "mcmcResultBlocked", "mcmcResult", "list")
  out
}
