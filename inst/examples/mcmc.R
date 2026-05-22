\dontrun{
  ## Bayesian posterior sampling examples for the unified mcmc() entry.
  ## Demonstrates three problem types:
  ##   (1) flat target (likelihood + Gaussian prior) with Langevin moves
  ##   (2) sequential SMC tempering on the same posterior
  ##   (3) hierarchical NLME: marginal (Pfad A) and joint (Pfad B)


  ## ----- Problem (1) and (2): ODE decay model ---------------------------
  reactions <- addReaction(eqnlist(), from = "A", to = "",
                           rate = "k*A", description = "decay")
  m   <- odemodel(reactions, modelname = "mcmc_decay_m",
                  solver = "CppODE", deriv2 = TRUE, compile = TRUE)
  xfn <- Xs(m)
  gfn <- Y(c(y = "A"), f = xfn, condition = NULL, attach.input = FALSE,
           modelname = "mcmc_decay_obs", compile = TRUE, deriv2 = TRUE)
  pfn <- Pexpl(c(A = "A", k = "k"), parameters = NULL,
               modelname = "mcmc_decay_p", compile = TRUE, deriv2 = TRUE,
               condition = "C1", derivMode = "symbolic")
  prd <- gfn * xfn * pfn

  times <- seq(0, 8, by = 1)
  truth <- prd(times, c(A = 1.0, k = 0.5))[[1]]
  yobs  <- truth[truth[, "time"] %in% times, "y"] +
    rnorm(length(times), 0, 0.08)
  data  <- datalist(C1 = data.frame(name = "y", time = times,
                                    value = yobs, sigma = 0.08))

  likObj   <- normL2(data, prd)
  priorObj <- constraintL2(c(A = 1.0, k = 0.5), sigma = 5)
  obj      <- likObj + priorObj

  ## MAP via trust gives a starting point + Hessian for the fixed-G chain
  mapfit <- trust(obj, c(A = 1.2, k = 0.4), rinit = 0.1, rmax = 10)


  ## ----- (1a) Single-chain Langevin with the MAP Hessian as fixed metric
  chainFixed <- mcmc(target = obj,
                     sequenceType = "single",
                     moveType     = "langevin",
                     metric       = "fixed",
                     moveControl  = langevinControl(GFixed = mapfit$hessian),
                     parinit      = mapfit$argument,
                     nIter        = 2000L, warmup = 800L)

  ## ----- (1b) Same problem, local-Fisher metric (recomputed each step)
  chainLocal <- mcmc(target = obj,
                     sequenceType = "single",
                     moveType     = "langevin",
                     metric       = "riemannFisher",
                     parinit      = mapfit$argument,
                     nIter        = 2000L, warmup = 800L)

  ## ----- (1c) HMC with Fisher mass matrix frozen at parinit
  chainHMC <- mcmc(target = obj,
                   sequenceType = "single",
                   moveType     = "hmc",
                   metric       = "riemannFisher",
                   moveControl  = hmcControl(leapfrogSteps = 8L),
                   parinit      = mapfit$argument,
                   nIter        = 1500L, warmup = 800L)

  ## ----- (1d) NUTS with Fisher mass matrix
  chainNUTS <- mcmc(target = obj,
                    sequenceType = "single",
                    moveType     = "nuts",
                    metric       = "riemannFisher",
                    parinit      = mapfit$argument,
                    nIter        = 1500L, warmup = 800L)

  rbind(MAP        = mapfit$argument,
        Fixed      = colMeans(chainFixed$samples),
        LocalMALA  = colMeans(chainLocal$samples),
        HMC        = colMeans(chainHMC$samples),
        NUTS       = colMeans(chainNUTS$samples))


  ## ----- (2) Sequential SMC over the same posterior, with log-evidence
  priorSampleFn <- function(n) {
    cbind(A = rnorm(n, 1.0, 5),
          k = rnorm(n, 0.5, 5))
  }
  smcTarget <- flatTarget(likObj      = likObj,
                          priorObj    = priorObj,
                          priorSample = priorSampleFn)
  chainSMC <- mcmc(target           = smcTarget,
                   sequenceType     = "sequential",
                   sequenceSchedule = "adaptiveEss",
                   populationSize   = 500L,
                   reSampling       = "systematic",
                   moveType         = "langevin",
                   sequenceControl  = smcControl(malaSteps = 5L, verbose = TRUE))
  chainSMC$logEvidence


  ## ----- (3) Hierarchical NLME: marginal (Pfad A) and joint (Pfad B) ----
  ## Same one-eta fixture as test-bayesNLME-marginal.R / test-bayesNLME-joint.R.
  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE, modelname = "mcmc_bnlme_obs")
  x <- Xt()
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2", "s3", "s4")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = "mcmc_bnlme_p")
  prd_nlme <- g * x * p

  true_mu  <- 2.0
  true_eta <- rnorm(length(subjects), 0, 0.3)
  y_nlme   <- true_mu * exp(true_eta) + rnorm(length(subjects), 0, 0.2)
  data_nlme <- as.datalist(data.frame(
    name = "y", time = 0, sigma = 0.2, value = y_nlme,
    condition = subjects, stringsAsFactors = FALSE))
  om  <- omega(eta = "eta", subjects = subjects)
  obj_nlme <- normL2(data_nlme, prd_nlme) + constraintL2(mu = 0, Omega = om)

  priorTheta    <- constraintL2(c(mu_pop = 2.0), sigma = 5.0)
  priorOmegaObj <- priorOmega(om, kind = "LKJHalfNormal", scaleSD = 1.0)
  priorSampleNLME <- function(n) {
    cbind(mu_pop        = rnorm(n, 2.0, 1.0),
          omega_eta_eta = rnorm(n, log(0.3), 0.3))
  }

  ## Pfad A: Laplace-marginalised etas, SMC over (theta_struct, omega_chol)
  tgtA <- bayesNLMEMarginal(obj = obj_nlme, omegaSpec = om,
                            prdfn = prd_nlme, data = data_nlme,
                            priorTheta    = priorTheta,
                            priorOmegaObj = priorOmegaObj,
                            priorSample   = priorSampleNLME)
  chainA <- mcmc(target           = tgtA,
                 sequenceType     = "sequential",
                 sequenceSchedule = "adaptiveEss",
                 populationSize   = 200L,
                 moveType         = "langevin",
                 sequenceControl  = smcControl(malaSteps = 5L))

  ## Pfad B: full joint Particle-Gibbs over (theta, omega, etas)
  tgtB <- bayesNLMEJoint(obj = obj_nlme, omegaSpec = om,
                         prdfn = prd_nlme, data = data_nlme,
                         priorTheta    = priorTheta,
                         priorOmegaObj = priorOmegaObj,
                         priorSample   = priorSampleNLME)
  chainB <- mcmc(target          = tgtB,
                 sequenceType    = "single",
                 moveType        = "langevin",
                 nIter           = 400L, warmup = 100L,
                 sequenceControl = smcControl(malaSteps = 5L))

  ## Cross-check against the FOCEI MAP
  mapNLME <- nlmeFit(obj_nlme, om, c(mu_pop = 2.0, omega_eta_eta = log(0.3)),
                     prdfn = prd_nlme, data = data_nlme, method = "focei")
  rbind(MAP   = mapNLME$argument,
        PfadA = colMeans(chainA$samples),
        PfadB = colMeans(chainB$samples))

  ## Pfad A logEvidence is available for Bayesian model comparison
  chainA$logEvidence

  ## Diagnostic + posterior plots
  plot(chainA);       plotTrace(chainA);  plotPairs(chainA)
  plot(chainB);       plotTrace(chainB)


  ## ----- (4) Multi-chain run with Gelman-Rubin R-hat --------------------
  chainsMulti <- mcmc(target           = tgtA,
                      sequenceType     = "sequential",
                      sequenceSchedule = "adaptiveEss",
                      populationSize   = 100L,
                      moveType         = "langevin",
                      sequenceControl  = smcControl(malaSteps = 3L),
                      chains           = 3L,
                      chainCores       = 1L)
  chainsMulti$rHat
  plot(chainsMulti)
}
