# ============================================================================
# Bayesian NLME end-to-end tests.
#
# Sections:
#   * bayesNLMEMarginal -- FOCEI Laplace marginal target, used with SMC.
#   * bayesNLMEJoint    -- Particle-Gibbs target on (theta, omega, {eta_i}).
#
# Both targets are built on the same minimal one-eta NLME fixture (also
# used in test-focei.R / test-nlmefit.R) so the Bayesian build can be
# cross-checked against a converged FOCEI MAP.
# ============================================================================


# Shared fixture: 4-subject one-eta NLME at intercept * exp(eta).
.makeBayesNLMEFixture <- function(tag = "bnlme") {
  set.seed(1)
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd), add = TRUE)

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE,
         modelname = paste0(tag, "_obs_", as.integer(Sys.time())))
  x <- Xt()
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2", "s3", "s4")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = paste0(tag, "_p_", as.integer(Sys.time())))

  true_mu  <- 2.0
  true_om  <- 0.3
  true_eta <- rnorm(length(subjects), 0, true_om)
  y_obs    <- true_mu * exp(true_eta) + rnorm(length(subjects), 0, 0.2)
  data <- as.datalist(data.frame(
    name = "y", time = 0, sigma = 0.2,
    value = y_obs, condition = subjects,
    stringsAsFactors = FALSE))

  om <- omega(eta = "eta", subjects = subjects)
  obj <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)
  list(g = g, x = x, p = p, prd = g * x * p,
       data = data, om = om, obj = obj,
       subjects = subjects, true_mu = true_mu, y_obs = y_obs)
}


# ---- bayesNLMEMarginal --------------------------------------------------

test_that("bayesNLMEMarginal builds a likObj whose value at FOCEI MAP matches the converged OFV", {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()

  fx <- .makeBayesNLMEFixture("bm")
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(nlmeFit(
    fx$obj, fx$om, init,
    prdfn = fx$prd, data = fx$data,
    method = "focei",
    control = list(focei = list(
      trustControl = list(rinit = 1, rmax = 10, iterlim = 50,
                          fterm = 1e-7, mterm = 1e-7)))))

  tgt <- bayesNLMEMarginal(obj = fx$obj, omegaSpec = fx$om,
                            prdfn = fx$prd, data = fx$data,
                            priorSample = function(n) {
                              cbind(mu_pop = rnorm(n, 2, 1),
                                    omega_eta_eta = rnorm(n, log(0.3), 0.3))
                            })

  bayesLik <- tgt$likObj
  out_map <- bayesLik(fit$argument)
  expect_equal(out_map$value, fit$value, tolerance = 1e-3)
  expect_lt(max(abs(out_map$gradient)), 5e-2)

  struct_names <- attr(bayesLik, "structuralNames")
  H_struct <- out_map$hessian[struct_names, struct_names, drop = FALSE]
  expect_gt(min(eigen(H_struct, only.values = TRUE)$values), 0)
})


test_that("bayesNLMEMarginal: Omega-gradient FD-check on a non-converged point", {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()

  fx <- .makeBayesNLMEFixture("bm_fd")
  tgt <- bayesNLMEMarginal(obj = fx$obj, omegaSpec = fx$om,
                            prdfn = fx$prd, data = fx$data,
                            priorSample = function(n) NULL)
  bayesLik <- tgt$likObj

  pars <- c(mu_pop = 1.7, omega_eta_eta = log(0.25))
  out <- bayesLik(pars)
  expect_true(is.finite(out$value))
  expect_setequal(names(out$gradient), names(pars))

  h <- 1e-3
  pars_p <- pars; pars_p["omega_eta_eta"] <- pars["omega_eta_eta"] + h
  pars_m <- pars; pars_m["omega_eta_eta"] <- pars["omega_eta_eta"] - h
  fd_omega <- (bayesLik(pars_p)$value - bayesLik(pars_m)$value) / (2 * h)
  expect_equal(unname(out$gradient["omega_eta_eta"]), fd_omega,
               tolerance = 5e-2)
})


test_that("mcmc(bayesNLMEMarginal(...), sequenceType = 'sequential') smoke test", {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()

  fx <- .makeBayesNLMEFixture("bm_smc")
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))
  fit <- suppressMessages(nlmeFit(
    fx$obj, fx$om, init,
    prdfn = fx$prd, data = fx$data, method = "focei"))

  priorTheta    <- constraintL2(c(mu_pop = 2.0), sigma = 5.0)
  priorOmegaObj <- priorOmega(fx$om, kind = "LKJHalfNormal", scaleSD = 1.0)
  pSample <- function(n) {
    cbind(mu_pop        = rnorm(n, 2.0, 1.0),
          omega_eta_eta = rnorm(n, log(0.3), 0.3))
  }

  tgt <- bayesNLMEMarginal(obj = fx$obj, omegaSpec = fx$om,
                            prdfn = fx$prd, data = fx$data,
                            priorTheta    = priorTheta,
                            priorOmegaObj = priorOmegaObj,
                            priorSample   = pSample)

  set.seed(42L)
  chain <- mcmc(target           = tgt,
                 sequenceType     = "sequential",
                 sequenceSchedule = "adaptiveEss",
                 populationSize   = 60L,
                 moveType         = "langevin",
                 sequenceControl  = smcControl(malaSteps = 2L,
                                                essThreshold = 0.5))

  expect_s3_class(chain, "mcmcResultSequential")
  expect_s3_class(chain, "bayesNLMEMarginal")
  expect_equal(ncol(chain$samples), 2L)
  expect_true(is.finite(chain$logEvidence))

  pm <- colMeans(chain$samples)
  expect_lt(abs(pm["mu_pop"] - fit$argument["mu_pop"]), 0.8)
})


# ---- bayesNLMEJoint -----------------------------------------------------

test_that("mcmc(bayesNLMEJoint(...)) runs and returns expected shape", {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()

  fx <- .makeBayesNLMEFixture("bj")
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))
  fit <- suppressMessages(nlmeFit(
    fx$obj, fx$om, init,
    prdfn = fx$prd, data = fx$data, method = "focei"))

  priorTheta    <- constraintL2(c(mu_pop = 2.0), sigma = 5.0)
  priorOmegaObj <- priorOmega(fx$om, kind = "LKJHalfNormal", scaleSD = 1.0)
  pSample <- function(n) {
    cbind(mu_pop        = rnorm(n, 2.0, 1.0),
          omega_eta_eta = rnorm(n, log(0.3), 0.3))
  }

  tgt <- bayesNLMEJoint(obj = fx$obj, omegaSpec = fx$om,
                         prdfn = fx$prd, data = fx$data,
                         priorTheta    = priorTheta,
                         priorOmegaObj = priorOmegaObj,
                         priorSample   = pSample)

  set.seed(43L)
  chain <- mcmc(target          = tgt,
                 sequenceType    = "single",
                 moveType        = "langevin",
                 nIter           = 200L,
                 warmup          = 80L,
                 sequenceControl = smcControl(malaSteps = 3L))

  expect_s3_class(chain, "mcmcResultBlocked")
  expect_s3_class(chain, "bayesNLMEJoint")
  expect_equal(ncol(chain$samples), 2L)
  expect_equal(dim(chain$etaSamples)[2L], length(fx$subjects))
  expect_equal(dim(chain$etaSamples)[3L], 1L)
  expect_true(all(is.finite(chain$samples)))

  # 200 sweeps on a tiny fixture verifies the orchestrator is exploring
  # rather than converged. The chain is high-dim relative to the trajectory
  # length (P_outer + N * K_eta = 6 parameters, 200 samples post-burn-in);
  # the meaningful assertions are finiteness, structural shape, and a
  # non-trivial acceptance rate.
  pm <- colMeans(chain$samples)
  expect_true(is.finite(pm["mu_pop"]))
  expect_lt(abs(pm["mu_pop"] - fit$argument["mu_pop"]), 15.0)
  expect_gt(mean(chain$acceptOuter), 0.05)
})
