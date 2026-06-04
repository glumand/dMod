# ============================================================================
# Tests for the unified mcmc() sampler API.
#
# Sections:
#   * Linear-Gaussian fixtures (no ODE)        - exact closed-form posterior
#   * Langevin moves                            - single + ODE end-to-end
#   * Sequential / SMC                          - logZ + bimodality
#   * Multi-chain (chains = N)                   - R-hat shape check
#   * Plot helpers                              - smoke tests on fake outputs
#
# The Bayesian NLME target constructors (bayesNLMEMarginal / Joint) are
# tested in test-bayesNLME.R; here we exercise the sampler kernels and
# generic plot dispatch.
# ============================================================================


# ---- Linear-Gaussian fixtures -------------------------------------------

# Quadratic -2 log L on a Gaussian: closed-form gradient/Hessian, used to
# verify the sampler moves are well-calibrated.
.makeLinGaussObj <- function(mu, Sigma) {
  Prec <- solve(Sigma)
  K    <- length(mu)
  par_names <- names(mu)
  if (is.null(par_names)) par_names <- paste0("p", seq_len(K))
  names(mu) <- par_names
  fn <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                 conditions = NULL, env = NULL) {
    p <- pars[par_names]
    d <- p - mu
    val  <- as.numeric(t(d) %*% Prec %*% d)
    grad <- as.numeric(2 * Prec %*% d); names(grad) <- par_names
    H    <- 2 * Prec; dimnames(H) <- list(par_names, par_names)
    structure(list(value = val, gradient = grad, hessian = H),
              class = c("objlist", "list"))
  }
  class(fn) <- c("objfn", "fn")
  attr(fn, "parameters") <- par_names
  fn
}


# Linear-regression SMC fixture: conjugate Gaussian posterior with known
# log-evidence. Used to verify both posterior moments and logZ recovery.
.makeLinGaussRegression <- function(seed = 1L) {
  set.seed(seed)
  K  <- 2L
  n  <- 30L
  X  <- cbind(1, rnorm(n))
  colnames(X) <- c("a", "b")
  theta_true <- c(a = 0.5, b = -1.2)
  sigma_y    <- 0.7
  y <- as.numeric(X %*% theta_true) + rnorm(n, 0, sigma_y)

  mu0   <- c(a = 0, b = 0)
  s0    <- c(a = 2, b = 2)
  Sigma0 <- diag(s0^2); dimnames(Sigma0) <- list(names(mu0), names(mu0))

  Lambda <- t(X) %*% X / sigma_y^2 + diag(1 / s0^2)
  m_post <- solve(Lambda, t(X) %*% y / sigma_y^2 + mu0 / s0^2)
  Sigma_post <- solve(Lambda)
  dimnames(Sigma_post) <- list(names(mu0), names(mu0))
  m_post <- setNames(as.numeric(m_post), names(mu0))

  logZ <- -0.5 * (
    n * log(2 * pi * sigma_y^2)
    + sum(log(diag(Sigma0))) - sum(log(diag(Sigma_post)))
    + sum(y * y) / sigma_y^2
    + sum(mu0 * mu0 / s0^2)
    - drop(t(m_post) %*% Lambda %*% m_post))

  likObj <- local({
    par_names <- names(mu0)
    fn <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                   conditions = NULL, env = NULL, ...) {
      th <- pars[par_names]
      r  <- as.numeric(X %*% th) - y
      val <- sum(r * r) / sigma_y^2 + n * log(2 * pi * sigma_y^2)
      g <- 2 * as.numeric(t(X) %*% r) / sigma_y^2; names(g) <- par_names
      H <- 2 * (t(X) %*% X) / sigma_y^2
      dimnames(H) <- list(par_names, par_names)
      structure(list(value = val, gradient = g, hessian = H),
                class = c("objlist", "list"))
    }
    class(fn) <- c("objfn", "fn")
    attr(fn, "parameters") <- par_names
    fn
  })

  priorObj <- local({
    par_names <- names(mu0)
    fn <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                   conditions = NULL, env = NULL, ...) {
      th <- pars[par_names]
      r  <- th - mu0
      val <- sum(r * r / s0^2) + sum(log(2 * pi * s0^2))
      g <- 2 * r / s0^2; names(g) <- par_names
      H <- diag(2 / s0^2); dimnames(H) <- list(par_names, par_names)
      structure(list(value = val, gradient = g, hessian = H),
                class = c("objlist", "list"))
    }
    class(fn) <- c("objfn", "fn")
    attr(fn, "parameters") <- par_names
    fn
  })

  priorSample <- function(npart) {
    cbind(a = rnorm(npart, mu0["a"], s0["a"]),
          b = rnorm(npart, mu0["b"], s0["b"]))
  }

  list(likObj = likObj, priorObj = priorObj, priorSample = priorSample,
       mPost = m_post, SigmaPost = Sigma_post, logZ = logZ,
       thetaTrue = theta_true)
}


# ---- Langevin (moveType = "langevin") -----------------------------------

test_that("mcmc/langevin on linear-Gaussian toy recovers analytic posterior moments", {
  set.seed(101)
  K   <- 3L
  mu  <- c(a = 0.5, b = -0.3, c = 1.2)
  Sig <- matrix(c(1.0, 0.3, 0.0,
                  0.3, 0.8, -0.1,
                  0.0, -0.1, 1.2), 3, 3,
                dimnames = list(names(mu), names(mu)))
  obj <- .makeLinGaussObj(mu, Sig)

  chain <- mcmc(target = obj, sequenceType = "single",
                 moveType = "langevin", metric = "riemannFisher",
                 parinit  = mu + c(1, 1, 1),
                 nIter    = 4000L, warmup = 1500L)
  expect_s3_class(chain, "mcmcResultSingle")
  expect_s3_class(chain, "mcmcResult")

  emp_mean <- colMeans(chain$samples)
  emp_cov  <- cov(chain$samples)
  ess <- mean(chain$ess)
  mc_sd_mean <- sqrt(diag(Sig) / ess)
  for (j in seq_len(K))
    expect_lt(abs(emp_mean[j] - mu[j]), 4 * mc_sd_mean[j])
  for (j in seq_len(K))
    expect_lt(abs(emp_cov[j, j] / Sig[j, j] - 1), 0.25)
  expect_true(chain$acceptRate > 0.30 && chain$acceptRate < 0.85)
})


test_that("mcmc/langevin fixed-metric matches local on linear-Gaussian toy", {
  set.seed(202)
  mu  <- c(a = 0.0, b = 0.0)
  Sig <- matrix(c(1.0, 0.4, 0.4, 1.0), 2, 2, dimnames = list(names(mu), names(mu)))
  obj <- .makeLinGaussObj(mu, Sig)

  GFixed <- 2 * solve(Sig)  # -2 log L units; the kernel divides by 2 internally

  c_local <- mcmc(target = obj, sequenceType = "single",
                   moveType = "langevin", metric = "riemannFisher",
                   parinit = mu + 1, nIter = 3000L, warmup = 1000L)
  c_fixed <- mcmc(target = obj, sequenceType = "single",
                   moveType = "langevin", metric = "fixed",
                   moveControl = langevinControl(GFixed = GFixed),
                   parinit = mu + 1, nIter = 3000L, warmup = 1000L)

  m_local <- colMeans(c_local$samples)
  m_fixed <- colMeans(c_fixed$samples)
  expect_lt(max(abs(m_local - m_fixed)), 0.15)
})


test_that("mcmc/langevin is reproducible under set.seed()", {
  mu  <- c(a = 0, b = 0)
  Sig <- diag(c(1, 1)); dimnames(Sig) <- list(c("a", "b"), c("a", "b"))
  obj <- .makeLinGaussObj(mu, Sig)

  set.seed(7)
  c1 <- mcmc(target = obj, sequenceType = "single", moveType = "langevin",
              parinit = c(a = 0.1, b = -0.1), nIter = 200L, warmup = 100L)
  set.seed(7)
  c2 <- mcmc(target = obj, sequenceType = "single", moveType = "langevin",
              parinit = c(a = 0.1, b = -0.1), nIter = 200L, warmup = 100L)

  expect_equal(c1$samples, c2$samples, tolerance = 1e-12)
  expect_equal(c1$acceptRate, c2$acceptRate)
})


test_that("mcmc respects parlower / parupper", {
  mu  <- c(a = 0, b = 0)
  Sig <- diag(c(1, 1)); dimnames(Sig) <- list(c("a", "b"), c("a", "b"))
  obj <- .makeLinGaussObj(mu, Sig)

  lo <- c(a = -0.5, b = -0.5)
  hi <- c(a =  0.5, b =  0.5)

  set.seed(1)
  chain <- mcmc(target = obj, sequenceType = "single", moveType = "langevin",
                 parinit = c(a = 0, b = 0), nIter = 1000L, warmup = 200L,
                 parlower = lo, parupper = hi)

  expect_true(all(chain$samples[, "a"] >= lo["a"] - 1e-12))
  expect_true(all(chain$samples[, "a"] <= hi["a"] + 1e-12))
  expect_true(all(chain$samples[, "b"] >= lo["b"] - 1e-12))
  expect_true(all(chain$samples[, "b"] <= hi["b"] + 1e-12))
})


test_that("fixed-vs-local-Fisher agree on a constant-G linear-Gauss toy", {
  # For linear-Gauss, G is constant in theta. Fixed and local should produce
  # the same chain up to MC noise.
  mu  <- c(a = 0.0, b = 0.0, c = 0.0)
  Sig <- diag(rep(1, 3)); dimnames(Sig) <- list(names(mu), names(mu))
  obj <- .makeLinGaussObj(mu, Sig)

  set.seed(303)
  c_a <- mcmc(target = obj, sequenceType = "single", moveType = "langevin",
               metric = "fixed",
               moveControl = langevinControl(GFixed = 2 * solve(Sig)),
               parinit = mu, nIter = 4000L, warmup = 1000L)
  c_b <- mcmc(target = obj, sequenceType = "single", moveType = "langevin",
               metric = "riemannFisher",
               parinit = mu, nIter = 4000L, warmup = 1000L)

  expect_lt(max(abs(colMeans(c_a$samples) - colMeans(c_b$samples))), 0.1)
})


test_that("mcmc/langevin on ODE-based decay model runs end-to-end", {
  # Structural smoke test: the C++ chain runner must integrate cleanly with
  # a real ODE-backed objfn (normL2 -> CppODE solver -> back to objfn).
  # Quantitative posterior recovery on this fixture is fragile because dual
  # averaging can overshoot when the metric eigenvalues are O(100-1000);
  # quantitative checks are done on the linear-Gauss fixtures above.
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  obj_data  <- normL2(data, bench$prd_id)
  obj_prior <- constraintL2(c(A = 1.0, k = 0.5), sigma = 5)
  obj <- obj_data + obj_prior

  init <- c(A = 1.2, k = 0.4)
  mapfit <- trust(obj, init, rinit = 0.1, rmax = 10, iterlim = 50)
  expect_true(mapfit$converged)

  set.seed(11)
  chain <- mcmc(target = obj, sequenceType = "single", moveType = "langevin",
                 metric = "fixed",
                 moveControl = langevinControl(GFixed = mapfit$hessian),
                 parinit = mapfit$argument, nIter = 200L, warmup = 200L)
  expect_s3_class(chain, "mcmcResultSingle")
  expect_equal(dim(chain$samples), c(200L, 2L))
  expect_true(all(is.finite(chain$samples)))
  expect_true(is.finite(chain$acceptRate))
})


test_that("detailed-balance: linear-Gauss 2D marginal KS test", {
  set.seed(404)
  mu  <- c(a = 0.5, b = -0.3)
  Sig <- matrix(c(1.0, 0.2, 0.2, 1.5), 2, 2,
                dimnames = list(names(mu), names(mu)))
  obj <- .makeLinGaussObj(mu, Sig)

  chain <- mcmc(target = obj, sequenceType = "single", moveType = "langevin",
                 metric = "riemannFisher",
                 parinit = mu + c(1.5, -1), nIter = 12000L, warmup = 2000L)

  thinned <- chain$samples[seq(1L, nrow(chain$samples), by = 10L), ]
  ks <- suppressWarnings(stats::ks.test(
    thinned[, "a"], "pnorm",
    mean = mu[["a"]], sd = sqrt(Sig["a", "a"])))
  expect_gt(ks$p.value, 0.01)
})


# ---- SMC numerics --------------------------------------------------------

test_that("smcLogSumExp / smcESS / smcSystematicResample / mcmcSmcReweight: numerics", {
  set.seed(1L)
  x <- rnorm(100L)
  expect_equal(dMod:::smcLogSumExp(x),
               log(sum(exp(x - max(x)))) + max(x),
               tolerance = 1e-12)

  logw <- rnorm(50L)
  w    <- exp(logw - dMod:::smcLogSumExp(logw))
  expect_equal(dMod:::smcESS(logw), sum(w)^2 / sum(w * w), tolerance = 1e-10)

  w <- c(rep(0, 9L), 1.0)
  idx <- dMod:::smcSystematicResample(weights = w, u = 0.1)
  expect_true(all(idx == 10L))

  N <- 50L
  logL <- rnorm(N)
  logwPrev <- rep(-log(N), N)
  rw <- dMod:::mcmcSmcReweight(logL, logwPrev, betaOld = 0.0, betaNew = 1.0)
  ref <- (logwPrev + logL) - dMod:::smcLogSumExp(logwPrev + logL)
  expect_equal(rw$logw, ref, tolerance = 1e-12)
})


test_that("stratified / residual / multinomial resamplers concentrate on the one heavy particle", {
  set.seed(2L)
  w <- c(rep(0, 9L), 1.0)
  for (fn in list(dMod:::mcmcStratifiedResample,
                  dMod:::mcmcResidualResample,
                  dMod:::mcmcMultinomialResample)) {
    idx <- fn(weights = w)
    expect_length(idx, length(w))
    expect_true(all(idx == 10L))
  }
})


# ---- Sequential / SMC ---------------------------------------------------

test_that("mcmc/sequential recovers linear-Gaussian posterior moments and evidence", {
  set.seed(11L)
  prob <- .makeLinGaussRegression(seed = 17L)

  tgt <- flatTarget(likObj = prob$likObj, priorObj = prob$priorObj,
                    priorSample = prob$priorSample)
  chain <- mcmc(target           = tgt,
                 sequenceType     = "sequential",
                 sequenceSchedule = "adaptiveEss",
                 populationSize   = 800L,
                 moveType         = "langevin",
                 sequenceControl  = smcControl(malaSteps = 6L))

  expect_s3_class(chain, "mcmcResultSequential")
  expect_s3_class(chain, "mcmcResult")
  expect_equal(ncol(chain$samples), 2L)
  expect_equal(nrow(chain$samples), 800L)

  pm <- colMeans(chain$samples)
  expect_lt(max(abs(pm - prob$mPost)), 0.1)

  pv <- apply(chain$samples, 2L, var)
  for (j in seq_along(pv))
    expect_lt(abs(pv[j] / prob$SigmaPost[j, j] - 1), 0.40)

  expect_lt(abs(chain$logEvidence - prob$logZ), 0.4)
})


test_that("mcmc/sequential on a bimodal mixture finds both modes", {
  set.seed(42L)
  likObj <- local({
    par_names <- "x"
    fn <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                   conditions = NULL, env = NULL, ...) {
      x <- pars[par_names]
      lp <- log(0.5 * exp(-0.5 * (x + 3)^2) + 0.5 * exp(-0.5 * (x - 3)^2)) -
            0.5 * log(2 * pi)
      val <- -2 * lp
      l_neg <- exp(-0.5 * (x + 3)^2)
      l_pos <- exp(-0.5 * (x - 3)^2)
      w_neg <- l_neg / (l_neg + l_pos)
      w_pos <- 1 - w_neg
      dlp <- -w_neg * (x + 3) - w_pos * (x - 3)
      g   <- -2 * dlp; names(g) <- par_names
      H   <- matrix(2, 1, 1, dimnames = list(par_names, par_names))
      structure(list(value = as.numeric(val), gradient = g, hessian = H),
                class = c("objlist", "list"))
    }
    class(fn) <- c("objfn", "fn")
    attr(fn, "parameters") <- par_names
    fn
  })
  priorSample <- function(n) matrix(rnorm(n, 0, 10), n, 1L,
                                    dimnames = list(NULL, "x"))

  tgt <- flatTarget(likObj = likObj, priorSample = priorSample)
  chain <- mcmc(target           = tgt,
                 sequenceType     = "sequential",
                 sequenceSchedule = "adaptiveEss",
                 populationSize   = 1500L,
                 moveType         = "langevin",
                 sequenceControl  = smcControl(malaSteps = 5L))

  n_neg <- sum(chain$samples[, 1L] < 0)
  n_pos <- sum(chain$samples[, 1L] > 0)
  expect_gt(n_neg, 1500L * 0.25)
  expect_gt(n_pos, 1500L * 0.25)
})


test_that("mcmc/sequential is reproducible under set.seed()", {
  prob <- .makeLinGaussRegression(seed = 3L)
  tgt <- flatTarget(likObj = prob$likObj, priorObj = prob$priorObj,
                    priorSample = prob$priorSample)

  set.seed(99L)
  c1 <- mcmc(target = tgt, sequenceType = "sequential",
              populationSize = 200L, moveType = "langevin",
              sequenceControl = smcControl(malaSteps = 2L))
  set.seed(99L)
  c2 <- mcmc(target = tgt, sequenceType = "sequential",
              populationSize = 200L, moveType = "langevin",
              sequenceControl = smcControl(malaSteps = 2L))

  expect_equal(c1$samples,     c2$samples, tolerance = 1e-10)
  expect_equal(c1$logEvidence, c2$logEvidence)
  expect_equal(c1$betaPath,    c2$betaPath)
})


# ---- Multi-chain (chains = N) -------------------------------------------

test_that("mcmc(target, chains = N) returns expected shape and R-hat", {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()

  set.seed(1)
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)
  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE,
         modelname = paste0("bnlme_mc_obs_", as.integer(Sys.time())))
  x <- Xt()
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2", "s3", "s4")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = paste0("bnlme_mc_p_", as.integer(Sys.time())))
  true_mu  <- 2.0
  true_eta <- rnorm(length(subjects), 0, 0.3)
  y_obs    <- true_mu * exp(true_eta) + rnorm(length(subjects), 0, 0.2)
  data <- as.datalist(data.frame(
    name = "y", time = 0, sigma = 0.2,
    value = y_obs, condition = subjects,
    stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)
  obj <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)

  priorTheta    <- constraintL2(c(mu_pop = 2.0), sigma = 5.0)
  priorOmegaObj <- priorOmega(om, kind = "LKJHalfNormal", scaleSD = 1.0)
  pSample <- function(n) {
    cbind(mu_pop        = rnorm(n, 2.0, 1.0),
          omega_eta_eta = rnorm(n, log(0.3), 0.3))
  }

  tgt <- bayesNLMEMarginal(obj = obj, omegaSpec = om,
                            prdfn = g * x * p, data = data,
                            priorTheta    = priorTheta,
                            priorOmegaObj = priorOmegaObj,
                            priorSample   = pSample)

  chains <- mcmc(target           = tgt,
                  sequenceType     = "sequential",
                  sequenceSchedule = "adaptiveEss",
                  populationSize   = 40L,
                  moveType         = "langevin",
                  sequenceControl  = smcControl(malaSteps = 2L),
                  chains           = 3L,
                  chainCores       = 1L)

  expect_s3_class(chains, "mcmcResultMulti")
  expect_s3_class(chains, "mcmcResultSequential")
  expect_equal(chains$nChains, 3L)
  expect_equal(ncol(chains$samples), 2L)
  expect_equal(length(unique(chains$chainId)), 3L)
  expect_true(all(is.finite(chains$rHat)))
  expect_lt(max(chains$rHat), 3.0)
  expect_gt(min(chains$rHat), 0.5)
})


# ---- Plot helpers (smoke tests on fake fixtures) ------------------------

.makeFakeSmc <- function(n = 200L, K = 3L) {
  par_names <- paste0("p", seq_len(K))
  S <- matrix(rnorm(n * K), n, K, dimnames = list(NULL, par_names))
  out <- list(
    samples      = S,
    weights      = rep(1 / n, n),
    logEvidence  = -42.5,
    betaPath     = c(0, 0.1, 0.3, 0.6, 0.9, 1.0),
    ESSPath      = c(n, 100, 100, 100, 100, 100),
    acceptRates  = c(0.55, 0.50, 0.48, 0.46, 0.45),
    stepsizePath = c(0.1, 0.12, 0.15, 0.13, 0.11),
    nLevels      = 5L,
    call         = sys.call())
  class(out) <- c("mcmcResultSequential", "mcmcResult", "list")
  out
}


.makeFakeSingle <- function(n = 300L, K = 2L) {
  par_names <- paste0("q", seq_len(K))
  S <- matrix(rnorm(n * K), n, K, dimnames = list(NULL, par_names))
  out <- list(samples = S, logp = rnorm(n), accept = rep(TRUE, n),
              acceptRate = 0.6, stepsize = rep(0.1, n),
              ess = setNames(rep(80, K), par_names),
              parinit = setNames(numeric(K), par_names),
              moveType = "langevin", metric = "riemannFisher",
              sequenceType = "single",
              call = sys.call())
  class(out) <- c("mcmcResultSingle", "mcmcResult", "list")
  out
}


test_that("plot.mcmcResult renders without error (SMC fixture)", {
  smc <- .makeFakeSmc()
  p <- plot(smc)
  expect_s3_class(p, "ggplot")

  p2 <- plot(smc, parameters = c("p1", "p2"))
  expect_s3_class(p2, "ggplot")
})


test_that("plotTrace.mcmcResultSequential and plotTrace.mcmcResultSingle render", {
  smc <- .makeFakeSmc()
  expect_s3_class(plotTrace(smc), "ggplot")

  mc <- .makeFakeSingle()
  expect_s3_class(plotTrace(mc), "ggplot")
})


test_that("plotPairs renders without error", {
  smc <- .makeFakeSmc(K = 3L)
  expect_s3_class(plotPairs(smc), "ggplot")
})


test_that("plot.mcmcResult errors on unknown parameter name", {
  smc <- .makeFakeSmc()
  expect_error(plot(smc, parameters = "nonexistent"),
               "Unknown parameter name")
})


test_that("plot.mcmcResultMulti returns a ggplot", {
  S <- matrix(rnorm(60 * 2), 60, 2L, dimnames = list(NULL, c("a", "b")))
  chains <- list(
    samples     = S,
    chainId     = rep(1:3, each = 20L),
    rHat        = c(a = 1.01, b = 1.02),
    logEvidence = c(-12.1, -12.3, -12.0),
    seeds       = 1:3,
    nChains     = 3L,
    method      = "sequential",
    chains      = list(),
    call        = sys.call())
  class(chains) <- c("mcmcResultMulti", "mcmcResultSequential",
                      "mcmcResult", "list")
  expect_s3_class(plot(chains), "ggplot")
})


# ---- SMC init retry -----------------------------------------------------

# Flaky likelihood: first `fail_first` calls return Inf (logL = -Inf), then
# revert to a standard Gaussian objective.
make_flaky_lik <- function(fail_first = 5L) {
  counter <- 0L
  fn <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, ...) {
    counter <<- counter + 1L
    if (counter <= fail_first) {
      list(value = Inf, gradient = c(a = 0),
           hessian = matrix(1, 1, 1, dimnames = list("a", "a")))
    } else {
      v <- as.numeric(pars["a"])^2
      list(value = v, gradient = c(a = 2 * pars["a"]),
           hessian = matrix(2, 1, 1, dimnames = list("a", "a")))
    }
  }
  class(fn) <- c("objfn", "fn")
  attr(fn, "parameters") <- "a"
  fn
}


test_that("SMC init retries -Inf particles by redrawing from the prior", {
  set.seed(42)
  pop <- 6L
  likObj  <- make_flaky_lik(fail_first = 3L)
  pSample <- function(n) matrix(rnorm(n, 0, 1), n, 1L,
                                dimnames = list(NULL, "a"))
  tgt <- flatTarget(likObj = likObj, priorSample = pSample)

  chain <- mcmc(target           = tgt,
                sequenceType     = "sequential",
                sequenceSchedule = "fixed",
                populationSize   = pop,
                moveType         = "langevin",
                metric           = "euclidean",
                sequenceControl  = smcControl(schedule = c(0, 1),
                                              malaSteps = 1L,
                                              continuousAdaption = FALSE,
                                              verbose = FALSE),
                retry = TRUE, nTries = 10L)

  expect_s3_class(chain, "mcmcResultSequential")
  expect_equal(nrow(chain$samples), pop)
  expect_true(is.finite(chain$logEvidence))
})


test_that("SMC init with retry = FALSE flags every initial -Inf as failed", {
  set.seed(7)
  pop <- 4L
  likObj  <- make_flaky_lik(fail_first = 10000L)
  pSample <- function(n) matrix(rnorm(n, 0, 1), n, 1L,
                                dimnames = list(NULL, "a"))
  tgt <- flatTarget(likObj = likObj, priorSample = pSample)

  expect_error(
    mcmc(target           = tgt,
         sequenceType     = "sequential",
         sequenceSchedule = "fixed",
         populationSize   = pop,
         moveType         = "langevin",
         metric           = "euclidean",
         sequenceControl  = smcControl(schedule = c(0, 1),
                                       malaSteps = 1L,
                                       continuousAdaption = FALSE,
                                       verbose = FALSE),
         retry = FALSE),
    "initial likelihood non-finite"
  )
})
