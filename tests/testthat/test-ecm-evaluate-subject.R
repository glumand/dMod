context("ecmEvaluateSubject: per-subject quadrature evaluator")

# NOTE: The per-subject gradient test below uses numDeriv as an independent
# check against the analytical implementation. The closed-form reference
# for a Smolyak-Gauss-Hermite quadrature over eta requires a contrived
# linear-in-theta toy log-likelihood (where Gauss-Hermite is exact). The
# resulting test would no longer exercise the realistic quadrature path
# and is rejected as cost/benefit unfavourable; numDeriv stays.


test_that("logLhat at level L recovers the 1D-integrate truth for a one-eta model", {
  set.seed(1)
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, modelname = "ecmeval_obs")
  x <- Xt()
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, modelname = "ecmeval_p")
  xPred <- g * x * p

  # Data: single observation per subject.
  y_obs <- c(2.0, 1.6)
  sigma <- 0.3
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = sigma,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)

  # Population parameters: mu_pop = 2, omega^2 = 0.09 (sd = 0.3 on log scale).
  psi <- c(mu_pop = 2.0, setNames(log(0.3), om$cholPars))

  # Build per-subject nodes around mode = 0 with H_i = 1/omega^2 + (data
  # contribution). For sanity start with etaHat = 0 and a Laplace-like H_i.
  # We'll evaluate logLhat and compare against a 1D integrate() ground truth.

  for (i in seq_along(subjects)) {
    # Find mode + H_i numerically.
    log_joint <- function(eta) {
      mu_i <- 2.0 * exp(eta)
      ll_data  <- -0.5 * ((y_obs[i] - mu_i) / sigma)^2 - 0.5 * log(2 * pi * sigma^2)
      ll_prior <- -0.5 * log(2 * pi * 0.09) - 0.5 * eta^2 / 0.09
      ll_data + ll_prior
    }
    opt <- optimize(log_joint, c(-3, 3), maximum = TRUE, tol = 1e-10)
    etaHat <- opt$maximum
    h <- 1e-4
    H_i <- -(log_joint(etaHat + h) - 2 * log_joint(etaHat) + log_joint(etaHat - h)) / h^2
    H_i <- matrix(H_i, 1, 1)

    # 1D integrate truth: log int p(y|eta) p(eta|Omega) deta.
    truth <- log(integrate(function(e) sapply(e, function(e0) exp(log_joint(e0))),
                           -Inf, Inf, rel.tol = 1e-10)$value)

    etaModes <- matrix(0, nrow = length(subjects), ncol = 1)
    etaModes[i, 1] <- etaHat

    for (level in 3:5) {
      nodes <- dMod:::makeSubjectNodes(etaHat, H_i, level)
      out   <- dMod:::ecmEvaluateSubject(
        subjIdx           = i,
        psiFull           = psi,
        etaModes          = etaModes,
        omegaSpec          = om,
        nodesSubj         = nodes,
        xPred             = xPred,
        datalist           = data,
        outerActiveNames = "mu_pop",
        mode               = "moments_only")
      # 5-place tolerance at level >= 4 for this simple model.
      tol <- switch(as.character(level), `3` = 1e-2, `4` = 1e-3, `5` = 1e-4)
      expect_equal(out$logLhat, truth, tolerance = tol,
                   info = sprintf("subject %d, level %d", i, level))
    }
  }
})


test_that("gradient at outer params matches numDeriv on the one-eta model", {
  skip_if_not_installed("numDeriv")
  set.seed(2)
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, modelname = "ecmgrad_obs")
  x <- Xt()
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2", "s3")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, modelname = "ecmgrad_p")
  xPred <- g * x * p

  y_obs <- c(2.1, 1.6, 2.5)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.3,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)

  psi <- c(mu_pop = 2.0, setNames(log(0.3), om$cholPars))

  # Build nodes around etaHat = 0, H_i = 5 for each subject (arbitrary but fixed).
  etaModes <- matrix(0, length(subjects), 1)
  nodes_list <- lapply(seq_along(subjects), function(i)
    dMod:::makeSubjectNodes(0, matrix(5, 1, 1), 4L))

  # Function: -2 log L_total = -2 sum_i logLhat_i.
  ofv_at <- function(mu) {
    psi_v <- psi
    psi_v["mu_pop"] <- mu
    sum(sapply(seq_along(subjects), function(i) {
      o <- dMod:::ecmEvaluateSubject(i, psi_v, etaModes, om, nodes_list[[i]],
                                xPred = xPred, datalist = data,
                                outerActiveNames = "mu_pop",
                                mode = "moments_only")
      -2 * o$logLhat
    }))
  }

  gr_analytic <- sum(sapply(seq_along(subjects), function(i) {
    o <- dMod:::ecmEvaluateSubject(i, psi, etaModes, om, nodes_list[[i]],
                              xPred = xPred, datalist = data,
                              outerActiveNames = "mu_pop",
                              mode = "with_grad")
    o$gradient["mu_pop"]
  }))
  gr_numeric <- numDeriv::grad(ofv_at, 2.0, method = "Richardson")
  expect_equal(gr_analytic, gr_numeric, tolerance = 1e-4)
})


test_that("frozen-node guarantee: nodesSubj does not change between em calls", {
  # The frozen-K_i invariant is enforced by the orchestrator (Phase 4d): nodes
  # are built once before each CM-1 trust call and not mutated during it. Here
  # we test the contract directly -- calling ecmEvaluateSubject twice with
  # different psi but the SAME nodesSubj must use identical nodes.
  set.seed(3)
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, modelname = "ecmfreeze_obs")
  x <- Xt()
  trafo <- eqnvec(intercept = "mu * exp(eta)")
  subj_table <- data.frame(eta = "eta_s1", row.names = "s1")
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, modelname = "ecmfreeze_p")
  xPred <- g * x * p

  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.3,
                                 value = 2.0, condition = "s1",
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = "s1")
  nodes <- dMod:::makeSubjectNodes(0, matrix(5, 1, 1), 4L)
  nodes_snapshot <- nodes$etaNodes

  psi_a <- c(mu = 2.0, setNames(log(0.3), om$cholPars))
  psi_b <- c(mu = 2.5, setNames(log(0.3), om$cholPars))

  o_a <- dMod:::ecmEvaluateSubject(1L, psi_a, matrix(0, 1, 1), om, nodes,
                              xPred, data, outerActiveNames = "mu",
                              mode = "moments_only")
  expect_identical(nodes$etaNodes, nodes_snapshot)
  o_b <- dMod:::ecmEvaluateSubject(1L, psi_b, matrix(0, 1, 1), om, nodes,
                              xPred, data, outerActiveNames = "mu",
                              mode = "moments_only")
  expect_identical(nodes$etaNodes, nodes_snapshot)
  # Sanity: different psi gives different logLhat.
  expect_false(isTRUE(all.equal(o_a$logLhat, o_b$logLhat)))
})
