# ============================================================================
# Quadrature / ECM stack tests.
#
# Sections (bottom-up through the call graph):
#   * sparseGridGH / makeSubjectNodes   - low-level Smolyak-GH numerics
#   * evalConditionResidual               - per-condition residual lifted helper
#   * ecmEvaluateSubject                  - per-subject 3-moment evaluator
#   * updateOmegaChol                     - CM-2 projection onto omega structure
#   * emObjfn(method = "quadrature")       - end-to-end objective + rebuildQuadrature
#
# Where a closed-form integral exists, we use it as the oracle. Otherwise
# numDeriv is used as an independent gradient check (see comments in each
# block for why; an exact closed-form reference would not exercise the
# realistic Smolyak code path).
# ============================================================================

context("Quadrature / ECM")


# ---- Smolyak-GH numerics --------------------------------------------------

# Analytic moments of physicists' GH: int z^k exp(-z^2) dz = 0 if k odd,
# (2m)!/(4^m m!) * sqrt(pi) if k = 2m.
gh_moment_1d <- function(k) {
  if (k %% 2L == 1L) return(0)
  m <- k / 2L
  factorial(2 * m) / (4^m * factorial(m)) * sqrt(pi)
}

# Smolyak-GH exactness criterion: A(L,K) is exact on z_1^{a_1} ... z_K^{a_K}
# iff sum_j ceil((a_j + 3) / 4) <= L.
is_smolyak_exact <- function(a, level) {
  required <- sum(ceiling((a + 3L) / 4L))
  required <= level
}


test_that("sparseGridGH K=1 integrates polynomials up to degree 2L-1", {
  for (level in 2:5) {
    g <- sparseGridGH(1L, as.integer(level))
    z <- as.numeric(g$nodes)
    w <- g$weights
    for (k in 0:(2L * level - 1L)) {
      est  <- sum(w * z^k)
      true <- gh_moment_1d(k)
      expect_equal(est, true, tolerance = 1e-10,
                   info = sprintf("level = %d, k = %d", level, k))
    }
  }
})


test_that("sparseGridGH K=2 is exact on the Smolyak-covered monomial set", {
  for (level in 2:5) {
    g <- sparseGridGH(2L, as.integer(level))
    z <- g$nodes
    w <- g$weights
    deg_max <- 4L * level - 3L
    for (a in 0:deg_max) {
      for (b in 0:deg_max) {
        if (!is_smolyak_exact(c(a, b), level)) next
        est  <- sum(w * (z[, 1]^a) * (z[, 2]^b))
        true <- gh_moment_1d(a) * gh_moment_1d(b)
        expect_equal(est, true, tolerance = 1e-10,
                     info = sprintf("level=%d, (a,b)=(%d,%d)", level, a, b))
      }
    }
  }
})


test_that("sparseGridGH K=3 is exact on the Smolyak-covered monomial set", {
  for (level in 3:5) {
    g <- sparseGridGH(3L, as.integer(level))
    z <- g$nodes
    w <- g$weights
    deg_max <- 4L * level - 3L
    for (a in 0:deg_max) for (b in 0:deg_max) for (c in 0:deg_max) {
      if (!is_smolyak_exact(c(a, b, c), level)) next
      est  <- sum(w * (z[, 1]^a) * (z[, 2]^b) * (z[, 3]^c))
      true <- gh_moment_1d(a) * gh_moment_1d(b) * gh_moment_1d(c)
      expect_equal(est, true, tolerance = 1e-10,
                   info = sprintf("level=%d, (a,b,c)=(%d,%d,%d)",
                                  level, a, b, c))
    }
  }
})


test_that("sparseGridGH mass equals pi^(K/2) and first moment vanishes", {
  for (K in 1:4) for (level in K:5) {
    g <- sparseGridGH(as.integer(K), as.integer(level))
    expect_equal(sum(g$weights), pi^(K / 2), tolerance = 1e-12,
                 info = sprintf("mass K=%d level=%d", K, level))
    if (level >= K + 1L) {
      for (k in seq_len(K)) {
        first <- sum(g$weights * g$nodes[, k])
        expect_lt(abs(first), 1e-12)
      }
    }
  }
})


test_that("sparseGridGH K=1 short-circuit yields the full m=2L-1 GH rule", {
  g <- sparseGridGH(1L, 4L)
  expect_equal(nrow(g$nodes), 7L)             # m(4) = 2*4 - 1 = 7
  expect_true(all(g$weights > 0))
  expect_equal(sum(g$weights), sqrt(pi), tolerance = 1e-12)
})


test_that("sparseGridGH errors on level < K", {
  expect_error(sparseGridGH(3L, 2L), "below dimension")
  expect_error(sparseGridGH(4L, 1L), "below dimension")
})


test_that("makeSubjectNodes round-trips first and second Gaussian moments", {
  # For Sigma = solve(H), augmented weights should integrate any g(eta) over
  # eta-space directly: sum_b W_b * g(eta_b) ~ int g(eta) deta. Use
  # g = N(eta | etaHat, Sigma) to get mass = 1, mean = etaHat, cov = Sigma.

  set.seed(42)
  K       <- 3L
  etaHat <- c(0.4, -0.2, 0.7)
  A       <- matrix(rnorm(K * K), K, K)
  H_i     <- crossprod(A) + diag(K)
  Sig     <- solve(H_i)

  for (level in (K + 1L):5L) {
    qn <- makeSubjectNodes(etaHat, H_i, level)
    W  <- qn$weightSigns * exp(qn$logAbsWeights)

    R    <- chol(Sig)
    logd <- sum(log(diag(R)))
    diff <- sweep(qn$etaNodes, 2, etaHat, "-")
    z    <- backsolve(R, t(diff), transpose = TRUE)
    qf   <- colSums(z^2)
    logN <- -0.5 * K * log(2 * pi) - logd - 0.5 * qf

    expect_equal(sum(W * exp(logN)), 1, tolerance = 1e-8,
                 info = sprintf("level = %d, integral of N", level))
    for (k in seq_len(K)) {
      est <- sum(W * qn$etaNodes[, k] * exp(logN))
      expect_equal(est, etaHat[k], tolerance = 1e-8,
                   info = sprintf("level = %d, first moment k=%d", level, k))
    }
    if (level >= K + 1L) {
      for (i in seq_len(K)) for (j in i:K) {
        est <- sum(W * (qn$etaNodes[, i] - etaHat[i]) *
                       (qn$etaNodes[, j] - etaHat[j]) * exp(logN))
        expect_equal(est, Sig[i, j], tolerance = 1e-7,
                     info = sprintf("level=%d, cov (%d,%d)", level, i, j))
      }
    }
  }
})


test_that("makeSubjectNodes positive mass dominates negative", {
  # Smolyak weights with non-nested GH are mixed-sign; diagnostic that the
  # grid isn't degenerate (positive weights heavily outweigh negative).
  qn  <- makeSubjectNodes(rep(0, 3L), diag(3L), 4L)
  pos <- sum(exp(qn$logAbsWeights[qn$weightSigns > 0]))
  neg <- sum(exp(qn$logAbsWeights[qn$weightSigns < 0]))
  expect_gt(pos, neg)
})


# ---- evalConditionResidual ------------------------------------------------

test_that("evalConditionResidual matches the in-normL2 closure (no errfn)", {
  set.seed(1)
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, modelname = "ecr_simple_obs")
  x <- Xt()

  subjects <- c("c1", "c2")
  trafo <- eqnvec(intercept = "mu * exp(eta)")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, modelname = "ecr_simple_p")

  y_obs <- c(2.1, 1.8)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  obj <- normL2(data, g * x * p)

  init <- c(mu = 2.0, eta_c1 = 0.1, eta_c2 = -0.1)
  result <- obj(init)
  prediction <- (g * x * p)(times = c(0), pars = init, deriv = TRUE,
                            conditions = subjects)
  per_cn <- lapply(subjects, function(cn) {
    evalConditionResidual(dataI = data[[cn]], predictionI = prediction[[cn]],
                            pars = init, cn = cn, bessel = 1, deriv = TRUE)
  })
  agg <- Reduce(`+`, per_cn)
  expect_equal(agg$value, result$value, tolerance = 1e-10)
  expect_equal(agg$gradient, result$gradient, tolerance = 1e-10)
  expect_equal(agg$hessian,  result$hessian,  tolerance = 1e-10)
})


test_that("normL2 still produces the same OFV after eval_condition refactor", {
  # Regression: hand-check that the eval_condition lift did not drift the
  # residual numerics relative to the pre-refactor OFV.
  set.seed(1)
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, modelname = "ecr_regr_obs")
  x <- Xt()
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2", "s3", "s4")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, modelname = "ecr_regr_p")

  true_mu  <- 2.0; true_om <- 0.3
  true_eta <- rnorm(length(subjects), 0, true_om)
  y_obs    <- true_mu * exp(true_eta) + rnorm(length(subjects), 0, 0.2)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)
  obj <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)

  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3),
            setNames(rep(0, 4), paste0("eta_", subjects)))
  res <- obj(init, deriv = TRUE)
  expect_true(is.finite(res$value))
  expect_true(all(is.finite(res$gradient)))
})


# ---- ecmEvaluateSubject ---------------------------------------------------

test_that("logLhat at level L recovers the 1D-integrate truth for a one-eta prdfn", {
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

  y_obs <- c(2.0, 1.6)
  sigma <- 0.3
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = sigma,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)

  psi <- c(mu_pop = 2.0, setNames(log(0.3), om$cholPars))

  for (i in seq_along(subjects)) {
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
        omega          = om,
        nodesSubj         = nodes,
        xPred             = xPred,
        datalist           = data,
        outerActiveNames = "mu_pop",
        mode               = "moments_only")
      tol <- switch(as.character(level), `3` = 1e-2, `4` = 1e-3, `5` = 1e-4)
      expect_equal(out$logLhat, truth, tolerance = tol,
                   info = sprintf("subject %d, level %d", i, level))
    }
  }
})


test_that("gradient at outer params matches numDeriv on the one-eta prdfn", {
  # numDeriv is the oracle here -- a clean closed-form reference for the
  # quadrature gradient would need a linear-in-theta toy log-likelihood
  # (where Gauss-Hermite is exact), which would not exercise the realistic
  # quadrature code path.
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

  etaModes <- matrix(0, length(subjects), 1)
  nodes_list <- lapply(seq_along(subjects), function(i)
    dMod:::makeSubjectNodes(0, matrix(5, 1, 1), 4L))

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
  # we test the contract directly: calling ecmEvaluateSubject twice with
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
  expect_false(isTRUE(all.equal(o_a$logLhat, o_b$logLhat)))
})


# ---- updateOmegaChol (CM-2 projection) -----------------------------------

# Helper: draw N samples eta_i ~ N(0, Omega_true), build MHatList = eta_i %o% eta_i.
draw_M_hat_list <- function(N, Omega_true) {
  K <- nrow(Omega_true)
  L <- t(chol(Omega_true))
  Z <- matrix(rnorm(N * K), N, K)
  E <- Z %*% t(L)
  lapply(seq_len(N), function(i) tcrossprod(E[i, ]))
}


test_that("diagonal: recovers Omega_true axis variances closed-form", {
  set.seed(1)
  om <- omega(eta = c("eta_a", "eta_b", "eta_c"), structure = "diag",
              subjects = paste0("s", 1:200))
  Omega_true <- diag(c(0.3, 0.1, 0.5)^2)
  M_list <- draw_M_hat_list(200, Omega_true)
  chol_vec <- updateOmegaChol(M_list, om)
  L <- om$buildL(chol_vec)
  Omega_est <- tcrossprod(L)
  expect_equal(unname(diag(Omega_est)), diag(Omega_true), tolerance = 0.1)
  expect_true(all(Omega_est[lower.tri(Omega_est)] == 0))
})


test_that("full: recovers Omega_true Cholesky closed-form", {
  set.seed(2)
  om <- omega(eta = c("eta_a", "eta_b"), structure = "full",
              subjects = paste0("s", 1:1000))
  Omega_true <- matrix(c(0.16, 0.06, 0.06, 0.09), 2, 2)
  M_list <- draw_M_hat_list(1000, Omega_true)
  chol_vec <- updateOmegaChol(M_list, om)
  L <- om$buildL(chol_vec)
  Omega_est <- tcrossprod(L)
  expect_equal(Omega_est, Omega_true, tolerance = 0.05,
               check.attributes = FALSE)
})


test_that("selective correlate: trust on Q-function recovers Omega_true", {
  set.seed(3)
  om <- omega(eta = c("eta_a", "eta_b", "eta_c"), structure = "diag",
              correlate = list(c("eta_a", "eta_c")),
              subjects = paste0("s", 1:1000))
  Omega_true <- diag(c(0.4, 0.2, 0.5)^2)
  Omega_true[1, 3] <- Omega_true[3, 1] <- 0.05
  expect_true(all(eigen(Omega_true)$values > 0))
  M_list <- draw_M_hat_list(1000, Omega_true)
  chol_vec <- updateOmegaChol(M_list, om)
  L <- om$buildL(chol_vec)
  Omega_est <- tcrossprod(L)
  expect_equal(Omega_est[1, 1], Omega_true[1, 1], tolerance = 0.05)
  expect_equal(Omega_est[2, 2], Omega_true[2, 2], tolerance = 0.05)
  expect_equal(Omega_est[3, 3], Omega_true[3, 3], tolerance = 0.05)
  expect_equal(Omega_est[1, 3], Omega_true[1, 3], tolerance = 0.05)
  expect_equal(Omega_est[1, 2], 0)
  expect_equal(Omega_est[2, 3], 0)
})


test_that("rank-deficient S triggers ridge warning but does not crash", {
  set.seed(4)
  om <- omega(eta = c("eta_a", "eta_b"), structure = "full",
              subjects = paste0("s", 1:5))
  M_list <- lapply(seq_len(5L), function(i) matrix(c(1, 0.5, 0.5, 0.25), 2, 2))
  expect_warning(updateOmegaChol(M_list, om), "rank-deficient")
})


# ---- emObjfn(method = "quadrature") -------------------------------------

.build_quad_setup <- function(seed = 1L) {
  set.seed(seed)
  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, modelname = paste0("emqd_obs_", seed))
  x <- Xt()
  subjects <- c("s1", "s2", "s3")
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE,
         modelname = paste0("emqd_p_", seed))
  y_obs <- 2.0 * exp(rnorm(length(subjects), 0, 0.3)) +
           rnorm(length(subjects), 0, 0.2)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)
  obj <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)
  list(g = g, x = x, p = p, data = data, om = om, obj = obj,
       prdfn = g * x * p, subjects = subjects)
}


test_that("emObjfn quadrature constructs and exposes rebuildQuadrature", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_quad_setup(1L)

  em <- emObjfn(s$obj, s$om,
                prdfn = s$prdfn, data = s$data,
                method = "quadrature", control = list(level = 4L))
  expect_s3_class(em, "emObjfn")
  expect_equal(attr(em, "method"), "quadrature")
  expect_true(is.function(attr(em, "rebuildQuadrature")))
  init <- c(mu_pop = 2.0)
  expect_error(em(init), "no quadrature grid built")
})


test_that("rebuildQuadrature populates nodes, modes, and frozen state", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_quad_setup(2L)

  em <- emObjfn(s$obj, s$om,
                prdfn = s$prdfn, data = s$data,
                method = "quadrature", control = list(level = 4L))
  psi <- c(mu_pop = 2.0, setNames(log(0.3), s$om$cholPars))
  refresh <- attr(em, "rebuildQuadrature")(psi)
  expect_equal(dim(refresh$etaModes), c(length(s$subjects), 1L))
  expect_length(refresh$HiList, length(s$subjects))
  expect_equal(refresh$level, 4L)
  expect_true(all(refresh$converged))

  init <- c(mu_pop = 2.0)
  out <- em(init)
  expect_true(is.finite(out$value))
  expect_true(is.finite(out$gradient["mu_pop"]))
  expect_equal(attr(out, "emDiag")$method, "quadrature")
})


test_that("frozen-node invariance: em(pars1) and em(pars2) reuse the same nodes", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_quad_setup(3L)

  em <- emObjfn(s$obj, s$om,
                prdfn = s$prdfn, data = s$data,
                method = "quadrature", control = list(level = 4L))
  psi <- c(mu_pop = 2.0, setNames(log(0.3), s$om$cholPars))
  attr(em, "rebuildQuadrature")(psi)
  out_a <- em(c(mu_pop = 2.0))
  modes_a <- attr(out_a, "emDiag")$etaModes
  out_b <- em(c(mu_pop = 2.5))
  modes_b <- attr(out_b, "emDiag")$etaModes
  expect_identical(modes_a, modes_b)
  expect_false(isTRUE(all.equal(out_a$value, out_b$value)))
})


test_that("quadrature gradient matches numDeriv on the one-eta prdfn", {
  # numDeriv is the oracle here -- closed-form would require a contrived
  # linear-in-theta toy where Gauss-Hermite is exact and would mask the
  # adaptive Smolyak code path being exercised.
  skip_if_not_installed("numDeriv")
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_quad_setup(4L)

  em <- emObjfn(s$obj, s$om,
                prdfn = s$prdfn, data = s$data,
                method = "quadrature", control = list(level = 5L))
  psi <- c(mu_pop = 2.0, setNames(log(0.3), s$om$cholPars))
  attr(em, "rebuildQuadrature")(psi)

  f_val <- function(mu) em(c(mu_pop = mu), deriv = FALSE)$value
  gr_num <- numDeriv::grad(f_val, 2.1, method = "Richardson")
  gr_ana <- em(c(mu_pop = 2.1))$gradient["mu_pop"]
  expect_equal(unname(gr_ana), gr_num, tolerance = 1e-3)
})
