context("sparse_grid_gh + makeSubjectNodes")


# Analytic moments of physicists' GH: int z^k exp(-z^2) dz = 0 if k odd,
# (2m)!/(4^m m!) * sqrt(pi) if k = 2m.
gh_moment_1d <- function(k) {
  if (k %% 2L == 1L) return(0)
  m <- k / 2L
  factorial(2 * m) / (4^m * factorial(m)) * sqrt(pi)
}

# Smolyak-GH exactness criterion: A(L,K) is exact on monomial z_1^{a_1} ... z_K^{a_K}
# iff sum_j ceil((a_j + 3) / 4) <= L, since 1D level l has degree of exactness
# 2 * (2l - 1) - 1 = 4l - 3, requiring level l_j = ceil((a_j + 3)/4) per axis.
is_smolyak_exact <- function(a, level) {
  required <- sum(ceiling((a + 3L) / 4L))
  required <= level
}


test_that("sparse_grid_gh K=1 integrates polynomials up to degree 2L-1", {
  for (level in 2:5) {
    g <- sparse_grid_gh(1L, as.integer(level))
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


test_that("sparse_grid_gh K=2 is exact on the Smolyak-covered monomial set", {
  for (level in 2:5) {
    g <- sparse_grid_gh(2L, as.integer(level))
    z <- g$nodes
    w <- g$weights
    # Sweep all (a, b) within the per-axis exactness ceiling 4*level - 3.
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


test_that("sparse_grid_gh K=3 is exact on the Smolyak-covered monomial set", {
  for (level in 3:5) {
    g <- sparse_grid_gh(3L, as.integer(level))
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


test_that("sparse_grid_gh mass equals pi^(K/2) and first moment vanishes", {
  for (K in 1:4) for (level in K:5) {
    g <- sparse_grid_gh(as.integer(K), as.integer(level))
    expect_equal(sum(g$weights), pi^(K / 2), tolerance = 1e-12,
                 info = sprintf("mass K=%d level=%d", K, level))
    if (level >= K + 1L) {
      # second moments require one axis at degree 2 + others at degree 0
      for (k in seq_len(K)) {
        first <- sum(g$weights * g$nodes[, k])
        expect_lt(abs(first), 1e-12)
      }
    }
  }
})


test_that("sparse_grid_gh K=1 short-circuit yields the full m=2L-1 GH rule", {
  g <- sparse_grid_gh(1L, 4L)
  expect_equal(nrow(g$nodes), 7L)             # m(4) = 2*4 - 1 = 7
  expect_true(all(g$weights > 0))
  expect_equal(sum(g$weights), sqrt(pi), tolerance = 1e-12)
})


test_that("sparse_grid_gh errors on level < K", {
  expect_error(sparse_grid_gh(3L, 2L), "below dimension")
  expect_error(sparse_grid_gh(4L, 1L), "below dimension")
})


test_that("makeSubjectNodes round-trips first and second Gaussian moments", {
  # For Sigma = solve(H), augmented weights should integrate any g(eta) over
  # eta-space with no implicit weight: sum_b W_b * g(eta_b) approx int g(eta) deta.
  # Use g = N(eta | etaHat, Sigma) to get mass = 1, mean = etaHat, cov = Sigma.

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

    # mass
    expect_equal(sum(W * exp(logN)), 1, tolerance = 1e-8,
                 info = sprintf("level = %d, integral of N", level))
    # mean
    for (k in seq_len(K)) {
      est <- sum(W * qn$etaNodes[, k] * exp(logN))
      expect_equal(est, etaHat[k], tolerance = 1e-8,
                   info = sprintf("level = %d, first moment k=%d", level, k))
    }
    # covariance: needs the (2,0,0)/(0,2,0)/(0,0,2)/(1,1,0)/... set covered
    # , that's the same Smolyak-exactness requirement as level >= K + 1.
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
  # Smolyak weights with non-nested GH are mixed-sign; diagnostic check that
  # the grid isn't degenerate (positive weights heavily outweigh negative).
  qn  <- makeSubjectNodes(rep(0, 3L), diag(3L), 4L)
  pos <- sum(exp(qn$logAbsWeights[qn$weightSigns > 0]))
  neg <- sum(exp(qn$logAbsWeights[qn$weightSigns < 0]))
  expect_gt(pos, neg)
})
