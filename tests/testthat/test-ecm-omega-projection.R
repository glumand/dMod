context("updateOmegaChol projection onto omegaSpec structure")


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
  # Off-diagonals must remain structurally zero.
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
  # Free the (eta_a, eta_c) correlation only; eta_b stays diagonal.
  om <- omega(eta = c("eta_a", "eta_b", "eta_c"), structure = "diag",
              correlate = list(c("eta_a", "eta_c")),
              subjects = paste0("s", 1:1000))
  # Omega_true with the structural-zero pattern.
  Omega_true <- diag(c(0.4, 0.2, 0.5)^2)
  Omega_true[1, 3] <- Omega_true[3, 1] <- 0.05
  expect_true(all(eigen(Omega_true)$values > 0))
  M_list <- draw_M_hat_list(1000, Omega_true)
  chol_vec <- updateOmegaChol(M_list, om)
  L <- om$buildL(chol_vec)
  Omega_est <- tcrossprod(L)
  # Tolerances are wider than full case because the constrained ML is slower.
  expect_equal(Omega_est[1, 1], Omega_true[1, 1], tolerance = 0.05)
  expect_equal(Omega_est[2, 2], Omega_true[2, 2], tolerance = 0.05)
  expect_equal(Omega_est[3, 3], Omega_true[3, 3], tolerance = 0.05)
  expect_equal(Omega_est[1, 3], Omega_true[1, 3], tolerance = 0.05)
  # Structurally-zero off-diag should remain zero.
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
