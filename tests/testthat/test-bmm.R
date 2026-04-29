context("bmm (batched matrix multiplication, batch-first)")

test_that("bmm_lb matches per-slice reference", {
  set.seed(1)
  B <- 5; M <- 3; K <- 4; N <- 2
  A <- array(rnorm(B * M * K), dim = c(B, M, K))
  Bmat <- matrix(rnorm(K * N), K, N)
  got <- A %bmm% Bmat
  expect_equal(dim(got), c(B, M, N))
  ref <- array(0, c(B, M, N))
  for (b in seq_len(B)) ref[b,,] <- A[b,,] %*% Bmat
  expect_lt(max(abs(got - ref)), 1e-10)
})

test_that("bmm_rb matches per-slice reference", {
  set.seed(2)
  B <- 6; M <- 4; K <- 3; N <- 5
  Amat <- matrix(rnorm(M * K), M, K)
  Barr <- array(rnorm(B * K * N), dim = c(B, K, N))
  got <- Amat %bmm% Barr
  expect_equal(dim(got), c(B, M, N))
  ref <- array(0, c(B, M, N))
  for (b in seq_len(B)) ref[b,,] <- Amat %*% Barr[b,,]
  expect_lt(max(abs(got - ref)), 1e-10)
})

test_that("bmm_bb matches per-slice reference", {
  set.seed(3)
  B <- 7; M <- 3; K <- 2; N <- 4
  A <- array(rnorm(B * M * K), dim = c(B, M, K))
  Barr <- array(rnorm(B * K * N), dim = c(B, K, N))
  got <- A %bmm% Barr
  expect_equal(dim(got), c(B, M, N))
  ref <- array(0, c(B, M, N))
  for (b in seq_len(B)) ref[b,,] <- A[b,,] %*% Barr[b,,]
  expect_lt(max(abs(got - ref)), 1e-10)
})

test_that("Bn = 1 edge case", {
  A <- array(rnorm(1 * 3 * 4), dim = c(1, 3, 4))
  X <- matrix(rnorm(4 * 2), 4, 2)
  got <- A %bmm% X
  expect_identical(dim(got), c(1L, 3L, 2L))
  expect_lt(max(abs(got[1,,] - A[1,,] %*% X)), 1e-10)
})

test_that("dimension-mismatch errors are raised", {
  A <- array(0, c(3, 2, 4))
  Bmat <- matrix(0, 3, 2)  # K mismatch: should be 4
  expect_error(A %bmm% Bmat)
})
