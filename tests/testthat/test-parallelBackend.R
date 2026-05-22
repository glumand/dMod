# Smoke tests for the cross-platform parallel-apply helper. Mirrors the
# pattern used by mstrust / profile() (foreach + doParallel under the
# hood) so we can be confident the Windows path works without an actual
# Windows runner.


test_that(".parallelLapply with cores=1 falls back to lapply", {
  out <- dMod:::.parallelLapply(1:5, function(i) i * 2, cores = 1L)
  expect_identical(out, lapply(1:5, function(i) i * 2))
})


test_that(".parallelLapply with cores>1 returns correct results", {
  # Either Unix fork via doParallel or PSOCK via makeCluster -- both go
  # through the same %dopar% loop and must produce identical results to
  # the serial path.
  out_par <- dMod:::.parallelLapply(1:8, function(i) i^2, cores = 2L)
  out_ser <- lapply(1:8, function(i) i^2)
  expect_identical(out_par, out_ser)
})


test_that(".parallelLapply preserves order across workers", {
  # R CMD check caps cores at 2 via _R_CHECK_LIMIT_CORES_, so honour that
  # ceiling here. 2 workers across 20 items still exercises the order-
  # preservation property of the foreach %dopar% backend.
  n_cores <- if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) 2L else 4L
  set.seed(7L)
  X <- as.list(rnorm(20))
  out <- dMod:::.parallelLapply(X, function(x) x + 1, cores = n_cores)
  expect_identical(out, lapply(X, function(x) x + 1))
})


test_that(".parallelLapply propagates worker errors", {
  # foreach surfaces worker errors as a regular R error -- the inner
  # stop("kaboom") should reach the caller via tryCatch.
  err <- tryCatch(dMod:::.parallelLapply(1:3, function(i) {
    if (i == 2L) stop("kaboom") else i
  }, cores = 2L), error = function(e) e)
  expect_s3_class(err, "error")
})
