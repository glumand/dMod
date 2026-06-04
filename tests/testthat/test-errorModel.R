# Behavioral tests for reduceReplicates() and fitErrorModel().
#
# reduceReplicates aggregates multi-replicate (name, time, condition) data
# down to a per-group mean and standard error. fitErrorModel fits a
# parametric variance model to the reduced data.
#
# Tests:
#   reduceReplicates: synthetic n-replicate samples with known mean and
#     standard deviation; the reduction should produce mean(value) and
#     standard error (= sd / sqrt(n)) up to a tolerance set by N(0,1)
#     sampling noise.
#   fitErrorModel: feed the reduced data to a constant-variance error
#     model and check that the recovered sigma equals the simulated
#     sigma_true within the n-replicates sampling envelope.


## ---- reduceReplicates: mean and standard error -----------------------

test_that("reduceReplicates returns the sample mean and standard error per group", {
  set.seed(42)
  sigma_true <- 0.1
  n_rep <- 5
  # 3 time points, single condition, 5 replicates each.
  raw <- do.call(rbind, lapply(c(1, 2, 5), function(tt) {
    data.frame(name = "A", time = tt,
               value = 1.0 + rnorm(n_rep, sd = sigma_true),
               condition = "C1",
               stringsAsFactors = FALSE)
  }))
  red <- reduceReplicates(raw)

  # Each (time, condition) group is a separate row.
  expect_equal(nrow(red), 3)
  expect_setequal(red$time, c(1, 2, 5))

  for (tt in c(1, 2, 5)) {
    row <- red[red$time == tt, ]
    grp <- raw[raw$time == tt, "value"]
    # `value` column is the mean.
    expect_equal(row$value, mean(grp), tolerance = 1e-12)
    # `sigma` is the standard error (sd / sqrt(n)).
    expect_equal(row$sigma, sd(grp) / sqrt(n_rep), tolerance = 1e-10)
    expect_equal(row$n, n_rep)
  }
})


## ---- fitErrorModel: recovers sigma_true on constant-variance data ----

test_that("fitErrorModel recovers exp(s0) ~ sigma_true^2 for constant variance", {
  testthat::skip_if_not_installed("optimx")
  set.seed(123)
  sigma_true <- 0.15
  n_rep <- 20
  # Many time points, many replicates, single condition: lots of evidence
  # to pin down the variance.
  raw <- do.call(rbind, lapply(seq(0.5, 5, by = 0.5), function(tt) {
    data.frame(name = "A", time = tt,
               value = 1.0 + rnorm(n_rep, sd = sigma_true),
               condition = "C1",
               stringsAsFactors = FALSE)
  }))
  red <- reduceReplicates(raw)

  # Constant variance error model: sigma^2 = exp(s0). Optimum: s0 ~ log(sigma_true^2).
  fit <- fitErrorModel(red, factors = "condition",
                       errorModel = "exp(s0)",
                       par = c(s0 = log(0.01)),
                       plotting = FALSE, blather = TRUE)
  s0_hat <- unique(fit$s0)[1]
  sigma_hat_sq <- exp(s0_hat)
  # Variance of sample variance with N-1 dof ~= 2 * sigma^4 / (N - 1).
  # Across 10 time points * 20 replicates we expect ~5% accuracy.
  expect_lt(abs(sqrt(sigma_hat_sq) - sigma_true) / sigma_true, 0.20)
})
