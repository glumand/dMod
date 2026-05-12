# Behavioral tests for res() (data <-> prediction residual operator).
#
# Verifies:
#   * weighted.residual = (prediction - value) / sigma  per row
#   * weighted.0       = prediction / sigma
#   * BLOQ handling: value = pmax(value_raw, lloq); bloq flag matches
#   * errmodel sigma fills rows where data$sigma is NA
#   * deriv attribute layout: [n_rows, n_params]
#
# These tests build a prediction matrix by hand to keep the math
# transparent. res() is otherwise exercised end-to-end by every
# normL2() / nll() test.


# Build a small prediction matrix for A(t) = exp(-0.5 * t) at t in {0,1,2}.
# The "deriv" attribute carries dA/dA_par and dA/dk_par populated from the
# linear-decay closed form.
.make_prdframe <- function() {
  prdf <- matrix(c(0, 1, 2,
                   1.0, exp(-0.5), exp(-1.0)),
                 nrow = 3, ncol = 2,
                 dimnames = list(NULL, c("time", "A")))
  d_arr <- array(0, c(3, 1, 2),
                 dimnames = list(NULL, "A", c("A_par", "k_par")))
  d_arr[, "A", "A_par"] <- c(1.0, exp(-0.5), exp(-1.0))
  d_arr[, "A", "k_par"] <- c(0, -1 * 1.0 * exp(-0.5), -2 * 1.0 * exp(-1.0))
  attr(prdf, "deriv") <- d_arr
  prdf
}


## ---- Value: weighted residual --------------------------------------------

test_that("res computes weighted.residual = (pred - value) / sigma per row", {
  prdf <- .make_prdframe()
  dat <- data.frame(time = c(1, 2), name = "A",
                    value = c(0.55, 0.30), sigma = c(0.1, 0.2),
                    lloq = c(-Inf, -Inf),
                    stringsAsFactors = FALSE)
  r <- res(dat, prdf)
  pred <- c(exp(-0.5), exp(-1.0))
  expect_equal(r$prediction, pred, tolerance = 1e-12)
  expect_equal(r$weighted.residual, (pred - dat$value) / dat$sigma,
               tolerance = 1e-12)
  expect_equal(r$weighted.0, pred / dat$sigma, tolerance = 1e-12)
  expect_true(all(r$bloq == FALSE))
})


## ---- BLOQ handling -------------------------------------------------------

test_that("res sets value = pmax(value_raw, lloq) and bloq mask correctly", {
  prdf <- .make_prdframe()
  dat <- data.frame(time = c(1, 2), name = "A",
                    value = c(0.55, 0.05),  # 0.05 < lloq -> censored
                    sigma = c(0.1, 0.1),
                    lloq  = c(-Inf, 0.10),
                    stringsAsFactors = FALSE)
  r <- res(dat, prdf)
  # Row 1: above LLOQ -> value unchanged; row 2: below -> bumped to LLOQ.
  expect_equal(r$value, c(0.55, 0.10), tolerance = 1e-12)
  expect_equal(r$bloq, c(FALSE, TRUE))
  # weighted.residual uses the LLOQ-substituted value.
  expect_equal(r$weighted.residual[2], (exp(-1) - 0.10) / 0.1,
               tolerance = 1e-12)
})


## ---- Errmodel-supplied sigma --------------------------------------------

test_that("res fills sigma from errmodel matrix where data$sigma is NA", {
  prdf <- .make_prdframe()
  err <- matrix(c(0, 1, 2,
                  0.2, 0.2, 0.2),
                nrow = 3, ncol = 2,
                dimnames = list(NULL, c("time", "A")))

  dat <- data.frame(time = c(1, 2), name = "A",
                    value = c(0.55, 0.30),
                    sigma = c(NA, 0.1),
                    lloq  = c(-Inf, -Inf),
                    stringsAsFactors = FALSE)
  r <- res(dat, prdf, err = err)
  expect_equal(r$sigma, c(0.2, 0.1), tolerance = 1e-12)
  pred <- c(exp(-0.5), exp(-1.0))
  expect_equal(r$weighted.residual,
               (pred - dat$value) / r$sigma, tolerance = 1e-12)
})


## ---- Deriv attribute layout --------------------------------------------

test_that("res 'deriv' attribute has shape [n_rows, n_params] with expected names", {
  prdf <- .make_prdframe()
  dat <- data.frame(time = c(1, 2), name = "A",
                    value = c(0.55, 0.30), sigma = c(0.1, 0.1),
                    lloq = c(-Inf, -Inf),
                    stringsAsFactors = FALSE)
  r <- res(dat, prdf)
  d <- attr(r, "deriv")
  expect_equal(dim(d), c(2, 2))
  expect_equal(colnames(d), c("A_par", "k_par"))
  expect_equal(d[, "A_par"], c(exp(-0.5), exp(-1.0)), tolerance = 1e-12)
  expect_equal(d[, "k_par"], c(-1 * exp(-0.5), -2 * exp(-1.0)),
               tolerance = 1e-12)
})


## ---- Observable mismatch error ----------------------------------------

test_that("res errors if an observable in data is missing from the prediction", {
  prdf <- .make_prdframe()
  dat <- data.frame(time = 1, name = "B",  # B does not exist in prdf
                    value = 0.5, sigma = 0.1, lloq = -Inf,
                    stringsAsFactors = FALSE)
  expect_error(res(dat, prdf), regexp = "Observable not found")
})
