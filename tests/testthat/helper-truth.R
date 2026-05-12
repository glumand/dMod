# Analytical ground-truth fixtures for behavioral tests.
#
# Pure-R closed-form solutions, no compilation. Tests use these to validate
# dMod's compiled / numerical machinery against the mathematics, not against
# another implementation.
#
# Conventions:
#   * Each truth_* function returns either a numeric vector (state values)
#     or a list with $value, $gradient, $hessian where useful.
#   * make_noisy_data() builds a data.frame in dMod's expected layout
#     (columns: name, time, value, sigma[, lloq], condition).
#   * numderiv_grad / numderiv_hess wrap numDeriv with project tolerances
#     and consistent argument order for use inside tests.


## ---- Closed-form solutions ----------------------------------------------

# Linear decay  A(t) = x0 * exp(-k * t).
# Returns a list with the state and its first / second derivatives w.r.t.
# the parameters (x0, k). Vectorised over t.
truth_decay <- function(t, x0, k) {
  E <- exp(-k * t)
  value <- x0 * E
  grad  <- cbind(x0 = E, k = -t * x0 * E)
  # d2A / d(x0)^2 = 0, d2A / dx0 dk = -t * E, d2A / dk^2 = t^2 * x0 * E
  hess <- array(0, c(length(t), 2L, 2L),
                dimnames = list(NULL, c("x0", "k"), c("x0", "k")))
  hess[, "x0", "k"] <- -t * E
  hess[, "k", "x0"] <- -t * E
  hess[, "k", "k"]  <- t^2 * x0 * E
  list(value = value, gradient = grad, hessian = hess)
}

# Two-step linear cascade A -> B -> C with rates (k1, k2), initial
# A(0) = A0, B(0) = C(0) = 0. Returns A, B, C at times t.
# Branch on k1 == k2 to keep the formula well-defined.
truth_two_step <- function(t, A0, k1, k2) {
  A <- A0 * exp(-k1 * t)
  if (isTRUE(all.equal(k1, k2))) {
    B <- A0 * k1 * t * exp(-k1 * t)
  } else {
    B <- A0 * k1 / (k2 - k1) * (exp(-k1 * t) - exp(-k2 * t))
  }
  C <- A0 - A - B
  data.frame(time = t, A = A, B = B, C = C)
}


## ---- Data simulation -----------------------------------------------------

# Simulate a noisy data frame in dMod's expected layout. truth_fn is a
# function (times, pars) -> numeric vector of length(times).
make_noisy_data <- function(truth_fn, pars, times,
                            name = "y",
                            sigma = 0.05,
                            condition = "C1",
                            lloq = NULL,
                            seed = 1L) {
  set.seed(seed)
  vals_true <- truth_fn(times, pars)
  noise <- rnorm(length(times), mean = 0, sd = sigma)
  out <- data.frame(
    name      = name,
    time      = times,
    value     = vals_true + noise,
    sigma     = sigma,
    condition = condition,
    stringsAsFactors = FALSE
  )
  if (!is.null(lloq)) out$lloq <- lloq
  out
}


## ---- numDeriv wrappers ---------------------------------------------------

# Thin wrappers around numDeriv with a single source-of-truth tolerance for
# behavioral tests. Pass-through ... to numDeriv for method tweaks.
numderiv_grad <- function(fn, x, ...) {
  if (!requireNamespace("numDeriv", quietly = TRUE))
    testthat::skip("numDeriv not installed")
  numDeriv::grad(fn, x, method = "Richardson", ...)
}

numderiv_hess <- function(fn, x, ...) {
  if (!requireNamespace("numDeriv", quietly = TRUE))
    testthat::skip("numDeriv not installed")
  numDeriv::hessian(fn, x, method = "Richardson", ...)
}


## ---- Closed-form objective values ---------------------------------------

# Closed-form value of the (negative) log-likelihood that nll_ALOQ computes
# for an above-LOQ Gaussian dataset:
#   value = sum(((pred - obs) / sigma)^2) + sum(log(2 * pi * sigma^2))
# pred, obs, sigma are length-n numeric. Used as an independent reference
# for normL2 / nll behavioral tests.
truth_nll_aloq <- function(pred, obs, sigma) {
  wr <- (pred - obs) / sigma
  sum(wr^2) + sum(log(2 * pi * sigma^2))
}

# Closed-form value of nll_BLOQ under M3:
#   value = -2 * sum(log(Phi(-(pred - lloq) / sigma)))
# Where pred is the model prediction at BLOQ time points and lloq is the
# substituted observation (data$value is pmax-ed to lloq before residual).
truth_nll_bloq_m3 <- function(pred, lloq, sigma) {
  wr <- (pred - lloq) / sigma
  -2 * sum(stats::pnorm(-wr, log.p = TRUE))
}

# Closed-form M4NM / M4BEAL BLOQ contribution:
#   value = -2 * sum(log(1 - Phi(wr) / Phi(w0)))
# with wr = (pred - lloq) / sigma and w0 = pred / sigma. Used as a closed
# reference even though the dMod implementation adds stability fallbacks for
# extreme arguments (which we test separately under nominal conditions).
truth_nll_bloq_m4 <- function(pred, lloq, sigma) {
  wr <- (pred - lloq) / sigma
  w0 <- pred / sigma
  -2 * sum(log(1 - stats::pnorm(wr) / stats::pnorm(w0)))
}
