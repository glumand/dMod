# ============================================================================
# Residual / per-data-point likelihood primitives.
#
# Sections:
#   * res()         - data <-> prediction residual operator
#   * datapointL2() - validation-point L2 constraint that reuses
#                     env$prediction populated by an upstream normL2
#   * C++ residual kernel (accumulate_aloq / accumulate_bloq) - FD-validated
#
# These tests build prediction matrices by hand to keep the math
# transparent. End-to-end coverage is in test-normL2.R.
# ============================================================================


# Shared helper: prediction matrix for A(t) = exp(-0.5 * t) with analytical
# sensitivities dA/dA_par, dA/dk_par.
.make_prdframe <- function(times = c(0, 1, 2)) {
  prdf <- matrix(c(times, exp(-0.5 * times)),
                 nrow = length(times), ncol = 2,
                 dimnames = list(NULL, c("time", "A")))
  d_arr <- array(0, c(length(times), 1, 2),
                 dimnames = list(NULL, "A", c("A_par", "k_par")))
  d_arr[, "A", "A_par"] <- exp(-0.5 * times)
  d_arr[, "A", "k_par"] <- -times * 1.0 * exp(-0.5 * times)
  attr(prdf, "deriv") <- d_arr
  prdf
}


# ---- res ----------------------------------------------------------------

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


test_that("res sets value = pmax(value_raw, lloq) and bloq mask correctly", {
  prdf <- .make_prdframe()
  dat <- data.frame(time = c(1, 2), name = "A",
                    value = c(0.55, 0.05),  # 0.05 < lloq -> censored
                    sigma = c(0.1, 0.1),
                    lloq  = c(-Inf, 0.10),
                    stringsAsFactors = FALSE)
  r <- res(dat, prdf)
  expect_equal(r$value, c(0.55, 0.10), tolerance = 1e-12)
  expect_equal(r$bloq, c(FALSE, TRUE))
  expect_equal(r$weighted.residual[2], (exp(-1) - 0.10) / 0.1,
               tolerance = 1e-12)
})


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


test_that("res errors if an observable in data is missing from the prediction", {
  prdf <- .make_prdframe()
  dat <- data.frame(time = 1, name = "B",  # B does not exist in prdf
                    value = 0.5, sigma = 0.1, lloq = -Inf,
                    stringsAsFactors = FALSE)
  expect_error(res(dat, prdf), regexp = "Observable not found")
})


# ---- datapointL2 --------------------------------------------------------

# datapointL2 penalises (prediction[name, t] - target_par) / sigma in L2.
# It reads prediction from the shared `env` set by an upstream normL2.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


test_that("datapointL2 value equals ((pred - target) / sigma)^2 at a known point", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  obj_main <- normL2(data, bench$prd_id)
  obj_val  <- datapointL2(name = "y", time = 2.0, value = "newpoint",
                          sigma = 0.05, condition = "C1")
  obj <- obj_main + obj_val

  pars <- c(bench$outerpars_id, newpoint = 0.5)
  o    <- obj(pars)
  o_main <- obj_main(bench$outerpars_id)

  pred_at_t <- unname(bench$prd_id(times = c(0, 2), pars = bench$outerpars_id,
                                   deriv = FALSE)$C1[2, "y"])
  expected_val <- ((pred_at_t - pars[["newpoint"]]) / 0.05)^2
  expect_equal(unname(o$value - o_main$value), expected_val, tolerance = 1e-3)
})


test_that("normL2 + datapointL2 value equals the sum of the parts", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  obj_main <- normL2(data, bench$prd_id)
  obj_val  <- datapointL2(name = "y", time = 5.0, value = "vtarget",
                          sigma = 0.1, condition = "C1")
  obj <- obj_main + obj_val
  pars <- c(bench$outerpars_id, vtarget = 0.08)

  # +.objfn evaluates main first, then datapointL2 reuses env$prediction.
  v_main <- obj_main(bench$outerpars_id)
  env_main <- attr(v_main, "env")
  v_pt   <- obj_val(pars = pars, env = env_main)

  v_total <- obj(pars)
  expect_equal(v_total$value, v_main$value + v_pt$value, tolerance = 1e-10)
})


test_that("datapointL2 contribution to the combined gradient follows the closed form", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  obj_main <- normL2(data, bench$prd_id)
  obj_val  <- datapointL2(name = "y", time = 3.0, value = "target",
                          sigma = 0.05, condition = "C1")
  obj <- obj_main + obj_val
  pars <- c(A = 1.0, k = 0.5, target = 0.2)

  # val_pt = ((pred(3) - target) / sigma_pt)^2 with pred(t) = A * exp(-k*t).
  # dpt/dA      =  2 (pred - target) * exp(-k*t) / sigma_pt^2
  # dpt/dk      =  2 (pred - target) * (-t * A * exp(-k*t)) / sigma_pt^2
  # dpt/dtarget = -2 (pred - target) / sigma_pt^2
  t_pt <- 3.0; sigma_pt <- 0.05
  pred_pt <- pars[["A"]] * exp(-pars[["k"]] * t_pt)
  r <- pred_pt - pars[["target"]]
  contrib <- c(A      = 2 * r * exp(-pars[["k"]] * t_pt) / sigma_pt^2,
               k      = 2 * r * (-t_pt * pars[["A"]] *
                                   exp(-pars[["k"]] * t_pt)) / sigma_pt^2,
               target = -2 * r / sigma_pt^2)

  g_combined <- obj(pars)$gradient
  g_main     <- obj_main(pars[c("A", "k")])$gradient

  expect_equal(unname(g_combined[["A"]]),
               unname(g_main[["A"]] + contrib[["A"]]), tolerance = 1e-3)
  expect_equal(unname(g_combined[["k"]]),
               unname(g_main[["k"]] + contrib[["k"]]), tolerance = 1e-3)
  expect_equal(unname(g_combined[["target"]]), unname(contrib[["target"]]),
               tolerance = 1e-3)
})


# ============================================================================
# C++ residual kernel parity (src/residual_kernel.{h,cpp})
# ============================================================================

# Helper: build a randomized residual-kernel test problem with quadratic
# pred/sigma dependence so d2pred is non-trivial and FD-Hessian tests make
# sense. theta0 is the evaluation point. n_obs ALOQ rows + n_bloq BLOQ rows.
make_setup <- function(n_obs, n_par, sigma_dep, n_bloq = 0, seed = 1L,
                       sigma_curved = FALSE) {
  set.seed(seed)
  n_total <- n_obs + n_bloq
  par_names <- paste0("p", seq_len(n_par))
  theta0    <- runif(n_par, -0.5, 0.5)
  names(theta0) <- par_names

  # Linear-plus-quadratic pred and sigma in theta, evaluated at theta0.
  pred0  <- runif(n_total, 1.0, 3.0)
  dpred0 <- matrix(runif(n_total * n_par, -0.4, 0.4), n_total, n_par,
                   dimnames = list(NULL, par_names))
  Q <- array(0.0, dim = c(n_total, n_par, n_par),
             dimnames = list(NULL, par_names, par_names))
  for (i in seq_len(n_total)) {
    M <- matrix(runif(n_par * n_par, -0.2, 0.2), n_par, n_par)
    Q[i, , ] <- 0.5 * (M + t(M))  # symmetric mixed partials
  }
  sigma0 <- runif(n_total, 0.4, 1.0)
  dsigma0 <- if (sigma_dep)
    matrix(runif(n_total * n_par, -0.05, 0.05), n_total, n_par,
           dimnames = list(NULL, par_names))
  else
    matrix(0.0, n_total, n_par, dimnames = list(NULL, par_names))
  # Optional symmetric quadratic sigma term so d2sigma is non-zero. Kept
  # small so sigma(theta) stays strictly positive in a neighbourhood of
  # theta0 (Richardson FD perturbs theta by ~1e-2).
  Q_sigma <- array(0.0, dim = c(n_total, n_par, n_par),
                   dimnames = list(NULL, par_names, par_names))
  if (sigma_curved && sigma_dep) {
    for (i in seq_len(n_total)) {
      M <- matrix(runif(n_par * n_par, -0.02, 0.02), n_par, n_par)
      Q_sigma[i, , ] <- 0.5 * (M + t(M))
    }
  }

  # y_data: ALOQ rows random with some signal; BLOQ rows < lloq.
  lloq <- rep(-Inf, n_total)
  y_data <- runif(n_total, 0.5, 2.5)
  if (n_bloq > 0) {
    bloq_rows <- seq_len(n_bloq) + n_obs
    lloq[bloq_rows]  <- 2.5
    y_data[bloq_rows] <- 2.5  # res(): val = pmax(value, lloq)
  }

  list(
    n_obs       = n_total,
    n_aloq      = n_obs,
    n_bloq      = n_bloq,
    n_par       = n_par,
    par_names   = par_names,
    theta0      = theta0,
    pred0       = pred0,
    dpred0      = dpred0,
    sigma0      = sigma0,
    dsigma0     = dsigma0,
    Q           = Q,
    Q_sigma     = Q_sigma,
    y_data      = y_data,
    lloq        = lloq
  )
}

# Predicts pred(theta), sigma(theta) under the quadratic model used by
# make_setup. Used for FD-Hessian checks. Both pred and sigma can have a
# symmetric quadratic theta-dependence (Q, Q_sigma).
eval_model <- function(setup, theta) {
  dt <- theta - setup$theta0
  pred  <- setup$pred0  + setup$dpred0  %*% dt
  sigma <- setup$sigma0 + setup$dsigma0 %*% dt
  # Quadratic contribution: 0.5 * dt^T Q[i] dt per row (pred and sigma).
  for (i in seq_len(setup$n_obs)) {
    pred[i, 1]  <- pred[i, 1]  + 0.5 * sum(dt * (setup$Q[i, , ]       %*% dt))
    sigma[i, 1] <- sigma[i, 1] + 0.5 * sum(dt * (setup$Q_sigma[i, , ] %*% dt))
  }
  list(pred = as.numeric(pred), sigma = as.numeric(sigma))
}

default_opts <- function() {
  list(
    use_deriv2_exact     = FALSE,
    bloq_mode            = "NONE",
    sigma_depends_on_par = FALSE,
    d2sigma_present      = FALSE,
    bessel               = 1.0,
    aloq_part1 = TRUE, aloq_part2 = TRUE, aloq_part3 = TRUE,
    bloq_part1 = TRUE, bloq_part2 = TRUE, bloq_part3 = TRUE
  )
}

# ---- BLOQ-deriv2-exact: FD validation (new functionality) ----
# With BOTH d2pred and d2sigma propagated, the kernel produces the analytical
# Hessian and matches FD to Richardson order for both the sigma-independent
# and sigma-dependent (curved sigma) cases.
test_that("BLOQ M3 + use_deriv2_exact Hessian matches finite-difference Hessian", {
  skip_if_not_installed("numDeriv")
  for (sigma_dep in c(FALSE, TRUE)) {
    setup <- make_setup(n_obs = 0, n_par = 3, sigma_dep = sigma_dep,
                        n_bloq = 8, seed = 7L,
                        sigma_curved = sigma_dep)
    obj_fn <- function(theta) {
      m  <- eval_model(setup, theta)
      wr <- (m$pred - setup$y_data) / m$sigma
      sum(-2 * stats::pnorm(-wr, log.p = TRUE))
    }
    H_fd <- numDeriv::hessian(obj_fn, setup$theta0, method = "Richardson")

    opts <- default_opts()
    opts$use_deriv2_exact     <- TRUE
    opts$bloq_mode            <- "M3"
    opts$sigma_depends_on_par <- sigma_dep
    cpp <- dMod:::residual_kernel_bloq(
      pred = setup$pred0, dpred = setup$dpred0, d2pred = setup$Q,
      y_data = setup$y_data, sigma = setup$sigma0,
      dsigma  = if (sigma_dep) setup$dsigma0 else NULL,
      d2sigma = if (sigma_dep) setup$Q_sigma else NULL,
      lloq = setup$lloq, opts = opts
    )
    label <- sprintf("M3 deriv2 FD-Hess sigma_dep=%s", sigma_dep)
    expect_lt(max(abs(cpp$hessian - H_fd)), 1e-4, label = label)
    expect_lt(max(abs(cpp$hessian - t(cpp$hessian))), 1e-12,
              label = paste0(label, " symmetry"))
  }
})


# The M4* BLOQ contribution -2 log(1 - Phi(wr)/Phi(w0)) has an exact analytic
# Hessian (not a Gauss-Newton surrogate). Validate it against finite differences
# for both M4NM and M4BEAL (which share the BLOQ formula) and for fixed and
# parameter-dependent (curved) sigma. This also exercises the d2pred / d2sigma
# exact-Hessian path.
test_that("BLOQ M4NM/M4BEAL + use_deriv2_exact Hessian matches finite-difference Hessian", {
  skip_if_not_installed("numDeriv")
  for (mode in c("M4NM", "M4BEAL")) {
    for (sigma_dep in c(FALSE, TRUE)) {
      setup <- make_setup(n_obs = 0, n_par = 3, sigma_dep = sigma_dep,
                          n_bloq = 8, seed = 9L, sigma_curved = sigma_dep)
      obj_fn <- function(theta) {
        m  <- eval_model(setup, theta)
        wr <- (m$pred - setup$y_data) / m$sigma
        w0 <- m$pred / m$sigma
        sum(-2 * log(1 - stats::pnorm(wr) / stats::pnorm(w0)))
      }
      H_fd <- numDeriv::hessian(obj_fn, setup$theta0, method = "Richardson")

      opts <- default_opts()
      opts$use_deriv2_exact     <- TRUE
      opts$bloq_mode            <- mode
      opts$sigma_depends_on_par <- sigma_dep
      opts$d2sigma_present      <- sigma_dep
      cpp <- dMod:::residual_kernel_bloq(
        pred = setup$pred0, dpred = setup$dpred0, d2pred = setup$Q,
        y_data = setup$y_data, sigma = setup$sigma0,
        dsigma  = if (sigma_dep) setup$dsigma0 else NULL,
        d2sigma = if (sigma_dep) setup$Q_sigma else NULL,
        lloq = setup$lloq, opts = opts
      )
      label <- sprintf("%s BLOQ deriv2 FD-Hess sigma_dep=%s", mode, sigma_dep)
      expect_lt(max(abs(cpp$hessian - H_fd)), 1e-4, label = label)
      expect_lt(max(abs(cpp$hessian - t(cpp$hessian))), 1e-12,
                label = paste0(label, " symmetry"))
    }
  }
})


# ---- ALOQ d2sigma-exact: FD validation ----
# Verify against the true Hessian of the ALOQ objective when both pred and
# sigma are quadratic in theta.
test_that("ALOQ d2sigma + d2pred Hessian matches finite-difference Hessian", {
  skip_if_not_installed("numDeriv")
  setup <- make_setup(n_obs = 10, n_par = 3, sigma_dep = TRUE,
                      n_bloq = 0, seed = 13L, sigma_curved = TRUE)
  obj_fn <- function(theta) {
    m  <- eval_model(setup, theta)
    wr <- (m$pred - setup$y_data) / m$sigma
    sum(wr^2 + log(2 * pi * m$sigma^2))
  }
  H_fd <- numDeriv::hessian(obj_fn, setup$theta0, method = "Richardson")

  opts <- default_opts()
  opts$use_deriv2_exact     <- TRUE
  opts$bloq_mode            <- "M3"
  opts$sigma_depends_on_par <- TRUE
  opts$d2sigma_present      <- TRUE
  cpp <- dMod:::residual_kernel_aloq(
    pred = setup$pred0, dpred = setup$dpred0, d2pred = setup$Q,
    y_data = setup$y_data, sigma = setup$sigma0,
    dsigma = setup$dsigma0, d2sigma = setup$Q_sigma,
    lloq = NULL, opts = opts
  )
  expect_lt(max(abs(cpp$hessian - H_fd)), 1e-4, label = "ALOQ d2 FD-Hess")
  expect_lt(max(abs(cpp$hessian - t(cpp$hessian))), 1e-12,
            label = "ALOQ d2 symmetry")
})


# ---- M4BEAL ALOQ correction: exact second-order FD validation ----
# M4BEAL adds +2 log Phi(w0) to every ALOQ row (truncation of the normal at 0).
# Its gradient and full Hessian, including the d2pred / d2sigma exact terms, are
# analytic and must match the finite-difference Hessian of the truncated
# objective wr^2 + log(2 pi sigma^2) + 2 log Phi(pred/sigma).
test_that("ALOQ M4BEAL + use_deriv2_exact Hessian matches finite-difference Hessian", {
  skip_if_not_installed("numDeriv")
  setup <- make_setup(n_obs = 8, n_par = 3, sigma_dep = TRUE,
                      n_bloq = 0, seed = 23L, sigma_curved = TRUE)
  obj_fn <- function(theta) {
    m  <- eval_model(setup, theta)
    wr <- (m$pred - setup$y_data) / m$sigma
    w0 <- m$pred / m$sigma
    sum(wr^2 + log(2 * pi * m$sigma^2) + 2 * stats::pnorm(w0, log.p = TRUE))
  }
  H_fd <- numDeriv::hessian(obj_fn, setup$theta0, method = "Richardson")

  opts <- default_opts()
  opts$use_deriv2_exact     <- TRUE
  opts$bloq_mode            <- "M4BEAL"
  opts$sigma_depends_on_par <- TRUE
  opts$d2sigma_present      <- TRUE
  cpp <- dMod:::residual_kernel_aloq(
    pred = setup$pred0, dpred = setup$dpred0, d2pred = setup$Q,
    y_data = setup$y_data, sigma = setup$sigma0,
    dsigma = setup$dsigma0, d2sigma = setup$Q_sigma,
    lloq = NULL, opts = opts
  )
  expect_lt(max(abs(cpp$hessian - H_fd)), 1e-4, label = "M4BEAL ALOQ d2 FD-Hess")
  expect_lt(max(abs(cpp$hessian - t(cpp$hessian))), 1e-12,
            label = "M4BEAL ALOQ d2 symmetry")
})


# ---- Edge case: zero rows ----
test_that("kernels are no-ops on zero-row inputs", {
  opts <- default_opts()
  opts$bloq_mode <- "M3"
  cpp_aloq <- dMod:::residual_kernel_aloq(
    pred = numeric(0), dpred = matrix(0, 0, 3),
    d2pred = NULL, y_data = numeric(0), sigma = numeric(0),
    dsigma = NULL, d2sigma = NULL, lloq = NULL, opts = default_opts()
  )
  cpp_bloq <- dMod:::residual_kernel_bloq(
    pred = numeric(0), dpred = matrix(0, 0, 3),
    d2pred = NULL, y_data = numeric(0), sigma = numeric(0),
    dsigma = NULL, d2sigma = NULL, lloq = NULL, opts = opts
  )
  expect_equal(cpp_aloq$value, 0)
  expect_equal(sum(abs(cpp_aloq$gradient)), 0)
  expect_equal(sum(abs(cpp_aloq$hessian)), 0)
  expect_equal(cpp_bloq$value, 0)
})
