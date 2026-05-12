# Parity tests for the shared C++ residual kernel
# (src/residual_kernel.{h,cpp}, exposed via residual_kernel_aloq /
# residual_kernel_bloq in src/residual_kernel_rwrap.cpp).
#
# For each (bloq_mode, sigma_depends_on_par, use_deriv2_exact) combination
# we compare the C++ result against R's nll_ALOQ / nll_BLOQ on randomized
# inputs. For the new BLOQ-deriv2-exact paths (which have no R reference),
# we additionally verify the Hessian against a finite-difference Hessian of
# the underlying objective.

# Helper: build a randomized residual-kernel test problem with quadratic
# pred/sigma dependence so d2pred is non-trivial and FD-Hessian tests make
# sense. theta0 is the evaluation point. n_obs ALOQ rows + n_bloq BLOQ rows.
make_setup <- function(n_obs, n_par, sigma_dep, n_bloq = 0, seed = 1L) {
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
    y_data      = y_data,
    lloq        = lloq
  )
}

# Predicts pred(theta), sigma(theta) under the quadratic model used by
# make_setup. Used for FD-Hessian checks.
eval_model <- function(setup, theta) {
  dt <- theta - setup$theta0
  pred  <- setup$pred0  + setup$dpred0  %*% dt
  sigma <- setup$sigma0 + setup$dsigma0 %*% dt
  # Quadratic contribution: 0.5 * dt^T Q[i] dt per row
  for (i in seq_len(setup$n_obs)) {
    pred[i, 1] <- pred[i, 1] + 0.5 * sum(dt * (setup$Q[i, , ] %*% dt))
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

# ---- ALOQ: parity vs dMod:::nll_ALOQ ----
test_that("ALOQ kernel matches R nll_ALOQ across (sigma_dep, deriv2, M4BEAL)", {
  for (sigma_dep in c(FALSE, TRUE)) {
    for (use_d2 in c(FALSE, TRUE)) {
      for (bloq_mode in c("NONE", "M3", "M4BEAL")) {
        # M4BEAL needs lloq > 0; we just feed valid inputs with no BLOQ rows.
        setup <- make_setup(n_obs = 12, n_par = 4, sigma_dep = sigma_dep,
                            n_bloq = 0, seed = 1L)

        # Subset to ALOQ rows.
        pred  <- setup$pred0
        sigma <- setup$sigma0
        dpred <- setup$dpred0
        dsigma <- if (sigma_dep) setup$dsigma0 else NULL
        d2pred <- if (use_d2) setup$Q else NULL
        y_data <- setup$y_data

        # R reference via nll_ALOQ.
        nout <- data.table::data.table(
          weighted.residual = (pred - y_data) / sigma,
          weighted.0        = pred / sigma,
          sigma             = sigma
        )
        R_ref <- dMod:::nll_ALOQ(
          nout            = nout,
          derivs          = dpred,
          derivs.err      = if (sigma_dep) dsigma else NULL,
          derivs2         = d2pred,
          par_names       = setup$par_names,
          opt.BLOQ        = bloq_mode,
          opt.hessian     = c(ALOQ_part1 = TRUE, ALOQ_part2 = TRUE,
                              ALOQ_part3 = TRUE,
                              BLOQ_part1 = TRUE, BLOQ_part2 = TRUE,
                              BLOQ_part3 = TRUE),
          bessel.correction = 1
        )

        opts <- default_opts()
        opts$use_deriv2_exact     <- use_d2
        opts$bloq_mode            <- bloq_mode
        opts$sigma_depends_on_par <- sigma_dep
        cpp <- dMod:::residual_kernel_aloq(
          pred = pred, dpred = dpred, d2pred = d2pred,
          y_data = y_data, sigma = sigma,
          dsigma = dsigma, lloq = NULL, opts = opts
        )

        label <- sprintf("ALOQ sigma_dep=%s use_d2=%s bloq=%s",
                         sigma_dep, use_d2, bloq_mode)
        expect_equal(cpp$value,    R_ref$value,    tolerance = 1e-10, info = label)
        expect_equal(unname(cpp$gradient), unname(R_ref$gradient),
                     tolerance = 1e-10, info = label)
        # Hessian: cpp returns unnamed, R returns named. Compare numeric.
        expect_equal(unname(cpp$hessian), unname(R_ref$hessian),
                     tolerance = 1e-10, info = label)
        # Symmetry
        expect_lt(max(abs(cpp$hessian - t(cpp$hessian))), 1e-12,
                  label = paste0(label, " symmetry"))
      }
    }
  }
})


# ---- BLOQ: parity vs dMod:::nll_BLOQ (no deriv2_exact, matches R math) ----
test_that("BLOQ kernel matches R nll_BLOQ across modes and sigma dependence", {
  for (sigma_dep in c(FALSE, TRUE)) {
    for (bloq_mode in c("M3", "M4NM", "M4BEAL")) {
      setup <- make_setup(n_obs = 0, n_par = 4, sigma_dep = sigma_dep,
                          n_bloq = 10, seed = 2L)
      pred  <- setup$pred0
      sigma <- setup$sigma0
      dpred <- setup$dpred0
      dsigma <- if (sigma_dep) setup$dsigma0 else NULL
      y_data <- setup$y_data
      lloq   <- setup$lloq

      nout_bloq <- data.table::data.table(
        weighted.residual = (pred - y_data) / sigma,
        weighted.0        = pred / sigma,
        sigma             = sigma,
        value             = y_data   # used by M4 NM/BEAL guard against lloq < 0
      )
      R_ref <- dMod:::nll_BLOQ(
        nout.bloq       = nout_bloq,
        derivs.bloq     = dpred,
        derivs.err.bloq = if (sigma_dep) dsigma else NULL,
        par_names       = setup$par_names,
        opt.BLOQ        = bloq_mode,
        opt.hessian     = c(BLOQ_part1 = TRUE, BLOQ_part2 = TRUE,
                            BLOQ_part3 = TRUE)
      )

      opts <- default_opts()
      opts$use_deriv2_exact     <- FALSE
      opts$bloq_mode            <- bloq_mode
      opts$sigma_depends_on_par <- sigma_dep
      cpp <- dMod:::residual_kernel_bloq(
        pred = pred, dpred = dpred, d2pred = NULL,
        y_data = y_data, sigma = sigma,
        dsigma = dsigma, lloq = lloq, opts = opts
      )

      label <- sprintf("BLOQ sigma_dep=%s bloq=%s", sigma_dep, bloq_mode)
      expect_equal(cpp$value, R_ref$value, tolerance = 1e-10, info = label)
      expect_equal(unname(cpp$gradient), unname(R_ref$gradient),
                   tolerance = 1e-10, info = label)
      expect_equal(unname(cpp$hessian), unname(R_ref$hessian),
                   tolerance = 1e-10, info = label)
      expect_lt(max(abs(cpp$hessian - t(cpp$hessian))), 1e-12,
                label = paste0(label, " symmetry"))
    }
  }
})


# ---- ALOQ + bessel parity ----
test_that("ALOQ kernel matches R nll_ALOQ under bessel correction", {
  setup <- make_setup(n_obs = 8, n_par = 3, sigma_dep = TRUE, seed = 5L)
  bessel <- 1.13
  nout <- data.table::data.table(
    weighted.residual = (setup$pred0 - setup$y_data) / setup$sigma0,
    weighted.0        = setup$pred0 / setup$sigma0,
    sigma             = setup$sigma0
  )
  R_ref <- dMod:::nll_ALOQ(
    nout            = nout,
    derivs          = setup$dpred0,
    derivs.err      = setup$dsigma0,
    derivs2         = NULL,
    par_names       = setup$par_names,
    opt.BLOQ        = "M3",
    opt.hessian     = c(ALOQ_part1 = TRUE, ALOQ_part2 = TRUE,
                        ALOQ_part3 = TRUE,
                        BLOQ_part1 = TRUE, BLOQ_part2 = TRUE,
                        BLOQ_part3 = TRUE),
    bessel.correction = bessel
  )
  opts <- default_opts()
  opts$bessel <- bessel
  opts$sigma_depends_on_par <- TRUE
  cpp <- dMod:::residual_kernel_aloq(
    pred = setup$pred0, dpred = setup$dpred0, d2pred = NULL,
    y_data = setup$y_data, sigma = setup$sigma0,
    dsigma = setup$dsigma0, lloq = NULL, opts = opts
  )
  expect_equal(cpp$value, R_ref$value, tolerance = 1e-10)
  expect_equal(unname(cpp$gradient), unname(R_ref$gradient),
               tolerance = 1e-10)
  expect_equal(unname(cpp$hessian), unname(R_ref$hessian),
               tolerance = 1e-10)
})


# ---- BLOQ-deriv2-exact: FD validation (new functionality) ----
test_that("BLOQ M3 + use_deriv2_exact Hessian matches finite-difference Hessian", {
  skip_if_not_installed("numDeriv")
  for (sigma_dep in c(FALSE, TRUE)) {
    setup <- make_setup(n_obs = 0, n_par = 3, sigma_dep = sigma_dep,
                        n_bloq = 8, seed = 7L)
    # M3 BLOQ objective as a function of theta.
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
      dsigma = if (sigma_dep) setup$dsigma0 else NULL,
      lloq = setup$lloq, opts = opts
    )

    # Without sigma_dep, full Hessian = FD Hessian to ~1e-5 (Richardson order).
    # With sigma_dep, kernel keeps GN approximation on the sigma cross terms,
    # so we only check that the d2pred contribution is in the right ballpark
    # (and that the symmetry is preserved); the bulk should still be close.
    tol <- if (sigma_dep) 5e-2 else 1e-4
    label <- sprintf("M3 deriv2 FD-Hess sigma_dep=%s", sigma_dep)
    expect_lt(max(abs(cpp$hessian - H_fd)), tol, label = label)
    expect_lt(max(abs(cpp$hessian - t(cpp$hessian))), 1e-12,
              label = paste0(label, " symmetry"))
  }
})


test_that("BLOQ M4NM + use_deriv2_exact: kernel adds d2pred contribution", {
  skip_if_not_installed("numDeriv")
  # NOTE: R's nll_BLOQ Hessian for M4 modes is not a clean
  # "Newton-minus-d2-wr" Gauss-Newton form (the GN parts use specific A1..A6
  # coefficients whose sum does not equal the f''-outer-product of the true
  # Newton Hessian). So even with use_deriv2_exact = TRUE we cannot expect
  # the C++ kernel to match the FD Hessian exactly under M4 modes; what we
  # DO get back from the d2 path is the f' * d^2 pred contribution that R's
  # GN form drops. We verify that contribution is non-trivial and that the
  # Hessian remains symmetric.
  setup <- make_setup(n_obs = 0, n_par = 3, sigma_dep = FALSE,
                      n_bloq = 8, seed = 9L)
  opts <- default_opts()
  opts$bloq_mode <- "M4NM"

  opts$use_deriv2_exact <- FALSE
  cpp_gn <- dMod:::residual_kernel_bloq(
    pred = setup$pred0, dpred = setup$dpred0, d2pred = NULL,
    y_data = setup$y_data, sigma = setup$sigma0,
    dsigma = NULL, lloq = setup$lloq, opts = opts
  )
  opts$use_deriv2_exact <- TRUE
  cpp_ex <- dMod:::residual_kernel_bloq(
    pred = setup$pred0, dpred = setup$dpred0, d2pred = setup$Q,
    y_data = setup$y_data, sigma = setup$sigma0,
    dsigma = NULL, lloq = setup$lloq, opts = opts
  )
  # The d2 path must change the Hessian (else the code path is dead).
  expect_gt(max(abs(cpp_ex$hessian - cpp_gn$hessian)), 1e-6,
            label = "M4NM deriv2 changes Hessian")
  # Symmetry preserved.
  expect_lt(max(abs(cpp_ex$hessian - t(cpp_ex$hessian))), 1e-12,
            label = "M4NM deriv2 symmetry")
})


# ---- Edge case: zero rows ----
test_that("kernels are no-ops on zero-row inputs", {
  opts <- default_opts()
  opts$bloq_mode <- "M3"
  cpp_aloq <- dMod:::residual_kernel_aloq(
    pred = numeric(0), dpred = matrix(0, 0, 3),
    d2pred = NULL, y_data = numeric(0), sigma = numeric(0),
    dsigma = NULL, lloq = NULL, opts = default_opts()
  )
  cpp_bloq <- dMod:::residual_kernel_bloq(
    pred = numeric(0), dpred = matrix(0, 0, 3),
    d2pred = NULL, y_data = numeric(0), sigma = numeric(0),
    dsigma = NULL, lloq = NULL, opts = opts
  )
  expect_equal(cpp_aloq$value, 0)
  expect_equal(sum(abs(cpp_aloq$gradient)), 0)
  expect_equal(sum(abs(cpp_aloq$hessian)), 0)
  expect_equal(cpp_bloq$value, 0)
})
