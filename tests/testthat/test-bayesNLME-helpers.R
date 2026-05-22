# ============================================================================
# Helper functions used by the Bayesian NLME pipeline.
#
# Sections:
#   * priorOmega           - LKJ + half-Normal / half-Cauchy hyperprior on
#                            omega-Cholesky parameters
#   * foceiOmegaGradient   - analytical d OFV_FOCEI / d omega_chol via the
#                            envelope theorem along eta
#   * .computeMetricDerivative - dG/dtheta tensor for the full-RMALA
#                            Christoffel correction
#   * nonCenteredOmega     - z <-> eta = L z bijection
#
# All four are non-orchestrator helpers; the end-to-end Bayesian NLME flow
# is exercised in test-bayesNLME.R, and the sampler kernels in test-mcmc.R.
# ============================================================================


# Shared FD helpers (priorOmega + foceiOmegaGradient use these).
.fd_grad <- function(f, x, h = 1e-5) {
  K <- length(x)
  g <- numeric(K)
  for (i in seq_len(K)) {
    xp <- x; xm <- x
    xp[i] <- xp[i] + h; xm[i] <- xm[i] - h
    g[i] <- (f(xp) - f(xm)) / (2 * h)
  }
  setNames(g, names(x))
}

.fd_hess <- function(f, x, h = 1e-4) {
  K <- length(x)
  H <- matrix(0, K, K, dimnames = list(names(x), names(x)))
  for (i in seq_len(K)) for (j in seq_len(i)) {
    xpp <- x; xpm <- x; xmp <- x; xmm <- x
    xpp[i] <- xpp[i] + h; xpp[j] <- xpp[j] + h
    xpm[i] <- xpm[i] + h; xpm[j] <- xpm[j] - h
    xmp[i] <- xmp[i] - h; xmp[j] <- xmp[j] + h
    xmm[i] <- xmm[i] - h; xmm[j] <- xmm[j] - h
    H[i, j] <- (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4 * h * h)
    if (i != j) H[j, i] <- H[i, j]
  }
  H
}


# ---- priorOmega ---------------------------------------------------------

test_that("priorOmega LKJHalfNormal diag-only: value/grad/hess match FD", {
  om <- omega(eta = c("eta_A", "eta_B", "eta_C"), structure = "diag")
  pr <- priorOmega(om, kind = "LKJHalfNormal", lkjEta = 2.0, scaleSD = 0.7)
  pars <- c(omega_A_A = log(0.3), omega_B_B = log(0.5), omega_C_C = log(0.2))
  out <- pr(pars)

  f_val <- function(p) pr(p, deriv = FALSE)$value
  expect_equal(out$value, f_val(pars), tolerance = 1e-12)

  g_fd <- .fd_grad(f_val, pars)
  expect_equal(out$gradient[names(g_fd)], g_fd, tolerance = 1e-4)

  H_fd <- .fd_hess(f_val, pars)
  expect_equal(out$hessian[names(g_fd), names(g_fd)], H_fd, tolerance = 1e-3)
})


test_that("priorOmega LKJHalfNormal full structure: value/grad/hess match FD", {
  om <- omega(eta = c("eta_A", "eta_B"), structure = "full")
  pr <- priorOmega(om, kind = "LKJHalfNormal", lkjEta = 2.0, scaleSD = 1.0)
  pars <- c(omega_A_A = log(0.5),
            omega_B_B = log(0.3),
            omega_B_A = 0.15)
  out <- pr(pars)

  f_val <- function(p) pr(p, deriv = FALSE)$value
  g_fd  <- .fd_grad(f_val, pars)
  expect_equal(out$gradient[names(g_fd)], g_fd, tolerance = 1e-4)

  H_fd <- .fd_hess(f_val, pars)
  expect_equal(out$hessian[names(g_fd), names(g_fd)], H_fd, tolerance = 1e-3)
})


test_that("priorOmega LKJHalfCauchy: value/grad/hess match FD on full structure", {
  om <- omega(eta = c("eta_A", "eta_B", "eta_C"), structure = "full")
  pr <- priorOmega(om, kind = "LKJHalfCauchy", lkjEta = 1.5, scaleSD = 0.5)
  pars <- c(omega_A_A = log(0.4),
            omega_B_B = log(0.6),
            omega_C_C = log(0.3),
            omega_B_A = 0.2,
            omega_C_A = -0.1,
            omega_C_B = 0.05)
  out <- pr(pars)

  f_val <- function(p) pr(p, deriv = FALSE)$value
  g_fd  <- .fd_grad(f_val, pars)
  expect_equal(out$gradient[names(g_fd)], g_fd, tolerance = 1e-3)

  H_fd <- .fd_hess(f_val, pars)
  expect_equal(out$hessian[names(g_fd), names(g_fd)], H_fd, tolerance = 5e-3)
})


test_that("priorOmega: Hessian is symmetric and block-diagonal by row of L", {
  om <- omega(eta = c("eta_A", "eta_B", "eta_C"), structure = "full")
  pr <- priorOmega(om, kind = "LKJHalfNormal", lkjEta = 2.0, scaleSD = 1.0)
  pars <- c(omega_A_A = log(0.4),
            omega_B_B = log(0.6),
            omega_C_C = log(0.3),
            omega_B_A = 0.2,
            omega_C_A = -0.1,
            omega_C_B = 0.05)
  out <- pr(pars)
  H <- out$hessian
  expect_equal(H, t(H), tolerance = 1e-12)

  # Block structure: pairs (m, n) from different rows of L must be 0.
  cholLoc <- om$cholLoc
  for (m in seq_len(nrow(cholLoc))) for (n in seq_len(nrow(cholLoc))) {
    if (cholLoc[m, "row"] != cholLoc[n, "row"]) {
      expect_equal(H[om$cholPars[m], om$cholPars[n]], 0, tolerance = 1e-12)
    }
  }
})


test_that("priorOmega: gradient at the prior mode is small", {
  # For half-Normal(0, 1) on each sigma_k with diag-only omega, the gradient
  # of -2 log p at omega_kk = 0 (-> sigma_k = 1) is
  #   d/d omega_kk = -2 * (1 - L_kk^2 / tau^2) = -2 * (1 - 1) = 0.
  om <- omega(eta = c("eta_A", "eta_B"), structure = "diag")
  pr <- priorOmega(om, kind = "LKJHalfNormal", lkjEta = 2.0, scaleSD = 1.0)
  pars <- c(omega_A_A = 0, omega_B_B = 0)
  out <- pr(pars)
  expect_lt(max(abs(out$gradient)), 1e-12)
})


# ---- foceiOmegaGradient -------------------------------------------------

# Setup: invent N subject mode vectors eta_hat_i and per-subject data
# Hessian contributions H_GN_data_i. The OFV-as-a-function-of-omega is
#   f(omega) = sum_i [ eta_i^T Omega(omega)^-1 eta_i
#                    + log|Omega(omega)|
#                    + log|H_GN_data_i + 2 Omega(omega)^-1| ]
.make_toy_omega_problem <- function(K, N, seed = 1L) {
  set.seed(seed)
  eta <- paste0("eta_", letters[seq_len(K)])
  om  <- omega(eta = eta, structure = "full")

  etaModes <- matrix(rnorm(N * K, sd = 0.3), N, K,
                     dimnames = list(NULL, eta))

  H_data <- lapply(seq_len(N), function(i) {
    A <- matrix(rnorm(K * K), K, K)
    crossprod(A) + diag(K)
  })

  list(K = K, N = N, om = om, etaModes = etaModes, H_data = H_data)
}


# OFV-as-function-of-omega given fixed eta_modes and H_data per subject.
.ofv_omega <- function(omegaChol, om, etaModes, H_data) {
  L <- om$buildL(omegaChol[om$cholPars])
  Omega <- tcrossprod(L)
  Omega_inv <- chol2inv(t(L))
  K <- om$K
  log_det_O <- 2 * sum(log(diag(L)))

  s <- 0
  for (i in seq_along(H_data)) {
    eta_i <- as.numeric(etaModes[i, ])
    quad  <- drop(crossprod(eta_i, Omega_inv %*% eta_i))
    H_i   <- H_data[[i]] + 2 * Omega_inv
    log_det_H <- determinant(H_i, logarithm = TRUE)$modulus
    s <- s + quad + log_det_O + as.numeric(log_det_H)
  }
  s
}

# Per-subject H_inv at the supplied omegaChol value.
.inv_list_at <- function(omegaChol, om, H_data) {
  L <- om$buildL(omegaChol[om$cholPars])
  Omega_inv <- chol2inv(t(L))
  lapply(H_data, function(H_d) solve(H_d + 2 * Omega_inv))
}


test_that("foceiOmegaGradient matches FD on diag Omega (K=3, N=8)", {
  pr <- .make_toy_omega_problem(K = 3, N = 8, seed = 11L)
  om_diag <- omega(eta = paste0("eta_", letters[1:3]), structure = "diag")
  pr$om <- om_diag

  oc <- setNames(log(c(0.4, 0.3, 0.5)), om_diag$cholPars)

  H_inv <- .inv_list_at(oc, om_diag, pr$H_data)
  g_ana <- foceiOmegaGradient(oc, om_diag, pr$etaModes, H_inv)
  g_fd  <- .fd_grad(function(x) .ofv_omega(x, om_diag, pr$etaModes, pr$H_data), oc)
  expect_equal(unname(g_ana), unname(g_fd), tolerance = 1e-4)
})


test_that("foceiOmegaGradient matches FD on full Omega (K=3, N=8)", {
  pr <- .make_toy_omega_problem(K = 3, N = 8, seed = 22L)
  om <- pr$om
  oc <- setNames(c(log(0.5), 0.15, log(0.4),
                   -0.10, 0.05, log(0.6)),
                 om$cholPars)

  H_inv <- .inv_list_at(oc, om, pr$H_data)
  g_ana <- foceiOmegaGradient(oc, om, pr$etaModes, H_inv)
  g_fd  <- .fd_grad(function(x) .ofv_omega(x, om, pr$etaModes, pr$H_data), oc)
  expect_equal(unname(g_ana), unname(g_fd), tolerance = 1e-3)
})


test_that("foceiOmegaGradient: K=2 N=4 full structure agrees with FD", {
  pr <- .make_toy_omega_problem(K = 2, N = 4, seed = 33L)
  om <- pr$om
  oc <- setNames(c(log(0.3), 0.1, log(0.5)), om$cholPars)

  H_inv <- .inv_list_at(oc, om, pr$H_data)
  g_ana <- foceiOmegaGradient(oc, om, pr$etaModes, H_inv)
  g_fd  <- .fd_grad(function(x) .ofv_omega(x, om, pr$etaModes, pr$H_data), oc)
  expect_equal(unname(g_ana), unname(g_fd), tolerance = 1e-4)
})


test_that("foceiOmegaGradient: K=1 reduces to closed-form", {
  set.seed(44L)
  K <- 1; N <- 6
  om <- omega(eta = "eta_a", structure = "diag")
  etaModes <- matrix(rnorm(N, sd = 0.3), N, 1, dimnames = list(NULL, "eta_a"))
  H_data <- lapply(seq_len(N), function(i) matrix(2 + abs(rnorm(1)), 1, 1))

  oc <- c(omega_a_a = log(0.4))
  Omega_inv <- exp(-2 * oc)
  # Analytical d OFV / d omega_aa (full hand-derivation, K=1 case):
  #   q_i = eta_i^2 Omega_inv,  d q_i / d omega = -2 q_i
  #   log|Omega| = 2 omega,     d / d omega = 2
  #   H_i = H_data_i + 2 Omega_inv, d H_i / d omega = -4 Omega_inv
  #   d log|H_i|/d omega = -4 Omega_inv / H_i
  H_inv <- lapply(H_data, function(H_d) 1 / (H_d + 2 * Omega_inv))
  g_ana <- foceiOmegaGradient(oc, om, etaModes, H_inv)

  q_sum <- sum(etaModes^2) * Omega_inv
  log_det_term <- 2 * N
  log_h_term <- sum(-4 * Omega_inv * sapply(H_inv, function(x) as.numeric(x)))
  g_closed <- -2 * q_sum + log_det_term + log_h_term

  expect_equal(unname(g_ana), unname(g_closed), tolerance = 1e-12)
})


# ---- .computeMetricDerivative -------------------------------------------

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


.make_decay_deriv2 <- function() {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  withr::local_dir(tempdir())
  f <- c(A = "-k*A")
  tag <- as.integer(Sys.time())
  m <- odemodel(f, modelname = paste0("md_decay_m_", tag),
                solver = "CppODE", deriv2 = TRUE, nStack = 4L, verbose = FALSE)
  ode_opts <- list(atol = 1e-12, rtol = 1e-12)
  xfn <- Xs(m, condition = "C1",
            optionsOde = ode_opts, optionsSens = ode_opts)
  gfn <- Y(c(y = "A"), f = f, parameters = "A",
           modelname = paste0("md_decay_obs_", tag),
           compile = TRUE, deriv2 = TRUE, attach.input = FALSE,
           condition = "C1")
  pfn <- Pexpl(c(A = "A", k = "k"), parameters = NULL,
               modelname = paste0("md_decay_p_", tag),
               compile = TRUE, deriv2 = TRUE, derivMode = "symbolic",
               condition = "C1")
  list(prd = gfn * xfn * pfn, m = m)
}


.fd_dG <- function(obj, theta, h = 1e-4) {
  K <- length(theta)
  Gp <- function(t) obj(t, deriv = TRUE, deriv2 = FALSE)$hessian / 2
  dG <- array(0, c(K, K, K), dimnames = list(names(theta), names(theta), names(theta)))
  for (c_idx in seq_len(K)) {
    tp <- theta; tp[c_idx] <- tp[c_idx] + h
    tm <- theta; tm[c_idx] <- tm[c_idx] - h
    Gp_v <- Gp(tp); Gm_v <- Gp(tm)
    dG[, , c_idx] <- (Gp_v - Gm_v) / (2 * h)
  }
  dG
}


test_that(".computeMetricDerivative matches FD on a decay model", {
  skip_if_no_compile()
  ds <- .make_decay_deriv2()
  times <- seq(0.5, 5, by = 0.5)
  pars  <- c(A = 1.0, k = 0.5)
  truth <- ds$prd(times, pars, deriv = FALSE)[[1]]
  ydata <- truth[truth[, "time"] %in% times, "y"] +
    stats::rnorm(length(times), sd = 0.05)
  data  <- datalist(C1 = data.frame(name = "y", time = times,
                                    value = ydata, sigma = 0.05))

  obj <- normL2(data, ds$prd)
  mctx <- metricContext(data, ds$prd)

  th <- c(A = 1.3, k = 0.4)
  ana <- dMod:::.computeMetricDerivative(mctx, th)
  fd  <- .fd_dG(obj, th, h = 1e-4)

  # Loose tolerance because FD is O(h^2) ~ 1e-8 plus ODE integration noise.
  expect_equal(ana, fd, tolerance = 5e-4, ignore_attr = TRUE)
})


test_that("metricContext requires a prdfn", {
  expect_error(metricContext(data.frame(), list()), "prdfn")
})


# ---- nonCenteredOmega ---------------------------------------------------

test_that("nonCenteredOmega: z -> eta -> z is the identity", {
  set.seed(1L)
  subjects <- c("s1", "s2", "s3")
  om <- omega(eta = c("eta_a", "eta_b"), structure = "full",
              subjects = subjects)
  nc <- nonCenteredOmega(om)

  omegaChol <- c(omega_a_a = log(0.5),
                 omega_b_b = log(0.3),
                 omega_b_a = 0.2)

  z0 <- setNames(rnorm(prod(dim(nc$zNames))), as.vector(nc$zNames))
  eta <- nc$mapZtoEta(z0, omegaChol)
  z1  <- nc$mapEtaToZ(eta, omegaChol)
  expect_equal(z1, z0, tolerance = 1e-12)
})


test_that("nonCenteredOmega: eta = L z element-wise check", {
  set.seed(2L)
  subjects <- c("s1", "s2")
  om <- omega(eta = c("eta_a", "eta_b"), structure = "full",
              subjects = subjects)
  nc <- nonCenteredOmega(om)
  omegaChol <- c(omega_a_a = log(1.0),
                 omega_b_b = log(2.0),
                 omega_b_a = 0.5)
  L <- om$buildL(omegaChol)

  z0 <- c(z_a_s1 = 1.0, z_b_s1 = -0.5,
          z_a_s2 = 0.3, z_b_s2 = 0.7)
  eta <- nc$mapZtoEta(z0, omegaChol)
  expected_s1 <- as.numeric(L %*% c(1.0, -0.5))
  expected_s2 <- as.numeric(L %*% c(0.3,  0.7))
  expect_equal(unname(eta[c("eta_a_s1", "eta_b_s1")]),
               expected_s1, tolerance = 1e-12)
  expect_equal(unname(eta[c("eta_a_s2", "eta_b_s2")]),
               expected_s2, tolerance = 1e-12)
})


test_that("nonCenteredOmega: standard normal z gives eta with empirical cov ~ Omega", {
  set.seed(3L)
  N <- 5000L
  subjects <- paste0("s", seq_len(N))
  om <- omega(eta = c("eta_a", "eta_b"), structure = "full",
              subjects = subjects)
  nc <- nonCenteredOmega(om)
  omegaChol <- c(omega_a_a = log(1.0),
                 omega_b_b = log(1.5),
                 omega_b_a = 0.3)
  L <- om$buildL(omegaChol)
  Omega_truth <- L %*% t(L)

  z <- setNames(rnorm(prod(dim(nc$zNames))), as.vector(nc$zNames))
  eta <- nc$mapZtoEta(z, omegaChol)
  E <- matrix(eta[as.vector(om$subjectEtas)], N, 2L,
              dimnames = list(NULL, om$eta))

  C_emp <- cov(E)
  expect_lt(abs(C_emp[1, 1] / Omega_truth[1, 1] - 1), 0.10)
  expect_lt(abs(C_emp[2, 2] / Omega_truth[2, 2] - 1), 0.10)
  expect_lt(abs(C_emp[1, 2] - Omega_truth[1, 2]),     0.10)
})
