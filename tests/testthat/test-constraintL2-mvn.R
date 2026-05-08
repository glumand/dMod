context("constraintL2 MVN path (FOCEI prior on random effects)")


test_that("backwards compatibility: existing scalar/diagonal path still works", {
  prior <- structure(rep(0, 3), names = letters[1:3])
  obj <- constraintL2(mu = prior)
  res <- obj(c(a = 1, b = -1, c = 0.5))
  expect_equal(res$value, 1 + 1 + 0.25)
})


test_that("MVN value matches direct computation for diagonal Omega", {
  set.seed(2)
  subjects <- paste0("s", 1:5)
  om <- omega(eta = c("eta_a", "eta_b"), subjects = subjects)
  obj <- constraintL2(mu = 0, Omega = om)

  eta_vals <- runif(10, -0.5, 0.5)
  names(eta_vals) <- as.vector(om$subject_etas)
  chol_vals <- c(omega_a_a = log(0.3), omega_b_b = log(0.4))
  pars <- c(eta_vals, chol_vals)

  L <- om$build_L(chol_vals)
  Omega_mat <- L %*% t(L)
  Omega_inv <- solve(Omega_mat)
  N <- length(subjects); K <- 2L
  eta_mat <- matrix(eta_vals, nrow = N, ncol = K)
  expected_quad <- sum(vapply(seq_len(N), function(i)
    drop(eta_mat[i, ] %*% Omega_inv %*% eta_mat[i, ]), 0.0))
  expected_logdet <- N * log(det(Omega_mat))
  expected_val <- expected_quad + expected_logdet

  expect_equal(obj(pars)$value, expected_val, tolerance = 1e-10)
})


test_that("MVN gradient agrees with numDeriv::grad (diagonal Omega)", {
  skip_if_not_installed("numDeriv")
  set.seed(3)

  subjects <- paste0("s", 1:4)
  om <- omega(eta = c("eta_a", "eta_b"), subjects = subjects)
  obj <- constraintL2(mu = 0, Omega = om)

  eta_vals <- rnorm(8, sd = 0.3)
  names(eta_vals) <- as.vector(om$subject_etas)
  chol_vals <- c(omega_a_a = log(0.5), omega_b_b = log(0.4))
  pars <- c(eta_vals, chol_vals)

  analytic <- obj(pars)$gradient
  numeric_grad  <- numDeriv::grad(function(x) obj(x)$value, pars)
  names(numeric_grad) <- names(pars)
  expect_equal(analytic[names(numeric_grad)], numeric_grad, tolerance = 1e-6)
})


test_that("MVN gradient agrees with numDeriv::grad (full Omega with correlations)", {
  skip_if_not_installed("numDeriv")
  set.seed(4)

  subjects <- paste0("s", 1:6)
  om <- omega(eta = c("eta_Cl", "eta_V", "eta_Ka"),
              structure = "full", subjects = subjects)
  obj <- constraintL2(mu = 0, Omega = om)

  eta_vals <- rnorm(18, sd = 0.4)
  names(eta_vals) <- as.vector(om$subject_etas)
  chol_vals <- c(omega_Cl_Cl = log(0.3), omega_V_V = log(0.2), omega_Ka_Ka = log(0.5),
                 omega_V_Cl  = 0.05,     omega_Ka_Cl = -0.02,  omega_Ka_V  = 0.07)
  pars <- c(eta_vals, chol_vals)

  analytic <- obj(pars)$gradient
  numeric_grad <- numDeriv::grad(function(x) obj(x)$value, pars)
  names(numeric_grad) <- names(pars)
  expect_equal(analytic[names(numeric_grad)], numeric_grad, tolerance = 1e-5)
})


test_that("MVN gradient agrees with numDeriv::grad (selective correlation)", {
  skip_if_not_installed("numDeriv")
  set.seed(5)

  subjects <- paste0("s", 1:5)
  om <- omega(eta = c("eta_Cl", "eta_V", "eta_Ka"),
              correlate = list(c("eta_Cl", "eta_V")),
              subjects = subjects)
  obj <- constraintL2(mu = 0, Omega = om)

  eta_vals <- rnorm(15, sd = 0.3)
  names(eta_vals) <- as.vector(om$subject_etas)
  chol_vals <- c(omega_Cl_Cl = log(0.3), omega_V_V = log(0.2), omega_Ka_Ka = log(0.4),
                 omega_V_Cl  = 0.04)
  pars <- c(eta_vals, chol_vals)

  analytic <- obj(pars)$gradient
  numeric_grad <- numDeriv::grad(function(x) obj(x)$value, pars)
  names(numeric_grad) <- names(pars)
  expect_equal(analytic[names(numeric_grad)], numeric_grad, tolerance = 1e-5)
})


test_that("MVN Hessian agrees with numDeriv::hessian on quadratic-only term", {
  # The Hessian of N*log|Omega| w.r.t. log-Cholesky parameters is zero
  # (linear in omega_kk). The Hessian of the quadratic term equals the
  # Gauss-Newton outer product of dz/dp. So full numerical Hessian should
  # match the analytic Hessian when residuals z_i are small (linearization
  # is exact at the mode). Here we test at zero etas, where the GN H equals
  # the true Hessian.
  skip_if_not_installed("numDeriv")
  set.seed(6)

  subjects <- paste0("s", 1:3)
  om <- omega(eta = c("eta_a", "eta_b"), structure = "full", subjects = subjects)
  obj <- constraintL2(mu = 0, Omega = om)

  eta_vals <- rep(0, 6)
  names(eta_vals) <- as.vector(om$subject_etas)
  chol_vals <- c(omega_a_a = log(0.3), omega_b_b = log(0.4),
                 omega_b_a = 0.05)
  pars <- c(eta_vals, chol_vals)

  analytic_h <- obj(pars)$hessian
  numeric_h  <- numDeriv::hessian(function(x) obj(x)$value, pars)
  rownames(numeric_h) <- colnames(numeric_h) <- names(pars)

  # eta-eta block should match exactly: it's the GN-exact 2*Omega^{-1} block diag
  eta_idx <- 1:6
  expect_equal(analytic_h[eta_idx, eta_idx],
               numeric_h[eta_idx, eta_idx], tolerance = 1e-4)
})


test_that("MVN value uses mu correctly when mu != 0", {
  set.seed(7)
  subjects <- paste0("s", 1:3)
  om <- omega(eta = c("eta_a", "eta_b"), subjects = subjects)
  obj <- constraintL2(mu = c(eta_a = 0.1, eta_b = -0.05), Omega = om)

  eta_vals <- c(eta_a_s1 = 0.1, eta_b_s1 = -0.05,
                eta_a_s2 = 0.1, eta_b_s2 = -0.05,
                eta_a_s3 = 0.1, eta_b_s3 = -0.05)
  chol_vals <- c(omega_a_a = log(0.3), omega_b_b = log(0.4))
  pars <- c(eta_vals, chol_vals)

  # All etas exactly at their priors -> quadratic part is zero
  result <- obj(pars)
  expected_val <- 3 * log(det(om$build_L(chol_vals) %*% t(om$build_L(chol_vals))))
  expect_equal(result$value, expected_val, tolerance = 1e-10)
})


test_that("MVN errors helpfully when subject expansion is missing", {
  om <- omega(eta = c("eta_a", "eta_b"))   # no subjects
  expect_error(constraintL2(mu = 0, Omega = om), "subject expansion")
})


test_that("MVN summable with normL2-style objfn via +.objfn", {
  # Sanity check: sum two constraintL2 objfns (one diagonal-style on theta,
  # one MVN on etas). +.objfn must merge their parameter sets cleanly.
  subjects <- paste0("s", 1:3)
  om <- obj_mvn <- NULL
  om <- omega(eta = c("eta_a"), subjects = subjects)

  obj_theta <- constraintL2(mu = c(theta_a = 0))
  obj_mvn   <- constraintL2(mu = 0, Omega = om)
  obj_sum   <- obj_theta + obj_mvn

  pars <- c(theta_a = 1.5,
            eta_a_s1 = 0.1, eta_a_s2 = -0.2, eta_a_s3 = 0.3,
            omega_a_a = log(0.4))

  vsum  <- obj_sum(pars)$value
  vsep  <- obj_theta(pars)$value + obj_mvn(pars)$value
  expect_equal(vsum, vsep, tolerance = 1e-10)
})
