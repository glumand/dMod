context("focei() FOCEI wrapper end-to-end smoke")


test_that("focei runs end-to-end on a minimal one-eta NLME model", {
  set.seed(1)

  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  # Observation y = intercept (no states), composed with Xt for the time axis
  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, modelname = "focei_smoke_obs")
  x <- Xt()

  # Per-subject parameter trafos via branch + insert
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2", "s3", "s4")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE,
         modelname = "focei_smoke_p")

  # Synthetic data: one observation per subject at time = 0
  true_mu  <- 2.0
  true_om  <- 0.3
  true_eta <- rnorm(length(subjects), 0, true_om)
  y_obs    <- true_mu * exp(true_eta) + rnorm(length(subjects), 0, 0.2)

  data <- data.frame(name = "y", time = 0, sigma = 0.2,
                     value = y_obs, condition = subjects,
                     stringsAsFactors = FALSE)
  data <- as.datalist(data)

  # Random-effect spec
  om <- omega(eta = "eta", subjects = subjects)

  # Joint = data likelihood + MVN prior on etas
  joint <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)

  focei_obj <- focei(joint, om, innerControl = list(rtol = 1e-7, maxit = 50))

  # focei's outer parameters
  outer_init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  # Evaluate once: structure checks ----------------------------------------
  res <- focei_obj(outer_init)
  expect_true(is.finite(res$value))
  expect_true(all(is.finite(res$gradient)))
  expect_true(all(is.finite(res$hessian)))
  expect_setequal(names(res$gradient), names(outer_init))
  expect_equal(dim(res$hessian),
               c(length(outer_init), length(outer_init)))

  diag <- attr(res, "focei_diag")
  expect_equal(dim(diag$eta_star), c(length(subjects), 1L))
  expect_true(all(diag$converged))
  expect_true(all(diag$iter <= 50L))
  expect_true(all(diag$logdet > -Inf & diag$logdet < Inf))
  expect_true(all(diag$n_floored == 0L))

  # Inner-solve correctness: at returned eta_i*, the joint gradient wrt eta_i
  # must vanish (within tolerance) ------------------------------------------
  full_pars <- c(outer_init, setNames(as.numeric(diag$eta_star),
                                      as.vector(om$subject_etas)))
  joint_full <- joint(full_pars, deriv = TRUE)
  for (nm in om$subject_etas) {
    expect_lt(abs(joint_full$gradient[nm]), 1e-3)
  }

  # Outer trust step converges to a sensible point --------------------------
  fit <- suppressMessages(trust(focei_obj, outer_init,
                                rinit = 1, rmax = 10, iterlim = 50,
                                fterm = 1e-7, mterm = 1e-7))
  expect_true(fit$converged)
  expect_true(is.finite(fit$value))
  # mu_pop should land somewhere near the data mean (loose check)
  expect_true(abs(fit$argument["mu_pop"] - mean(y_obs)) < 0.5)
})


test_that("trust on_step hook fires on every step decision", {
  qobj <- function(x, ...) {
    nm <- names(x) %||% paste0("p", seq_along(x))
    names(x) <- nm
    list(value    = sum(x^2),
         gradient = setNames(2 * x, nm),
         hessian  = matrix(2 * diag(length(x)),
                           length(x), length(x),
                           dimnames = list(nm, nm)))
  }
  events <- list()
  recorder <- function(rho, accepted, iter, r) {
    events[[length(events) + 1L]] <<- list(rho = rho, accepted = accepted,
                                           iter = iter, r = r)
  }

  fit <- trust(qobj, parinit = c(p1 = 5, p2 = 5), rinit = 1, rmax = 10,
               iterlim = 30, fterm = 1e-10, mterm = 1e-10,
               on_step = recorder)

  expect_gt(length(events), 0L)
  expect_true(all(vapply(events, function(e) is.logical(e$accepted), TRUE)))
  # Every accepted non-terminating step must have rho >= 1/4 by construction
  # (terminating steps short-circuit the rho check and may carry NaN/-Inf rho).
  rhos <- vapply(events, function(e) e$rho, 0.0)
  acc  <- vapply(events, function(e) isTRUE(e$accepted), TRUE)
  finite_acc <- acc & is.finite(rhos)
  if (any(finite_acc)) expect_true(all(rhos[finite_acc] >= 1/4 - 1e-12))
})


test_that("focei correction = 'lagged' matches FD on the OFV within ~10%", {
  set.seed(1)

  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE, modelname = "focei_corr_obs")
  x <- Xt()
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2", "s3", "s4")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = "focei_corr_p")

  true_mu  <- 2.0; true_om  <- 0.3
  true_eta <- rnorm(length(subjects), 0, true_om)
  y_obs    <- true_mu * exp(true_eta) + rnorm(length(subjects), 0, 0.2)
  data <- data.frame(name = "y", time = 0, sigma = 0.2,
                     value = y_obs, condition = subjects,
                     stringsAsFactors = FALSE)
  data <- as.datalist(data)
  om <- omega(eta = "eta", subjects = subjects)

  joint <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)
  focei_lag  <- focei(joint, om, model = g * x * p, data = data,
                      correction = "lagged",
                      innerControl = list(rtol = 1e-8, maxit = 50))
  focei_eag  <- focei(joint, om, model = g * x * p, data = data,
                      correction = "eager",
                      innerControl = list(rtol = 1e-8, maxit = 50))
  focei_none <- focei(joint, om,
                      innerControl = list(rtol = 1e-8, maxit = 50))

  outer_init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  res_none <- focei_none(outer_init)
  res_lag  <- focei_lag (outer_init)
  res_eag  <- focei_eag (outer_init)
  diag_lag <- attr(res_lag, "focei_diag")
  diag_eag <- attr(res_eag, "focei_diag")

  # Eager and lagged refresh on cold start, so first-call gradients agree
  expect_equal(res_lag$gradient, res_eag$gradient, tolerance = 1e-12)
  expect_equal(diag_lag$refresh_count, 1L)
  expect_equal(diag_eag$refresh_count, 1L)
  expect_true(!is.null(diag_lag$correction_value))

  # FD ground truth via correction = "none" focei (re-solves inner each call)
  fd_grad <- numeric(length(outer_init))
  h <- 1e-5
  for (k in seq_along(outer_init)) {
    pp <- outer_init; pp[k] <- pp[k] + h
    pm <- outer_init; pm[k] <- pm[k] - h
    fd_grad[k] <- (focei_none(pp, deriv = FALSE)$value -
                   focei_none(pm, deriv = FALSE)$value) / (2 * h)
  }
  fd_correction <- fd_grad - res_none$gradient

  expect_equal(as.numeric(diag_lag$correction_value),
               as.numeric(fd_correction), tolerance = 0.1)
  # And the corrected gradient should be much closer to FD total than the
  # uncorrected envelope-only one.
  err_corr   <- max(abs(as.numeric(res_lag$gradient)  - fd_grad))
  err_uncorr <- max(abs(as.numeric(res_none$gradient) - fd_grad))
  expect_lt(err_corr, 0.3 * err_uncorr)
})


test_that("focei lagged cache holds across small steps and refreshes on big ones", {
  set.seed(2)

  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE, modelname = "focei_lag5_obs")
  x <- Xt()
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2", "s3")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = "focei_lag5_p")

  y_obs <- 2.0 * exp(rnorm(length(subjects), 0, 0.3)) +
           rnorm(length(subjects), 0, 0.2)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)
  joint <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)

  # Big tau so the anchor-distance trigger does not fire on small steps;
  # M = 100 to disable periodic refresh for this test.
  fobj <- focei(joint, om, model = g * x * p, data = data,
                correction = "lagged",
                correctionControl = list(tau = 1.0, M = 100L),
                innerControl = list(rtol = 1e-8, maxit = 50))

  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  r0 <- fobj(init)               # cold start -> refresh
  expect_equal(attr(r0, "focei_diag")$refresh_count, 1L)

  init_small <- init + 0.001     # tiny step
  r1 <- fobj(init_small)
  expect_equal(attr(r1, "focei_diag")$refresh_count, 1L)  # NO refresh

  init_big <- c(mu_pop = 2.0 + 5, omega_eta_eta = log(0.3))  # |dmu|/|2|>1
  r2 <- fobj(init_big)
  expect_equal(attr(r2, "focei_diag")$refresh_count, 2L)  # refresh
})


test_that("focei lagged correction errors clearly when model/data are missing", {
  om <- omega(eta = "eta", subjects = c("a", "b"))
  joint <- structure(function(...) NULL, class = c("objfn", "fn"))
  attr(joint, "parameters") <- c("mu_pop", "omega_eta_eta", "eta_a", "eta_b")
  attr(joint, "conditions") <- c("a", "b")
  expect_error(focei(joint, om, correction = "lagged"),
               "model.*required")
  expect_error(focei(joint, om, model = identity, correction = "lagged"),
               "data.*required")
})


test_that("focei attr 'on_step' is the closure that records rejected steps", {
  set.seed(3)

  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE, modelname = "focei_on_step_obs")
  x <- Xt()
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = "focei_on_step_p")

  y_obs <- c(2.1, 1.9)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)
  joint <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)

  fobj <- focei(joint, om, model = g * x * p, data = data,
                correction = "lagged",
                innerControl = list(rtol = 1e-8, maxit = 50))

  expect_true(is.function(attr(fobj, "on_step")))
  expect_equal(attr(fobj, "correction"), "lagged")

  # Manually wire the closure: simulate one accepted + one rejected step.
  hk <- attr(fobj, "on_step")
  hk(rho = 0.9, accepted = TRUE,  iter = 1L, r = 1)
  diag1 <- attr(fobj(c(mu_pop = 2.0, omega_eta_eta = log(0.3))), "focei_diag")
  expect_equal(diag1$last_rho, 0.9)

  hk(rho = -0.1, accepted = FALSE, iter = 2L, r = 0.25)
  # Next call should refresh due to last_step_rejected = TRUE
  init2 <- c(mu_pop = 2.0001, omega_eta_eta = log(0.3))
  diag2 <- attr(fobj(init2), "focei_diag")
  expect_equal(diag2$last_rho, -0.1)
  expect_gt(diag2$refresh_count, diag1$refresh_count)
})
