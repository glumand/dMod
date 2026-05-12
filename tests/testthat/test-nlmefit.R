context("nlmeFit orchestrator (focei, quadrature, foceiQuadrature)")


build_one_eta_for_ecm <- function(seed = 1L, N = 4L) {
  set.seed(seed)
  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE, modelname = paste0("nlmefit_obs_", seed))
  x <- Xt()
  subjects <- paste0("s", seq_len(N))
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = paste0("nlmefit_p_", seed))
  true_mu  <- 2.0; true_om <- 0.3
  true_eta <- rnorm(N, 0, true_om)
  y_obs    <- true_mu * exp(true_eta) + rnorm(N, 0, 0.2)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)
  joint <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)
  list(joint = joint, om = om, model = g * x * p, data = data,
       subjects = subjects, true_mu = true_mu, true_om = true_om,
       y_obs = y_obs)
}


test_that("nlmeFit(method='foceiQuadrature') runs end-to-end on one-eta model", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- build_one_eta_for_ecm(1L)
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(nlmeFit(
    s$joint, s$om, init, model = s$model, data = s$data,
    method = "foceiQuadrature",
    control = list(quadrature = list(level = 4L,
                                     epsQuadLevels = c(3L, 4L),
                                     maxEcmPerStage = 3L)),
    verbose = FALSE))

  expect_s3_class(fit, "nlmeFit")
  expect_equal(fit$method, "foceiQuadrature")
  expect_true(!is.null(fit$foceiStart))
  expect_true(!is.null(fit$stageTrace))
  expect_true(nrow(fit$stageTrace) >= 2L)
  expect_true(is.finite(fit$value))
  expect_true(abs(fit$argument["mu_pop"] - mean(s$y_obs)) < 0.5)
})


test_that("nlmeFit(method='quadrature') runs cold without a FOCEI prelude", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- build_one_eta_for_ecm(2L)
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(nlmeFit(
    s$joint, s$om, init, model = s$model, data = s$data,
    method = "quadrature",
    control = list(quadrature = list(level = 4L,
                                     epsQuadLevels = c(3L, 4L),
                                     maxEcmPerStage = 3L)),
    verbose = FALSE))
  expect_s3_class(fit, "nlmeFit")
  expect_equal(fit$method, "quadrature")
  expect_null(fit$foceiStart)
  expect_true(is.finite(fit$value))
})


test_that("nlmeFit emits etaSE + shrinkage from Laplace inverse Hessian", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- build_one_eta_for_ecm(7L, N = 5L)
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(nlmeFit(
    s$joint, s$om, init, model = s$model, data = s$data,
    method = "focei", verbose = FALSE))

  expect_s3_class(fit, "nlmeFit")
  expect_true(!is.null(fit$etaSE))
  expect_equal(dim(fit$etaSE), c(length(s$subjects), 1L))
  expect_true(all(is.finite(fit$etaSE)))
  expect_true(all(fit$etaSE > 0))

  expect_true(!is.null(fit$shrinkage))
  expect_equal(dim(fit$shrinkage), c(length(s$subjects), 1L))
  expect_true(all(fit$shrinkage <= 1 + 1e-8))
  expect_true(all(fit$shrinkage > -1))

  printed <- capture.output(print(fit))
  expect_true(any(grepl("eta", printed, ignore.case = TRUE)))
  expect_true(any(grepl("shrink", printed, ignore.case = TRUE)))
})


test_that("nlmeFit(method='foceiQuadrature') polishes a FOCEI fit without OFV blow-up", {
  # Looser invariant: starting from a FOCEI warmstart, the ECM final OFV
  # should sit within a small distance of the FOCEI OFV (polish either
  # improves or leaves it roughly unchanged on a well-resolved problem).
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- build_one_eta_for_ecm(3L, N = 6L)
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(nlmeFit(
    s$joint, s$om, init, model = s$model, data = s$data,
    method = "foceiQuadrature",
    control = list(quadrature = list(level = 4L,
                                     epsQuadLevels = c(4L, 5L),
                                     maxEcmPerStage = 3L)),
    verbose = FALSE))
  expect_true(is.finite(fit$foceiStart$value))
  expect_lt(abs(fit$value - fit$foceiStart$value), 10)
})
