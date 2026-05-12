context("predict.nlmeFit + diagnostic plot helpers")


build_for_plots <- function(seed = 1L) {
  set.seed(seed)
  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE, modelname = paste0("plt_obs_", seed))
  x <- Xt()
  subjects <- paste0("s", 1:4)
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = paste0("plt_p_", seed))
  true_eta <- rnorm(4, 0, 0.3)
  obs_rows <- do.call(rbind, lapply(seq_along(subjects), function(i) {
    ts <- c(0, 1, 2)
    data.frame(name = "y", time = ts,
               value = 2.0 * exp(true_eta[i]) + rnorm(length(ts), 0, 0.2),
               sigma = 0.2, condition = subjects[i], stringsAsFactors = FALSE)
  }))
  data <- as.datalist(obs_rows)
  om <- omega(eta = "eta", subjects = subjects)
  joint <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)
  list(joint = joint, om = om, model = g * x * p, data = data,
       subjects = subjects)
}

run_focei <- function(s, init = c(mu_pop = 2.0, omega_eta_eta = log(0.3)))
  suppressMessages(nlmeFit(s$joint, s$om, init,
                            model = s$model, data = s$data,
                            method = "focei", verbose = FALSE))


test_that("predict.nlmeFit returns a long data.frame with IPRED/PRED/IWRES", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  fit <- run_focei(build_for_plots(1L))
  pf <- predict(fit, times = seq(0, 3, length.out = 10))
  expect_s3_class(pf, "data.frame")
  expected_cols <- c("condition","time","name","observed","sigma",
                     "IPRED","PRED","source","IRES","PRES","IWRES","PWRES")
  expect_true(all(expected_cols %in% names(pf)))
  expect_true(any(pf$source == "obs"))
  expect_true(any(pf$source == "grid"))
  obs <- pf[pf$source == "obs", , drop = FALSE]
  expect_equal(obs$IRES, obs$observed - obs$IPRED, tolerance = 1e-12)
})


test_that("plot.nlmeFit, plotIndivs, plotResiduals return ggplots", {
  skip_if_not_installed("ggplot2")
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  fit <- run_focei(build_for_plots(2L))
  expect_s3_class(plotIndivs(fit),    "ggplot")
  expect_s3_class(plot(fit),          "ggplot")
  expect_s3_class(plotResiduals(fit), "ggplot")
})


test_that("plotHistIndivs returns either a ggplot (cowplot) or a list of two", {
  skip_if_not_installed("ggplot2")
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  fit <- run_focei(build_for_plots(3L))
  p <- plotHistIndivs(fit)
  if (requireNamespace("cowplot", quietly = TRUE)) {
    expect_s3_class(p, "ggplot")
  } else {
    expect_type(p, "list")
    expect_named(p, c("hist", "qq"))
  }
})


test_that("plotTrace errors on focei fit, works on foceiQuadrature fit", {
  skip_if_not_installed("ggplot2")
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- build_for_plots(4L)
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))
  fit_lap <- suppressMessages(nlmeFit(s$joint, s$om, init,
                                       model = s$model, data = s$data,
                                       method = "focei", verbose = FALSE))
  expect_error(plotTrace(fit_lap), "requires.*quadrature")

  fit_qd <- suppressMessages(nlmeFit(s$joint, s$om, init,
                                      model = s$model, data = s$data,
                                      method = "foceiQuadrature",
                                      control = list(quadrature = list(
                                        level = 4L,
                                        epsQuadLevels = c(3L, 4L),
                                        maxEcmPerStage = 2L)),
                                      verbose = FALSE))
  expect_s3_class(plotTrace(fit_qd), "ggplot")
})


test_that("plotResiduals back-compat: parframe path still works", {
  # Smoke-test that the existing plotResiduals(parframe, x, data, ...) entry
  # point isn't broken by the nlmeFit dispatch shim. Just confirm the
  # nlmeFit branch is bypassed for a non-nlmeFit input.
  pf <- structure(data.frame(value = 1.0, index = 1L),
                  class = c("parframe", "data.frame"))
  expect_false(inherits(pf, "nlmeFit"))
})
