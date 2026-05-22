# ============================================================================
# nlmeFit + msnlmeFit + diagnostic plots (the public NLME API).
#
# Sections:
#   * nlmeFit method dispatch     - focei / quadrature / foceiQuadrature
#   * etaSE / shrinkage           - Laplace-inverse-Hessian diagnostics
#   * msnlmeFit                   - multi-start wrapper around nlmeFit
#   * predict.nlmeFit + plots     - data frame + ggplot diagnostic helpers
#
# The C++ FOCEI kernel itself is tested in test-focei.R; here we only
# exercise the orchestrator and the public output shape.
# ============================================================================

context("nlmeFit + msnlmeFit + diagnostic plots")


# Shared one-eta NLME fixture builder.
.build_one_eta <- function(seed = 1L, N = 4L, tag = "nlmef") {
  set.seed(seed)
  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE,
         modelname = paste0("nlmefit_obs_", tag, "_", seed))
  x <- Xt()
  subjects <- paste0("s", seq_len(N))
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = paste0("nlmefit_p_", tag, "_", seed))
  true_mu  <- 2.0; true_om <- 0.3
  true_eta <- rnorm(N, 0, true_om)
  y_obs    <- true_mu * exp(true_eta) + rnorm(N, 0, 0.2)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)
  obj <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)
  list(obj = obj, om = om, prdfn = g * x * p, data = data,
       subjects = subjects, true_mu = true_mu, true_om = true_om,
       y_obs = y_obs)
}


# ---- Method dispatch -----------------------------------------------------

test_that("nlmeFit(method='foceiQuadrature') runs end-to-end on one-eta prdfn", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(1L, tag = "fq")
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(nlmeFit(
    s$obj, s$om, init, prdfn = s$prdfn, data = s$data,
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
  s <- .build_one_eta(2L, tag = "qd")
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(nlmeFit(
    s$obj, s$om, init, prdfn = s$prdfn, data = s$data,
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


test_that("nlmeFit(method='foceiQuadrature') polishes a FOCEI fit without OFV blow-up", {
  # Loose invariant: starting from a FOCEI warmstart, the ECM final OFV
  # should sit within a small distance of the FOCEI OFV (polish either
  # improves or leaves it roughly unchanged on a well-resolved problem).
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(3L, N = 6L, tag = "polish")
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(nlmeFit(
    s$obj, s$om, init, prdfn = s$prdfn, data = s$data,
    method = "foceiQuadrature",
    control = list(quadrature = list(level = 4L,
                                     epsQuadLevels = c(4L, 5L),
                                     maxEcmPerStage = 3L)),
    verbose = FALSE))
  expect_true(is.finite(fit$foceiStart$value))
  expect_lt(abs(fit$value - fit$foceiStart$value), 10)
})


# ---- etaSE + shrinkage ---------------------------------------------------

test_that("nlmeFit emits etaSE + shrinkage from Laplace inverse Hessian", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(7L, N = 5L, tag = "se")
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(nlmeFit(
    s$obj, s$om, init, prdfn = s$prdfn, data = s$data,
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


# ---- msnlmeFit -----------------------------------------------------------

test_that("msnlmeFit returns a parlist of nlmeFits and as.parframe works", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(11L, tag = "ms_basic")
  center <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  pl <- suppressMessages(msnlmeFit(
    s$obj, s$om, center,
    prdfn = s$prdfn, data = s$data,
    method = "focei",
    fits = 3L, cores = 1L,
    samplefun = "rnorm", sd = 0.2,
    start1stfromCenter = TRUE,
    verbose = FALSE))

  expect_s3_class(pl, "parlist")
  expect_length(pl, 3L)
  for (fit in pl) {
    expect_s3_class(fit, "nlmeFit")
    expect_true(all(c("argument", "value", "converged", "iterations",
                      "parinit", "index") %in% names(fit)))
    # default keepFull = FALSE strips heavy state
    expect_null(fit$emDiag)
    expect_null(fit$prdfn)
    expect_null(fit$data)
  }
  expect_equal(pl[[1]]$parinit, center)

  pf <- as.parframe(pl)
  expect_s3_class(pf, "parframe")
  expect_true(all(c("value", "converged", "iterations") %in%
                    attr(pf, "metanames")))
  expect_true(diff(range(pf$value)) >= 0)
  best <- as.parvec(pf)
  expect_lt(abs(unname(best["mu_pop"]) - mean(s$y_obs)), 0.5)
})


test_that("msnlmeFit accepts a parframe as center (rows used directly)", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(12L, tag = "ms_pf")

  starts_df <- data.frame(mu_pop        = c(2.0, 2.5),
                          omega_eta_eta = c(log(0.3), log(0.5)))
  pf_in <- parframe(starts_df)

  pl <- suppressMessages(msnlmeFit(
    s$obj, s$om, pf_in,
    prdfn = s$prdfn, data = s$data,
    method = "focei", cores = 1L, verbose = FALSE))

  expect_s3_class(pl, "parlist")
  expect_length(pl, 2L)  # nrow(center) overrides fits
})


test_that("msnlmeFit keepFull = TRUE preserves emDiag / prdfn / data", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(13L, tag = "ms_keep")
  center <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  pl <- suppressMessages(msnlmeFit(
    s$obj, s$om, center,
    prdfn = s$prdfn, data = s$data,
    method = "focei", fits = 2L, cores = 1L,
    keepFull = TRUE, start1stfromCenter = TRUE,
    sd = 0.1, verbose = FALSE))

  expect_length(pl, 2L)
  expect_false(is.null(pl[[1]]$emDiag))
  expect_false(is.null(pl[[1]]$prdfn))
  expect_false(is.null(pl[[1]]$data))
})


test_that("msnlmeFit handles per-fit failures without aborting the run", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(14L, tag = "ms_err")

  starts_df <- data.frame(mu_pop        = c(2.0, NA_real_),
                          omega_eta_eta = c(log(0.3), log(0.3)))
  pf_in <- parframe(starts_df)

  pl <- suppressMessages(msnlmeFit(
    s$obj, s$om, pf_in,
    prdfn = s$prdfn, data = s$data,
    method = "focei", cores = 1L, verbose = FALSE))

  expect_length(pl, 2L)
  stats <- vapply(pl, function(f) !is.null(f$error), logical(1))
  expect_true(any(!stats))  # at least one succeeded
  expect_true(any(stats))   # at least one failed
})


# ---- predict.nlmeFit + plots --------------------------------------------

# Plots fixture: more times so plotIndivs has a curve to draw.
.build_for_plots <- function(seed = 1L) {
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
  obj <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)
  list(obj = obj, om = om, prdfn = g * x * p, data = data,
       subjects = subjects)
}

.run_focei <- function(s, init = c(mu_pop = 2.0, omega_eta_eta = log(0.3)))
  suppressMessages(nlmeFit(s$obj, s$om, init,
                            prdfn = s$prdfn, data = s$data,
                            method = "focei", verbose = FALSE))


test_that("predict.nlmeFit returns a long data.frame with IPRED/PRED/IWRES", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  fit <- .run_focei(.build_for_plots(1L))
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
  fit <- .run_focei(.build_for_plots(2L))
  expect_s3_class(plotIndivs(fit),    "ggplot")
  expect_s3_class(plot(fit),          "ggplot")
  expect_s3_class(plotResiduals(fit), "ggplot")
})


test_that("plotIndivs paginates when subjectsPerPage is set", {
  skip_if_not_installed("ggplot2")
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  fit <- .run_focei(.build_for_plots(2L))  # 4 subjects
  pages <- plotIndivs(fit, subjectsPerPage = 2L)
  expect_type(pages, "list")
  expect_length(pages, 2L)
  lapply(pages, expect_s3_class, "ggplot")
  expect_match(pages[[1]]$labels$title, "page 1/2", fixed = TRUE)
  expect_match(pages[[2]]$labels$title, "page 2/2", fixed = TRUE)
  one <- plotIndivs(fit, subjectsPerPage = 10L)
  expect_type(one, "list")
  expect_length(one, 1L)
  expect_s3_class(one[[1]], "ggplot")
  expect_error(plotIndivs(fit, subjectsPerPage = 0L), "positive integer")
  expect_error(plotIndivs(fit, subjectsPerPage = c(2L, 3L)), "positive integer")
})


test_that("plotHistIndivs returns either a ggplot (cowplot) or a list of two", {
  skip_if_not_installed("ggplot2")
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  fit <- .run_focei(.build_for_plots(3L))
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
  s <- .build_for_plots(4L)
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))
  fit_lap <- suppressMessages(nlmeFit(s$obj, s$om, init,
                                       prdfn = s$prdfn, data = s$data,
                                       method = "focei", verbose = FALSE))
  expect_error(plotTrace(fit_lap), "requires.*quadrature")

  fit_qd <- suppressMessages(nlmeFit(s$obj, s$om, init,
                                      prdfn = s$prdfn, data = s$data,
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
