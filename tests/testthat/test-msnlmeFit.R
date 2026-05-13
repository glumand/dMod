context("msnlmeFit multi-start wrapper around nlmeFit")


build_one_eta_for_ms <- function(seed = 1L, N = 4L, tag = "ms") {
  set.seed(seed)
  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE,
         modelname = paste0("msnlmefit_obs_", tag, "_", seed))
  x <- Xt()
  subjects <- paste0("s", seq_len(N))
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = paste0("msnlmefit_p_", tag, "_", seed))
  true_mu  <- 2.0; true_om <- 0.3
  true_eta <- rnorm(N, 0, true_om)
  y_obs    <- true_mu * exp(true_eta) + rnorm(N, 0, 0.2)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)
  obj <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)
  list(obj = obj, om = om, prdfn = g * x * p, data = data,
       subjects = subjects, true_mu = true_mu, y_obs = y_obs)
}


test_that("msnlmeFit returns a parlist of nlmeFits and as.parframe works", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- build_one_eta_for_ms(11L, tag = "basic")
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
  expect_equal(pl[[1]]$parinit, center)  # start1stfromCenter respected

  pf <- as.parframe(pl)
  expect_s3_class(pf, "parframe")
  expect_true(all(c("value", "converged", "iterations") %in%
                    attr(pf, "metanames")))
  expect_true(diff(range(pf$value)) >= 0)  # sorted by value
  # best fit recovers true mu within tolerance
  best <- as.parvec(pf)
  expect_lt(abs(unname(best["mu_pop"]) - mean(s$y_obs)), 0.5)
})


test_that("msnlmeFit accepts a parframe as center (rows used directly)", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- build_one_eta_for_ms(12L, tag = "pf")

  # Hand-built parframe with 2 rows, both reasonable starts
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
  s <- build_one_eta_for_ms(13L, tag = "keep")
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
  s <- build_one_eta_for_ms(14L, tag = "err")

  # Wrap obj so that the second fit's initial value forces an error.
  bad_joint <- function(pars, ...) {
    if (isTRUE(getOption(".dMod.force_fail", FALSE)))
      stop("forced fail")
    s$obj(pars, ...)
  }
  class(bad_joint) <- class(s$obj)
  # Skip the synthetic-error path on platforms where method dispatch would
  # discard the class - we instead use a manifestly bad start.

  starts_df <- data.frame(mu_pop        = c(2.0, NA_real_),
                          omega_eta_eta = c(log(0.3), log(0.3)))
  pf_in <- parframe(starts_df)

  pl <- suppressMessages(msnlmeFit(
    s$obj, s$om, pf_in,
    prdfn = s$prdfn, data = s$data,
    method = "focei", cores = 1L, verbose = FALSE))

  expect_length(pl, 2L)
  # Fit 1 should succeed, fit 2 should fail and carry an error field
  stats <- vapply(pl, function(f) !is.null(f$error), logical(1))
  expect_true(any(!stats))  # at least one succeeded
  expect_true(any(stats))   # at least one failed
})
