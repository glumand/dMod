context("nlmeFit(method='focei') orchestrator end-to-end")


test_that("nlmeFit(method='focei') runs on a minimal one-eta NLME model", {
  set.seed(1)

  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE, modelname = "focei_smoke_obs")
  x <- Xt()

  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2", "s3", "s4")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = "focei_smoke_p")

  true_mu  <- 2.0
  true_om  <- 0.3
  true_eta <- rnorm(length(subjects), 0, true_om)
  y_obs    <- true_mu * exp(true_eta) + rnorm(length(subjects), 0, 0.2)

  data <- as.datalist(data.frame(
    name = "y", time = 0, sigma = 0.2,
    value = y_obs, condition = subjects,
    stringsAsFactors = FALSE))

  om <- omega(eta = "eta", subjects = subjects)
  joint <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)

  outer_init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(nlmeFit(
    joint, om, outer_init,
    model = g * x * p, data = data,
    method = "focei",
    control = list(focei = list(
      trustControl = list(rinit = 1, rmax = 10, iterlim = 50,
                          fterm = 1e-7, mterm = 1e-7)))))
  expect_s3_class(fit, "nlmeFit")
  expect_equal(fit$method, "focei")
  expect_true(fit$converged)
  expect_true(is.finite(fit$value))
  expect_true(abs(fit$argument["mu_pop"] - mean(y_obs)) < 0.5)
  expect_true(!is.null(fit$omega))
  expect_equal(dim(fit$omega), c(1L, 1L))
  expect_equal(dim(fit$etaModes), c(length(subjects), 1L))
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
  rhos <- vapply(events, function(e) e$rho, 0.0)
  acc  <- vapply(events, function(e) isTRUE(e$accepted), TRUE)
  finite_acc <- acc & is.finite(rhos)
  if (any(finite_acc)) expect_true(all(rhos[finite_acc] >= 1/4 - 1e-12))
})


test_that("nlmeFit() rejects unknown method via match.arg", {
  expect_error(nlmeFit(joint = NULL, omegaSpec = NULL, init = c(p = 1),
                       method = "doesNotExist"),
               "'arg' should be one of")
})
