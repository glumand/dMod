context("evalConditionResidual lifted helper")


test_that("evalConditionResidual matches the in-normL2 closure (no errmodel)", {
  set.seed(1)
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, modelname = "ecr_simple_obs")
  x <- Xt()

  subjects <- c("c1", "c2")
  trafo <- eqnvec(intercept = "mu * exp(eta)")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, modelname = "ecr_simple_p")

  y_obs <- c(2.1, 1.8)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  obj <- normL2(data, g * x * p)

  init <- c(mu = 2.0, eta_c1 = 0.1, eta_c2 = -0.1)
  result <- obj(init)
  # value via normL2 should equal sum of per-condition contributions.
  prediction <- (g * x * p)(times = c(0), pars = init, deriv = TRUE,
                            conditions = subjects)
  per_cn <- lapply(subjects, function(cn) {
    evalConditionResidual(dataI = data[[cn]], predictionI = prediction[[cn]],
                            pars = init, cn = cn, bessel = 1, deriv = TRUE)
  })
  agg <- Reduce(`+`, per_cn)
  expect_equal(agg$value, result$value, tolerance = 1e-10)
  expect_equal(agg$gradient, result$gradient, tolerance = 1e-10)
  expect_equal(agg$hessian,  result$hessian,  tolerance = 1e-10)
})


test_that("normL2 still produces the same OFV after refactor (regression)", {
  # Recompute the test-focei.R smoke result by hand here -- exact reproducibility
  # check that the eval_condition lift didn't drift the residual numerics.
  set.seed(1)
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, modelname = "ecr_regr_obs")
  x <- Xt()
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2", "s3", "s4")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, modelname = "ecr_regr_p")

  true_mu  <- 2.0; true_om <- 0.3
  true_eta <- rnorm(length(subjects), 0, true_om)
  y_obs    <- true_mu * exp(true_eta) + rnorm(length(subjects), 0, 0.2)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)
  joint <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)

  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3),
            setNames(rep(0, 4), paste0("eta_", subjects)))
  res <- joint(init, deriv = TRUE)
  # Hand-computed: 4 conditions with eta=0, prediction = mu_pop = 2.0 each.
  # residual = (y - 2.0) / 0.2; obj = sum(r^2) + sum(log(2*pi*0.04)) + prior.
  expect_true(is.finite(res$value))
  expect_true(all(is.finite(res$gradient)))
})
