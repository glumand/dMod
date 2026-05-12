context("emObjfn(method='quadrature') + rebuildQuadrature")

# NOTE: The quadrature gradient test below uses numDeriv as an independent
# check on the analytic quadrature-gradient formula. A closed-form
# reference would need a hand-derived toy where Gauss-Hermite is exact
# (linear-in-theta log-likelihood), which would mask the adaptive
# Smolyak code path being exercised; numDeriv stays.


build_one_eta_setup <- function(seed = 1L) {
  set.seed(seed)
  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, modelname = paste0("emqd_obs_", seed))
  x <- Xt()
  subjects <- c("s1", "s2", "s3")
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE,
         modelname = paste0("emqd_p_", seed))
  y_obs <- 2.0 * exp(rnorm(length(subjects), 0, 0.3)) +
           rnorm(length(subjects), 0, 0.2)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)
  joint <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)
  list(g = g, x = x, p = p, data = data, om = om, joint = joint,
       model = g * x * p, subjects = subjects)
}


test_that("emObjfn quadrature constructs and exposes rebuildQuadrature", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- build_one_eta_setup(1L)

  em <- emObjfn(s$joint, s$om,
                model = s$model, data = s$data,
                method = "quadrature", control = list(level = 4L))
  expect_s3_class(em, "emObjfn")
  expect_equal(attr(em, "method"), "quadrature")
  expect_true(is.function(attr(em, "rebuildQuadrature")))
  # Calling em without rebuildQuadrature errors clearly.
  init <- c(mu_pop = 2.0)
  expect_error(em(init), "no quadrature grid built")
})


test_that("rebuildQuadrature populates nodes, modes, and frozen state", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- build_one_eta_setup(2L)

  em <- emObjfn(s$joint, s$om,
                model = s$model, data = s$data,
                method = "quadrature", control = list(level = 4L))
  psi <- c(mu_pop = 2.0, setNames(log(0.3), s$om$cholPars))
  refresh <- attr(em, "rebuildQuadrature")(psi)
  expect_equal(dim(refresh$etaModes), c(length(s$subjects), 1L))
  expect_length(refresh$HiList, length(s$subjects))
  expect_equal(refresh$level, 4L)
  expect_true(all(refresh$converged))

  init <- c(mu_pop = 2.0)
  out <- em(init)
  expect_true(is.finite(out$value))
  expect_true(is.finite(out$gradient["mu_pop"]))
  expect_equal(attr(out, "emDiag")$method, "quadrature")
})


test_that("frozen-node invariance: em(pars1) and em(pars2) reuse the same nodes", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- build_one_eta_setup(3L)

  em <- emObjfn(s$joint, s$om,
                model = s$model, data = s$data,
                method = "quadrature", control = list(level = 4L))
  psi <- c(mu_pop = 2.0, setNames(log(0.3), s$om$cholPars))
  attr(em, "rebuildQuadrature")(psi)
  out_a <- em(c(mu_pop = 2.0))
  modes_a <- attr(out_a, "emDiag")$etaModes
  out_b <- em(c(mu_pop = 2.5))
  modes_b <- attr(out_b, "emDiag")$etaModes
  # The frozen modes (used for node generation) must NOT have changed.
  expect_identical(modes_a, modes_b)
  # But the OFV should differ because the structural param moved.
  expect_false(isTRUE(all.equal(out_a$value, out_b$value)))
})


test_that("quadrature gradient matches numDeriv on the one-eta model", {
  skip_if_not_installed("numDeriv")
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- build_one_eta_setup(4L)

  em <- emObjfn(s$joint, s$om,
                model = s$model, data = s$data,
                method = "quadrature", control = list(level = 5L))
  psi <- c(mu_pop = 2.0, setNames(log(0.3), s$om$cholPars))
  attr(em, "rebuildQuadrature")(psi)

  f_val <- function(mu) em(c(mu_pop = mu), deriv = FALSE)$value
  gr_num <- numDeriv::grad(f_val, 2.1, method = "Richardson")
  gr_ana <- em(c(mu_pop = 2.1))$gradient["mu_pop"]
  expect_equal(unname(gr_ana), gr_num, tolerance = 1e-3)
})
