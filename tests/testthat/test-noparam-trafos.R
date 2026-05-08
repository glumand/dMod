context("Parameter transformations and observation functions with empty outer parameter sets")

# Regression for the case where an inner trafo / observable has no outer
# parameters: e.g. Pexpl(c(A = "1.0", B = "2.5")) where the right-hand sides
# are pure numerics. CppODE's funCpp handles parameters = character(0), and
# the dMod wrappers must propagate that through evaluation and composition.

test_that("Pexpl with pure-numeric trafo evaluates (symbolic and dual)", {

  setwd(tempdir())

  trafo <- c(A = "1.0", B = "2.5")

  p_sym <- Pexpl(trafo, derivMode = "symbolic", compile = FALSE,
                 modelname = "noparam_pexpl_sym")
  out_sym <- p_sym(c(dummy = 1.0))
  expect_equal(unclass(out_sym[[1]])[c("A", "B")], c(A = 1.0, B = 2.5))

  p_dual <- Pexpl(trafo, derivMode = "dual", compile = TRUE,
                  modelname = "noparam_pexpl_dual")
  out_dual <- p_dual(c(dummy = 1.0))
  expect_equal(unclass(out_dual[[1]])[c("A", "B")], c(A = 1.0, B = 2.5))
})


test_that("Pimpl with no outer parameters does not crash in build_jacobian", {

  setwd(tempdir())

  # One dependent state, no parameters: unique root A = 1.
  trafo <- c(A = "A - 1.0")
  p <- Pimpl(trafo, parameters = NULL, compile = TRUE,
             modelname = "noparam_pimpl")

  out <- p(c(dummy = 1.0))
  expect_true(is.numeric(unclass(out[[1]])["A"]))
  expect_equal(unname(unclass(out[[1]])["A"]), 1.0, tolerance = 1e-3)
})


test_that("Y with pure-numeric observable composes with an Xs prediction", {

  setwd(tempdir())

  f <- as.eqnvec(c(A = "-k*A"))
  m <- odemodel(f, modelname = "noparam_y_ode", compile = TRUE,
                solver = "CppODE")
  x <- Xs(m)

  g <- Y(c(y1 = "1.0"), f = NULL, states = c("A"),
         parameters = character(0),
         derivMode = "symbolic", compile = FALSE,
         modelname = "noparam_y_obs")

  out <- (g * x)(seq(0, 5, length.out = 3), c(A = 1.0, k = 0.1))
  pred <- out[[1]]
  expect_true(all(pred[, "y1"] == 1.0))
  expect_equal(pred[, "time"], c(0, 2.5, 5))
})


test_that("Full g*x*p chain with constant-only Pexpl evaluates", {

  setwd(tempdir())

  f <- as.eqnvec(c(A = "-k*A"))
  m <- odemodel(f, modelname = "noparam_full_ode", compile = TRUE,
                solver = "CppODE")
  x <- Xs(m)

  trafo <- c(A = "1.0", k = "0.5")
  p <- Pexpl(trafo, derivMode = "symbolic", compile = FALSE,
             modelname = "noparam_full_p")
  g <- Y(c(y1 = "A"), f = NULL, states = c("A"),
         parameters = character(0),
         derivMode = "symbolic", compile = FALSE,
         modelname = "noparam_full_g")

  out <- (g * x * p)(seq(0, 5, length.out = 3), c(dummy = 1.0))
  pred <- out[[1]]
  expect_equal(unname(pred[, "y1"]), unname(pred[, "A"]))
  expect_equal(unname(pred[1, "A"]), 1.0, tolerance = 1e-8)
  expect_equal(unname(pred[3, "A"]), exp(-0.5 * 5), tolerance = 1e-4)
})
