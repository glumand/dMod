# Behavioral tests for Xs() with forcing inputs.
#
# Verifies that a constant forcing input u(t) = u_const driving the linear
# system  dA/dt = -k*A + u(t)  produces the closed-form solution
#
#   A(t) = (A0 - u_const/k) * exp(-k*t) + u_const/k
#
# At t -> infinity this asymptotes to u_const/k, which we also check.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


test_that("Xs with constant forcing input matches the closed-form linear ODE solution", {
  skip_if_no_compile()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)

  # dA/dt = F - k*A.  F is a forcing input; k a parameter; A the state.
  reactions <- eqnlist() |>
    addReaction("",  "A", "F",     "production by forcing") |>
    addReaction("A", "",  "k * A", "decay")

  m <- odemodel(reactions, forcings = "F",
                modelname = paste0("xs_forc_", as.integer(Sys.time())),
                compile = TRUE)

  # Constant forcing input u_const at all times.
  u_const <- 0.6
  forc <- data.frame(name = "F",
                     time = seq(0, 20, by = 0.5),
                     value = u_const)
  xf <- Xs(m, forcings = forc, condition = "C1")
  pf <- P(eqnvec(A = "A", k = "k"), condition = "C1",
          modelname = paste0("xs_forc_p_", as.integer(Sys.time())),
          compile = TRUE)
  prd <- xf * pf

  A0 <- 0.1; k <- 0.3
  pars <- c(A = A0, k = k)
  times <- c(0, 1, 2, 5, 10, 20)
  out <- prd(times = times, pars = pars, deriv = FALSE)$C1

  expected <- (A0 - u_const / k) * exp(-k * times) + u_const / k
  expect_equal(out[match(times, out[, "time"]), "A"], expected,
               tolerance = 1e-4)

  # Asymptote sanity: A(t = 20) should be close to u_const / k. We use the
  # closed-form residual exp(-k*20) bound rather than a hard absolute tol.
  expect_lt(abs(out[match(20, out[, "time"]), "A"] - u_const / k),
            abs(A0 - u_const / k) * exp(-k * 20) * 2)
})
