# Behavioral tests for Xd() (data-grid prediction) and Xf() (no-sensitivity
# ODE prediction).
#
# Xd: linear interpolation through a fixed (time, value)-grid; verify both
# exact grid recovery and midpoint linear-interpolation values.
# Xf: same ODE as Xs but without sensitivities; verify the analytical
# decay trajectory and the absence of the "deriv" attribute.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


## ---- Xd: exact recovery at grid points ---------------------------------

test_that("Xd returns the grid values at the grid times", {
  # Data grid: state A at four times. Rownames are the parameter names
  # whose values define the corresponding A-values.
  grid <- data.frame(
    name = "A",
    time = c(0, 1, 2, 3),
    row.names = c("p0", "p1", "p2", "p3"))
  pars <- c(p0 = 1.0, p1 = 0.5, p2 = 0.25, p3 = 0.125)

  xfn <- Xd(grid, condition = "C1")
  out <- xfn(times = c(0, 1, 2, 3), pars = pars)$C1
  expect_equal(out[match(c(0, 1, 2, 3), out[, "time"]), "A"],
               c(1.0, 0.5, 0.25, 0.125), tolerance = 1e-12)
})


## ---- Xd: linear interpolation between grid points ---------------------

test_that("Xd linearly interpolates between grid points", {
  grid <- data.frame(
    name = "A",
    time = c(0, 1, 2, 3),
    row.names = c("p0", "p1", "p2", "p3"))
  pars <- c(p0 = 1.0, p1 = 0.5, p2 = 0.25, p3 = 0.125)

  xfn <- Xd(grid, condition = "C1")
  mid_times <- c(0.5, 1.5, 2.5)
  out <- xfn(times = mid_times, pars = pars)$C1
  expected <- c((1.0 + 0.5) / 2, (0.5 + 0.25) / 2, (0.25 + 0.125) / 2)
  expect_equal(out[match(mid_times, out[, "time"]), "A"],
               expected, tolerance = 1e-12)
})


## ---- Xf: closed-form decay value, no sensitivities -------------------

test_that("Xf reproduces the linear-decay closed form and emits no deriv attribute", {
  skip_if_no_compile()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)

  reactions <- eqnlist() |>
    addReaction("A", "", "k * A", "decay")
  m <- odemodel(reactions,
                modelname = paste0("xf_decay_", as.integer(Sys.time())),
                compile = TRUE)
  xfn <- Xf(m, condition = "C1")

  times <- c(0, 1, 2, 5)
  pars <- c(A = 1.3, k = 0.42)
  out <- xfn(times = times, pars = pars)$C1

  expect_equal(out[match(times, out[, "time"]), "A"],
               pars[["A"]] * exp(-pars[["k"]] * times),
               tolerance = 1e-5)
  # Xf is the no-sensitivity path: no "deriv" attribute on the prdframe.
  expect_null(attr(out, "deriv"))
})
