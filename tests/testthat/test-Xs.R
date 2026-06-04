# ============================================================================
# Behavioral tests for the prediction functions Xs / Xd / Xf.
#
# Closed-form analytical references throughout (linear decay, two-step
# cascade, Xd grid recovery, forced linear ODE). Hessian and second-order
# chain-rule semantics live in test-deriv2.R.
# ============================================================================

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


# ---- Xs: state trajectories (closed-form) -------------------------------

test_that("Xs on linear decay matches A0 * exp(-k * t)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  prd_states <- bench$xfn * bench$pfn_id

  times <- c(0, 0.5, 1, 2, 5, 10)
  pars  <- c(A = 1.7, k = 0.42)
  out <- prd_states(times = times, pars = pars, deriv = FALSE)

  closed <- pars[["A"]] * exp(-pars[["k"]] * times)
  expect_equal(out$C1[, "A"], closed, tolerance = 1e-5)
})


test_that("Xs on two-step cascade matches the closed-form A(t), B(t), C(t)", {
  skip_if_no_compile()
  testthat::skip_if_not_installed("CppODE")
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  reactions <- eqnlist() |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "C", "k2 * B")
  m  <- odemodel(reactions, modelname = "test_xs_cascade", compile = FALSE)
  xf <- Xs(m)
  pf <- P(eqnvec(A = "A", B = "0", C = "0", k1 = "k1", k2 = "k2"),
          condition = "C1",
          modelname = "test_xs_cascade_p", compile = FALSE)
  compile(xf, pf)
  prd <- xf * pf

  times <- c(0, 0.5, 1, 2, 4)
  pars <- c(A = 1.0, k1 = 0.7, k2 = 0.3)
  out <- prd(times = times, pars = pars, deriv = FALSE)

  closed <- truth_two_step(times, pars["A"], pars["k1"], pars["k2"])
  expect_equal(out$C1[, "A"], closed$A, tolerance = 1e-5)
  expect_equal(out$C1[, "B"], closed$B, tolerance = 1e-5)
  expect_equal(out$C1[, "C"], closed$C, tolerance = 1e-5)
})


# ---- Xs: first-order sensitivities --------------------------------------

test_that("Xs sensitivities on linear decay match analytical d/d(A0,k)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  prd_states <- bench$xfn * bench$pfn_id

  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.0, k = 0.5)
  out <- prd_states(times = times, pars = pars, deriv = TRUE)
  d <- attr(out$C1, "deriv")  # shape [time, var, par]

  expect_equal(d[, "A", "A"], exp(-pars[["k"]] * times),     tolerance = 1e-5)
  expect_equal(d[, "A", "k"], -times * pars[["A"]] * exp(-pars[["k"]] * times),
               tolerance = 1e-5)
})


test_that("Xs predictions across conditions are independent and parameter-local", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  .dmod_with_fx_workdir({
    pfn_C2 <- P(eqnvec(A = "A_C2", k = "k_C2"), condition = "C2",
                modelname = "test_xs_p_C2", compile = TRUE)
  })
  prd_multi <- bench$xfn * (bench$pfn_id + pfn_C2)
  times <- c(0, 1, 2, 5)
  pars <- c(A = 1.0, k = 0.5, A_C2 = 2.0, k_C2 = 1.0)
  out <- prd_multi(times = times, pars = pars, deriv = TRUE)

  expect_equal(out$C1[, "A"], 1.0 * exp(-0.5 * times), tolerance = 1e-5)
  expect_equal(out$C2[, "A"], 2.0 * exp(-1.0 * times), tolerance = 1e-5)

  dC1 <- attr(out$C1, "deriv")
  dC2 <- attr(out$C2, "deriv")
  expect_setequal(dimnames(dC1)[[3]], c("A", "k"))
  expect_setequal(dimnames(dC2)[[3]], c("A_C2", "k_C2"))
})


test_that("Xs sensitivities match the analytical decay formulas at non-integer times", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  prd_states <- bench$xfn * bench$pfn_id
  times <- c(0, 0.7, 1.4, 3.5)
  pars  <- c(A = 1.3, k = 0.65)

  out <- prd_states(times = times, pars = pars, deriv = TRUE)
  d <- attr(out$C1, "deriv")

  expect_equal(d[, "A", "A"], exp(-pars[["k"]] * times), tolerance = 1e-5)
  expect_equal(d[, "A", "k"],
               -times * pars[["A"]] * exp(-pars[["k"]] * times),
               tolerance = 1e-5)
})


# ---- Xs: events ---------------------------------------------------------

test_that("Xs with an 'add' event reproduces the analytical post-event trajectory", {
  # Linear decay with an additive jump at t0:
  #   pre   A(t)        = A0 * exp(-k * t)
  #   at t0 A(t0)       = A0 * exp(-k * t0) + Delta
  #   post  A(t)        = (A0 * exp(-k * t0) + Delta) * exp(-k * (t - t0))
  skip_if_no_compile()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)

  reactions <- eqnlist() |>
    addReaction("A", "", "k * A", "decay")
  ev <- eventlist(var = "A", time = 5, value = "A_add", method = "add")

  m  <- odemodel(reactions, events = ev,
                 modelname = paste0("xs_event_", as.integer(Sys.time())),
                 compile = FALSE)
  xf <- Xs(m)
  pf <- P(eqnvec(A = "A", k = "k", A_add = "A_add"),
          condition = "C1",
          modelname = paste0("xs_event_p_", as.integer(Sys.time())),
          compile = FALSE)
  compile(xf, pf)
  prd <- xf * pf

  A0 <- 1.0; k <- 0.4; Delta <- 0.5; t0 <- 5
  pars <- c(A = A0, k = k, A_add = Delta)
  times <- c(0, 1, 3, 4.99, 5.01, 6, 8, 10)
  out <- prd(times = times, pars = pars, deriv = FALSE)$C1

  A_t0_pre <- A0 * exp(-k * t0)
  A_t0_post <- A_t0_pre + Delta
  expected <- ifelse(times < t0,
                     A0 * exp(-k * times),
                     A_t0_post * exp(-k * (times - t0)))

  expect_equal(out[match(times, out[, "time"]), "A"], expected,
               tolerance = 1e-4)
})


# ---- Xs: forcings -------------------------------------------------------

test_that("Xs with constant forcing input matches the closed-form linear ODE solution", {
  # dA/dt = F - k*A with constant F gives A(t) = (A0 - F/k)*exp(-k*t) + F/k.
  skip_if_no_compile()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)

  reactions <- eqnlist() |>
    addReaction("",  "A", "F",     "production by forcing") |>
    addReaction("A", "",  "k * A", "decay")

  m <- odemodel(reactions, forcings = "F",
                modelname = paste0("xs_forc_", as.integer(Sys.time())),
                compile = FALSE)

  u_const <- 0.6
  forc <- data.frame(name = "F",
                     time = seq(0, 20, by = 0.5),
                     value = u_const)
  xf <- Xs(m, forcings = forc, condition = "C1")
  pf <- P(eqnvec(A = "A", k = "k"), condition = "C1",
          modelname = paste0("xs_forc_p_", as.integer(Sys.time())),
          compile = FALSE)
  compile(xf, pf)
  prd <- xf * pf

  A0 <- 0.1; k <- 0.3
  pars <- c(A = A0, k = k)
  times <- c(0, 1, 2, 5, 10, 20)
  out <- prd(times = times, pars = pars, deriv = FALSE)$C1

  expected <- (A0 - u_const / k) * exp(-k * times) + u_const / k
  expect_equal(out[match(times, out[, "time"]), "A"], expected,
               tolerance = 1e-4)

  expect_lt(abs(out[match(20, out[, "time"]), "A"] - u_const / k),
            abs(A0 - u_const / k) * exp(-k * 20) * 2)
})


# ---- Xd: linear-interpolation grid prediction ---------------------------

test_that("Xd returns the grid values at the grid times", {
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


# ---- Xf: no-sensitivity ODE prediction ----------------------------------

test_that("Xf reproduces the linear-decay closed form and emits no deriv attribute", {
  skip_if_no_compile()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)

  reactions <- eqnlist() |>
    addReaction("A", "", "k * A", "decay")
  m <- odemodel(reactions,
                modelname = paste0("xf_decay_", as.integer(Sys.time())),
                compile = FALSE)
  xfn <- Xf(m, condition = "C1"); compile(xfn)

  times <- c(0, 1, 2, 5)
  pars <- c(A = 1.3, k = 0.42)
  out <- xfn(times = times, pars = pars)$C1

  expect_equal(out[match(times, out[, "time"]), "A"],
               pars[["A"]] * exp(-pars[["k"]] * times),
               tolerance = 1e-5)
  expect_null(attr(out, "deriv"))
})


# ============================================================================
# Xs.CppODE theta-sensitivity path: heap vs stack AD slab parity
# (Phi'(theta) as sens1ini; per-condition varying theta counts via two-condition setup)
# ============================================================================

test_that("Heap and stack AD slabs match on a single-condition linear model", {

  withr::local_dir(tempdir())
  f <- c(A = "-k1*A + k2*B",
         B =  "k1*A - k2*B")

  # Default heap slab vs explicit stack slab (nStack = 4 covers {A,B,log_k1,log_k2}).
  mod_v1 <- odemodel(f, modelname = "rep_v1", solver = "CppODE")
  mod_v2 <- odemodel(f, modelname = "rep_v2", solver = "CppODE", nStack = 4L)

  # Same parameter transformation for both.
  trafo <- c(A = "A", B = "B", k1 = "exp(log_k1)", k2 = "exp(log_k2)")
  p1 <- P(trafo, modelname = "rep_trafo_v1", compile = TRUE)
  p2 <- P(trafo, modelname = "rep_trafo_v2", compile = TRUE)

  tight <- list(atol = 1e-10, rtol = 1e-10)
  x1 <- Xs(mod_v1, optionsSens = tight) * p1
  x2 <- Xs(mod_v2, optionsSens = tight) * p2

  theta <- c(A = 1.0, B = 0.2, log_k1 = log(0.5), log_k2 = log(0.3))
  times <- seq(0, 3, length.out = 7)

  pred1 <- x1(times, theta,
              conditions = NULL,
              deriv = TRUE)
  pred2 <- x2(times, theta,
              conditions = NULL,
              deriv = TRUE)

  d1 <- getDerivs(pred1)[[1]]
  d2 <- getDerivs(pred2)[[1]]

  arr1 <- attr(d1, "deriv")
  arr2 <- attr(d2, "deriv")
  expect_equal(dim(arr1), dim(arr2))
  expect_equal(dimnames(arr1)[[2]], dimnames(arr2)[[2]])
  expect_equal(dimnames(arr1)[[3]], dimnames(arr2)[[3]])
  expect_equal(as.numeric(arr1), as.numeric(arr2), tolerance = 1e-6)
})


test_that("Heap/stack parity holds with per-condition varying theta subsets", {

  withr::local_dir(tempdir())
  f <- c(A = "-k1*A + k2*B",
         B =  "k1*A - k2*B")

  mod_v1 <- odemodel(f, modelname = "repmulti_v1", solver = "CppODE")
  # Stack upper bound: any condition may activate up to 4 thetas.
  mod_v2 <- odemodel(f, modelname = "repmulti_v2", solver = "CppODE", nStack = 4L)

  # Condition "closed" uses log_k1; condition "open" uses log_k_open instead.
  # Global theta set has 5 elements; each condition activates 4.
  trafo_closed <- c(A = "A", B = "B", k1 = "exp(log_k1)",      k2 = "exp(log_k2)")
  trafo_open   <- c(A = "A", B = "B", k1 = "exp(log_k_open)",  k2 = "exp(log_k2)")

  p1 <-
    P(trafo_closed, condition = "closed",
      modelname = "repmulti_trafo_cl_v1", compile = TRUE) +
    P(trafo_open,   condition = "open",
      modelname = "repmulti_trafo_op_v1", compile = TRUE)

  p2 <-
    P(trafo_closed, condition = "closed",
      modelname = "repmulti_trafo_cl_v2", compile = TRUE) +
    P(trafo_open,   condition = "open",
      modelname = "repmulti_trafo_op_v2", compile = TRUE)

  tight <- list(atol = 1e-10, rtol = 1e-10)
  x1 <- Xs(mod_v1, optionsSens = tight) * p1
  x2 <- Xs(mod_v2, optionsSens = tight) * p2

  theta <- c(A = 1.0, B = 0.2,
             log_k1 = log(0.5), log_k_open = log(0.8), log_k2 = log(0.3))
  times <- seq(0, 2, length.out = 5)

  pred1 <- x1(times, theta, deriv = TRUE)
  pred2 <- x2(times, theta, deriv = TRUE)

  for (cond in c("closed", "open")) {
    arr1 <- attr(pred1[[cond]], "deriv")
    arr2 <- attr(pred2[[cond]], "deriv")
    expect_equal(dim(arr1), dim(arr2),
                 info = paste("shape mismatch for condition", cond))
    expect_equal(dim(arr2)[3], 4L,
                 info = paste("ncol mismatch for condition", cond))
    expect_equal(dimnames(arr1)[[3]], dimnames(arr2)[[3]],
                 info = paste("theta names mismatch for condition", cond))
    expect_equal(as.numeric(arr1), as.numeric(arr2), tolerance = 1e-6,
                 info = paste("values mismatch for condition", cond))
  }
})
