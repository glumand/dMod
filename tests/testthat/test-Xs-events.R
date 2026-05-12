# Behavioral tests for Xs() with event handling.
#
# Verifies that an event "add" at time t0 with value Delta produces the
# analytical post-event trajectory for a linear decay system:
#
#   pre-event:   A(t) = A0 * exp(-k * t)
#   at t0:       A(t0) = A0 * exp(-k * t0) + Delta
#   post-event:  A(t) = (A0 * exp(-k * t0) + Delta) * exp(-k * (t - t0))

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


test_that("Xs with an 'add' event reproduces the analytical post-event trajectory", {
  skip_if_no_compile()
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)

  reactions <- eqnlist() |>
    addReaction("A", "", "k * A", "decay")
  ev <- eventlist(var = "A", time = 5, value = "A_add", method = "add")

  m  <- odemodel(reactions, events = ev,
                 modelname = paste0("xs_event_", as.integer(Sys.time())),
                 compile = TRUE)
  xf <- Xs(m)
  pf <- P(eqnvec(A = "A", k = "k", A_add = "A_add"),
          condition = "C1",
          modelname = paste0("xs_event_p_", as.integer(Sys.time())),
          compile = TRUE)
  prd <- xf * pf

  A0 <- 1.0; k <- 0.4; Delta <- 0.5; t0 <- 5
  pars <- c(A = A0, k = k, A_add = Delta)
  # Avoid the event time exactly; pre/post sides are unambiguous.
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
