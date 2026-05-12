# Behavioral tests for nll() (the per-condition -2 log L kernel).
#
# nll consumes the output of res() and an external parameter vector. It
# routes ALOQ rows to nll_ALOQ and BLOQ rows to nll_BLOQ depending on the
# opt.BLOQ argument. These tests build res() outputs from a hand-made
# prediction matrix so each assertion is a direct test against the
# documented formula rather than against another implementation.


# Shared helper: a 4-time prediction for A(t) = exp(-0.5 * t) plus its
# analytical sensitivities dA/dA_par, dA/dk_par.
.make_prdframe_nll <- function(times = c(0, 1, 2, 3)) {
  prdf <- matrix(c(times, exp(-0.5 * times)), nrow = length(times), ncol = 2,
                 dimnames = list(NULL, c("time", "A")))
  d_arr <- array(0, c(length(times), 1, 2),
                 dimnames = list(NULL, "A", c("A_par", "k_par")))
  d_arr[, "A", "A_par"] <- exp(-0.5 * times)
  d_arr[, "A", "k_par"] <- -times * 1.0 * exp(-0.5 * times)
  attr(prdf, "deriv") <- d_arr
  prdf
}


## ---- ALOQ value formula ------------------------------------------------

test_that("nll on a pure-ALOQ dataset matches sum(wr^2) + sum(log(2*pi*sigma^2))", {
  prdf <- .make_prdframe_nll()
  dat <- data.frame(time = c(1, 2, 3), name = "A",
                    value = c(0.55, 0.30, 0.18),
                    sigma = c(0.1, 0.1, 0.1),
                    lloq  = -Inf,
                    stringsAsFactors = FALSE)
  nout <- res(dat, prdf)
  pars <- c(A_par = 1.0, k_par = 0.5)

  o <- nll(nout, pars = pars, deriv = TRUE)
  pred <- prdf[match(dat$time, prdf[, "time"]), "A"]
  expected <- truth_nll_aloq(pred, dat$value, dat$sigma)
  expect_equal(unname(o$value), expected, tolerance = 1e-12)
})


## ---- Bessel correction --------------------------------------------------

test_that("nll bessel.correction inflates wr^2 by exactly factor^2", {
  prdf <- .make_prdframe_nll()
  dat <- data.frame(time = c(1, 2, 3), name = "A",
                    value = c(0.55, 0.30, 0.18),
                    sigma = 0.1, lloq = -Inf,
                    stringsAsFactors = FALSE)
  nout <- res(dat, prdf)
  pars <- c(A_par = 1.0, k_par = 0.5)

  o_no <- nll(nout, pars = pars, deriv = TRUE, bessel.correction = 1)
  o_bs <- nll(nout, pars = pars, deriv = TRUE, bessel.correction = 1.2)

  log_sigma_term <- sum(log(2 * pi * dat$sigma^2))
  chi_no <- o_no$value - log_sigma_term
  chi_bs <- o_bs$value - log_sigma_term
  expect_equal(unname(chi_bs), unname(chi_no * 1.2^2), tolerance = 1e-12)
})


## ---- BLOQ M1: drop BLOQ rows ------------------------------------------

test_that("opt.BLOQ='M1' drops BLOQ rows entirely (value equals ALOQ-only nll)", {
  prdf <- .make_prdframe_nll()
  dat <- data.frame(time = c(1, 2, 3), name = "A",
                    value = c(0.55, 0.30, 0.05),
                    sigma = 0.1,
                    lloq  = c(-Inf, -Inf, 0.10),
                    stringsAsFactors = FALSE)
  nout <- res(dat, prdf)
  pars <- c(A_par = 1.0, k_par = 0.5)

  # M1: BLOQ row dropped -> only rows 1, 2 contribute.
  o_m1 <- nll(nout, pars = pars, deriv = TRUE, opt.BLOQ = "M1")

  # Reference: ALOQ-only res via subsetting raw data.
  nout_aloq <- res(dat[1:2, , drop = FALSE], prdf)
  o_ref <- nll(nout_aloq, pars = pars, deriv = TRUE)

  expect_equal(unname(o_m1$value), unname(o_ref$value), tolerance = 1e-12)
})


## ---- BLOQ M3: closed-form additive term -------------------------------

test_that("opt.BLOQ='M3' adds -2 * sum(log(Phi(-wr_bloq))) on top of ALOQ value", {
  prdf <- .make_prdframe_nll()
  dat <- data.frame(time = c(1, 2, 3), name = "A",
                    value = c(0.55, 0.30, 0.05),
                    sigma = 0.1,
                    lloq  = c(-Inf, -Inf, 0.10),
                    stringsAsFactors = FALSE)
  nout <- res(dat, prdf)
  pars <- c(A_par = 1.0, k_par = 0.5)

  o_m1 <- nll(nout, pars = pars, deriv = TRUE, opt.BLOQ = "M1")
  o_m3 <- nll(nout, pars = pars, deriv = TRUE, opt.BLOQ = "M3")

  # BLOQ row: pred at t=3, value substituted to lloq.
  pred_bloq <- prdf[match(3, prdf[, "time"]), "A"]
  sigma_bloq <- 0.1
  lloq_val <- 0.10
  m3_term <- truth_nll_bloq_m3(pred_bloq, lloq_val, sigma_bloq)
  expect_equal(unname(o_m3$value - o_m1$value), m3_term, tolerance = 1e-12)
})


## ---- BLOQ M4NM: closed-form ratio formula -----------------------------

test_that("opt.BLOQ='M4NM' adds -2*sum(log(1 - Phi(wr)/Phi(w0))) on top of ALOQ value", {
  prdf <- .make_prdframe_nll()
  dat <- data.frame(time = c(1, 2, 3), name = "A",
                    value = c(0.55, 0.30, 0.05),
                    sigma = 0.1,
                    lloq  = c(-Inf, -Inf, 0.10),
                    stringsAsFactors = FALSE)
  nout <- res(dat, prdf)
  pars <- c(A_par = 1.0, k_par = 0.5)

  o_m1   <- nll(nout, pars = pars, deriv = TRUE, opt.BLOQ = "M1")
  o_m4nm <- nll(nout, pars = pars, deriv = TRUE, opt.BLOQ = "M4NM")

  pred_bloq <- prdf[match(3, prdf[, "time"]), "A"]
  m4_term <- truth_nll_bloq_m4(pred_bloq, 0.10, 0.1)
  expect_equal(unname(o_m4nm$value - o_m1$value), m4_term, tolerance = 1e-12)
})


## ---- Gradient: closed form ---------------------------------------------

test_that("nll gradient on ALOQ equals 2 * sum_i wr_i * dpred_i/dtheta / sigma_i", {
  prdf <- .make_prdframe_nll()
  dat <- data.frame(time = c(1, 2, 3), name = "A",
                    value = c(0.55, 0.30, 0.18),
                    sigma = 0.1, lloq = -Inf,
                    stringsAsFactors = FALSE)
  nout <- res(dat, prdf)
  pars <- c(A_par = 1.0, k_par = 0.5)

  o <- nll(nout, pars = pars, deriv = TRUE)

  # Closed form for A(t) = A0 * exp(-k * t), value = sum(wr^2) + log_term:
  #   dvalue/dA0 = 2 * sum_i wr_i * exp(-k * t_i) / sigma_i
  #   dvalue/dk  = 2 * sum_i wr_i * (-t_i * A0 * exp(-k * t_i)) / sigma_i
  t <- dat$time
  pred <- prdf[match(t, prdf[, "time"]), "A"]
  wr <- (pred - dat$value) / dat$sigma
  J_A <- exp(-pars[["k_par"]] * t)
  J_k <- -t * pars[["A_par"]] * exp(-pars[["k_par"]] * t)
  g_ref <- c(A_par = 2 * sum(wr * J_A / dat$sigma),
             k_par = 2 * sum(wr * J_k / dat$sigma))
  expect_equal(unname(o$gradient[names(g_ref)]), unname(g_ref),
               tolerance = 1e-12)
})
