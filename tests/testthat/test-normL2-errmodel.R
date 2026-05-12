# Behavioral tests for normL2() with an error model supplying sigma.
#
# When data$sigma is NA, res() fills it from the errmodel matrix. This file
# verifies:
#   * proportional error: sigma = srel * y yields the closed-form
#     proportional-error log-likelihood
#   * gradient picks up the d log(sigma) / d theta contribution
#     (analytic via quotient rule)
#   * mixed-source case: rows with explicit sigma keep it, NA rows fall
#     through to the errmodel
#   * getParameters() on the errmodel-augmented objective includes the
#     errmodel parameters

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


# Build a prediction chain plus a proportional-error errmodel:
#   sigma(y) = srel * y    where y = A (the observable)
# srel passes through the trafo identically so it appears in the inner-par
# set of both prd and errmodel.
.build_prop_chain <- function(mn_suffix) {
  bench <- fx_decay_compiled()
  .dmod_with_fx_workdir({
    e_prop <- Y(c(y = "srel * y"), f = bench$gfn, attach.input = FALSE,
                condition = "C1",
                modelname = paste0("fx_decay_err_prop_", mn_suffix),
                compile = TRUE)
    pfn_prop <- P(eqnvec(A = "A", k = "k", srel = "srel"),
                  condition = "C1",
                  modelname = paste0("fx_decay_p_prop_", mn_suffix),
                  compile = TRUE)
    prd_prop <- bench$gfn * bench$xfn * pfn_prop
  })
  list(prd = prd_prop, e = e_prop)
}


## ---- Proportional error: closed-form value -----------------------------

test_that("normL2 with sigma = srel*y matches the proportional-error log-likelihood", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .build_prop_chain("val")

  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data(pars = pars[c("A", "k")], sigma = 0.05)
  data$C1$sigma <- NA_real_

  # The errmodel returns sigma_i = srel * y_i = srel * pred_i. Closed-form
  # -2 log L: chi^2 + sum log(2 pi sigma_i^2).
  prd_at <- ec$prd(times = data$C1$time, pars = pars, deriv = FALSE)$C1
  pred <- prd_at[, "y"]
  sigma_pred <- pars[["srel"]] * pred
  expected <- truth_nll_aloq(pred, data$C1$value, sigma_pred)

  for_each_backend(function(cpp) {
    obj <- normL2(data, ec$prd, errmodel = ec$e, use.bessel = FALSE)
    o   <- obj(pars)
    expect_equal(o$value, expected, tolerance = 1e-4,
                 info = paste0("cpp=", cpp))
  })
})


## ---- Sigma-dependent gradient: closed form -----------------------------

test_that("normL2 gradient with proportional errmodel follows the analytic closed form", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .build_prop_chain("grad")

  data <- fx_decay_data(pars = c(A = 1.0, k = 0.5), sigma = 0.05)
  data$C1$sigma <- NA_real_
  pars <- c(A = 1.2, k = 0.45, srel = 0.08)  # perturbed

  # Model: pred = A * exp(-k*t),  sigma = srel * pred
  # wr = (pred - obs) / sigma. value = sum(wr^2) + sum(log(2 pi sigma^2)).
  # dpred/dA = exp(-kt),  dpred/dk = -t * A * exp(-kt),  dpred/dsrel = 0
  # dsigma/dA = srel * dpred/dA,  dsigma/dk = srel * dpred/dk,  dsigma/dsrel = pred
  # dwr/dtheta = (dpred/dtheta * sigma - (pred - obs) * dsigma/dtheta) / sigma^2
  # d value / d theta = 2 sum wr * dwr/dtheta  +  sum 2 * dlog(sigma)/dtheta
  t <- data$C1$time; obs <- data$C1$value
  A <- pars[["A"]]; k <- pars[["k"]]; srel <- pars[["srel"]]
  pred <- A * exp(-k * t)
  sigma_vec <- srel * pred
  wr <- (pred - obs) / sigma_vec
  dpred <- cbind(A = exp(-k * t), k = -t * A * exp(-k * t), srel = 0)
  dsigma <- srel * dpred
  dsigma[, "srel"] <- pred
  dlog_sigma <- dsigma / sigma_vec
  dwr <- (dpred * sigma_vec - (pred - obs) * dsigma) / sigma_vec^2
  g_ref <- colSums(2 * wr * dwr) + colSums(2 * dlog_sigma)

  for_each_backend(function(cpp) {
    obj <- normL2(data, ec$prd, errmodel = ec$e, use.bessel = FALSE)
    g_ana <- obj(pars)$gradient
    expect_equal(unname(g_ana[names(g_ref)]), unname(g_ref),
                 tolerance = 1e-3, info = paste0("cpp=", cpp))
  })
})


## ---- Mixed sigma source: explicit + errmodel ---------------------------

test_that("rows with explicit sigma keep it; NA rows fall through to errmodel", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .build_prop_chain("mix")

  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data(pars = pars[c("A", "k")], sigma = 0.05)
  # Half NA (will fall through to errmodel), half explicit (keep 0.05).
  n <- nrow(data$C1)
  data$C1$sigma <- ifelse(seq_len(n) <= n %/% 2, NA_real_, 0.05)

  for_each_backend(function(cpp) {
    obj <- normL2(data, ec$prd, errmodel = ec$e, use.bessel = FALSE)
    o   <- obj(pars)
    expect_true(is.finite(o$value),
                label = paste0("cpp=", cpp, " mixed sigma finite"))
  })

  # Cross-check: split the dataset into pure-NA and pure-explicit halves,
  # build two separate objectives, and sum. With use.bessel = FALSE the
  # value is exactly additive across condition subsets.
  data_na <- data; data_na$C1 <- data_na$C1[is.na(data$C1$sigma), ]
  data_ex <- data; data_ex$C1 <- data_ex$C1[!is.na(data$C1$sigma), ]
  with_cpp_backend(FALSE, {
    v_full <- normL2(data,    ec$prd, errmodel = ec$e, use.bessel = FALSE)(pars)$value
    v_na   <- normL2(data_na, ec$prd, errmodel = ec$e, use.bessel = FALSE)(pars)$value
    v_ex   <- normL2(data_ex, ec$prd, errmodel = ec$e, use.bessel = FALSE)(pars)$value
  })
  expect_equal(v_full, v_na + v_ex, tolerance = 1e-9)
})


## ---- getParameters() exposes errmodel pars -----------------------------

test_that("getParameters(normL2(..., errmodel = ec$e)) includes errmodel pars", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .build_prop_chain("pars")
  data <- fx_decay_data()
  data$C1$sigma <- NA_real_
  obj <- normL2(data, ec$prd, errmodel = ec$e)
  expect_true("srel" %in% attr(obj, "parameters"))
})
