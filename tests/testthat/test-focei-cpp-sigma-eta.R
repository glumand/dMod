context("focei-cpp-sigma-eta")

# Generalized fast inner: when the error model depends on the prediction
# (proportional / combined errors), sigma is a function of eta and the kernel
# must include the dsigma/deta contributions in gradient and GN Hessian.
#
# Oracle: closed-form R replica of the kernel math, plus numDeriv on the OFV
# to validate the analytical gradient. The R replica is bit-exact since it
# evaluates the same formulas off the same prdframe sensitivities.
#
# NOTE: This file intentionally uses numDeriv as one of two independent
# checks on the analytic kernel gradient. Deriving a clean closed-form
# reference would require constructing a separate linear-in-theta toy
# model that exercises the FOCEI inner trust loop; the cost/benefit does
# not warrant rewriting. The R replica above is the primary oracle;
# numDeriv on the value is the independent backup.

test_that("fast inner: value/gradient/H_GN match R oracle for sigma(eta) (proportional)", {
  skip_on_cran()
  if (!requireNamespace("CppODE", quietly = TRUE))
    skip("CppODE not installed")
  if (!requireNamespace("numDeriv", quietly = TRUE))
    skip("numDeriv not installed")

  set.seed(7)
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
         force = TRUE)

  # Tiny synthetic PK: 1-cmt iv bolus, 3 subjects, K_eta = 1 (eta_V only).
  subjects <- c("A", "B", "C")
  times <- c(0.5, 1, 2, 4, 8)
  set.seed(123)
  V_true <- 20; Cl_true <- 5; dose <- 100
  pred_at <- function(t, V, Cl) (dose / V) * exp(-Cl / V * t)
  data_rows <- do.call(rbind, lapply(seq_along(subjects), function(i) {
    s <- subjects[i]
    eta <- c(-0.2, 0.1, 0.05)[i]
    V_i <- V_true * exp(eta)
    pred <- pred_at(times, V_i, Cl_true)
    sigma_i <- 0.1 * pred
    data.frame(name = "y", time = times,
               value = pred + rnorm(length(times), 0, sigma_i),
               sigma = NA_real_, condition = s, stringsAsFactors = FALSE)
  }))
  dlist <- as.datalist(data_rows)

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "Cc", "", "Cl/V * Cc", "elimination")
  m <- odemodel(reactions, modelname = "sigeta_ode", compile = TRUE,
                solver = "CppODE", deriv2 = FALSE)
  x <- Xs(m)
  g <- Y(c(y = "Cc"), x, modelname = "sigeta_obs", compile = TRUE,
         deriv2 = FALSE)
  # Proportional error: sigma_y = sigma_prop * y. sigma is a function of the
  # prediction y, hence of eta_V via the chain V -> y -> sigma.
  err <- Y(eqnvec(y = "sigma_prop * y"), g, attach.input = FALSE,
           compile = TRUE, modelname = "sigeta_err")

  trafo <- eqnvec(V  = "exp(tv + eta_V)",
                  Cl = "exp(tcl)",
                  Cc = "dose / exp(tv + eta_V)",
                  sigma_prop = "exp(log_sigma_prop)")
  subj_table <- data.frame(eta_V = paste0("eta_V_", subjects),
                           dose  = rep(dose, length(subjects)),
                           row.names = subjects, stringsAsFactors = FALSE)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE,
         modelname = "sigeta_p", deriv2 = FALSE)
  model <- g * x * p
  om <- omega(eta = "eta_V", subjects = subjects)

  outer_pars <- c(tv = log(V_true), tcl = log(Cl_true),
                  log_sigma_prop = log(0.1))
  outer_pars[om$cholPars] <- log(0.3)
  eta_vals <- c(eta_V_A = -0.15, eta_V_B = 0.05, eta_V_C = 0.1)
  pars_full <- c(outer_pars, eta_vals)

  K <- om$K
  pars_full_names <- names(pars_full)
  eta_names_list <- lapply(seq_along(subjects),
                           function(i) paste0("eta_V_", subjects[i]))
  pars_probe <- pars_full
  pars_probe[grep("^eta_V_", names(pars_probe))] <- 0
  fast_meta <- dMod:::.buildFastMeta(model, err, dlist, subjects,
                                     eta_names_list, pars_full_names,
                                     pars_probe)
  subject_meta <- list(
    subjects        = subjects,
    eta_idx_global  = matrix(match(unlist(eta_names_list), pars_full_names),
                             nrow = length(subjects), ncol = K, byrow = TRUE),
    eta_names       = eta_names_list,
    K               = K,
    fast_meta       = fast_meta)

  L <- om$buildL(pars_full[om$cholPars])
  Omega <- tcrossprod(L)
  Omega_inv <- solve(Omega)
  Omega_log_det <- as.numeric(determinant(Omega, logarithm = TRUE)$modulus)

  eta_warmstart <- matrix(eta_vals, nrow = length(subjects), ncol = K)

  res <- dMod:::focei_inner_trust(
    model_cb = model, err_cb = err,
    pars_full = pars_full, eta_warmstart = eta_warmstart,
    subject_meta = subject_meta,
    Omega_inv_mat = Omega_inv, Omega_log_det = Omega_log_det,
    fixed = NULL,
    control = list(iterlim = 0))
  # focei_inner_trust returns iterations = 1 when iterlim = 0 (the
  # trust loop's for-init bumps the counter once before the bound check).
  # The point is that fast_eval is evaluated exactly at the warmstart.
  expect_lte(unname(res$iterations[[1]]), 1L)

  # R replica of fast_eval_one_subject (bit-exact match to the kernel).
  fast_eval_R <- function(eta_block, subject) {
    pars <- pars_full
    pars[paste0("eta_V_", subject)] <- eta_block
    pred_list <- model(times = times, pars = pars, deriv = TRUE,
                       conditions = subject)
    prdf <- pred_list[[1]]
    pinner <- attr(prdf, "parameters")
    err_list <- err(out = prdf, pars = pinner, conditions = subject)
    erm <- err_list[[1]]

    data_i <- data_rows[data_rows$condition == subject, ]
    ti_p <- match(data_i$time, prdf[, "time"])
    ti_e <- match(data_i$time, erm [, "time"])
    eta_nm <- paste0("eta_V_", subject)
    pred_vals  <- prdf[ti_p, "y"]
    sigma_vals <- erm [ti_e, "y"]
    Jp <- attr(prdf, "deriv")[ti_p, "y", eta_nm]
    Js <- attr(erm,  "deriv")[ti_e, "y", eta_nm]

    wr     <- (pred_vals - data_i$value) / sigma_vals
    inv_s  <- 1 / sigma_vals
    inv_s2 <- inv_s * inv_s
    dwr    <- inv_s * (Jp - wr * Js)
    dlogs  <- inv_s * Js

    value <- sum(wr^2 + log(2 * pi * sigma_vals^2))
    grad  <- 2 * sum(wr * dwr + dlogs)
    # Hessian (Part0..3, mirrors fast_eval_one_subject)
    H_GN <- 2 * sum(dwr^2)
    H_GN <- H_GN + sum(-2 * wr * inv_s2 * 2 * Jp * Js)
    H_GN <- H_GN + sum( 4 * wr^2 * inv_s2 * Js^2)
    H_GN <- H_GN + sum(-2 * inv_s2 * Js^2)

    # eta-prior
    eta_q <- sum(eta_block * (Omega_inv %*% eta_block))
    value <- value + eta_q + Omega_log_det
    grad  <- grad + 2 * (Omega_inv %*% eta_block)[1]
    H_GN  <- H_GN + 2 * Omega_inv[1, 1]
    list(value = value, grad = grad, hess = H_GN)
  }

  for (i in seq_along(subjects)) {
    s <- subjects[i]
    eta_i <- eta_vals[paste0("eta_V_", s)]
    R <- fast_eval_R(eta_i, s)

    # Value: R replica matches kernel bit-exactly (modulo float-add order).
    expect_equal(unname(res$value[[s]]), R$value, tolerance = 1e-12,
                 info = paste("value, subject", s))
    # Gradient: same.
    expect_equal(as.numeric(res$gradient[s, ]), R$grad, tolerance = 1e-12,
                 info = paste("gradient, subject", s))
    # Hessian: same.
    expect_equal(as.numeric(res$H_GN[[s]]), R$hess, tolerance = 1e-12,
                 info = paste("H_GN, subject", s))

    # Cross-check: analytical gradient must match numDeriv on the same OFV.
    # numDeriv ~ 1e-6 absolute on a function of magnitude O(1e2) with grad O(1e1).
    f_val <- function(eta_block) {
      pars <- pars_full
      pars[paste0("eta_V_", s)] <- eta_block
      pred_list <- model(times = times, pars = pars, deriv = FALSE,
                         conditions = s)
      prdf <- pred_list[[1]]
      pinner <- attr(prdf, "parameters")
      err_list <- err(out = prdf, pars = pinner, conditions = s)
      erm <- err_list[[1]]
      data_i <- data_rows[data_rows$condition == s, ]
      ti_p <- match(data_i$time, prdf[, "time"])
      ti_e <- match(data_i$time, erm [, "time"])
      wr <- (prdf[ti_p, "y"] - data_i$value) / erm[ti_e, "y"]
      val <- sum(wr^2 + log(2 * pi * erm[ti_e, "y"]^2))
      eta_q <- sum(eta_block * (Omega_inv %*% eta_block))
      val + eta_q + Omega_log_det
    }
    g_num <- numDeriv::grad(f_val, eta_i)
    # Relative tolerance: Richardson finite diff caps out around 1e-4 here
    # because of ODE-solver tolerance amplification through the chain rule.
    rel_err <- abs(as.numeric(res$gradient[s, ]) - g_num) /
               max(abs(g_num), 1)
    expect_lt(rel_err, 5e-3,
              label = paste("rel err kernel grad vs numDeriv, subject", s))
  }
})
