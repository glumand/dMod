# ============================================================================
# FOCEI tests: end-to-end nlmeFit orchestrator + C++ kernel parity.
#
# Sections:
#   * End-to-end nlmeFit(method = "focei") on a minimal one-eta NLME prdfn.
#   * Pre-rewrite Theoph regression vs fixtures/focei_theoph_reference.rds.
#   * C++ kernel parity against R replica on a sigma(eta) (proportional) model.
#   * C++ kernel parity on a 2-output (parent/metabolite) model.
# ============================================================================

context("FOCEI orchestrator + C++ kernel")


# ---- End-to-end nlmeFit ---------------------------------------------------

test_that("nlmeFit(method='focei') runs on a minimal one-eta NLME prdfn", {
  set.seed(1)

  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))

  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE, modelname = "focei_smoke_obs")
  x <- Xt()

  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subjects <- c("s1", "s2", "s3", "s4")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = "focei_smoke_p")

  true_mu  <- 2.0
  true_om  <- 0.3
  true_eta <- rnorm(length(subjects), 0, true_om)
  y_obs    <- true_mu * exp(true_eta) + rnorm(length(subjects), 0, 0.2)

  data <- as.datalist(data.frame(
    name = "y", time = 0, sigma = 0.2,
    value = y_obs, condition = subjects,
    stringsAsFactors = FALSE))

  om <- omega(eta = "eta", subjects = subjects)
  obj <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)

  outer_init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(nlmeFit(
    obj, om, outer_init,
    prdfn = g * x * p, data = data,
    method = "focei",
    control = list(focei = list(
      trustControl = list(rinit = 1, rmax = 10, iterlim = 50,
                          fterm = 1e-7, mterm = 1e-7)))))
  expect_s3_class(fit, "nlmeFit")
  expect_equal(fit$method, "focei")
  expect_true(fit$converged)
  expect_true(is.finite(fit$value))
  expect_true(abs(fit$argument["mu_pop"] - mean(y_obs)) < 0.5)
  expect_true(!is.null(fit$Omega))
  expect_equal(dim(fit$Omega), c(1L, 1L))
  expect_equal(dim(fit$etaModes), c(length(subjects), 1L))
})


test_that("nlmeFit() rejects unknown method via match.arg", {
  expect_error(nlmeFit(obj = NULL, omega = NULL, init = c(p = 1),
                       method = "doesNotExist"),
               "'arg' should be one of")
})


# ---- Pre-rewrite Theoph regression ---------------------------------------

test_that("nlmeFit(method='focei') matches the pre-rewrite Theoph baseline", {
  # Anchor against the Phase 0 baseline recorded before the C++ kernel landed.
  # The fixture stores ($value, $argument) at convergence with eager Stage-2
  # correction; the consolidated kernel must reproduce them within
  # Schur/eigen tolerance.
  fixture_path <- "fixtures/focei_theoph_reference.rds"

  skip_on_cran()
  if (!requireNamespace("CppODE", quietly = TRUE))
    skip("CppODE not installed")
  if (!file.exists(fixture_path))
    skip(paste0("regression fixture not found at ", fixture_path))

  ref <- readRDS(fixture_path)

  set.seed(1)
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
         force = TRUE)

  data(Theoph, package = "datasets")
  Theoph$Subject <- as.character(Theoph$Subject)
  subjects <- sort(unique(Theoph$Subject))

  doses <- vapply(subjects, function(s) {
    rec <- Theoph[Theoph$Subject == s, ][1, ]
    rec$Dose * rec$Wt
  }, 0.0)
  dlist <- as.datalist(data.frame(
    name = "y", time = Theoph$Time, value = Theoph$conc,
    sigma = NA_real_, condition = Theoph$Subject,
    stringsAsFactors = FALSE))

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
  reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
  reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")
  m <- odemodel(reactions, modelname = "theoph_cppreg", compile = TRUE,
                solver = "CppODE", deriv2 = TRUE)
  x <- Xs(m)
  g <- Y(c(y = "Cc"), x, modelname = "theoph_cppreg_obs",
         compile = TRUE, deriv2 = TRUE)
  err <- Y(eqnvec(y = "sigma_add"), g, attach.input = FALSE,
           compile = TRUE, modelname = "theoph_cppreg_err")

  trafo <- eqnvec(Ka        = "exp(tka + eta_Ka)",
                  V         = "exp(tv  + eta_V)",
                  Cl        = "exp(tcl + eta_Cl)",
                  Ag        = "Ag_init",
                  Cc        = "0",
                  sigma_add = "exp(log_sigma_add)")
  subj_table <- data.frame(
    eta_Ka  = paste0("eta_Ka_", subjects),
    eta_V   = paste0("eta_V_",  subjects),
    eta_Cl  = paste0("eta_Cl_", subjects),
    Ag_init = doses, row.names = subjects, stringsAsFactors = FALSE)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE,
         modelname = "theoph_cppreg_p", deriv2 = TRUE)
  prdfn <- g * x * p

  om <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
  obj <- normL2(dlist, prdfn, errmodel = err, use.bessel = FALSE) +
           constraintL2(mu = 0, Omega = om)

  fit <- nlmeFit(obj, om, ref$init,
                 prdfn    = prdfn,
                 data     = dlist,
                 errfn = err,
                 method   = "focei",
                 control  = list(focei = list(
                   innerControl = list(rtol = 1e-7, maxit = 30),
                   trustControl = list(iterlim = 100))),
                 verbose = FALSE)

  expect_true(fit$converged)
  expect_lt(abs(fit$value - ref$value), 1e-2)
  expect_lt(max(abs(fit$argument[names(ref$argument)] - ref$argument)), 5e-3)
})


# ---- C++ kernel parity: sigma(eta) (proportional error) ------------------

test_that("fast inner: value/gradient/H_GN match R oracle for sigma(eta) (proportional)", {
  # Generalized fast inner: when the error model depends on the prediction
  # (proportional / combined errors), sigma is a function of eta and the
  # kernel must include the dsigma/deta contributions in gradient and GN
  # Hessian. Oracle: closed-form R replica of the kernel math, plus numDeriv
  # on the OFV to validate the analytical gradient.

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
    H_GN <- 2 * sum(dwr^2)
    H_GN <- H_GN + sum(-2 * wr * inv_s2 * 2 * Jp * Js)
    H_GN <- H_GN + sum( 4 * wr^2 * inv_s2 * Js^2)
    H_GN <- H_GN + sum(-2 * inv_s2 * Js^2)

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

    expect_equal(unname(res$value[[s]]), R$value, tolerance = 1e-12,
                 info = paste("value, subject", s))
    expect_equal(as.numeric(res$gradient[s, ]), R$grad, tolerance = 1e-12,
                 info = paste("gradient, subject", s))
    expect_equal(as.numeric(res$H_GN[[s]]), R$hess, tolerance = 1e-12,
                 info = paste("H_GN, subject", s))

    # Cross-check: analytical gradient must match numDeriv on the same OFV.
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
    rel_err <- abs(as.numeric(res$gradient[s, ]) - g_num) /
               max(abs(g_num), 1)
    expect_lt(rel_err, 5e-3,
              label = paste("rel err kernel grad vs numDeriv, subject", s))
  }
})


# ---- C++ kernel parity: 2-output (parent / metabolite) model -------------

test_that("fast inner: value/gradient/H_GN match R oracle for a 2-output model", {
  # Multi-output fast inner: data has multiple observables per subject, with
  # distinct sigma per observable. Long-format meta carries per-row indices
  # into the model and err deriv arrays; the kernel loops over rows, not
  # over time x observable.

  skip_on_cran()
  if (!requireNamespace("CppODE", quietly = TRUE))
    skip("CppODE not installed")

  set.seed(11)
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
         force = TRUE)

  # Tiny parent / metabolite model: A -> B -> 0. Observe both.
  subjects <- c("S1", "S2", "S3")
  times <- c(0.5, 1, 2, 4, 8)
  pred_at <- function(t, ka, ke) {
    # A(t) = A0 exp(-ka t), B(t) = (ka A0)/(ke-ka) (exp(-ka t) - exp(-ke t))
    A0 <- 10
    A <- A0 * exp(-ka * t)
    B <- (ka * A0 / (ke - ka)) * (exp(-ka * t) - exp(-ke * t))
    list(A = A, B = B)
  }
  ka_true <- 0.7; ke_true <- 0.3
  set.seed(321)
  data_rows <- do.call(rbind, lapply(seq_along(subjects), function(i) {
    s <- subjects[i]
    eta_ka <- c(-0.2, 0.15, 0.05)[i]
    eta_ke <- c(0.1, -0.05, 0.0)[i]
    pred <- pred_at(times, ka_true * exp(eta_ka), ke_true * exp(eta_ke))
    rows_A <- data.frame(name = "yA", time = times,
                         value = pred$A + rnorm(length(times), 0, 0.15),
                         sigma = NA_real_, condition = s,
                         stringsAsFactors = FALSE)
    rows_B <- data.frame(name = "yB", time = times,
                         value = pred$B + rnorm(length(times), 0, 0.1),
                         sigma = NA_real_, condition = s,
                         stringsAsFactors = FALSE)
    rbind(rows_A, rows_B)
  }))
  dlist <- as.datalist(data_rows)

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "A", "B", "ka * A", "absorption")
  reactions <- addReaction(reactions, "B", "",  "ke * B", "elimination")
  m <- odemodel(reactions, modelname = "mo_ode", compile = TRUE,
                solver = "CppODE", deriv2 = FALSE)
  x <- Xs(m)
  g <- Y(c(yA = "A", yB = "B"), x, modelname = "mo_obs", compile = TRUE,
         deriv2 = FALSE)
  err <- Y(eqnvec(yA = "sigA", yB = "sigB"), g, attach.input = FALSE,
           compile = TRUE, modelname = "mo_err")

  trafo <- eqnvec(ka = "exp(tka + eta_ka)",
                  ke = "exp(tke + eta_ke)",
                  A  = "A0",
                  B  = "0",
                  sigA = "exp(lsigA)",
                  sigB = "exp(lsigB)")
  subj_table <- data.frame(eta_ka = paste0("eta_ka_", subjects),
                           eta_ke = paste0("eta_ke_", subjects),
                           A0     = rep(10, length(subjects)),
                           row.names = subjects, stringsAsFactors = FALSE)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE,
         modelname = "mo_p", deriv2 = FALSE)
  model <- g * x * p
  om <- omega(eta = c("eta_ka", "eta_ke"), subjects = subjects)

  outer_pars <- c(tka = log(ka_true), tke = log(ke_true),
                  lsigA = log(0.15), lsigB = log(0.1))
  outer_pars[om$cholPars] <- log(0.3)
  eta_vals <- c(eta_ka_S1 = -0.15, eta_ke_S1 = 0.05,
                eta_ka_S2 = 0.10,  eta_ke_S2 = -0.02,
                eta_ka_S3 = 0.00,  eta_ke_S3 = 0.00)
  pars_full <- c(outer_pars, eta_vals)

  K <- om$K
  pars_full_names <- names(pars_full)
  eta_names_list <- lapply(seq_along(subjects), function(i) {
    c(paste0("eta_ka_", subjects[i]), paste0("eta_ke_", subjects[i]))
  })
  pars_probe <- pars_full
  pars_probe[grep("^eta_", names(pars_probe))] <- 0
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

  eta_warmstart <- matrix(eta_vals, nrow = length(subjects), ncol = K,
                          byrow = TRUE)

  res <- dMod:::focei_inner_trust(
    model_cb = model, err_cb = err,
    pars_full = pars_full, eta_warmstart = eta_warmstart,
    subject_meta = subject_meta,
    Omega_inv_mat = Omega_inv, Omega_log_det = Omega_log_det,
    fixed = NULL,
    control = list(iterlim = 0))

  # R replica: same per-row math, looping over all data rows (mixed yA/yB).
  fast_eval_R <- function(eta_block, subject) {
    pars <- pars_full
    eta_nms <- paste0(c("eta_ka_", "eta_ke_"), subject)
    pars[eta_nms] <- eta_block
    pred_list <- model(times = times, pars = pars, deriv = TRUE,
                       conditions = subject)
    prdf <- pred_list[[1]]
    pinner <- attr(prdf, "parameters")
    err_list <- err(out = prdf, pars = pinner, conditions = subject)
    erm <- err_list[[1]]

    data_i <- data_rows[data_rows$condition == subject, ]
    ti_p <- match(data_i$time, prdf[, "time"])
    ti_e <- match(data_i$time, erm [, "time"])
    n <- nrow(data_i)

    value <- 0
    grad <- numeric(K)
    H_GN <- matrix(0, K, K)
    for (j in seq_len(n)) {
      o <- data_i$name[j]
      pred_val <- prdf[ti_p[j], o]
      sigma_val <- erm [ti_e[j], o]
      Jp <- attr(prdf, "deriv")[ti_p[j], o, eta_nms]
      Js <- attr(erm,  "deriv")[ti_e[j], o, eta_nms]
      wr <- (pred_val - data_i$value[j]) / sigma_val
      inv_s <- 1 / sigma_val
      inv_s2 <- inv_s * inv_s
      dwr <- inv_s * (Jp - wr * Js)
      value <- value + wr^2 + log(2 * pi * sigma_val^2)
      grad <- grad + 2 * (wr * dwr + inv_s * Js)
      for (k1 in seq_len(K)) for (k2 in seq_len(K)) {
        H_GN[k1, k2] <- H_GN[k1, k2] +
          2 * dwr[k1] * dwr[k2] +
          (-2 * wr * inv_s2) * (Jp[k1] * Js[k2] + Js[k1] * Jp[k2]) +
          (4 * wr^2 * inv_s2 - 2 * inv_s2) * Js[k1] * Js[k2]
      }
    }
    eta_q <- as.numeric(t(eta_block) %*% Omega_inv %*% eta_block)
    value <- unname(value) + eta_q + Omega_log_det
    grad <- grad + 2 * (Omega_inv %*% eta_block)
    H_GN <- H_GN + 2 * Omega_inv
    list(value = value, grad = as.numeric(grad), hess = unname(H_GN))
  }

  for (i in seq_along(subjects)) {
    s <- subjects[i]
    eta_nms_i <- paste0(c("eta_ka_", "eta_ke_"), s)
    eta_i <- eta_vals[eta_nms_i]
    R <- fast_eval_R(eta_i, s)

    expect_equal(unname(res$value[[s]]), R$value, tolerance = 1e-10,
                 info = paste("value, subject", s))
    expect_equal(as.numeric(res$gradient[s, ]), R$grad, tolerance = 1e-10,
                 info = paste("gradient, subject", s))
    expect_equal(as.numeric(res$H_GN[[s]]),
                 as.numeric(R$hess), tolerance = 1e-10,
                 info = paste("H_GN, subject", s))
  }
})
