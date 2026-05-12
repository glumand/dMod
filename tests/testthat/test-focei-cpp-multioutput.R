context("focei-cpp-multioutput")

# Multi-output fast inner: data has multiple observables per subject, with
# distinct sigma per observable. Long-format meta carries per-row indices
# into the model and err deriv arrays; the kernel loops over rows, not over
# time x observable.

test_that("fast inner: value/gradient/H_GN match R oracle for a 2-output model", {
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
  # Additive errors, distinct sigma per observable.
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
