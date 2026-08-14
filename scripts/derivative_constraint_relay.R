## ---------------------------------------------------------------------------
## Derivative (curvature) regularization on a NON-OSCILLATING model that is
## deliberately (but moderately) underdetermined by its data -- companion to
## scripts/derivative_constraint.R, which used the (feedback-driven, mildly
## oscillation-prone) MAPK cascade.
##
## Model: a purely FEEDFORWARD 15-state phospho-relay (receptor pulse -> 4
## kinase-style tiers -> transcription factor -> nuclear translocation).
## There is no feedback loop anywhere -- every reaction only regulates
## something downstream of itself -- so this model structurally cannot
## oscillate or ring; each tier's response is a single decaying hump.
##
## BALANCED for stability/identifiability (revised from an earlier all-MM,
## single-observable version that had condition numbers up to 1e30-1e40 and
## made even `profile(method = "optimize")` impractically slow -- see the
## conversation record):
##  - Michaelis-Menten (saturating) kinetics are kept at exactly the two
##    "bookend" reactions -- S1's activation by RecA, and TF's activation by
##    S4PP -- where a genuine saturating nonlinearity is the point. Every
##    OTHER activation step (S2, S3, S4) uses plain bimolecular mass action
##    (no Km, never saturates), and every deactivation/dephosphorylation
##    step is plain first-order decay (no Km). This removes the compounding
##    near-singular denominators that made the all-MM cascade so
##    pathologically ill-conditioned, while keeping two real saturating
##    nonlinearities where a "sensor" and an "effector" stage plausibly
##    warrant one. Parameter count also drops from 27 to 21 (8 Km's removed).
##  - THREE observables now, not one: obsS2 = scale_S2*S2P (early-mid),
##    obsS4 = scale_S4*S4PP (late), and obsTF = scale_TF*TFPnuc (final).
##    This leaves only the S3 tier (S3/S3P/S3PP) genuinely unobserved,
##    sandwiched between two flanking measurements -- a smaller, more
##    tractable "blind spot" than the previous 4-tier one, but still a real
##    one: the sensitivity condition number at a reference fit is still
##    around 1e7-1e9 (vs. 1e30-1e40 before) -- much more tractable, but far
##    from perfectly identifiable, so there is still a genuine role for the
##    curvature penalty to play on the S3 tier specifically.
##
## Because this is simulated data, the TRUE trajectory is known, which lets
## this script directly quantify recovery error (RMSE vs. truth) and
## roughness (vs. truth) for baseline and curvature-regularized fits alike
## -- something you would not have with real data, but which is exactly
## what demonstrates the regularization is fixing a real problem, not just
## trading likelihood for an arbitrary-seeming penalty.
##
## Everything else -- the curvatureConstraint()/secondDiffTime() machinery,
## the dense-grid-fitting-via-g_full design, self-normalizing per-quantity
## weights -- is copied verbatim from derivative_constraint.R. See that
## file's header for the full derivation.
## ---------------------------------------------------------------------------

.dmod_find_repo_root <- function(start) {
  here <- normalizePath(start, mustWork = FALSE, winslash = "/")
  for (i in seq_len(8)) {
    desc <- file.path(here, "DESCRIPTION")
    if (file.exists(desc)) {
      pkg <- tryCatch(unname(read.dcf(desc, "Package")[1, 1]), error = function(e) NA_character_)
      if (identical(pkg, "dMod")) return(here)
    }
    parent <- dirname(here)
    if (parent == here) break
    here <- parent
  }
  NA_character_
}

KNOWN_REPO_ROOT <- "/home/mio/Phd/dMod"

repo_root <- local({
  full <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", full, value = TRUE)
  ofile <- if (!is.null(sys.frames()) && length(sys.frames()) >= 1L) sys.frame(1)$ofile else NULL
  candidates <- c(
    if (length(hit)) dirname(sub("^--file=", "", hit[1])),
    if (!is.null(ofile)) dirname(ofile),
    getwd(),
    KNOWN_REPO_ROOT
  )
  for (cand in candidates) {
    if (is.null(cand) || is.na(cand)) next
    found <- .dmod_find_repo_root(cand)
    if (!is.na(found)) return(found)
  }
  stop("Could not locate the dMod repo root (no DESCRIPTION with Package: dMod found ",
       "walking up from any of: ", paste(candidates, collapse = ", "), "). ",
       "Run this script from inside the dMod checkout, or edit KNOWN_REPO_ROOT above.")
})

repo_lib <- file.path(repo_root, ".Rlib")
dir.create(repo_lib, showWarnings = FALSE, recursive = TRUE)

desc_version <- unname(read.dcf(file.path(repo_root, "DESCRIPTION"), "Version")[1, 1])
installed_version <- tryCatch(as.character(packageVersion("dMod", lib.loc = repo_lib)),
                               error = function(e) NA_character_)
if (!identical(installed_version, desc_version)) {
  cat("Building dMod", desc_version, "from", repo_root, "into", repo_lib, "...\n")
  status <- system2(file.path(R.home("bin"), "R"),
                     c("CMD", "INSTALL", "--no-multiarch", "--with-keep.source",
                       "-l", shQuote(repo_lib), shQuote(repo_root)))
  if (status != 0)
    stop("R CMD INSTALL failed -- install dependencies first (see DESCRIPTION Imports/Remotes).")
}

if ("package:dMod" %in% search()) detach("package:dMod", unload = TRUE)

library(dMod, lib.loc = repo_lib)
stopifnot(identical(find.package("dMod"), file.path(repo_lib, "dMod")))
cat("Using dMod", as.character(packageVersion("dMod")),
    "from", find.package("dMod"), "\n\n")

set.seed(1)

plot_dir <- file.path(repo_root, "scripts", "plots", "relay_derivative_constraint")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

setwd(tempdir())

## --- 1. Model: 15-state feedforward phospho-relay, NO feedback -------------
## Receptor pulse (RecA, decaying) drives S1, which drives S2, which drives
## a double-phospho S3 tier, which drives a double-phospho S4 tier, which
## drives TF -> TFP -> nuclear TFPnuc. Every arrow points strictly
## downstream -- nothing here can create a loop, hence no oscillation.
## MM kinetics only at the two bookends (S1, TF activation); everything
## else is mass action / first-order (see header for the rationale).
f <- NULL
f <- addReaction(f, "RecA", "Rec", "kdeact_Rec*RecA", "receptor deactivation")

f <- addReaction(f, "S1",  "S1P", "kcat_S1*RecA*S1/(Km_S1+S1)", "S1 activation (MM)")
f <- addReaction(f, "S1P", "S1",  "kcat_dS1*S1P", "S1 deactivation (linear)")

f <- addReaction(f, "S2",  "S2P", "kcat_S2*S1P*S2", "S2 activation (mass action)")
f <- addReaction(f, "S2P", "S2",  "kcat_dS2*S2P", "S2 deactivation (linear)")

f <- addReaction(f, "S3",   "S3P",  "kcat_S3a*S2P*S3", "S3 phosphorylation 1 (mass action)")
f <- addReaction(f, "S3P",  "S3PP", "kcat_S3b*S2P*S3P", "S3 phosphorylation 2 (mass action)")
f <- addReaction(f, "S3PP", "S3P",  "Vmax_dS3*S3PP", "S3 dephosphorylation 2 (linear)")
f <- addReaction(f, "S3P",  "S3",   "Vmax_dS3*S3P", "S3 dephosphorylation 1 (linear)")

f <- addReaction(f, "S4",   "S4P",  "kcat_S4a*S3PP*S4", "S4 phosphorylation 1 (mass action)")
f <- addReaction(f, "S4P",  "S4PP", "kcat_S4b*S3PP*S4P", "S4 phosphorylation 2 (mass action)")
f <- addReaction(f, "S4PP", "S4P",  "Vmax_dS4*S4PP", "S4 dephosphorylation 2 (linear)")
f <- addReaction(f, "S4P",  "S4",   "Vmax_dS4*S4P", "S4 dephosphorylation 1 (linear)")

f <- addReaction(f, "TF",     "TFP",    "kcat_TF*S4PP*TF/(Km_TF+TF)", "TF phosphorylation (MM)")
f <- addReaction(f, "TFP",    "TF",     "Vmax_dTF*TFP", "TF dephosphorylation (linear)")
f <- addReaction(f, "TFP",    "TFPnuc", "ktrans*TFP", "TF nuclear translocation")
f <- addReaction(f, "TFPnuc", "TF",     "kdeg_nuc*TFPnuc", "nuclear TF turnover")

model <- odemodel(f, modelname = "toyRelayModel", solver = "CppODE")
x <- Xs(model)

## --- 2. Observation: three readouts, leaving only the S3 tier unobserved --
## obsS2 and obsS4 flank the S3 tier, which stays the genuinely
## under-constrained part of the model. `g` (attach.input = FALSE) generates
## the data; `g_full` (attach.input = TRUE) is used for fitting so the
## curvature penalty can see every state.
observables <- eqnvec(obsS2 = "scale_S2*S2P", obsS4 = "scale_S4*S4PP", obsTF = "scale_TF*TFPnuc")
g      <- Y(observables, x, compile = TRUE, modelname = "toyRelayObs",     attach.input = FALSE)
g_full <- Y(observables, x, compile = TRUE, modelname = "toyRelayObsFull", attach.input = TRUE)

innerpars <- getParameters(g, x)

## --- 3. Parameter transformation: log-parametrization ----------------------
states <- f$states
trafo <- repar("x~x", x = innerpars)
trafo <- repar("x~0", x = c("S1P","S2P","S3P","S3PP","S4P","S4PP","TFP","TFPnuc"), trafo)
trafo <- repar("x~1", x = c("S1","S2","S3","S4","TF"), trafo)
trafo <- repar("x~dose", x = "RecA", trafo)
trafo <- repar("x~1-dose", x = "Rec", trafo)
ratepars <- setdiff(innerpars, states)
trafo <- repar("x~exp(x)", x = c(ratepars, "dose"), trafo)

p <- P(trafo, condition = "sim", compile = TRUE, modelname = "toyRelayTrafo")

## --- 4. "True" parameters and simulated data -------------------------------
## Smooth/graded kinetics throughout (moderate Km's at the two remaining MM
## reactions, no near-zero Km that would make either tier switch-like by
## construction) -- any sharp edges seen in an unregularized FIT are
## therefore an estimation artifact, not a feature of the true dynamics.
true_pars_inner <- c(
  dose = 0.8,
  kdeact_Rec = 0.03,
  kcat_S1 = 2, Km_S1 = 0.3, kcat_dS1 = 2,
  kcat_S2 = 2, kcat_dS2 = 2,
  kcat_S3a = 2, kcat_S3b = 2, Vmax_dS3 = 1,
  kcat_S4a = 2, kcat_S4b = 2, Vmax_dS4 = 1,
  kcat_TF = 2, Km_TF = 0.3, Vmax_dTF = 0.5,
  ktrans = 0.3, kdeg_nuc = 0.1,
  ## Per-observable scales, tuned so all three peak at a comparable O(1)
  ## amplitude relative to sigma below (checked empirically -- see
  ## conversation record for the earlier scale_TF=1/sigma-swamps-signal bug
  ## this generalizes the fix for).
  scale_S2 = 1, scale_S4 = 3, scale_TF = 1
)
stopifnot(setequal(names(true_pars_inner), getParameters(p)))
pouter_true <- log(true_pars_inner)

## Sparse (every 5 time units) and moderately noisy -- still deliberately
## underdetermined for the unobserved S3 tier, even with three observables
## and a much better-conditioned model than the earlier all-MM version.
times_data <- seq(0, 80, by = 5)
prediction_true <- (g * x * p)(times_data, pouter_true)

sigma <- 0.03
sim <- as.data.frame(prediction_true)
sim$value <- sim$value + rnorm(nrow(sim), sd = sigma)
sim$sigma <- sigma
data <- as.datalist(sim, split.by = "condition")

## --- 5. Fitting objective on a dense grid, via g_full ----------------------
times_fit <- seq(0, 80, by = 0.1)
obj <- normL2(data, g_full * x * p, times = times_fit)

## Ground truth on the fitting grid -- ONLY used for the diagnostic
## RMSE/roughness-vs-truth comparison in Part 9; never used by either fit.
true_full <- (g_full * x * p)(times_fit, pouter_true)[[1]]
cols_all  <- setdiff(colnames(true_full), "time")
hstep     <- diff(times_fit)[1]

roughnessAll <- function(pars) {
  pred <- (g_full * x * p)(times_fit, pars)[[1]][, cols_all, drop = FALSE]
  n <- nrow(pred)
  curv <- (pred[3:n, , drop = FALSE] - 2*pred[2:(n-1), , drop = FALSE] + pred[1:(n-2), , drop = FALSE]) / hstep^2
  sqrt(colMeans(curv^2))
}
rmseVsTruth <- function(pars) {
  pred <- (g_full * x * p)(times_fit, pars)[[1]][, cols_all, drop = FALSE]
  sqrt(colMeans((pred - true_full[, cols_all, drop = FALSE])^2))
}

## Keep only the `n` best (lowest-value) fits before plotting a waterfall, so
## that a single very bad fit (non-convergence, blown-up objective, ...)
## doesn't stretch the y-axis and squash everything else into an unreadable
## sliver. Dropping the point entirely (rather than clipping the axis) keeps
## the plotted range driven only by the fits actually shown. See
## fim_l2_constraintL2.R, where this was introduced.
topByValue <- function(pf, n = 400) {
  n_keep <- min(n, nrow(pf))
  pf[order(pf$value)[seq_len(n_keep)], ]
}

## --- 6. Baseline: multistart of 100 ----------------------------------------
n_cores <- max(1, min(15, parallel::detectCores() - 1))
n_fits  <- 100

pouter_init <- setNames(rep(-1, length(pouter_true)), names(pouter_true))
starts <- msParframe(pouter_init, n = n_fits, seed = 21, samplefun = runif, min = -3, max = 3)

cat("=== Baseline: multistart run (", n_cores, "cores ) ===\n")
fits <- mstrust(obj, starts, studyname = "relayDerivBaseline", rinit = 1, rmax = 10,
                 cores = n_cores, iterlim = 400)
pf <- as.parframe(fits)
cat("Converged:", sum(pf$converged), "of", nrow(pf), "fits.\n")

waterfall_base <- plotValues(topByValue(pf), tol = 0.1) + ggplot2::ggtitle("baseline (unregularized)")
ggplot2::ggsave(file.path(plot_dir, "relay_derivative_constraint_baseline_waterfall.pdf"),
                waterfall_base, width = 7, height = 4, dpi = 120)

## Best among CONVERGED fits only (an unconverged run can report a lower
## `value` mid-trajectory without that value being a real optimum).
pf_conv <- pf[pf$converged, ]
best <- as.parvec(pf_conv, 1)
cat("Best (converged) objective: truth =", obj(pouter_true)$value, " fit =", pf_conv$value[1], "\n\n")

## --- 7. Curvature (second-derivative) constraint ---------------------------
## Verbatim from derivative_constraint.R.
secondDiffTime <- function(arr, h) {
  n <- dim(arr)[1]
  lo <- 1:(n - 2); mid <- 2:(n - 1); hi <- 3:n
  if (length(dim(arr)) == 2) {
    (arr[hi, , drop = FALSE] - 2 * arr[mid, , drop = FALSE] + arr[lo, , drop = FALSE]) / h^2
  } else {
    (arr[hi, , , drop = FALSE] - 2 * arr[mid, , , drop = FALSE] + arr[lo, , , drop = FALSE]) / h^2
  }
}

curvatureConstraint <- function(base_obj, best, condition,
                                 rel_strength = 0.1, eps_curv_frac = 1e-2) {

  out_best  <- base_obj(best, deriv = TRUE)
  pred_best <- attr(out_best, "env")$prediction[[condition]]

  time_grid <- pred_best[, "time"]
  hs <- diff(time_grid)
  stopifnot("curvatureConstraint() assumes a uniform time grid" =
              max(abs(hs - hs[1])) < 1e-8 * hs[1])
  hs <- hs[1]

  cols <- dimnames(attr(pred_best, "deriv"))[[2]]
  Y_best <- pred_best[, cols, drop = FALSE]
  curv_best <- secondDiffTime(Y_best, hs)

  ref <- colMeans(curv_best^2)
  ## eps_curv_frac = 1e-2 (not derivative_constraint.R's 1e-6): found
  ## empirically that with THIS model's mix of MM and mass-action tiers,
  ## max(ref) itself is small (~2e-3, since the model is smoother overall
  ## by design) and some quantities' ref is many orders of magnitude
  ## smaller still (e.g. S4PP ~1.5e-9) -- a proportional floor anchored to
  ## max(ref) with the smaller fraction was too small in absolute terms to
  ## matter, letting that quantity's weight reach ~2.7e8 (vs ~450 for a
  ## typical quantity) and dominate the whole penalty: rel_strength=0.1
  ## caused the raw likelihood to blow up from -235 to +687 before this
  ## was raised. The weight cap below is a secondary safety net.
  eps_curv <- eps_curv_frac * max(ref)
  weight <- 1 / (ref + eps_curv)
  weight <- pmin(weight, 100 * stats::median(weight))

  n_cols <- length(cols)
  n_t    <- nrow(curv_best)
  strength_curv <- rel_strength * abs(out_best$value)
  C <- strength_curv / (n_cols * n_t)

  cat("=== curvatureConstraint(): per-quantity reference curvature at `best` ===\n")
  print(data.frame(quantity = cols, ref_curvature = ref, weight = weight))
  cat("strength_curv =", strength_curv, " (", rel_strength * 100, "% of |obj(best)$value| =", abs(out_best$value), ")\n\n")

  wrapped <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, env = NULL) {

    base_out <- base_obj(..., fixed = fixed, deriv = deriv, deriv2 = deriv2, env = env)
    penv <- attr(base_out, "env")
    pred <- penv$prediction[[condition]]

    Y    <- pred[, cols, drop = FALSE]
    curv <- secondDiffTime(Y, hs)

    pen_value <- C * sum(vapply(seq_len(n_cols), function(j) weight[j] * sum(curv[, j]^2), numeric(1)))

    pen_gradient <- NULL
    pen_hessian  <- NULL

    if (deriv) {
      dY    <- attr(pred, "deriv")[, cols, , drop = FALSE]
      dcurv <- secondDiffTime(dY, hs)
      theta <- dimnames(dcurv)[[3]]

      pen_gradient <- setNames(numeric(length(theta)), theta)
      pen_hessian  <- matrix(0, length(theta), length(theta), dimnames = list(theta, theta))

      for (j in seq_len(n_cols)) {
        Jc <- dcurv[, j, ]
        dim(Jc) <- c(dim(dcurv)[1], length(theta))
        pen_gradient <- pen_gradient + 2 * C * weight[j] * as.vector(curv[, j] %*% Jc)
        pen_hessian  <- pen_hessian  + 2 * C * weight[j] * (t(Jc) %*% Jc)
      }

      stopifnot(all(names(base_out$gradient) %in% theta))
      pen_gradient <- pen_gradient[names(base_out$gradient)]
      pen_hessian  <- pen_hessian[names(base_out$gradient), names(base_out$gradient)]
    }

    combined <- objlist(
      value    = base_out$value + pen_value,
      gradient = if (deriv) base_out$gradient + pen_gradient else NULL,
      hessian  = if (deriv) base_out$hessian  + pen_hessian  else NULL
    )
    attr(combined, "env") <- penv
    attr(combined, "curvature_penalty") <- pen_value
    combined
  }

  class(wrapped) <- class(base_obj)
  attr(wrapped, "conditions") <- attr(base_obj, "conditions")
  attr(wrapped, "parameters") <- attr(base_obj, "parameters")
  wrapped
}

obj_curv <- curvatureConstraint(obj, best, condition = "sim", rel_strength = 0.1)

## --- 8. Curvature-regularized fit: multistart, same starts -----------------
cat("=== Curvature-regularized: multistart run (", n_cores, "cores ) ===\n")
fits_curv <- mstrust(obj_curv, starts, studyname = "relayDerivCurvature", rinit = 1, rmax = 10,
                      cores = n_cores, iterlim = 400)
pf_curv <- as.parframe(fits_curv)
cat("Converged:", sum(pf_curv$converged), "of", nrow(pf_curv), "fits.\n\n")

waterfall_curv <- plotValues(topByValue(pf_curv), tol = 0.1) + ggplot2::ggtitle("curvature-regularized")
ggplot2::ggsave(file.path(plot_dir, "relay_derivative_constraint_curvature_waterfall.pdf"),
                waterfall_curv, width = 7, height = 4, dpi = 120)

pf_curv_conv <- pf_curv[pf_curv$converged, ]
best_curv <- as.parvec(pf_curv_conv, 1)

## --- 9. Compare against ground truth ---------------------------------------
pen_at_best      <- attr(obj_curv(best), "curvature_penalty")
pen_at_best_curv <- attr(obj_curv(best_curv), "curvature_penalty")

cat("=== Summary: baseline vs. curvature-regularized (raw likelihood) ===\n")
cat(sprintf("%-22s raw likelihood = %10.3f   curvature penalty = %10.3f\n",
            "baseline (best)", obj(best)$value, pen_at_best))
cat(sprintf("%-22s raw likelihood = %10.3f   curvature penalty = %10.3f\n",
            "curvature-reg (best)", obj(best_curv)$value, pen_at_best_curv))

cat("\n=== Recovery of the TRUE trajectory (only possible because this is\n",
    " simulated data -- neither fit ever sees `pouter_true`) ===\n")
rmse_base <- rmseVsTruth(best);      rough_base <- roughnessAll(best)
rmse_curv <- rmseVsTruth(best_curv); rough_curv <- roughnessAll(best_curv)
rough_true <- roughnessAll(pouter_true)

cmp <- data.frame(
  quantity           = cols_all,
  rmse_baseline       = rmse_base,
  rmse_curvature_reg  = rmse_curv,
  roughness_baseline_over_truth      = rough_base / rough_true,
  roughness_curvature_reg_over_truth = rough_curv / rough_true
)
print(cmp, row.names = FALSE)
cat("\n(roughness ratios near 1 mean the fit's curvature matches the TRUE\n",
    " trajectory's; ratios >> 1 mean the fit is spuriously rougher than truth.)\n\n")

## --- 10. Plots --------------------------------------------------------------
pred_time <- seq(0, 80, by = 0.5)
pred_obs <- as.prdlist(list(
  truth     = (g * x * p)(pred_time, pouter_true)[[1]],
  baseline  = (g * x * p)(pred_time, best)[[1]],
  curvature = (g * x * p)(pred_time, best_curv)[[1]]
))
pl_obs <- plotCombined(pred_obs, data)
ggplot2::ggsave(file.path(plot_dir, "relay_derivative_constraint_observable.pdf"),
                pl_obs, width = 8, height = 4, dpi = 120)

pred_full <- as.prdlist(list(
  truth     = (g_full * x * p)(pred_time, pouter_true)[[1]],
  baseline  = (g_full * x * p)(pred_time, best)[[1]],
  curvature = (g_full * x * p)(pred_time, best_curv)[[1]]
))
pl_full <- plotCombined(pred_full, data)
ggplot2::ggsave(file.path(plot_dir, "relay_derivative_constraint_all_states.pdf"),
                pl_full, width = 12, height = 10, dpi = 120)

cat("Plots saved to", plot_dir, "\n")

## --- 11. Strength sweep: rel_strength across a log-grid --------------------
## `rel_strength = 0.1` (Part 7/8) was fixed by hand. Here the same
## construction is repeated across a log-grid, so the tradeoff between
## likelihood cost and curvature suppression becomes visible directly. 0.1
## stays in the grid so the Part 7/8 run above is reused, not refit.
cat("\n=== Strength sweep: curvature regularization across rel_strength ===\n")

rel_strengths <- c(0.01, 0.1, 1)

fitCurvatureRegularized <- function(rel_strength_i, file_stub) {
  obj_curv_i <- curvatureConstraint(obj, best, condition = "sim", rel_strength = rel_strength_i)
  fits_i <- mstrust(obj_curv_i, starts, studyname = paste0("relayDerivCurv_str", rel_strength_i),
                     rinit = 1, rmax = 10, cores = n_cores, iterlim = 400)
  pf_i <- as.parframe(fits_i)
  cat("  rel_strength =", rel_strength_i, ": converged", sum(pf_i$converged), "of", nrow(pf_i), "\n")

  waterfall_i <- plotValues(topByValue(pf_i), tol = 0.1) +
    ggplot2::ggtitle(paste0("curvature-regularized, rel_strength=", rel_strength_i))
  ggplot2::ggsave(file.path(plot_dir, paste0(file_stub, "_waterfall.pdf")),
                   waterfall_i, width = 7, height = 4, dpi = 120)

  list(obj_curv = obj_curv_i, pf = pf_i)
}

sweep_curv <- list()
for (rel_strength_i in rel_strengths) {
  key <- as.character(rel_strength_i)
  if (isTRUE(all.equal(rel_strength_i, 0.1))) {
    cat("--- rel_strength =", rel_strength_i, "(reusing Part 7/8 run) ---\n")
    sweep_curv[[key]] <- list(obj_curv = obj_curv, pf = pf_curv)
  } else {
    sweep_curv[[key]] <- fitCurvatureRegularized(
      rel_strength_i, sprintf("relay_derivative_constraint_str%g", rel_strength_i))
  }
  sweep_curv[[key]]$pf$rel_strength <- rel_strength_i
}

## --- 11a. Combined waterfall: baseline + all five strengths, one plot -----
combined <- rbind(
  data.frame(value = pf$value, converged = pf$converged, group = "baseline"),
  do.call(rbind, lapply(sweep_curv, function(s)
    data.frame(value = s$pf$value, converged = s$pf$converged, group = as.character(s$pf$rel_strength))))
)
combined$group <- factor(combined$group, levels = c("baseline", as.character(rel_strengths)))
pooled <- topByValue(combined, n = 150)
pooled$index <- seq_len(nrow(pooled))

waterfall_sweep <- ggplot2::ggplot(pooled, ggplot2::aes(x = index, y = value, color = group, pch = converged)) +
  ggplot2::geom_point() +
  ggplot2::xlab("index (baseline + all strengths pooled, best 400 of 600 by value)") + ggplot2::ylab("value") +
  ggplot2::ggtitle("Baseline vs. curvature regularization: strength sweep") +
  theme_dMod()
ggplot2::ggsave(file.path(plot_dir, "relay_derivative_constraint_strength_sweep_waterfall.pdf"),
                waterfall_sweep, width = 8, height = 4.5, dpi = 120)

## --- 11b. Combined fit plots: truth, baseline, best at each strength ------
bestByStrength <- function(sweep_list) {
  setNames(lapply(sweep_list, function(s) {
    pf_conv_i <- s$pf[s$pf$converged, ]
    as.parvec(pf_conv_i, 1)
  }), names(sweep_list))
}
best_by_strength <- bestByStrength(sweep_curv)

pred_obs_sweep <- as.prdlist(c(
  list(truth = (g * x * p)(pred_time, pouter_true)[[1]],
       baseline = (g * x * p)(pred_time, best)[[1]]),
  setNames(lapply(best_by_strength, function(b) (g * x * p)(pred_time, b)[[1]]),
           paste0("str", names(best_by_strength)))
))
pl_obs_sweep <- plotCombined(pred_obs_sweep, data)
ggplot2::ggsave(file.path(plot_dir, "relay_derivative_constraint_strength_sweep_observable.pdf"),
                pl_obs_sweep, width = 9, height = 4.5, dpi = 120)

pred_full_sweep <- as.prdlist(c(
  list(truth = (g_full * x * p)(pred_time, pouter_true)[[1]],
       baseline = (g_full * x * p)(pred_time, best)[[1]]),
  setNames(lapply(best_by_strength, function(b) (g_full * x * p)(pred_time, b)[[1]]),
           paste0("str", names(best_by_strength)))
))
pl_full_sweep <- plotCombined(pred_full_sweep, data)
ggplot2::ggsave(file.path(plot_dir, "relay_derivative_constraint_strength_sweep_all_states.pdf"),
                pl_full_sweep, width = 13, height = 11, dpi = 120)

## --- 11c. Summary table across strengths -----------------------------------
cat("\n=== Strength-sweep summary: raw likelihood, curvature penalty, and\n",
    " recovery of S3P (in the unobserved S3 tier, flanked by obsS2/obsS4) ===\n")
for (rel_strength_i in rel_strengths) {
  key  <- as.character(rel_strength_i)
  s    <- sweep_curv[[key]]
  b_i  <- best_by_strength[[key]]
  pen_i   <- attr(s$obj_curv(b_i), "curvature_penalty")
  rmse_i  <- rmseVsTruth(b_i)["S3P"]
  rough_i <- roughnessAll(b_i)["S3P"] / rough_true["S3P"]
  cat(sprintf("rel_strength=%-6g  raw likelihood=%9.3f  curvature penalty=%9.3f  rmse(S3P)=%7.4f  rough(S3P)/truth=%7.3f\n",
              rel_strength_i, obj(b_i)$value, pen_i, rmse_i, rough_i))
}

cat("\nStrength-sweep plots saved to", plot_dir, "\n")

## --- 12. Profile likelihood: baseline + all five strengths -----------------
## For each arm, profile EVERY parameter around that arm's OWN best fit,
## using that arm's OWN objective -- so a regularized profile's shape
## includes the curvature penalty's own contribution, showing directly
## whether/how the penalty tightens a parameter's profile, not just how the
## raw likelihood alone behaves at a different point.
##
## method = "optimize": re-optimizes all other parameters at every profile
## step, rather than Euler-integrating along the profile path
## (method = "integrate"). More expensive per step, but far more robust --
## "integrate" repeatedly collapsed its stepsize and threw solver errors
## against this model's extreme ill-conditioning (condition numbers up to
## 1e30-1e40, see the per-fit diagnostics), making baseline alone take
## ~41 minutes for 27 parameters with 14+ non-fatal solver failures.
## `cores` parallelizes across parameters (up to n_cores at a time).
##
## `profile()` attaches whatever numeric attributes the objective returns
## (e.g. normL2's own "data" attribute for the baseline; this script's
## "curvature_penalty" attribute for the regularized arms) as extra
## decomposed "mode" lines. Those two sets differ across arms, so plotting
## is filtered to `mode == "total"` -- the actual value each fit optimized
## -- for a clean, consistent comparison across all six.
cat("\n=== Profile likelihood: baseline + all five strengths ===\n")

profile_list <- list()

cat("--- profiling: baseline ---\n")
profile_list[["baseline"]] <- profile(obj, best, whichPar = names(best),
                                       method = "optimize", cores = n_cores)

for (rel_strength_i in rel_strengths) {
  key <- as.character(rel_strength_i)
  s   <- sweep_curv[[key]]
  b_i <- best_by_strength[[key]]
  cat("--- profiling: rel_strength =", rel_strength_i, "---\n")
  profile_list[[key]] <- profile(s$obj_curv, b_i, whichPar = names(b_i),
                                  method = "optimize", cores = n_cores)
}

p_profiles <- plotProfile(profile_list, mode == "total") +
  ggplot2::facet_wrap(~name, scales = "free_x", ncol = 5) +
  theme_dMod()
ggplot2::ggsave(file.path(plot_dir, "relay_derivative_constraint_profiles.pdf"),
                p_profiles, width = 20, height = 22, dpi = 120)

cat("Profile plot saved to", plot_dir, "\n")
