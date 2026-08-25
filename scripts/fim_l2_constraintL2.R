## ---------------------------------------------------------------------------
## Variant of scripts/example_toy_mapk_fim_fit.R: same model, same data, same
## baseline/multistart setup (Parts 0-9 below are unchanged from that file),
## but Part 10's two regularization strategies are rebuilt on top of dMod's
## OWN constraintL2() instead of the hand-written quadraticPenalty()/
## quadraticPenalty(mu, Prec) helper.
##
## constraintL2(mu, sigma) only supports a DIAGONAL penalty --
## sum_i (p_i - mu_i)^2 / sigma_i^2 -- one sigma per parameter, no
## off-diagonal terms. That's exactly what's needed for the isotropic L2
## comparison arm (Prec_l2 = w_iso * I was already diagonal), but the
## FIM-anisotropic arm's Prec_fim = V diag(w) V^T is a full, rotated matrix
## with cross-parameter terms baked in from the reference Hessian's
## eigenvectors -- constraintL2 has no argument for that shape of prior.
##
## The bridge used below: extract just the DIAGONAL of Prec_fim,
## diag(Prec_fim)_i = sum_k V_ik^2 * w_k, and feed that in as one sigma per
## parameter (sigma_i = 1/sqrt(diag(Prec_fim)_i)). This is not an
## approximation pulled out of thin air -- for a quadratic form d^T P d,
## the SECOND DERIVATIVE with respect to a single p_i, holding every other
## parameter fixed, is exactly 2*P_ii; constraintL2's penalty gives exactly
## 2/sigma_i^2 for that same derivative. Matching sigma_i to diag(Prec_fim)_i
## therefore reproduces the FIM penalty's exact per-parameter curvature at
## `best`, while necessarily dropping the cross-parameter (off-diagonal)
## coupling that made the original penalty direction-aware ACROSS parameters,
## not just per-parameter. In other words: this keeps "how strongly does the
## FIM want to pull THIS parameter back", loses "...pull it back along with
## THAT other parameter specifically".
##
## Everything else -- model, simulated data, baseline fit/Hessian (Part
## B/6/7), multistart starts (Part C), and the fitAndPlotRegularized() driver
## -- is copied verbatim from example_toy_mapk_fim_fit.R so the two scripts'
## numbers are directly comparable. Part 10b's self-bootstrapping
## continuation regularization is NOT reproduced here -- this file only
## covers the two constructions Part 10 there compares (FIM vs. plain L2).
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

set.seed(5)

## Separate plot folder from example_toy_mapk_fim_fit.R's, even though the
## model/data/baseline are identical -- keeps the constraintL2-based
## regularization plots from overwriting (or being confused with) that
## script's quadraticPenalty()-based ones.
plot_dir <- file.path(repo_root, "scripts", "plots", "mapk_fim_fit_constraintL2")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

setwd(tempdir())

## --- 0. Waterfall diagnostics: condition number at first vs. last step ----
## Identical to example_toy_mapk_fim_fit.R -- see that file for the full
## derivation of these helpers.
sensitivityConditionNumber <- function(objfun, pars) {
  out <- objfun(pars)
  env <- attr(out, "env")
  prediction <- if (!is.null(env)) env$prediction else NULL
  if (is.null(prediction)) return(NA_real_)
  cn <- vapply(prediction, function(pred) {
    S <- attr(pred, "deriv")
    if (is.null(S)) return(NA_real_)
    J <- S
    dim(J) <- c(dim(S)[1] * dim(S)[2], dim(S)[3])
    sv <- svd(J, nu = 0, nv = 0)$d
    max(sv) / min(sv)
  }, numeric(1))
  max(cn, na.rm = TRUE)
}

sensitivitySingularValues <- function(objfun, pars) {
  out <- objfun(pars)
  env <- attr(out, "env")
  prediction <- if (!is.null(env)) env$prediction else NULL
  if (is.null(prediction)) return(NA_real_)
  unlist(lapply(prediction, function(pred) {
    S <- attr(pred, "deriv")
    if (is.null(S)) return(NULL)
    J <- S
    dim(J) <- c(dim(S)[1] * dim(S)[2], dim(S)[3])
    svd(J, nu = 0, nv = 0)$d
  }), use.names = FALSE)
}

## Keep only the `n` best (lowest-value) fits before plotting a waterfall, so
## that a single very bad fit (non-convergence, blown-up objective, ...)
## doesn't stretch the y-axis and squash everything else into an unreadable
## sliver. Dropping the point entirely (rather than clipping the axis) keeps
## the plotted range driven only by the fits actually shown.
topByValue <- function(pf, n = 400) {
  n_keep <- min(n, nrow(pf))
  pf[order(pf$value)[seq_len(n_keep)], ]
}

plotWaterfallColored <- function(pf, color_by, color_label, tol = 1, reverse = FALSE) {
  jumps <- dMod:::stepDetect(pf$value, tol)
  y.range <- c(min(pf$value), max(max(pf$value), min(pf$value) + tol))
  y.jumps <- seq(y.range[2], y.range[1], length.out = length(jumps))

  ord <- order(pf$value)
  pars <- pf[ord, ]
  pars$index <- seq_len(nrow(pars))
  pars$.color <- color_by[ord]

  col_low  <- if (reverse) "orange" else "dodgerblue"
  col_high <- if (reverse) "dodgerblue" else "orange"

  ggplot2::ggplot(pars, ggplot2::aes(x = index, y = value, pch = converged, color = .color)) +
    ggplot2::geom_vline(xintercept = jumps, lty = 2) +
    ggplot2::geom_point() +
    ggplot2::annotate("text", x = jumps + 1, y = y.jumps, label = jumps, hjust = 0, color = "firebrick", size = 3) +
    ggplot2::xlab("index") + ggplot2::ylab("value") +
    ggplot2::scale_color_gradient(low = col_low, high = col_high, trans = "log10", name = color_label) +
    ggplot2::coord_cartesian(ylim = y.range) +
    theme_dMod()
}

plotWaterfallDiagnostics <- function(pf, cond_start, cond_end, min_sv_start, tol = 1) {
  p_iter  <- plotValues(pf, tol = tol) +
    ggplot2::ggtitle("colored by iterations")
  p_start <- plotWaterfallColored(pf, cond_start, "condition number\n(first step)", tol = tol) +
    ggplot2::ggtitle("colored by sensitivity at first step")
  p_end   <- plotWaterfallColored(pf, cond_end, "condition number\n(last step)", tol = tol) +
    ggplot2::ggtitle("colored by sensitivity at last step")
  p_minsv <- plotWaterfallColored(pf, min_sv_start, "smallest singular\nvalue (first step)", tol = tol, reverse = TRUE) +
    ggplot2::ggtitle("colored by smallest singular value at first step")
  cowplot::plot_grid(p_iter, p_start, p_end, p_minsv, nrow = 2)
}

## --- 1. Model: 15-state phosphorylation cascade with negative feedback ----
## Identical to example_toy_mapk_fim_fit.R -- see that file's header for the
## full topology/non-identifiability writeup.
f <- NULL

f <- addReaction(f, "RecA", "Rec", "kdeact_Rec*RecA", "receptor deactivation")

f <- addReaction(f, "Raf",  "RafP", "kcat_Raf*RecA*Raf/(Km_Raf+Raf)/(1+InhP/Ki_fb)", "Raf activation")
f <- addReaction(f, "RafP", "Raf",  "kcat_dRaf*Etot_dRaf*RafP/(Km_dRaf+RafP)", "Raf deactivation")

f <- addReaction(f, "MEK",   "MEKP",  "kcat_MEK1*RafP*MEK/(Km_MEK+MEK)", "MEK phosphorylation 1")
f <- addReaction(f, "MEKP",  "MEKPP", "kcat_MEK2*RafP*MEKP/(Km_MEK+MEKP)", "MEK phosphorylation 2")
f <- addReaction(f, "MEKPP", "MEKP",  "Vmax_dMEK*MEKPP/(Km_dMEK+MEKPP)/(1+RafP/Ki_ff)", "MEK dephosphorylation 2")
f <- addReaction(f, "MEKP",  "MEK",   "Vmax_dMEK*MEKP/(Km_dMEK+MEKP)", "MEK dephosphorylation 1")

f <- addReaction(f, "ERK",   "ERKP",  "kcat_ERK1*MEKPP*ERK/(Km_ERK+ERK)", "ERK phosphorylation 1")
f <- addReaction(f, "ERKP",  "ERKPP", "kcat_ERK2*MEKPP*ERKP/(Km_ERK+ERKP)", "ERK phosphorylation 2")
f <- addReaction(f, "ERKPP", "ERKP",  "Vmax_dERK*ERKPP/(Km_dERK+ERKPP)", "ERK dephosphorylation 2")
f <- addReaction(f, "ERKP",  "ERK",   "Vmax_dERK*ERKP/(Km_dERK+ERKP)", "ERK dephosphorylation 1")

f <- addReaction(f, "Inh",  "InhP", "kcat_fb*ERKPP*Inh/(Km_fb+Inh)", "feedback activation")
f <- addReaction(f, "InhP", "Inh",  "kdeact_fb*InhP", "feedback deactivation")

f <- addReaction(f, "TF",     "TFP",    "kcat_TF*ERKPP*TF/(Km_TF+TF)", "TF phosphorylation")
f <- addReaction(f, "TFP",    "TF",     "Vmax_dTF*TFP/(Km_dTF+TFP)", "TF dephosphorylation")
f <- addReaction(f, "TFP",    "TFPnuc", "ktrans*TFP", "TF nuclear translocation")
f <- addReaction(f, "TFPnuc", "TF",     "kdeg_nuc*TFPnuc", "nuclear TF turnover")

model <- odemodel(f, modelname = "toyMapkModel", solver = "CppODE")
x <- Xs(model)

## --- 2. Observation: two arbitrary-unit readouts ---------------------------
observables <- eqnvec(obsERK = "scale_ERK*ERKPP", obsTF = "scale_TF*TFPnuc")
g <- Y(observables, x, compile = TRUE, modelname = "toyMapkObs", attach.input = FALSE)

innerpars <- getParameters(g, x)

## --- 3. Parameter transformation: log-parametrization ----------------------
states <- f$states
trafo <- repar("x~x", x = innerpars)
trafo <- repar("x~0", x = c("RafP", "MEKP", "MEKPP", "ERKP", "ERKPP", "InhP", "TFP", "TFPnuc"), trafo)
trafo <- repar("x~1", x = c("Raf", "MEK", "ERK", "Inh", "TF"), trafo)
trafo <- repar("x~dose", x = "RecA", trafo)
trafo <- repar("x~1-dose", x = "Rec", trafo)
ratepars <- setdiff(innerpars, c(states, "Etot_dRaf"))
trafo <- repar("x~exp(x)", x = c(ratepars, "dose"), trafo)

trafo <- repar("x~0.5", x = "Etot_dRaf", trafo)

p <- P(trafo, condition = "sim", compile = TRUE, modelname = "toyMapkTrafo")

## --- 4. "True" parameters and simulated data -------------------------------
true_pars_inner <- c(
  dose = 0.8,
  kdeact_Rec = 0.03,
  kcat_Raf = 2, Km_Raf = 0.3, Ki_fb = 0.2, kcat_dRaf = 2, Km_dRaf = 0.3,
  kcat_MEK1 = 3, kcat_MEK2 = 3, Km_MEK = 0.3, Vmax_dMEK = 1, Km_dMEK = 0.3, Ki_ff = 0.3,
  kcat_ERK1 = 3, kcat_ERK2 = 3, Km_ERK = 0.3, Vmax_dERK = 0.8, Km_dERK = 0.3,
  kcat_fb = 1.5, Km_fb = 0.3, kdeact_fb = 0.3,
  kcat_TF = 2, Km_TF = 0.3, Vmax_dTF = 0.5, Km_dTF = 0.3,
  ktrans = 0.3, kdeg_nuc = 0.1,
  scale_ERK = 1, scale_TF = 1
)
stopifnot(setequal(names(true_pars_inner), getParameters(p)))
pouter_true <- log(true_pars_inner)

times <- rbind(seq(0, 10, by=2), seq(10, 30, by = 4), seq(30, 80, by=8))
prediction_true <- (g * x * p)(times, pouter_true)

sigma <- 0.015
sim <- as.data.frame(prediction_true)
sim$value <- sim$value + rnorm(nrow(sim), sd = sigma)
sim$sigma <- sigma

data <- as.datalist(sim, split.by = "condition")
obj <- normL2(data, g * x * p)

## --- 5. Fit the model (Part 6 of example_toy_mapk_fim_fit.R) ---------------
n_cores <- max(1, min(15, parallel::detectCores() - 1))

pouter_init <- setNames(rep(-1, length(pouter_true)), names(pouter_true))

myfit <- trust(obj, pouter_init, rinit = 1, rmax = 10, iterlim = 600)
cat("=== trust() fit ===\n")
cat("Converged :", myfit$converged, " iterations:", myfit$iterations, "\n")
cat("Objective : truth =", obj(pouter_true)$value, " fit =", myfit$value, "\n\n")

best <- myfit$argument
names(best) <- names(pouter_init)

## --- 6. Reference Hessian at the fit (Part B/7) -----------------------------
o_best <- obj(best)
H <- o_best$hessian[names(best), names(best)]
ev <- eigen(H, symmetric = TRUE)$values
zero_tol <- 1e-6 * max(abs(ev))
is_zero <- abs(ev) < zero_tol
cat("=== Hessian eigenvalue spectrum at the fit ===\n")
cat(sum(is_zero), "of", length(ev), "eigenvalues are numerically zero.\n\n")

## --- 7. Multistart run (Part C) --------------------------------------------
cat("=== Multistart run (", n_cores, "cores ) ===\n")
n_fits <- 50
starts <- msParframe(pouter_init, n = n_fits, seed = 21, samplefun = runif, min = -5, max = 5)

sv_start <- lapply(seq_len(n_fits), function(i)
  sensitivitySingularValues(obj, as.parvec(starts, i)))
cond_start <- vapply(sv_start, function(sv) max(sv) / min(sv), numeric(1))
min_sv_start <- vapply(sv_start, min, numeric(1))

fits <- mstrust_fim(obj, starts, studyname = "toyMapkFitCL2", rinit = 1, rmax = 10,
                 cores = n_cores, iterlim = 600)
pf <- as.parframe(fits)
best_ms <- as.parvec(pf, 1)
cat("Converged:", sum(pf$converged), "of", nrow(pf), "fits.\n")
cat("Best objective:", pf$value[1], " (truth:", obj(pouter_true)$value, ")\n\n")

pf$cond_start <- cond_start[pf$index]
pf$min_sv_start <- min_sv_start[pf$index]
pf$cond_end <- vapply(seq_len(nrow(pf)), function(j)
  sensitivityConditionNumber(obj, as.parvec(pf, j)), numeric(1))

waterfall <- plotValues(topByValue(pf), tol = 0.1)
waterfall_path <- file.path(plot_dir, "mapk_baseline_waterfall.pdf")
ggplot2::ggsave(waterfall_path, waterfall, width = 7, height = 4, dpi = 120)
cat("Waterfall plot saved to", waterfall_path, "\n\n")

## --- 8. constraintL2-based regularization -----------------------------------
## Same eigendecomposition of the reference Hessian H (from Part 6) and the
## same `strength`/`eps_reg` calibration as example_toy_mapk_fim_fit.R's
## Part 10, so the weights below are numerically the same quantities -- the
## only thing that changes is how they're packaged for the optimizer:
## constraintL2(mu, sigma) instead of the hand-written quadraticPenalty(mu,
## Prec).
cat("=== Building the FIM- and L2-based constraintL2() penalties ===\n")
eig <- eigen(H, symmetric = TRUE)
s <- eig$values
eps_reg <- min(abs(s[!is_zero]))
strength <- 0.1

## FIM-anisotropic precision, exactly as in example_toy_mapk_fim_fit.R Part
## 10 -- built only to read off its DIAGONAL below, never used as a full
## matrix penalty here (constraintL2 has no slot for the off-diagonal terms).
w_fim <- strength / (s + eps_reg)
Prec_fim <- eig$vectors %*% diag(w_fim) %*% t(eig$vectors)
Prec_fim <- (Prec_fim + t(Prec_fim)) / 2
dimnames(Prec_fim) <- list(names(best), names(best))

## Per-parameter sigma from the diagonal of Prec_fim: sigma_i =
## 1/sqrt(Prec_fim[i,i]) reproduces the FIM penalty's exact single-parameter
## curvature at `best` (see header comment). A tiny floor on the diagonal
## guards against a pathological negative/zero entry turning into a
## negative or infinite sigma.
diag_fim <- pmax(diag(Prec_fim), .Machine$double.eps)
sigma_fim <- setNames(1 / sqrt(diag_fim), names(best))

## Isotropic L2: same w_iso as Part 10's Prec_l2 = w_iso * I -- since that
## penalty was already diagonal, this reproduces it EXACTLY (not just its
## diagonal), by construction (see header derivation).
w_iso <- strength / (median(s) + eps_reg)
sigma_l2 <- setNames(rep(1 / sqrt(w_iso), length(best)), names(best))

cat("eigenvalue range at reference fit:", range(s), "\n")
cat("FIM-derived sigma range:", range(sigma_fim), " (small sigma -> tight/strong pull toward `best`)\n")
cat("L2 sigma (uniform):", sigma_l2[1], "\n\n")

obj_fim <- obj + constraintL2(mu = best, sigma = sigma_fim)
obj_l2  <- obj + constraintL2(mu = best, sigma = sigma_l2)

## --- 9. Shared multistart + diagnostics + plotting routine -----------------
## Identical to example_toy_mapk_fim_fit.R's fitAndPlotRegularized().
fitAndPlotRegularized <- function(obj_reg, label, file_stub) {
  cat("=== ", label, ": multistart run (", n_cores, "cores ) ===\n")
  fits_reg <- mstrust_fim(obj_reg, starts, studyname = paste0("toyMapkFitCL2_", file_stub),
                       rinit = 1, rmax = 10, cores = n_cores, iterlim = 600)
  pf_reg <- as.parframe(fits_reg)

  pf_reg$cond_start <- cond_start[pf_reg$index]
  pf_reg$min_sv_start <- min_sv_start[pf_reg$index]
  pf_reg$cond_end <- vapply(seq_len(nrow(pf_reg)), function(j)
    sensitivityConditionNumber(obj, as.parvec(pf_reg, j)), numeric(1))
  pf_reg$value <- vapply(seq_len(nrow(pf_reg)), function(j)
    obj(as.parvec(pf_reg, j))$value, numeric(1))

  cat("Converged:", sum(pf_reg$converged), "of", nrow(pf_reg), "fits.\n")
  best_reg <- as.parvec(pf_reg, 1)
  cat("Best RAW likelihood:", min(pf_reg$value), " (truth:", obj(pouter_true)$value, ")\n\n")

  waterfall_reg <- plotValues(topByValue(pf_reg), tol = 0.1)
  ggplot2::ggsave(file.path(plot_dir, paste0(file_stub, "_waterfall.pdf")),
                   waterfall_reg, width = 7, height = 4, dpi = 120)

  n_keep <- min(90, nrow(pf_reg))
  pf_reg_top <- pf_reg[order(pf_reg$value)[seq_len(n_keep)], ]
  waterfall_diag_reg <- plotWaterfallDiagnostics(pf_reg_top, pf_reg_top$cond_start,
                                                  pf_reg_top$cond_end, pf_reg_top$min_sv_start, tol = 1)
  ggplot2::ggsave(file.path(plot_dir, paste0(file_stub, "_waterfall_diagnostics.pdf")),
                   waterfall_diag_reg, width = 10, height = 9, dpi = 120)

  pred_fit_reg <- (g * x * p)(seq(0, 80, by = 0.1), best_reg)
  pl_reg <- plot(pred_fit_reg, data)
  ggplot2::ggsave(file.path(plot_dir, paste0(file_stub, "_fit_plot.pdf")),
                   pl_reg, width = 8, height = 4, dpi = 120)

  cat(label, "plots saved to", plot_dir, "with prefix '", file_stub, "'\n\n")
  pf_reg
}

pf_fim <- fitAndPlotRegularized(obj_fim, "FIM (constraintL2) regularization", "mapk_fim_cL2")
pf_l2  <- fitAndPlotRegularized(obj_l2,  "Simple L2 (constraintL2) regularization", "mapk_l2_cL2")

## --- 10. Summary: baseline vs. FIM-regularized vs. L2-regularized ---------
cat("=== Summary comparison (raw likelihood, all three methods) ===\n")
summarize <- function(pf_x, label) {
  cat(sprintf("%-12s converged %3d/%3d, best = %8.2f, range among converged = [%8.2f, %8.2f]\n",
              label, sum(pf_x$converged), nrow(pf_x),
              min(pf_x$value),
              min(pf_x$value[pf_x$converged]), max(pf_x$value[pf_x$converged])))
}
summarize(pf,     "baseline")
summarize(pf_fim, "FIM (cL2)")
summarize(pf_l2,  "L2 (cL2)")
cat("(lower is better; a TIGHTER range among converged fits means less residual\n",
    " practical non-identifiability was left for the multistart to wander into.)\n\n")

## --- 11. Combined prediction plot: truth, baseline, FIM, L2 ----------------
plot_all_path <- file.path(plot_dir, "mapk_all_methods_fit_plot.pdf")
pred_time <- seq(0, 80, by = 0.5)
best_fm <- as.parvec(pf_fim, 1)
best_l2 <- as.parvec(pf_l2, 1)
pred_all <- as.prdlist(list(
  truth    = (g*x*p)(pred_time, pouter_true)[[1]],
  baseline = (g*x*p)(pred_time, best_ms)[[1]],
  fim      = (g*x*p)(pred_time, best_fm)[[1]],
  l2       = (g*x*p)(pred_time, best_l2)[[1]]
))
pl_all <- plotCombined(pred_all, data)
pl_all
ggplot2::ggsave(plot_all_path, pl_all, width = 8, height = 4, dpi = 120)
cat("Combined truth/baseline/FIM/L2 prediction plot saved to", plot_all_path, "\n")

## --- 12. Strength sweep: FIM- and L2-regularization across `strength` -----
## Part 8 fixed `strength <- 0.1` by hand. Here the same construction (same
## `eig`/`s`/`eps_reg`/`best` from Parts 6/8) is repeated across a log-grid of
## `strength` values, so the two regularizations' sensitivity to that one
## hyperparameter becomes visible directly instead of being baked into a
## single hard-coded number. `0.1` is kept in the grid so the Part 8/9 run
## can be reused instead of refit.
cat("=== Strength sweep: FIM and L2 regularization across strengths ===\n")

strengths <- c(0.01, 0.1, 1, 10, 100)

## Same formulas as Part 8, just parameterized by `strength_i` instead of the
## hard-coded `strength`.
buildSigmas <- function(strength_i) {
  w_fim_i <- strength_i / (s + eps_reg)
  Prec_fim_i <- eig$vectors %*% diag(w_fim_i) %*% t(eig$vectors)
  Prec_fim_i <- (Prec_fim_i + t(Prec_fim_i)) / 2
  diag_fim_i <- pmax(diag(Prec_fim_i), .Machine$double.eps)
  sigma_fim_i <- setNames(1 / sqrt(diag_fim_i), names(best))

  w_iso_i <- strength_i / (median(s) + eps_reg)
  sigma_l2_i <- setNames(rep(1 / sqrt(w_iso_i), length(best)), names(best))

  list(sigma_fim = sigma_fim_i, sigma_l2 = sigma_l2_i)
}

sweep_fim <- list()
sweep_l2  <- list()

for (strength_i in strengths) {

  if (isTRUE(all.equal(strength_i, strength))) {
    ## Reuse the Part 8/9 run already computed at the default strength --
    ## no need to refit what's already sitting in pf_fim/pf_l2.
    cat("--- strength =", strength_i, "(reusing Part 8/9 run) ---\n")
    pf_fim_i <- pf_fim
    pf_l2_i  <- pf_l2
  } else {
    cat("--- strength =", strength_i, "---\n")
    sig <- buildSigmas(strength_i)
    obj_fim_i <- obj + constraintL2(mu = best, sigma = sig$sigma_fim)
    obj_l2_i  <- obj + constraintL2(mu = best, sigma = sig$sigma_l2)

    stub_fim <- sprintf("mapk_fim_cL2_str%g", strength_i)
    stub_l2  <- sprintf("mapk_l2_cL2_str%g",  strength_i)

    pf_fim_i <- fitAndPlotRegularized(obj_fim_i, paste0("FIM strength=", strength_i), stub_fim)
    pf_l2_i  <- fitAndPlotRegularized(obj_l2_i,  paste0("L2 strength=", strength_i),  stub_l2)
  }

  pf_fim_i$strength <- strength_i
  pf_l2_i$strength  <- strength_i
  sweep_fim[[as.character(strength_i)]] <- pf_fim_i
  sweep_l2[[as.character(strength_i)]]  <- pf_l2_i
}

## --- 12a. Combined waterfall: all five strengths, one plot per method -----
plotWaterfallByStrength <- function(sweep_list, title, n_plotfits = 50) {
  combined <- do.call(rbind, lapply(sweep_list, function(pf_i)
    data.frame(value = pf_i$value, converged = pf_i$converged, strength = pf_i$strength)))
  combined$strength <- factor(combined$strength, levels = sort(unique(combined$strength)))

  pars <- topByValue(combined, n = n_plotfits)
  pars$index <- seq_len(nrow(pars))

  ggplot2::ggplot(pars, ggplot2::aes(x = index, y = value, color = strength, pch = converged)) +
    ggplot2::geom_point() +
    ggplot2::xlab("index (all strengths pooled, best 400 of 500 by value)") + ggplot2::ylab("value") +
    ggplot2::ggtitle(title) +
    theme_dMod()
}

waterfall_fim_sweep <- plotWaterfallByStrength(sweep_fim, "FIM (constraintL2): strength sweep", n_plotfits = 200)
ggplot2::ggsave(file.path(plot_dir, "mapk_fim_strength_sweep_waterfall.pdf"),
                waterfall_fim_sweep, width = 8, height = 4.5, dpi = 120)

waterfall_l2_sweep <- plotWaterfallByStrength(sweep_l2, "Simple L2 (constraintL2): strength sweep", n_plotfits = 200)
ggplot2::ggsave(file.path(plot_dir, "mapk_l2_strength_sweep_waterfall.pdf"),
                waterfall_l2_sweep, width = 8, height = 4.5, dpi = 120)

cat("Strength-sweep waterfalls saved to", plot_dir, "\n\n")

## --- 12b. Combined fit plot: truth, baseline, best fit at each strength ---
bestByStrength <- function(sweep_list) {
  setNames(lapply(sweep_list, function(pf_i) as.parvec(pf_i, 1)), names(sweep_list))
}

best_fim_sweep <- bestByStrength(sweep_fim)
best_l2_sweep  <- bestByStrength(sweep_l2)

predAtStrengths <- function(best_list, label_prefix) {
  setNames(lapply(best_list, function(b) (g*x*p)(pred_time, b)[[1]]),
           paste0(label_prefix, "_str", names(best_list)))
}

pred_fim_sweep <- as.prdlist(c(
  list(truth = (g*x*p)(pred_time, pouter_true)[[1]],
       baseline = (g*x*p)(pred_time, best_ms)[[1]]),
  predAtStrengths(best_fim_sweep, "fim")
))
pl_fim_sweep <- plotCombined(pred_fim_sweep, data)
ggplot2::ggsave(file.path(plot_dir, "mapk_fim_strength_sweep_fit_plot.pdf"),
                pl_fim_sweep, width = 9, height = 4.5, dpi = 120)

pred_l2_sweep <- as.prdlist(c(
  list(truth = (g*x*p)(pred_time, pouter_true)[[1]],
       baseline = (g*x*p)(pred_time, best_ms)[[1]]),
  predAtStrengths(best_l2_sweep, "l2")
))
pl_l2_sweep <- plotCombined(pred_l2_sweep, data)
ggplot2::ggsave(file.path(plot_dir, "mapk_l2_strength_sweep_fit_plot.pdf"),
                pl_l2_sweep, width = 9, height = 4.5, dpi = 120)

cat("Strength-sweep fit plots saved to", plot_dir, "\n\n")

## --- 12c. Summary table across strengths -----------------------------------
cat("=== Strength-sweep summary (raw likelihood) ===\n")
for (strength_i in strengths) {
  key <- as.character(strength_i)
  summarize(sweep_fim[[key]], sprintf("FIM (s=%-6g)", strength_i))
  summarize(sweep_l2[[key]],  sprintf("L2  (s=%-6g)", strength_i))
}
