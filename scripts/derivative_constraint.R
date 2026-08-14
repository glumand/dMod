## ---------------------------------------------------------------------------
## Derivative (curvature) regularization: penalize the discretized SECOND
## time-derivative of every state and observable, to discourage biologically
## implausible sharp kinks/edges in a fit while leaving fast-but-smooth
## dynamics (e.g. genuine oscillation or ultrasensitive switching) alone.
##
## Same model, same simulated data, same baseline fitting setup as
## scripts/fim_l2_constraintL2.R -- but this file is standalone and does NOT
## build the FIM-anisotropic / isotropic-L2 constraintL2() penalties, does
## NOT eigendecompose a reference Hessian, and does NOT run the
## identifiability/condition-number/waterfall diagnostics from that script.
## It only compares: baseline fit vs. baseline fit + curvature penalty.
##
## Design (see conversation record for the full derivation):
##  - Penalized quantities: ALL states and ALL observables, not just the two
##    observables that are actually fit to data. States come along "for
##    free" via Y(..., attach.input = TRUE): the state trajectory is already
##    computed internally as an intermediate step of evaluating g*x*p, so
##    attaching it is a cheap array-append, not a second ODE solve.
##  - Curvature measure: a central second difference IN TIME,
##      d2Y/dt2[t] ~= (Y[t+h] - 2*Y[t] + Y[t-h]) / h^2,
##    applied both to the trajectory itself (for the penalty VALUE) and to
##    the model's EXACT parameter-sensitivity array dY/dtheta (for the
##    penalty's GRADIENT/Hessian). Differentiation w.r.t. theta and the fixed
##    linear time-difference operator commute, so the parameter gradient is
##    exact given the time discretization -- theta is never finite-differenced.
##  - Grid: fitting itself now happens on a dense, uniform grid
##    (seq(0, 80, by = 0.1)) instead of the sparse/uneven grid used only to
##    simulate the data. Because normL2()'s kernel only scores rows that
##    match actual data points, this does not change the fitted likelihood
##    at all -- it only means the dense state/observable/sensitivity array
##    needed for the curvature penalty is already sitting in `env$prediction`
##    after every ordinary objective evaluation.
##  - Weighting: per-quantity and self-normalizing. Each state/observable is
##    weighted by the inverse of its OWN reference curvature at the baseline
##    fit `best`, so no single fast-moving species dominates the penalty just
##    because of its numeric scale (mirrors how eps_reg/median(s) calibrate
##    the FIM/L2 constructions in fim_l2_constraintL2.R against the model's
##    own reference behavior, rather than an arbitrary absolute constant).
##  - Architecture: a WRAPPER around the base objective, not a `+`-composed
##    term. It calls the base objective once, reads its already-computed
##    `env$prediction`, and adds the curvature term on top -- one dense-grid
##    ODE+sensitivity integration per evaluation, not two. This matters more
##    as the model grows.
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

plot_dir <- file.path(repo_root, "scripts", "plots", "mapk_derivative_constraint")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

setwd(tempdir())

## --- 1. Model: 15-state phosphorylation cascade with negative feedback ----
## Identical to fim_l2_constraintL2.R / example_toy_mapk_fim_fit.R.
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
## Two variants of the SAME observation function:
##  - `g`      (attach.input = FALSE): only the two observables. Used for
##    simulating the data -- keeps `sim`/`data` exactly as in the other
##    scripts, i.e. no "measurements" of internal states sneak in.
##  - `g_full` (attach.input = TRUE): observables AND all raw states,
##    with exact sensitivities for both. Used to build the FITTING
##    objective, so the curvature penalty has states to work with. This
##    does not change the fitted likelihood at all (normL2 only scores
##    columns that match actual data names) -- it only adds unused-by-the-
##    likelihood columns to env$prediction for the curvature wrapper to read.
observables <- eqnvec(obsERK = "scale_ERK*ERKPP", obsTF = "scale_TF*TFPnuc")
g      <- Y(observables, x, compile = TRUE, modelname = "toyMapkObs",     attach.input = FALSE)
g_full <- Y(observables, x, compile = TRUE, modelname = "toyMapkObsFull", attach.input = TRUE)

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
## Simulated on the ORIGINAL sparse/uneven grid, via `g` (observables only)
## -- identical to fim_l2_constraintL2.R, so `data` only ever contains real
## "measurements" of obsERK/obsTF.
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

times_data <- rbind(seq(0, 10, by = 2), seq(10, 80, by = 4))
prediction_true <- (g * x * p)(times_data, pouter_true)

sigma <- 0.015
sim <- as.data.frame(prediction_true)
sim$value <- sim$value + rnorm(nrow(sim), sd = sigma)
sim$sigma <- sigma

data <- as.datalist(sim, split.by = "condition")

## --- 5. Fitting objective on a DENSE grid, via g_full ----------------------
## `times_fit` covers [0, 80] uniformly at step 0.1 -- every point in
## `times_data` is already an exact multiple of 0.1, so the dense grid
## strictly extends the actual data times rather than conflicting with them.
times_fit <- seq(0, 80, by = 0.1)
obj <- normL2(data, g_full * x * p, times = times_fit)

## --- 6. Baseline fit: multistart of 100 --------------------------------------
## A single trust() run from a fixed init can land on a mediocre local
## optimum in this sloppy model (it did, in an earlier single-start version
## of this script -- see conversation record). Multistart gives a more
## reliable reference point `best` for calibrating the curvature weights.
n_cores <- max(1, min(15, parallel::detectCores() - 1))
n_fits  <- 100

pouter_init <- setNames(rep(-1, length(pouter_true)), names(pouter_true))
starts <- msParframe(pouter_init, n = n_fits, seed = 21, samplefun = runif, min = -5, max = 5)

cat("=== Baseline: multistart run (", n_cores, "cores ) ===\n")
fits <- mstrust(obj, starts, studyname = "derivConstraintBaseline", rinit = 1, rmax = 10,
                 cores = n_cores, iterlim = 600)
pf <- as.parframe(fits)
cat("Converged:", sum(pf$converged), "of", nrow(pf), "fits.\n")
cat("Objective: truth =", obj(pouter_true)$value, " best fit =", pf$value[1], "\n\n")

best <- as.parvec(pf, 1)

## --- 7. Curvature (second-derivative) constraint ---------------------------
## Central second difference in time, applied to either a [time, quantity]
## matrix (for the penalty value) or a [time, quantity, theta] sensitivity
## array (for its exact gradient/Hessian) -- same linear operator either way.
secondDiffTime <- function(arr, h) {
  n <- dim(arr)[1]
  lo <- 1:(n - 2); mid <- 2:(n - 1); hi <- 3:n
  if (length(dim(arr)) == 2) {
    (arr[hi, , drop = FALSE] - 2 * arr[mid, , drop = FALSE] + arr[lo, , drop = FALSE]) / h^2
  } else {
    (arr[hi, , , drop = FALSE] - 2 * arr[mid, , , drop = FALSE] + arr[lo, , , drop = FALSE]) / h^2
  }
}

## Builds a WRAPPER around `base_obj`: calls it once, reads its cached
## env$prediction (already on the dense grid `times_fit`), and adds a
## curvature penalty on top. One dense-grid integration per evaluation.
##
## Weighting: per-quantity, self-normalizing. `ref_j` is quantity j's own
## mean squared curvature AT `best` (computed once, at construction time);
## `weight_j = 1/(ref_j + eps)` so that, evaluated exactly at `best`, every
## quantity contributes ~1 to the (unweighted-by-strength) penalty --
## deviations from that reference are what actually drive the penalty.
##
## `rel_strength` sets the penalty's OWN value at `best` as a fraction of
## the raw likelihood's value at `best` (e.g. 0.1 = "10% of the baseline
## fit's likelihood") -- a self-calibrated, interpretable stand-in for an
## arbitrary absolute constant, in the same spirit as how `eps_reg`/
## `median(s)` calibrate the FIM/L2 penalties against the model's own
## reference behavior in fim_l2_constraintL2.R.
curvatureConstraint <- function(base_obj, best, condition,
                                 rel_strength = 0.1, eps_curv_frac = 1e-6) {

  out_best  <- base_obj(best, deriv = TRUE)
  pred_best <- attr(out_best, "env")$prediction[[condition]]

  time_grid <- pred_best[, "time"]
  hstep <- diff(time_grid)
  stopifnot("curvatureConstraint() assumes a uniform time grid" =
              max(abs(hstep - hstep[1])) < 1e-8 * hstep[1])
  hstep <- hstep[1]

  cols <- dimnames(attr(pred_best, "deriv"))[[2]]
  Y_best <- pred_best[, cols, drop = FALSE]
  curv_best <- secondDiffTime(Y_best, hstep)

  ref <- colMeans(curv_best^2)
  eps_curv <- eps_curv_frac * max(ref)
  weight <- 1 / (ref + eps_curv)

  n_cols <- length(cols)
  n_t    <- nrow(curv_best)
  ## abs(): the baseline objective's value need not be positive (e.g. a
  ## Gaussian NLL's normalization constant can dominate and make the total
  ## negative for small sigma, as it does here) -- strength_curv must stay
  ## positive regardless of that sign convention, or the "penalty" would
  ## reward curvature instead of suppressing it.
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
    curv <- secondDiffTime(Y, hstep)

    pen_value <- C * sum(vapply(seq_len(n_cols), function(j) weight[j] * sum(curv[, j]^2), numeric(1)))

    pen_gradient <- NULL
    pen_hessian  <- NULL

    if (deriv) {
      dY    <- attr(pred, "deriv")[, cols, , drop = FALSE]
      dcurv <- secondDiffTime(dY, hstep)
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

## --- 8. Fit the curvature-regularized objective: multistart, same starts --
## Same `starts` as Part 6, so the comparison isn't confounded by different
## starting points landing in different basins.
cat("=== Curvature-regularized: multistart run (", n_cores, "cores ) ===\n")
fits_curv <- mstrust(obj_curv, starts, studyname = "derivConstraintCurvature", rinit = 1, rmax = 10,
                      cores = n_cores, iterlim = 600)
pf_curv <- as.parframe(fits_curv)
cat("Converged:", sum(pf_curv$converged), "of", nrow(pf_curv), "fits.\n\n")

best_curv <- as.parvec(pf_curv, 1)

## --- 9. Compare: raw likelihood vs. curvature penalty, before and after ---
pen_at_best      <- attr(obj_curv(best), "curvature_penalty")
pen_at_best_curv <- attr(obj_curv(best_curv), "curvature_penalty")

cat("=== Summary: baseline vs. curvature-regularized ===\n")
cat(sprintf("%-22s raw likelihood = %10.3f   curvature penalty = %10.3f\n",
            "baseline (best)", obj(best)$value, pen_at_best))
cat(sprintf("%-22s raw likelihood = %10.3f   curvature penalty = %10.3f\n",
            "curvature-reg (best)", obj(best_curv)$value, pen_at_best_curv))
cat("(raw likelihood evaluated via the UNregularized `obj`, in both rows, so the\n",
    " two are directly comparable; curvature penalty evaluated via `obj_curv`.)\n\n")

## --- 10. Plots --------------------------------------------------------------
## Observables only (what's actually fit to data).
pred_time <- seq(0, 80, by = 0.5)
pred_obs <- as.prdlist(list(
  truth     = (g * x * p)(pred_time, pouter_true)[[1]],
  baseline  = (g * x * p)(pred_time, best)[[1]],
  curvature = (g * x * p)(pred_time, best_curv)[[1]]
))
pl_obs <- plotCombined(pred_obs, data)
ggplot2::ggsave(file.path(plot_dir, "mapk_derivative_constraint_observables.pdf"),
                pl_obs, width = 8, height = 4, dpi = 120)

## All states + observables, to show the smoothing effect on unobserved
## species too. `data` only has rows for the two observables, so the extra
## state panels simply show model curves with no data points overlaid.
pred_full <- as.prdlist(list(
  truth     = (g_full * x * p)(pred_time, pouter_true)[[1]],
  baseline  = (g_full * x * p)(pred_time, best)[[1]],
  curvature = (g_full * x * p)(pred_time, best_curv)[[1]]
))
pl_full <- plotCombined(pred_full, data)
ggplot2::ggsave(file.path(plot_dir, "mapk_derivative_constraint_all_states.pdf"),
                pl_full, width = 12, height = 10, dpi = 120)

cat("Plots saved to", plot_dir, "\n")
