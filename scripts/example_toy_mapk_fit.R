## ---------------------------------------------------------------------------
## Toy example #5: a larger (15-state, 29-parameter) phosphorylation cascade
## with negative feedback, built as a bigger, more realistic cousin of the
## smaller toy models in this directory, in the style of the MAPK/ERK
## pathway models used throughout the identifiability-analysis literature
## (Kholodenko 2000; Raue et al. 2009).
##
## --- Topology ----------------------------------------------------------
##
##   stimulus -> RecA -| Raf -> RafP -> MEK(P,PP) -> ERK(P,PP) -+-> TF -> TFP -> TFPnuc
##                        ^                                     |
##                        |______________ InhP <-----------------
##                              (negative feedback)
##
## An initial pulse of active receptor (RecA(0) = dose, everything else at
## rest) drives a three-tier kinase cascade Raf -> MEK -> ERK, each tier a
## Michaelis-Menten phosphorylation/dephosphorylation cycle (MEK and ERK are
## DISTRIBUTIVE dual-phosphorylation cycles, i.e. two sequential single-site
## steps, as in the standard Huang-Ferrell MAPK cascade model). Active ERK
## (ERKPP) does two things: it phosphorylates a downstream transcription
## factor that translocates to the nucleus (TF -> TFP -> TFPnuc, the
## "reporter" readout), and it activates an inhibitor (Inh -> InhP) that
## feeds back to throttle Raf activation -- a delayed negative feedback loop
## that turns the single input pulse into a damped oscillatory response (see
## the plot saved at the bottom). RafP additionally shields MEKPP from its
## phosphatase (a local feedforward loop, Ki_ff below) -- so there are two
## distinct non-cascade regulatory interactions, not just the feedback.
##
## 15 states: Rec, RecA, Raf, RafP, MEK, MEKP, MEKPP, ERK, ERKP, ERKPP, Inh,
## InhP, TF, TFP, TFPnuc. Only ERKPP and TFPnuc are observed, each through an
## unknown-scale readout (obsERK = scale_ERK*ERKPP, obsTF = scale_TF*TFPnuc)
## -- exactly as real phospho-Western/reporter data comes in arbitrary
## units, not absolute concentrations.
##
## --- Non-identifiability --------------------------------------------------
##
##   Raf's deactivation step is written with kcat and enzyme total as
##   SEPARATE factors, "kcat_dRaf*Etot_dRaf*RafP/(...)", rather than lumped
##   into one Vmax as every other phosphatase step is. Neither symbol
##   appears anywhere else in the model, so as originally written only
##   their PRODUCT was dynamically meaningful: (kcat_dRaf, Etot_dRaf) and
##   (kcat_dRaf*c, Etot_dRaf/c) gave an EXACTLY identical prediction for any
##   c > 0 -- the textbook enzyme-kinetics non-identifiability (kcat and
##   total enzyme concentration are never separately measurable from
##   product-formation data alone). That has since been fixed by pinning
##   Etot_dRaf to a known constant in the trafo (Part 3) -- as if total
##   phosphatase abundance had been measured independently, the same way
##   the other pool totals are already fixed -- leaving kcat_dRaf as the
##   sole free parameter in that rate law. Part A below confirms the fix:
##   perturbing kcat_dRaf alone now visibly changes the prediction, where
##   before it didn't.
##
##   PRACTICAL (Part B): with the cascade's intermediate states (Raf, RafP,
##   MEK*, Inh, InhP, TF, TFP) never observed directly and only two
##   arbitrary-unit readouts available, the other parameters are still far
##   from all independently pinned down -- the Hessian at a good fit has
##   only a handful of genuinely curved (informative) directions out of 29;
##   the rest are flat to numerical precision. A small multistart run (Part
##   C) shows this directly: several fits land at noticeably different
##   parameter combinations while remaining statistically indistinguishable
##   in fit quality. Fixing Etot_dRaf removed exactly one flat direction --
##   it does not touch this practical non-identifiability at all.
##
## Same repo-root/build bootstrap as scripts/example_toy_fit.R; see that file
## for a detailed explanation of the fallback chain below.
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
setwd(tempdir())

## --- 0. Waterfall diagnostics: condition number at first vs. last step ----
## Same helpers as scripts/example_toy_slide_fit.R (see that file for the
## full derivation): mstrust_fim() already prints a per-fit "practical
## identifiability" summary at each fit's FIRST evaluation -- singular
## values of the flattened prediction Jacobian d(prediction)/d(parameter)
## and their condition number (max/min) -- see
## .mstrust_sensitivity_wrapper() in R/statistics.R (plain mstrust() does
## not print this). That's printed only, not returned, so it can't be
## attached to a parframe.
## sensitivityConditionNumber() below recomputes the same number for an
## arbitrary parameter vector, letting us capture it both BEFORE
## optimization (at each fit's starting point) and AFTER (at its converged
## estimate), then use both as waterfall color scales alongside the
## existing iteration-count one from plotValues() -- particularly relevant
## here since Part B already found this model's Hessian to be extremely
## ill-conditioned at the optimum, so seeing how that condition number
## builds up (or doesn't) from a fit's starting point is informative in its
## own right for a model this large.
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

## sensitivityConditionNumber() collapses the singular-value spectrum down
## to a single max/min ratio; this returns the full sorted vector instead
## (same SVD, just not reduced), for fits where the whole spectrum -- not
## just its condition number -- is worth keeping around. Used below to
## additionally track the SMALLEST singular value at each fit's starting
## point: the single worst-constrained direction right at the start,
## before any optimizer step has a chance to move away from it, which the
## condition number alone can't distinguish from "many moderately bad
## directions" (same ratio, very different practical implications).
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

## Same layout as plotValues.parframe() (R/parClass.R) -- sorted waterfall,
## vertical dashed lines at tol-detected plateaus -- but colored by an
## arbitrary per-fit numeric vector (aligned to pf's row order) on a log
## color scale, since condition numbers span many orders of magnitude.
## stepDetect() isn't exported; reached via ::: to reuse the exact plateau
## logic plotValues() itself uses rather than duplicating it.
##
## reverse = TRUE flips which end is orange vs. blue: for condition numbers
## (bad = large) the default low=blue/high=orange already reads "orange is
## concerning"; for a smallest-singular-value color scale bad instead means
## SMALL, so reverse = TRUE keeps "orange is concerning" true there too,
## rather than silently flipping the meaning of the color between panels.
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

## The waterfall four ways: colored by iteration count (plotValues()'s own
## default), by the sensitivity condition number at each fit's FIRST
## evaluation (before any optimizer step), by the same number re-evaluated
## at each fit's FINAL (converged) estimate, and by the SMALLEST singular
## value at the first evaluation (the worst-constrained single direction
## at the start, as opposed to cond_start's overall ratio).
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
f <- NULL

## Receptor: initial stimulus pulse (RecA(0) = dose), decaying back to rest.
f <- addReaction(f, "RecA", "Rec", "kdeact_Rec*RecA", "receptor deactivation")

## Raf: activated by RecA, throttled by the feedback inhibitor InhP
## (division by (1+InhP/Ki_fb) is a competitive-inhibition-style term);
## deactivated by a phosphatase written with kcat and enzyme total kept as
## separate factors (Etot_dRaf is pinned to a constant in Part 3 below,
## not left free -- see the header comment and Part A for why).
f <- addReaction(f, "Raf",  "RafP", "kcat_Raf*RecA*Raf/(Km_Raf+Raf)/(1+InhP/Ki_fb)", "Raf activation")
f <- addReaction(f, "RafP", "Raf",  "kcat_dRaf*Etot_dRaf*RafP/(Km_dRaf+RafP)", "Raf deactivation")

## MEK: distributive dual phosphorylation by RafP. The second
## dephosphorylation step is shielded by RafP (division by (1+RafP/Ki_ff))
## -- a local feedforward loop layered on top of the cascade.
f <- addReaction(f, "MEK",   "MEKP",  "kcat_MEK1*RafP*MEK/(Km_MEK+MEK)", "MEK phosphorylation 1")
f <- addReaction(f, "MEKP",  "MEKPP", "kcat_MEK2*RafP*MEKP/(Km_MEK+MEKP)", "MEK phosphorylation 2")
f <- addReaction(f, "MEKPP", "MEKP",  "Vmax_dMEK*MEKPP/(Km_dMEK+MEKPP)/(1+RafP/Ki_ff)", "MEK dephosphorylation 2")
f <- addReaction(f, "MEKP",  "MEK",   "Vmax_dMEK*MEKP/(Km_dMEK+MEKP)", "MEK dephosphorylation 1")

## ERK: distributive dual phosphorylation by MEKPP, constitutive phosphatase.
f <- addReaction(f, "ERK",   "ERKP",  "kcat_ERK1*MEKPP*ERK/(Km_ERK+ERK)", "ERK phosphorylation 1")
f <- addReaction(f, "ERKP",  "ERKPP", "kcat_ERK2*MEKPP*ERKP/(Km_ERK+ERKP)", "ERK phosphorylation 2")
f <- addReaction(f, "ERKPP", "ERKP",  "Vmax_dERK*ERKPP/(Km_dERK+ERKPP)", "ERK dephosphorylation 2")
f <- addReaction(f, "ERKP",  "ERK",   "Vmax_dERK*ERKP/(Km_dERK+ERKP)", "ERK dephosphorylation 1")

## Feedback inhibitor: activated by ERKPP, decays constitutively -- closes
## the loop back onto Raf activation above.
f <- addReaction(f, "Inh",  "InhP", "kcat_fb*ERKPP*Inh/(Km_fb+Inh)", "feedback activation")
f <- addReaction(f, "InhP", "Inh",  "kdeact_fb*InhP", "feedback deactivation")

## Downstream transcription factor: phosphorylated by ERKPP, dephosphorylated
## constitutively, and translocates to the nucleus, where it eventually
## turns over back to cytoplasmic TF -- the delay stage of the readout.
f <- addReaction(f, "TF",     "TFP",    "kcat_TF*ERKPP*TF/(Km_TF+TF)", "TF phosphorylation")
f <- addReaction(f, "TFP",    "TF",     "Vmax_dTF*TFP/(Km_dTF+TFP)", "TF dephosphorylation")
f <- addReaction(f, "TFP",    "TFPnuc", "ktrans*TFP", "TF nuclear translocation")
f <- addReaction(f, "TFPnuc", "TF",     "kdeg_nuc*TFPnuc", "nuclear TF turnover")

model <- odemodel(f, modelname = "toyMapkModel", solver = "CppODE")
x <- Xs(model)

## --- 2. Observation: two arbitrary-unit readouts ---------------------------
observables <- eqnvec(obsERK = "scale_ERK*ERKPP", obsTF = "scale_TF*TFPnuc")
g <- Y(observables, x, compile = TRUE, modelname = "toyMapkObs", attach.input = FALSE)

## getParameters(g, x) (union form), not getParameters(g * x): the composed
## prdfn's "parameters" attribute is inherited from x alone, so it misses
## parameters only g introduces (scale_ERK, scale_TF) -- see
## example_toy_sloppy_fit.R for the same note.
innerpars <- getParameters(g, x)

## --- 3. Parameter transformation: log-parametrization ----------------------
## Every protein pool starts fully in its "resting" form (unphosphorylated,
## cytoplasmic, inactive) except RecA, which carries the initial stimulus.
## Pool totals themselves (Rec+RecA, Raf+RafP, MEK+MEKP+MEKPP, ...) are
## fixed at 1 (arbitrary units) rather than free parameters -- the unknown
## SCALE of what's actually measured is carried entirely by scale_ERK/
## scale_TF instead, which is the realistic place for it to live (Western
## blot / reporter units, not absolute protein concentration).
states <- f$states
trafo <- repar("x~x", x = innerpars)
trafo <- repar("x~0", x = c("RafP", "MEKP", "MEKPP", "ERKP", "ERKPP", "InhP", "TFP", "TFPnuc"), trafo)
trafo <- repar("x~1", x = c("Raf", "MEK", "ERK", "Inh", "TF"), trafo)
trafo <- repar("x~dose", x = "RecA", trafo)
trafo <- repar("x~1-dose", x = "Rec", trafo)
ratepars <- setdiff(innerpars, states)
trafo <- repar("x~exp(x)", x = c(ratepars, "dose"), trafo)

# fix structural non-identifiability 
trafo <- repar("x~0.5", x = "Etot_dRaf", trafo)

p <- P(trafo, condition = "sim", compile = TRUE, modelname = "toyMapkTrafo")

## --- 4. "True" parameters and simulated data -------------------------------
## Etot_dRaf is NOT listed here -- it's fixed inside the trafo now (Part 3),
## not part of the outer/fitted parameter vector at all. Its value there
## (0.5) already matches what the data below is generated with, so nothing
## about the simulated dynamics changes from before the fix.
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
cat("Model summary: ", length(states), "states,", length(true_pars_inner) - 3,
    "rate parameters (Etot_dRaf now fixed, not among them), 2 observation scales, 1 dose =",
    length(true_pars_inner), "total fitted parameters\n\n")
pouter_true <- log(true_pars_inner)

## Sampled out to 80 time units: the feedback-driven damped oscillation
## (see the plot at the bottom) has largely relaxed by then.
times <- seq(0, 80, by = 2)
prediction_true <- (g * x * p)(times, pouter_true)

sigma <- 0.03
sim <- as.data.frame(prediction_true)
sim$value <- sim$value + rnorm(nrow(sim), sd = sigma)
sim$sigma <- sigma

data <- as.datalist(sim, split.by = "condition")
obj <- normL2(data, g * x * p)

## --- 5. Part A: confirm the fix -- kcat_dRaf is no longer redundant --------
## Before the fix (see git history), any (kcat_dRaf, Etot_dRaf) pair with
## the same product reproduced the data exactly -- moving kcat_dRaf alone
## while Etot_dRaf silently compensated left the prediction unchanged. Now
## that Etot_dRaf is pinned to a constant inside the trafo (Part 3),
## nothing is left to compensate: perturbing kcat_dRaf alone MUST change
## the prediction, which is the whole point of removing the redundant
## degree of freedom.
cat("=== Part A: confirming kcat_dRaf is no longer redundant ===\n")
perturbed <- true_pars_inner
perturbed["kcat_dRaf"] <- 0.7 * true_pars_inner["kcat_dRaf"]

pred_true      <- as.data.frame((g * x * p)(times, log(true_pars_inner)))
pred_perturbed <- as.data.frame((g * x * p)(times, log(perturbed)))
max_diff <- max(abs(pred_true$value - pred_perturbed$value))
cat("kcat_dRaf:", true_pars_inner["kcat_dRaf"], "->", perturbed["kcat_dRaf"],
    " (Etot_dRaf fixed at 0.5 throughout, nothing left free to compensate)\n")
cat("max |prediction(true) - prediction(perturbed)| over", length(times) * 2,
    "data points:", max_diff, "(clearly nonzero -- kcat_dRaf now genuinely changes the",
    "dynamics, so the exact structural symmetry is gone)\n\n")

## --- 6. Fit the model -------------------------------------------------------
## This is a much bigger optimization than the other toy scripts (29
## parameters, a stiff/oscillatory 15-state ODE) -- single fits take
## seconds, not milliseconds, so unlike the smaller scripts this one
## genuinely benefits from running multiple starts in parallel.
n_cores <- max(1, min(15, parallel::detectCores() - 1))
pouter_init <- pouter_true + rnorm(length(pouter_true), sd = 0.15)

myfit <- trust(obj, pouter_init, rinit = 1, rmax = 10, iterlim = 600)
cat("=== trust() fit ===\n")
cat("Converged :", myfit$converged, " iterations:", myfit$iterations, "\n")
cat("Objective : truth =", obj(pouter_true)$value, " fit =", myfit$value, "\n\n")

best <- myfit$argument
names(best) <- names(pouter_init)

## --- 7. Part B: practical non-identifiability -- Hessian eigenvalue spectrum
## As in example_toy_chain10_fit.R: separate genuinely informative curvature
## from numerical noise before reporting a spread, since several
## eigenvalues typically come back as exact zero (to float precision) at a
## point this degenerate.
o_best <- obj(best)
H <- o_best$hessian[names(best), names(best)]
ev <- eigen(H, symmetric = TRUE)$values
zero_tol <- 1e-6 * max(abs(ev))
is_zero <- abs(ev) < zero_tol
cat("=== Part B: Hessian eigenvalue spectrum at the fit ===\n")
cat(sum(is_zero), "of", length(ev),
    "eigenvalues are numerically zero (flat/practically non-identifiable directions).\n")
cat("Spread among the", sum(!is_zero), "non-zero eigenvalues (max/min |eigenvalue|):",
    max(abs(ev[!is_zero])) / min(abs(ev[!is_zero])), "\n")
cat("Only", sum(!is_zero), "of", length(ev),
    "parameter directions are genuinely constrained by this data at all.\n\n")

## --- 8. Part C: a small multistart run --------------------------------------
## With a 29-parameter, largely-flat objective, don't expect every start to
## fully converge in a modest iterlim -- that's itself part of the practical
## non-identifiability story. plotValues() below shows the resulting spread.
##
## Starting points are precomputed as a parframe (msParframe(), a
## reproducible-random-starts helper) instead of letting mstrust_fim() sample
## them internally -- passed as `center`, mstrust_fim() then uses each row's
## parameters directly as that fit's parinit (no further randomization; see
## the `is.parframe(center)` branch in mstrust_fim()/R/statistics.R), so
## cond_start[i], computed from starts row i BEFORE fitting, lines up
## index-for-index with fit i.
cat("=== Part C: multistart run (", n_cores, "cores ) ===\n")
n_fits <- 100
starts <- msParframe(pouter_init, n = n_fits, seed = 21)

## Full singular-value spectrum of the sensitivity Jacobian at each fit's
## STARTING point, saved as a list (sv_start[[i]] is fit i's full sorted
## vector) so both the overall condition number and the single smallest
## singular value can be pulled from the same evaluation -- computing them
## separately would mean calling obj() twice per start for no reason.
sv_start <- lapply(seq_len(n_fits), function(i)
  sensitivitySingularValues(obj, as.parvec(starts, i)))
cond_start <- vapply(sv_start, function(sv) max(sv) / min(sv), numeric(1))
min_sv_start <- vapply(sv_start, min, numeric(1))

fits <- mstrust_fim(obj, starts, studyname = "toyMapkFit", rinit = 1, rmax = 10,
                 cores = n_cores, iterlim = 600)
pf <- as.parframe(fits)
best_ms <- as.parvec(pf, 1)
cat("Converged:", sum(pf$converged), "of", nrow(pf), "fits.\n")
cat("Best objective:", pf$value[1], " (truth:", obj(pouter_true)$value, ")\n")

## as.parframe() preserves each row's original 1..n_fits start slot in
## $index even after sorting by value, so cond_start[pf$index] (and
## min_sv_start[pf$index]) still line up correctly. cond_end is evaluated
## fresh at each row's own converged estimate.
pf$cond_start <- cond_start[pf$index]
pf$min_sv_start <- min_sv_start[pf$index]
pf$cond_end <- vapply(seq_len(nrow(pf)), function(j)
  sensitivityConditionNumber(obj, as.parvec(pf, j)), numeric(1))

waterfall <- plotValues(pf, tol = 1)
waterfall
waterfall_path <- file.path(tempdir(), "example_toy_mapk_waterfall.png")
ggplot2::ggsave(waterfall_path, waterfall, width = 7, height = 4, dpi = 120)
cat("Waterfall plot saved to", waterfall_path, "\n")

waterfall_diag <- plotWaterfallDiagnostics(pf[1:90], pf$cond_start, pf$cond_end, pf$min_sv_start, tol = 1)
waterfall_diag
waterfall_diag_path <- file.path(tempdir(), "example_toy_mapk_waterfall_diagnostics.png")
ggplot2::ggsave(waterfall_diag_path, waterfall_diag, width = 10, height = 9, dpi = 120)
cat("Waterfall diagnostics plot saved to", waterfall_diag_path, "\n")

## --- 9. Prediction vs. data plot --------------------------------------------
pred_fit <- (g * x * p)(seq(0, 80, by = 0.5), best_ms)
pl <- plot(pred_fit, data)
plot_path <- file.path(tempdir(), "example_toy_mapk_fit_plot.png")
ggplot2::ggsave(plot_path, pl, width = 8, height = 4, dpi = 120)
cat("Prediction plot saved to", plot_path, "\n")
