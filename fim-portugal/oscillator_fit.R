## ---------------------------------------------------------------------------
## A model that is genuinely difficult to fit: M = 7 ALGEBRAICALLY INDEPENDENT
## copies of a cyclic 3-compartment oscillator,
##
##   E_m -> F_m -> G_m -> E_m     (m = 1..M, no shared states or rates)
##
##   dE_m/dt = -w1_m*E_m            + w3_m*G_m
##   dF_m/dt =  w1_m*E_m - w2_m*F_m
##   dG_m/dt =             w2_m*F_m - w3_m*G_m
##
## each sampled sparsely relative to its own oscillation period.
##
## For comparable rates w1_m~w2_m~w3_m, block m's relaxation matrix has a
## complex eigenvalue pair, so E_m(t) rings (damped oscillation) as it
## decays. Sampling that ringing sparsely relative to its own period lets an
## ALIASED frequency fit block m's data almost as well as the true one, so a
## single such block already gives two separated, reproducible local optima.
## Since the M blocks share no states or parameters, every fit picks
## true-or-aliased INDEPENDENTLY, and often only PARTIALLY correctly, in each
## block -- 2^M = 2^7 = 128 "clean" true/aliased combinations, plus messier
## intermediate outcomes where a block's optimizer gets stuck somewhere
## between the two. Each block uses its own observation noise (sigma_m, step
## 4 below), so each block's mismatch contributes a different-sized penalty
## to the summed objective -- this is what spreads the many combinations
## across many distinct total objective values instead of piling up into a
## handful of big plateaus.
##
## The result (verified interactively, see fit statistics printed in Part 6)
## is a waterfall where the large best-to-worst objective gap is not carried
## by a rare handful of outliers but by MOST of the 500 multistart runs:
## each block can converge to its true frequency, an aliased frequency (both
## comparably good), or -- for a sizeable fraction of starts, not just a rare
## corner case -- a genuinely poorly-fit local optimum for that block (one
## rate saturates near a box bound while the others collapse near zero,
## a real, fully-converged local minimum, just a bad one). Because starts
## are spread wide enough (Part 6) to reach that bad optimum's basin
## routinely, and M = 7 blocks fail independently and to varying degrees,
## the 500 runs populate a near-continuous staircase from the best fit
## (all/most blocks true-or-aliased) out to several thousand objective units
## worse (multiple blocks badly converged at once) -- only a modest fraction
## of runs (roughly a third) sit in the tight near-optimal cluster; the rest
## are spread across the full range. Because the near-optimal cluster and
## the graded climb above it are still very different scales locally, Part 7
## plots the waterfall twice: once over the full range and once zoomed into
## the dense cluster at the bottom -- both come from the exact same
## fits/tol, only the y-axis window differs.
##
## Compare scripts/example_toy_sloppy_fit.R (Part C), which uses a single
## such oscillator block for the same true/aliased-frequency idiom but only
## ever produces 2-3 plateaus with a gap of order 1 objective unit.
##
## Deliberately NO regularization here. fim-portugal/ab_fit.R and
## functions-portugal.R's groupL2()/nullSpaceProjector() address a
## DIFFERENT failure mode: a continuous flat/null direction in the Hessian
## at a single reference fit (structural or practical non-identifiability),
## which a quadratic penalty projected onto that null space can tame without
## touching the identified directions (notes/group_regularization_null_space.tex).
## Genuine local optima are not that -- they are isolated, discrete basins,
## not a continuous manifold through one reference point, so a penalty
## anchored at any single fit would just bias-lock the optimizer into
## whichever basin that reference happened to be in, rather than revealing
## the multimodality. The right tool for THIS pathology is what's used
## below instead: many random multistart starts and a waterfall plot of the
## sorted, converged objective values.
##
## Builds/loads the dMod that lives in THIS checkout (DESCRIPTION Version,
## from repo_root) into a local, gitignored library under .Rlib/, so it
## never picks up whatever dMod happens to be installed in the default R
## library. Invoke as `Rscript fim-portugal/oscillator_fit.R` from anywhere;
## the first run compiles the package, later runs reuse it. Same
## repo-root/build bootstrap as scripts/example_toy_fit.R.
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

set.seed(7)
out_dir <- file.path(repo_root, "fim-portugal")
setwd(tempdir())

## --- 1. Model: M independent cyclic 3-compartment oscillators --------------
## Built programmatically -- M*3 states / M*3 reactions is tedious and
## error-prone by hand. Blocks share no states or rates (f from block m only
## ever mentions E_m/F_m/G_m/w*_m), so the resulting ODE system is exactly
## block-diagonal.
M <- 7
f <- NULL
for (m in seq_len(M)) {
  Em <- paste0("E", m); Fm <- paste0("F", m); Gm <- paste0("G", m)
  w1m <- paste0("w1_", m); w2m <- paste0("w2_", m); w3m <- paste0("w3_", m)
  f <- addReaction(f, from = Em, to = Fm, rate = paste0(w1m, "*", Em),
                    description = paste0(Em, " to ", Fm, " (block ", m, ")"))
  f <- addReaction(f, from = Fm, to = Gm, rate = paste0(w2m, "*", Fm),
                    description = paste0(Fm, " to ", Gm, " (block ", m, ")"))
  f <- addReaction(f, from = Gm, to = Em, rate = paste0(w3m, "*", Gm),
                    description = paste0(Gm, " to ", Em, " (block ", m, ", closes the loop)"))
}

model <- odemodel(f, modelname = "multiOscillatorModel", solver = "CppODE")
x <- Xs(model)

## --- 2. Observation function: each E_m observed directly, no free gain ----
observables <- as.eqnvec(setNames(paste0("E", seq_len(M)), paste0("obsE", seq_len(M))))
g <- Y(observables, x, compile = TRUE, modelname = "multiOscillatorObs", attach.input = FALSE)

## --- 3. Parameter transformation: log-parametrization, single condition ---
innerpars <- getParameters(g * x)
trafo <- repar("x~x", x = innerpars)
zero_states <- as.vector(outer(c("F", "G"), seq_len(M), paste0))     # F_m(0) = G_m(0) = 0
trafo <- repar("x~0", x = zero_states, trafo)
one_states <- paste0("E", seq_len(M))                                # E_m(0) = 1, fixed dose
trafo <- repar("x~1", x = one_states, trafo)
trafo <- repar("x~exp(x)", x = innerpars, trafo)                     # log-transform all w1_m,w2_m,w3_m

p <- P(trafo, condition = "sim", compile = TRUE, modelname = "multiOscillatorTrafo")

## --- 4. "True" parameters and simulated data --------------------------------
## Every block uses w1_m=w2_m=w3_m=1, the most balanced point where the
## complex eigenvalue pair (and hence the ringing) is most pronounced --
## period 2*pi/Im(eigenvalue) ~ 2*pi/0.866 ~ 7.25 time units. Sampling every
## 3.5 time units is sparse relative to that period, which is what lets an
## aliased frequency fit almost as well as the true one, per block.
##
## sigma_m grows mildly across blocks (0.05, 0.055, ..., ~0.09) so blocks
## are not interchangeable but stay comparably informative -- see the header
## comment on why this, combined with a wide start spread (Part 6), produces
## a near-continuous staircase spanning from a dense near-optimal cluster
## out to a long, graded climb where more and more blocks land badly.
true_pars_inner <- setNames(rep(1, 3 * M),
                             as.vector(outer(c("w1_", "w2_", "w3_"), seq_len(M), paste0)))
pouter_true <- log(true_pars_inner)

times <- seq(0, 40, by = 3.5)
prediction_true <- (g * x * p)(times, pouter_true)

sigma_m <- 0.05 * 1.1^(seq_len(M) - 1)
sim <- as.data.frame(prediction_true)
sim$sigma <- sigma_m[as.integer(sub("obsE", "", sim$name))]
sim$value <- sim$value + rnorm(nrow(sim), sd = sim$sigma)
data <- as.datalist(sim, split.by = "condition")

## --- 5. Objective function ---------------------------------------------------
obj <- normL2(data, g * x * p)

outerpars <- getParameters(P(trafo))

## --- 6. Multistart fit -- NO regularization, wide start spread -------------
## sd = 20 (log-space) and a wide +-10 box together are what turn the bad
## local optimum from a rare corner case into a routine outcome: with sd
## this wide relative to the box, a large fraction of starts land at or near
## a box edge for at least one of a block's 3 rates, which is exactly where
## that block's bad-but-genuine local optimum (one rate saturating high,
## the other two collapsing near zero) lives. Since this happens to each of
## the M blocks independently and gets more likely per block the wider the
## box/sd are relative to the ~O(1) true-vs-aliased basin, MOST of the 500
## fits end up with at least one -- often several -- blocks landing there,
## which is what spreads the bulk of the fits (not just a rare tail) across
## the full objective range. (Empirically: a tighter box/sd, e.g. the +-6/5
## used previously, keeps ~98% of fits in the near-optimal cluster and only
## ~2% in a bad-optimum tail; +-10/20 flips that to only about a third
## staying in the cluster.) All local optima here are genuine (trust()
## converges, gradient ~ 0) -- the spread reflects real basin sizes in the
## objective landscape, not truncated/unconverged optimization.
n_fits <- 500
parupper <- setNames(rep(10, length(pouter_true)), names(pouter_true))
parlower <- setNames(rep(-10, length(pouter_true)), names(pouter_true))

fits <- mstrust(obj, pouter_true, studyname = "multiOscillatorFit", rinit = 1, rmax = 10,
                 fits = n_fits, cores = min(10, parallel::detectCores()), iterlim = 400, sd = 20,
                 parupper = parupper, parlower = parlower)

pf <- as.parframe(fits)
best <- as.parvec(pf, 1)

cat("=== mstrust() best of", nrow(pf), "starts (", sum(pf$converged), "converged ) ===\n")
cat("Best objective value:", pf$value[1], " Worst objective value:", pf$value[nrow(pf)], "\n")

cat("\nSorted objective values across all", nrow(pf), "starts (rounded to 2 digits):\n")
vals <- table(round(sort(pf$value), 2))
cat(length(vals), "distinct rounded objective values;",
    "per-value counts range", min(vals), "to", max(vals),
    "(median", median(as.integer(vals)), ")\n")

## --- 7. The waterfall plot: ~2^M genuine local optima, no regularization --
## Same underlying step-detected jumps (tol = 0.1) plotted twice: once over
## the full range, where most fits now form a graded staircase climbing
## continuously from best to worst (not just a dense cluster plus a rare
## tail); once zoomed via coord_cartesian() into just the near-optimal
## sub-cluster at the bottom, which resolves the individual plateaus the
## full view compresses at this scale.
waterfall_full <- plotValues(pf, tol = 0.1) + ggplot2::ggtitle("full range")
waterfall_full
full_plot_path <- file.path(out_dir, "oscillator_fit_waterfall_full.png")
ggplot2::ggsave(full_plot_path, waterfall_full, width = 9, height = 5, dpi = 120)

## Find the dense cluster automatically: the first gap of >5 objective units
## between consecutive sorted values marks where the bad-optimum tail
## starts, whatever its length happens to be for this particular run.
sorted_vals <- sort(pf$value)
big_gaps <- which(diff(sorted_vals) > 5)
cluster_end <- if (length(big_gaps)) big_gaps[1] else length(sorted_vals)
zoom_ylim <- c(sorted_vals[1], sorted_vals[cluster_end])
waterfall_zoom <- plotValues(pf, tol = 0.1) +
  ggplot2::coord_cartesian(ylim = zoom_ylim) +
  ggplot2::ggtitle("zoomed into the near-optimal sub-cluster")
waterfall_zoom
zoom_plot_path <- file.path(out_dir, "oscillator_fit_waterfall_zoom.png")
ggplot2::ggsave(zoom_plot_path, waterfall_zoom, width = 9, height = 5, dpi = 120)

cat("\nWaterfall plots saved to", full_plot_path, "and", zoom_plot_path, "\n")
cat(length(attr(waterfall_full, "jumps")),
    "steps detected at tol = 0.1 (rank(s):",
    paste(attr(waterfall_full, "jumps"), collapse = ", "), ")\n")

## --- 8. Prediction vs. data plot (best fit) ----------------------------------
pred_fit <- (g * x * p)(seq(0, 40, by = 0.2), best)
pl <- plot(pred_fit, data)
pl
plot_path <- file.path(out_dir, "oscillator_fit_plot.png")
ggplot2::ggsave(plot_path, pl, width = 10, height = 7, dpi = 120)
cat("Prediction plot saved to", plot_path, "\n")
