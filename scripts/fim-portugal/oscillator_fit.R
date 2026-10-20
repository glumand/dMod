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
## fim_functions.R's groupL2()/nullSpaceProjector() address a
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
## Uses the installed dMod2 (loaded by ../setup.R). Invoke as
## `Rscript scripts/fim-portugal/oscillator_fit.R` from anywhere.
## ---------------------------------------------------------------------------

# Locate this script's folder (Rscript or source()); pasting the code into
# a console leaves no file-based signal, so fall back to the fixed path.
# Edit it if this checkout ever moves.
script_dir <- local({
  hit   <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  ofile <- if (length(sys.frames())) sys.frame(1)$ofile else NULL
  if (length(hit)) dirname(normalizePath(sub("^--file=", "", hit[1])))
  else if (!is.null(ofile)) dirname(normalizePath(ofile))
  else "/home/mio/Phd/2026_dMod_fim-regularization/scripts/fim-portugal"
})
source(file.path(script_dir, "..", "setup.R"))

set.seed(7)
out_dir <- script_dir
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

model <- odemodel(f, modelname = "multiOscillatorModel")
x <- Xs(model)

## --- 2. Observation function: each E_m observed directly, no free gain ----
observables <- as.eqnvec(setNames(paste0("E", seq_len(M)), paste0("obsE", seq_len(M))))
g <- Y(observables, x, compile = TRUE, modelname = "multiOscillatorObs", attachInput = FALSE)

## --- 3. Parameter transformation: log-parametrization, single condition ---
innerpars <- getParameters(g * x)
trafo <- repar(expr = "x~x", x = innerpars)
zero_states <- as.vector(outer(c("F", "G"), seq_len(M), paste0))     # F_m(0) = G_m(0) = 0
trafo <- repar(trafo, "x~0", x = zero_states)
one_states <- paste0("E", seq_len(M))                                # E_m(0) = 1, fixed dose
trafo <- repar(trafo, "x~1", x = one_states)
trafo <- repar(trafo, "x~exp(x)", x = innerpars)                     # log-transform all w1_m,w2_m,w3_m

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

fits <- mstrust(obj, pouter_true, name = "multiOscillatorFit", rinit = 1, rmax = 10,
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
