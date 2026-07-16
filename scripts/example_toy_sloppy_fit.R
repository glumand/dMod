## ---------------------------------------------------------------------------
## Toy example #2: two independent sub-models glued into one fit, chosen so
## every classic multistart-fitting pathology shows up at once:
##
##   * a linear cascade A -> B -> C -> D -> (degraded), observed only through
##     an unknown-scale readout of the last species. Structurally
##     non-identifiable: since the ODE is linear and homogeneous in the
##     states, D(t) is proportional to A(0), so the observation scale and the
##     initial amount only ever appear as the product scale*A(0). This is an
##     exact symmetry, certified below with symmetryDetection().
##
##   * practically non-identifiable / "sloppy" rate constants k1..k4 in that
##     same cascade: some directions in rate-constant space barely move the
##     prediction while others move it a lot, giving a Hessian eigenvalue
##     spectrum spanning many orders of magnitude (the textbook definition of
##     a "sloppy" model, Gutenkunst et al. 2007) -- unlike the scale/A(0)
##     direction, these WOULD eventually be resolved by enough precise data.
##
##   * a damped oscillator E -> F -> G -> E (a cyclic 3-compartment network;
##     for similar rates w1≈w2≈w3 its relaxation matrix has a complex
##     eigenvalue pair, i.e. it rings as it decays), observed directly and
##     sampled sparsely relative to its own period. This is what the first
##     cut of this script was missing: the cascade above is "sloppy" but
##     essentially unimodal -- trust()/mstrust() navigate its flat valley
##     reliably and land on (statistically) the same objective value every
##     time, so a waterfall plot of it alone shows a single flat step. A
##     frequency that can alias against the sample spacing gives genuine
##     *separate* optima instead of one smeared-out valley: fits that lock
##     onto the true period and fits that lock onto an aliased one land on
##     two distinct, reproducible objective-value plateaus. That is what
##     plotValues()'s waterfall plot is for (see warfarin.R for the same
##     idiom with msnlmeFit()) -- and why part C below needs a small,
##     non-default tol to resolve the steps: the gap between the true-period
##     and aliased-period plateaus is real but only ~0.1-0.3 in objective
##     value here, well under plotValues()'s default tol = 1.
##
## E/F/G/w1/w2/w3 are algebraically independent of A/B/C/D/k1..k4/scale (no
## shared states or rates), so they don't interfere with the cascade's
## structural/sloppy story -- they just add a second, independent source of
## multistart pathology (genuine local optima) to the same combined fit.
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

set.seed(42)
setwd(tempdir())

## --- 1a. Sub-model 1: 4-step irreversible cascade A -> B -> C -> D -> (degr) -
##     dA/dt = -k1*A
##     dB/dt =  k1*A - k2*B
##     dC/dt =  k2*B - k3*C
##     dD/dt =  k3*C - k4*D
## Only D is observed, and only through an unknown gain: obsD = scale*D.
## Kept as its own eqnlist (f_chain) so Part A below can run symmetryDetection
## on exactly this sub-network, undiluted by the oscillator's states.

f_chain <- NULL
f_chain <- addReaction(f_chain, from = "A", to = "B", rate = "k1*A", description = "A to B")
f_chain <- addReaction(f_chain, from = "B", to = "C", rate = "k2*B", description = "B to C")
f_chain <- addReaction(f_chain, from = "C", to = "D", rate = "k3*C", description = "C to D")
f_chain <- addReaction(f_chain, from = "D", to = "",  rate = "k4*D", description = "D degrades")

## --- 1b. Sub-model 2: cyclic 3-compartment oscillator E -> F -> G -> E -----
##     dE/dt = -w1*E          + w3*G
##     dF/dt =  w1*E - w2*F
##     dG/dt =         w2*F   - w3*G
## A catenary chain (A->B->C->D above) only ever has real, negative
## eigenvalues -- pure relaxation. Closing the loop (G feeds back into E)
## changes that: for comparable w1≈w2≈w3 the relaxation matrix picks up a
## complex-conjugate eigenvalue pair, i.e. E(t) rings (damped oscillation) as
## it relaxes. E(0) = 1 is fixed (a known "dose"), F(0) = G(0) = 0, and E is
## observed directly with no free gain -- deliberately no extra structural
## non-identifiability here, this sub-model exists purely to demonstrate
## genuine local optima.

f <- f_chain
f <- addReaction(f, from = "E", to = "F", rate = "w1*E", description = "E to F")
f <- addReaction(f, from = "F", to = "G", rate = "w2*F", description = "F to G")
f <- addReaction(f, from = "G", to = "E", rate = "w3*G", description = "G to E (closes the loop)")

model <- odemodel(f, modelname = "sloppyChain", solver = "CppODE")
x <- Xs(model)

observables <- eqnvec(obsD = "scale*D", obsE = "E")
g <- Y(observables, x, compile = TRUE, modelname = "sloppyObs", attach.input = TRUE)

## --- 2. Parameter transformation: log-parametrization, single condition ---
## getParameters(g, x) (union form), not getParameters(g * x): the composed
## prdfn's "parameters" attribute is inherited from x alone (see "obsfn *
## prdfn" in R/classes.R), so it misses parameters that only g introduces,
## like the observation gain "scale" here.
innerpars <- getParameters(g, x)
trafo <- repar("x~x", x = innerpars)
trafo <- repar("x~0", x = c("B", "C", "D", "F", "G"), trafo)  # B=C=D=F=G = 0 at t=0
trafo <- repar("x~1", x = "E", trafo)                          # E(0) = 1, fixed dose
trafo <- repar("x~exp(x)", x = innerpars, trafo)                # log-transform the rest

p <- P(trafo, condition = "sim", compile = TRUE, modelname = "sloppyTrafo")

## --- 3. "True" parameters and simulated data -------------------------------
## k3, k4 are deliberately close together (a near-degenerate pair of decay
## time scales), which is what makes them practically hard to pin down
## individually from a single noisy readout of D. w1 = w2 = w3 puts the
## oscillator's relaxation matrix at a balanced point where the complex
## eigenvalue pair (and hence the ringing) is most pronounced.
true_pars_inner <- c(k1 = 1.2, k2 = 0.9, k3 = 0.15, k4 = 0.12, A = 4, scale = 2.5,
                     w1 = 1.0, w2 = 1.0, w3 = 1.0)
pouter_true <- log(true_pars_inner)

## obsD needs a wide, densely-sampled window to see the cascade's slow decay;
## obsE needs SPARSE sampling relative to its own oscillation period
## (2*pi/Im(eigenvalue) ~ 2*pi/0.866 ~ 7.25 time units at w1=w2=w3=1) for the
## aliasing that produces genuine local optima. Both observables are
## evaluated on their union, then each keeps only its own time subset.
times_D <- seq(0, 40, by = 2)
times_E <- seq(0, 40, by = 3.5)
times_common <- sort(unique(c(times_D, times_E)))

prediction_true <- as.data.frame((g * x * p)(times_common, pouter_true))

sigma_D <- 0.2
sigma_E <- 0.07
simD <- subset(prediction_true, name == "obsD" & time %in% times_D)
simE <- subset(prediction_true, name == "obsE" & time %in% times_E)
simD$value <- simD$value + rnorm(nrow(simD), sd = sigma_D); simD$sigma <- sigma_D
simE$value <- simE$value + rnorm(nrow(simE), sd = sigma_E); simE$sigma <- sigma_E
sim <- rbind(simD, simE)

data <- as.datalist(sim, split.by = "condition")

## --- 4. Objective function --------------------------------------------------
obj <- normL2(data, g * x * p)

pouter_init <- pouter_true + rnorm(length(pouter_true), sd = 0.5)

## --- 5. Structural non-identifiability: certify the scale/A(0) symmetry ---
## f_chain (just the A-D sub-network, built in step 1a) fed straight into
## symmetryDetection() together with its observation function. This is a
## formal, noise-free statement: it holds for ANY data, not just this
## realization -- unlike the "sloppy" rate directions or the oscillator's
## local optima below, which are properties of what THIS particular (noisy,
## finite-time, finitely-sampled) data set can and cannot resolve.
## (Running symmetryDetection on the full 7-state network f instead would
## also flag a second, messy polynomial direction mixing F, G, w1, w2, w3 --
## an artifact of leaving F(0), G(0) free, which they are NOT in the actual
## fitted model (trafo pins both to 0). Restricting to f_chain avoids that
## noise and keeps the one symmetry that actually matters for the fit.)
cat("=== Part A: structural non-identifiability (symmetryDetection) ===\n")
sym_out <- symmetryDetection(f_chain, eqnvec(obsD = "scale*D"),
                             method = "observability", reconstruct = TRUE)

## --- 6. Fit via trust()/mstrust() -------------------------------------------
myfit <- trust(obj, pouter_init, rinit = 1, rmax = 10)
cat("\n=== trust() fit ===\n")
cat("Converged :", myfit$converged, " iterations:", myfit$iterations, "\n")
cat("Objective : truth =", obj(pouter_true)$value, " fit =", myfit$value, "\n\n")

## A wide start spread (sd) matters here: it's what lets some mstrust() runs
## fall into the oscillator's aliased-frequency local optimum in Part C
## instead of only ever finding the true one.
fits <- mstrust(obj, pouter_init, studyname = "sloppyFit", rinit = 1, rmax = 10,
                 fits = 600, cores = 1, iterlim = 300, sd = 1.5)
pf <- as.parframe(fits)
best <- as.parvec(pf, 1)

mytimes <- c(seq(0, max(data$sim$time), 1))

prd <- Reduce("*", list(g, x, p))
plotCombined(prd(times_common, best), data)
plotPrediction(prd(mytimes, best))

waterfall <-  plotValues(pf)
waterfall
cat("=== mstrust() best of", nrow(pf), "starts ===\n")
cat("Best objective value:", pf$value[1], "\n\n")

## --- 7. Part B: practical non-identifiability among near-optimal fits ------
## Fits within a chi-square(df=1) deviance of the best one are all
## statistically indistinguishable. Looking at how much each PARAMETER
## varies across that set separates the exactly-flat structural direction
## (scale, A individually) from the merely-wide practical one (k3, k4)
## from the well-constrained one (k1, k2).
near <- pf[pf$value - pf$value[1] < qchisq(0.95, df = 1), ]
inner_near <- exp(as.matrix(near[, names(true_pars_inner)]))

cat("=== Part B: near-optimal fits (", nrow(near), "of", nrow(pf), "starts) ===\n")
cat("scale*A product (the identifiable combination):\n")
print(summary(inner_near[, "scale"] * inner_near[, "A"]))
cat("\nscale alone (structurally non-identifiable):\n")
print(summary(inner_near[, "scale"]))
cat("\nA alone (structurally non-identifiable):\n")
print(summary(inner_near[, "A"]))
cat("\nCoefficient of variation per rate constant (bigger = harder to pin down):\n")
print(apply(inner_near[, c("k1", "k2", "k3", "k4")], 2,
            function(v) sd(v) / mean(v)))

## --- 8. Sloppiness: Hessian eigenvalue spectrum at the best fit ------------
o_best <- obj(best)
H <- o_best$hessian[names(best), names(best)]
ev <- eigen(H, symmetric = TRUE)   # eigen() already returns values decreasing: stiffest first
eigenvalues <- ev$values
cat("\n=== Sloppiness: Hessian eigenvalue spectrum at best fit ===\n")
print(eigenvalues)
cat("Spread (max/min eigenvalue):", max(eigenvalues) / min(abs(eigenvalues)), "\n")

stiffest_dir <- ev$vectors[, 1]
sloppiest_dir <- ev$vectors[, ncol(ev$vectors)]
names(stiffest_dir) <- names(sloppiest_dir) <- names(best)
cat("\nStiffest direction (best-constrained parameter combination):\n")
print(round(stiffest_dir, 3))
cat("\nSloppiest direction (worst-constrained parameter combination --\n",
    "expect scale and A to dominate with opposite sign, the structural\n",
    "symmetry certified in Part A):\n")
print(round(sloppiest_dir, 3))

## --- 9. Part C: genuine local optima -- the waterfall plot -----------------
## plotValues() is dMod's canonical multistart diagnostic (see
## inst/examples/warfarin.R, section 5): sort the converged objective values
## and look for plateaus. A single flat plateau means every start found the
## same optimum (globally identified, up to the flat directions already
## characterized in Parts A/B). Separate plateaus mean separate basins of
## attraction -- genuinely different fits, not just noise on one fit. That
## second kind of structure is exactly what the cascade in Parts A/B does
## NOT produce (it's sloppy but unimodal); it's what the oscillator was added
## for. The plateau gap here is real but modest (tens of a log-likelihood
## unit) -- well under plotValues()'s default tol = 1, so it needs a smaller
## tol to resolve; check attr(p_waterfall, "jumps") / the printed value table
## against whatever tol you pick.
cat("\n=== Part C: genuine local optima (waterfall plot) ===\n")
cat("Sorted objective values across all", nrow(pf), "starts:\n")
print(table(round(sort(pf$value), 2)))

p_waterfall <- plotValues(pf, tol = 0.05)
p_waterfall
waterfall_plot_path <- file.path(tempdir(), "example_toy_sloppy_fit_waterfall.png")
ggplot2::ggsave(waterfall_plot_path, p_waterfall, width = 7, height = 4, dpi = 120)
cat("Waterfall plot saved to", waterfall_plot_path, "\n")
cat("(steps detected at rank(s):", paste(attr(p_waterfall, "jumps"), collapse = ", "), ")\n")

## --- 10. Profile likelihood: unbounded (structural) vs. wide-but-bounded ---
## (practical) confidence interval -----------------------------------------
cat("\n=== Profile likelihood: scale (structural) vs. k4 (practical) ===\n")
profs <- do.call(rbind, list(
  profile(obj, best, whichPar = "scale", limits = c(-2, 2), method = "optimize"),
  profile(obj, best, whichPar = "k1",    limits = c(-2, 2), method = "optimize"),
  profile(obj, best, whichPar = "k2",    limits = c(-2, 2), method = "optimize"),
  profile(obj, best, whichPar = "k3",    limits = c(-2, 2), method = "optimize"),
  profile(obj, best, whichPar = "k4",    limits = c(-2, 2), method = "optimize"),
  profile(obj, best, whichPar = "A",    limits = c(-2, 2), method = "optimize"),
  profile(obj, best, whichPar = "w1",    limits = c(-2, 2), method = "optimize"),
  profile(obj, best, whichPar = "w2",    limits = c(-2, 2), method = "optimize"),
  profile(obj, best, whichPar = "w3",    limits = c(-2, 2), method = "optimize")
))
print(confint(profs, val.column = "value"))

pl_prof <- plotProfile(profs)
pl_prof

path <- plotProfilesAndPaths(profs, c("k4", "scale"))
path
profile_plot_path <- file.path(tempdir(), "example_toy_sloppy_fit_profiles.png")
ggplot2::ggsave(profile_plot_path, pl_prof, width = 7, height = 4, dpi = 120)
cat("\nProfile plot saved to", profile_plot_path, "\n")

## --- 11. Prediction vs. data plot -------------------------------------------
pred_fit <- (g * x * p)(seq(0, 40, by = 0.2), best)
pl <- plot(pred_fit, data)
plot_path <- file.path(tempdir(), "example_toy_sloppy_fit_plot.png")
ggplot2::ggsave(plot_path, pl, width = 7, height = 4, dpi = 120)
cat("Prediction plot saved to", plot_path, "\n")
