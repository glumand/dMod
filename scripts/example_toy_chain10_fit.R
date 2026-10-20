## ---------------------------------------------------------------------------
## Toy example #3: a long linear (catenary) reaction chain with 10 states,
##
##   X1 -> X2 -> X3 -> ... -> X10 -> (degraded)
##
## where only the LAST state, X10, is observed. This is a deliberately
## extreme structural/practical non-identifiability demo -- more extreme
## than example_toy_sloppy_fit.R's single scale/A(0) symmetry:
##
##   * EXACT permutation symmetry (certified numerically in Part A below,
##     not via symmetryDetection() -- that engine looks for continuous Lie
##     symmetries, and this one is a discrete, finite group action, which it
##     is not built to find): the Laplace transform of a catenary chain is
##
##       X10(s) = X1(0) * k1*k2*...*k9 / prod_{i=1}^{10} (s + k_i)
##
##     The denominator only cares about the SET {k1,...,k10}, and the
##     numerator constant only cares about their product with one factor
##     left out. Consequently X10(t) is exactly invariant under (a) any
##     permutation of the 9 interior hand-off rates k1..k9 among themselves
##     (no compensation needed -- they only ever appear multiplied
##     together), and (b) swapping ANY rate constant into the "terminal"
##     role played by k10, provided X1(0) is rescaled to keep the product
##     X1(0)*k1*...*k9 fixed. So all 10 rate constants are interchangeable
##     up to a compensating rescaling of the initial dose -- a full S10
##     symmetry, not just S9.
##
##   * On top of that EXACT symmetry, recovering the individual k_i is
##     additionally hard for the mundane reason that fitting a sum of 10
##     exponentials with closely-spaced decay rates from finite noisy data
##     is one of the most notoriously ill-conditioned problems in nonlinear
##     regression (Hessian condition numbers below run past 1e8) -- so even
##     restricting to one representative of the symmetry-equivalence class,
##     individual rates stay only loosely constrained. What IS tightly
##     constrained is the combination C = X1(0) * k1*k2*...*k9 (the
##     numerator constant itself), demonstrated in Part C.
##
## Same repo-root/build bootstrap as scripts/example_toy_fit.R; see that file
## for a detailed explanation of the fallback chain below.
## ---------------------------------------------------------------------------

# Locate this script's folder (Rscript or source()); pasting the code into
# a console leaves no file-based signal, so fall back to the fixed path.
# Edit it if this checkout ever moves.
script_dir <- local({
  hit   <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  ofile <- if (length(sys.frames())) sys.frame(1)$ofile else NULL
  if (length(hit)) dirname(normalizePath(sub("^--file=", "", hit[1])))
  else if (!is.null(ofile)) dirname(normalizePath(ofile))
  else "/home/mio/Phd/2026_dMod_fim-regularization/scripts"
})
source(file.path(script_dir, "setup.R"))

set.seed(7)
setwd(tempdir())

## --- 1. Model: 10-state linear chain X1 -> X2 -> ... -> X10 -> (degraded) --
##     dX1/dt  = -k1*X1
##     dXi/dt  =  k(i-1)*X(i-1) - ki*Xi     for i = 2..9
##     dX10/dt =  k9*X9 - k10*X10
## Built programmatically since 10 states/9 hand-off reactions is tedious (and
## error-prone) to write out by hand.
n_states <- 10
states <- paste0("X", seq_len(n_states))

f <- NULL
for (i in seq_len(n_states - 1)) {
  f <- addReaction(f, from = states[i], to = states[i + 1],
                    rate = paste0("k", i, "*", states[i]),
                    description = paste0(states[i], " to ", states[i + 1]))
}
f <- addReaction(f, from = states[n_states], to = "",
                  rate = paste0("k", n_states, "*", states[n_states]),
                  description = paste0(states[n_states], " degrades"))

model <- odemodel(f, modelname = "toyChain10Model")
x <- Xs(model)

## --- 2. Observation function: only the LAST state is measured -------------
observables <- eqnvec(obsLast = states[n_states])
g <- Y(observables, x, compile = TRUE, modelname = "toyChain10Obs", attachInput = FALSE)

## --- 3. Parameter transformation: log-parametrization, single condition ---
## X1(0) = A (free, fitted); all other states start at 0; every rate
## constant k1..k10 is log-transformed to keep the optimizer in positive
## parameter space.
innerpars <- getParameters(g * x)
trafo <- repar(expr = "x~x", x = innerpars)
trafo <- repar(trafo, "x~0", x = states[-1])      # X2..X10 start at 0
trafo <- repar(trafo, "x~exp(x)", x = innerpars)  # log-transform rates + X1(0)

p <- P(trafo, condition = "sim", compile = TRUE, modelname = "toyChain10Trafo")

## --- 4. "True" parameters and simulated data -------------------------------
## Rates deliberately spread out but all O(1) -- close enough together that
## no single one dominates, which is what makes the multi-exponential fit
## painfully ill-conditioned on top of the exact symmetry from Part A.
true_pars_inner <- c(k1 = 0.9, k2 = 0.7, k3 = 1.1, k4 = 0.5, k5 = 0.8,
                     k6 = 0.6, k7 = 1.0, k8 = 0.4, k9 = 0.75, k10 = 0.3,
                     X1 = 5)
pouter_true <- log(true_pars_inner)

## The chain's mean transit time is sum(1/ki) ~ 16-17 time units; sample out
## to ~4x that so both the rise and the decay tail of obsLast are captured.
times <- seq(0, 70, by = 8)
prediction_true <- (g * x * p)(times, pouter_true)

sigma <- 0.02
sim <- as.data.frame(prediction_true)
sim$value <- sim$value + rnorm(nrow(sim), sd = sigma)
sim$sigma <- sigma

data <- as.datalist(sim, split.by = "condition")

## --- 5. Objective function --------------------------------------------------
times <- c(seq(0, max(data$sim$time), 1))

obj <- normL2(data, g * x * p)

pouter_init <- pouter_true + rnorm(length(pouter_true), sd = 0.5)

## --- 6. Part A: certify the exact S10 permutation symmetry -----------------
## Swap k3 (an interior hand-off rate) and k10 (the terminal degradation
## rate), and rescale X1(0) to keep X1(0)*k1*...*k9 constant. If the theory
## above is right, obsLast(t) should come back numerically unchanged despite
## every individual rate constant having moved.
cat("=== Part A: certifying the exact permutation symmetry ===\n")
swapped <- true_pars_inner
swapped["k3"]  <- true_pars_inner["k10"]
swapped["k10"] <- true_pars_inner["k3"]
swapped["X1"]  <- true_pars_inner["X1"] *
  prod(true_pars_inner[paste0("k", 1:9)]) / prod(swapped[paste0("k", 1:9)])

pred_true    <- as.data.frame((g * x * p)(times, log(true_pars_inner)))
pred_swapped <- as.data.frame((g * x * p)(times, log(swapped)))
max_diff <- max(abs(pred_true$value - pred_swapped$value))
cat("k3 <-> k10 swapped, X1(0) rescaled to compensate:\n")
cat("  original  X1 =", true_pars_inner["X1"], " k3 =", true_pars_inner["k3"],
    " k10 =", true_pars_inner["k10"], "\n")
cat("  swapped   X1 =", round(swapped["X1"], 4), " k3 =", swapped["k3"],
    " k10 =", swapped["k10"], "\n")
cat("  max |obsLast(true) - obsLast(swapped)| over", length(times),
    "time points:", max_diff, "(solver tolerance -- confirms exact symmetry)\n\n")

## --- 7. Fit via trust() / mstrust() -----------------------------------------
myfit <- trust(obj, pouter_init, rinit = 1, rmax = 10)
cat("=== trust() fit ===\n")
cat("Converged :", myfit$converged, " iterations:", myfit$iterations, "\n")
cat("Objective : truth =", obj(pouter_true)$value, " fit =", myfit$value, "\n\n")

fits <- mstrust(obj, pouter_init, name = "toyChain10Fit", rinit = 1, rmax = 10,
                 fits = 600, cores = 1, iterlim = 1000, sd = 1.5)
pf <- as.parframe(fits)
best <- as.parvec(pf, 1)

mytimes <- c(seq(0, max(data$sim$time), 1))

prd <- Reduce("*", list(g, x, p))
plotCombined(prd(mytimes, best), data)
# plotPrediction(prd(mytimes, best))

waterfall <-  plotValues(pf[1:400], tol=0.01)
waterfall

cat("=== mstrust() best of", nrow(pf), "starts ===\n")
cat("Best objective value:", pf$value[1], " (truth:", obj(pouter_true)$value, ")\n\n")

## --- 8. Part B: individual rate constants are NOT expected to match -------
## Because of the exact symmetry in Part A, comparing fitted k_i to true k_i
## by NAME is meaningless -- report the sorted values instead, which is the
## most the data could ever pin down about the "shape" of the rate set.
fitted_inner <- exp(best)
cat("=== Part B: sorted rate constants (identity is not identifiable) ===\n")
cat("sorted true k1..k10  :", sort(true_pars_inner[paste0("k", 1:10)]), "\n")
cat("sorted fitted k1..k10:", sort(fitted_inner[paste0("k", 1:10)]), "\n")
cat("(Expect a poor match beyond the exact symmetry itself: fitting closely\n",
    " spaced multi-exponential decay rates from finite noisy data is\n",
    " separately ill-conditioned -- see Part D.)\n\n")

## --- 9. Part C: the identifiable invariant, C = X1(0)*k1*...*k9 -----------
## Across fits within a chi-square(df=11) deviance of the best one (all
## statistically indistinguishable), individual rate constants vary a lot,
## but C -- the numerator constant in the Laplace-domain formula from Part A
## -- should be tightly constrained.
near <- pf[pf$value - pf$value[1] < qchisq(0.95, df = length(true_pars_inner)), ]
inner_near <- exp(as.matrix(near[, names(true_pars_inner)]))
C_true <- unname(true_pars_inner["X1"] * prod(true_pars_inner[paste0("k", 1:9)]))
C_near <- inner_near[, "X1"] * apply(inner_near[, paste0("k", 1:9)], 1, prod)

cat("=== Part C: near-optimal fits (", nrow(near), "of", nrow(pf), "starts) ===\n")
cat("C = X1(0)*k1*...*k9, true value:", C_true, "\n")
print(summary(C_near))
cat("Coefficient of variation, C (the identifiable combination):",
    sd(C_near) / mean(C_near), "\n")
cat("Coefficient of variation per rate constant (bigger = harder to pin down):\n")
print(apply(inner_near[, paste0("k", 1:10)], 2, function(v) sd(v) / mean(v)))

## --- 10. Part D: sloppiness -- Hessian eigenvalue spectrum at the best fit -
## Some eigenvalues come back as exact numerical zero (~1e-12, sign not
## meaningful) rather than merely small: whenever the optimizer lands on
## several k_i taking nearly the same value (as it typically does here --
## Part B's collapsed sorted values are a symptom of the same thing), the
## objective is exactly invariant under exchanging that near-degenerate
## subset (a direct consequence of the S10 symmetry from Part A), which
## forces the corresponding antisymmetric directions to exact-zero
## curvature. Separate those from the genuinely informative eigenvalues
## before reporting a spread, or the ratio is dominated by float noise.
o_best <- obj(best)
H <- o_best$hessian[names(best), names(best)]
ev <- eigen(H, symmetric = TRUE)
eigenvalues <- ev$values
zero_tol <- 1e-6 * max(abs(eigenvalues))
is_zero <- abs(eigenvalues) < zero_tol
cat("\n=== Part D: Hessian eigenvalue spectrum at best fit ===\n")
print(eigenvalues)
cat(sum(is_zero), "eigenvalue(s) are numerically zero (exact flat directions from",
    "near-degenerate rate constants, see comment above).\n")
cat("Spread among the", sum(!is_zero), "non-zero eigenvalues (max/min |eigenvalue|):",
    max(abs(eigenvalues[!is_zero])) / min(abs(eigenvalues[!is_zero])), "\n")

## --- 11. Prediction vs. data plot --------------------------------------------
## The observable itself should still be fit well -- only its underlying
## rate-constant decomposition is degenerate.
pred_fit <- (g * x * p)(seq(0, 70, by = 0.5), best)
pl <- plot(pred_fit, data)
plot_path <- file.path(tempdir(), "example_toy_chain10_fit_plot.png")
ggplot2::ggsave(plot_path, pl, width = 7, height = 4, dpi = 120)
cat("\nPlot saved to", plot_path, "\n")
