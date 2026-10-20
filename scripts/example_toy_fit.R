## ---------------------------------------------------------------------------
## Small example: simulate data from a toy ODE model, then re-estimate the
## parameters with dMod2 (the installed package, loaded by setup.R).
## Invoke as `Rscript scripts/example_toy_fit.R` from anywhere.
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

set.seed(1)
setwd(tempdir())

## --- 1. Model: simple irreversible reaction chain A -> B -> (degraded) ----
##     dA/dt = -k1*A
##     dB/dt =  k1*A - k2*B

f <- NULL
f <- addReaction(f, from = "A", to = "B", rate = "k1*q1*A",
                  description = "conversion of A to B")
f <- addReaction(f, from = "B", to = "",  rate = "k2*B",
                  description = "degradation of B")

model <- odemodel(f, modelname = "toyModel")
x <- Xs(model)

## --- 2. Observation function: observe A and B directly --------------------
observables <- eqnvec(obsA = "A", obsB = "B")
g <- Y(observables, x, compile = TRUE, modelname = "toyObs", attachInput = FALSE)

## --- 3. Parameter transformation: log-parametrization, single condition ---
innerpars <- getParameters(g * x)
trafo <- repar(expr = "x~x", x = innerpars)
trafo <- repar(trafo, "x~0", x = "B")             # B(0) = 0
trafo <- repar(trafo, "x~exp(x)", x = innerpars)  # log-transform k1, k2, A0

p <- P(trafo, condition = "sim", compile = TRUE, modelname = "toyTrafo")

## --- 4. "True" parameters and simulated data -------------------------------
true_pars_inner <- c(k1 = 0.8, q1 = 1.2, k2 = 0.3, A = 5)
pouter_true <- log(true_pars_inner)

times <- seq(0, 12, by = 1)
prediction_true <- (g * x * p)(times, pouter_true)

sigma <- 0.15
sim <- as.data.frame(prediction_true)
sim$value <- sim$value + rnorm(nrow(sim), sd = sigma)
sim$sigma <- sigma

data <- as.datalist(sim, split.by = "condition")

## --- 5. Objective function --------------------------------------------------
obj <- normL2(data, g * x * p)

pouter_init <- pouter_true + rnorm(length(pouter_true), sd = 0.5)

## --- 6. Fit via trust() / mstrust() -----------------------------------------
## In the old dMod fork, trust()/mstrust() stalled here: obj()'s
## gradient/hessian came back in normL2()'s internal parameter order, which
## trust() consumed positionally. dMod2 returns them in pouter_init's order
## (see the diagnostic in section 8), so no fix is needed any more.
myfit <- trust(obj, pouter_init, rinit = 1, rmax = 10)
cat("=== trust() fit ===\n")
cat("Converged :", myfit$converged, " iterations:", myfit$iterations, "\n")
cat("Objective : truth =", obj(pouter_true)$value, " fit =", myfit$value, "\n\n")

fits <- mstrust(obj, pouter_init, name = "toyFit", rinit = 1, rmax = 10,
                 fits = 20, cores = 1, iterlim = 200)
pf <- as.parframe(fits)
best <- as.parvec(pf, 1)

cat("=== mstrust() best of", nrow(pf), "starts ===\n")
cat("Best objective value:", pf$value[1], "\n\n")

fitted_inner <- exp(best)
comparison <- data.frame(
  parameter    = names(true_pars_inner),
  true         = true_pars_inner,
  fitted       = fitted_inner[names(true_pars_inner)],
  rel_error_pc = 100 * (fitted_inner[names(true_pars_inner)] - true_pars_inner) / true_pars_inner
)
print(comparison, row.names = FALSE)

## --- 7. Prediction vs. data plot --------------------------------------------
pred_fit <- (g * x * p)(seq(0, 12, by = 0.1), best)
pl <- plot(pred_fit, data)
plot_path <- file.path(tempdir(), "example_toy_fit_plot.png")
ggplot2::ggsave(plot_path, pl, width = 7, height = 4, dpi = 120)
cat("\nPlot saved to", plot_path, "\n")

## --- 8. Diagnostic: name order mismatch that used to break trust() ---------
## Kept for the record -- shows the root cause directly.
o0 <- obj(pouter_init)
cat("\n=== Diagnostic: parameter name order ===\n")
cat("pouter_init names :", names(pouter_init), "\n")
cat("obj() gradient names:", names(o0$gradient), "\n")
cat("(dMod2 returns the gradient in pouter_init's order; the old dMod fork\n")
cat(" did not, which broke its positional trust().)\n")

