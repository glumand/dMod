## ---------------------------------------------------------------------------
## Small example: simulate data from a toy ODE model, then re-estimate the
## parameters with dMod. Builds/loads the dMod that lives in THIS checkout
## (DESCRIPTION Version, from repo_root) into a local, gitignored library
## under .Rlib/, so it never picks up whatever dMod happens to be installed
## in the default R library. Invoke as `Rscript scripts/example_toy_fit.R`
## from anywhere; the first run compiles the package (needs the toolchain +
## dependencies from DESCRIPTION Imports/Remotes), later runs reuse it.
## ---------------------------------------------------------------------------

repo_root <- {
  full <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", full, value = TRUE)
  script_path <- if (length(hit)) sub("^--file=", "", hit[1]) else NA_character_
  if (is.na(script_path)) normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  else normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
}

repo_lib <- file.path(repo_root, ".Rlib")
dir.create(repo_lib, showWarnings = FALSE, recursive = TRUE)

desc_version <- read.dcf(file.path(repo_root, "DESCRIPTION"), "Version")[1, 1]
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

# Make sure no other dMod is already attached from a different library.
if ("package:dMod" %in% search()) detach("package:dMod", unload = TRUE)

library(dMod, lib.loc = repo_lib)
stopifnot(identical(find.package("dMod"), file.path(repo_lib, "dMod")))
cat("Using dMod", as.character(packageVersion("dMod")),
    "from", find.package("dMod"), "\n\n")

set.seed(1)
setwd(tempdir())

## --- 1. Model: simple irreversible reaction chain A -> B -> (degraded) ----
##     dA/dt = -k1*A
##     dB/dt =  k1*A - k2*B

f <- NULL
f <- addReaction(f, from = "A", to = "B", rate = "k1*A",
                  description = "conversion of A to B")
f <- addReaction(f, from = "B", to = "",  rate = "k2*B",
                  description = "degradation of B")

model <- odemodel(f, modelname = "toyModel", solver = "CppODE")
x <- Xs(model)

## --- 2. Observation function: observe A and B directly --------------------
observables <- eqnvec(obsA = "A", obsB = "B")
g <- Y(observables, x, compile = TRUE, modelname = "toyObs", attach.input = FALSE)

## --- 3. Parameter transformation: log-parametrization, single condition ---
innerpars <- getParameters(g * x)
trafo <- repar("x~x", x = innerpars)
trafo <- repar("x~0", x = "B", trafo)             # B(0) = 0
trafo <- repar("x~exp(x)", x = innerpars, trafo)  # log-transform k1, k2, A0

p <- P(trafo, condition = "sim", compile = TRUE, modelname = "toyTrafo")

## --- 4. "True" parameters and simulated data -------------------------------
true_pars_inner <- c(k1 = 0.8, k2 = 0.3, A = 5)
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
## Previously trust()/mstrust() rejected every trial step on this model and
## stalled at the initial point (while still reporting converged = TRUE).
## Root cause: obj()'s gradient/hessian come back named in normL2()'s own
## internal parameter order (e.g. "A","k1","k2"), which need not match the
## order pouter_init happens to be in ("k1","k2","A"). trust_impl() (src/
## trust_kernel.cpp) copied them into its internal arrays *positionally*,
## silently scrambling the trust-region model while still taking real steps
## in the correctly-labelled parameter vector -- hence every step looked
## like a bad direction, at every radius. Fixed by aligning objfun's
## returned gradient/hessian to parinit's name order before use (see
## align_grad()/align_hess() in trust_kernel.cpp). trust()/mstrust() now
## work directly again.
myfit <- trust(obj, pouter_init, rinit = 1, rmax = 10)
cat("=== trust() fit ===\n")
cat("Converged :", myfit$converged, " iterations:", myfit$iterations, "\n")
cat("Objective : truth =", obj(pouter_true)$value, " fit =", myfit$value, "\n\n")

fits <- mstrust(obj, pouter_init, studyname = "toyFit", rinit = 1, rmax = 10,
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
cat("(trust_kernel.cpp now aligns these by name before use; previously it\n")
cat(" assumed they lined up positionally, which they don't here.)\n")
