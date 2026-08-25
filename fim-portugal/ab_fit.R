## ---------------------------------------------------------------------------
## Simple AB model: A -> B with rate k*q, only A observed.
##
## Note: k and q enter the ODE only as the product k*q, so individually they
## are structurally non-identifiable from data on A alone -- only the
## product k*q (and the initial value A0) is identifiable. This is
## deliberate: it's a minimal example of the non-identifiability the
## FIM-null-space regularization work in this repo is about (see
## notes/group_regularization_null_space.tex). Expect the profiles for k and
## q individually to be flat/unbounded to one side, while the profile for
## the derived quantity k*q (and for A0) is well-identified.
##
## Builds/loads the dMod that lives in THIS checkout (DESCRIPTION Version,
## from repo_root) into a local, gitignored library under .Rlib/, so it
## never picks up whatever dMod happens to be installed in the default R
## library. Invoke as `Rscript fim-portugal/ab_fit.R` from anywhere; the
## first run compiles the package, later runs reuse it.
## ---------------------------------------------------------------------------

# Locate the repo root robustly across invocation styles -- same pattern as
# scripts/example_toy_fit.R.
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
out_dir <- file.path(repo_root, "fim-portugal")
setwd(tempdir())

library(dplyr)

## --- 1. Model: A -> B, rate k*q ---------------------------------------------
##     dA/dt = -k*q*A
##     dB/dt =  k*q*A

f <- NULL
f <- addReaction(f, from = "A", to = "B", rate = "k*q*A", description = "A to B by rate k*q")

model <- odemodel(f, modelname = "abModel", solver = "CppODE")
x <- Xs(model)

## --- 2. Observation function: only A is observed ---------------------------
observables <- eqnvec(obsA = "A")
g <- Y(observables, x, compile = TRUE, modelname = "abObs", attach.input = TRUE)

## --- 3. Parameter transformation: log-parametrization, single condition ----
innerpars <- getParameters(g * x)
trafo <- repar("x~x", x = innerpars)
trafo <- repar("x~0", x = "B", trafo)             # B(0) = 0
trafo <- repar("x~exp(x)", x = innerpars, trafo)  # log-transform k, q, A0

p <- P(trafo, condition = "sim", compile = TRUE, modelname = "abTrafo")

## --- 4. "True" parameters and simulated data --------------------------------
true_pars_inner <- c(k = 0.5, q = 0.6, A = 5)
pouter_true <- log(true_pars_inner)

times <- seq(0, 12, by = 1)
prediction_true <- (g * x * p)(times, pouter_true)

sigma <- 0.15
sim <- as.data.frame(prediction_true)
sim$value <- sim$value + rnorm(nrow(sim), sd = sigma)
sim$sigma <- sigma
sim <- sim[sim$name == "obsA", ]
data <- as.datalist(sim, split.by = "condition")

## --- 5. Objective function ---------------------------------------------------
obj <- normL2(data, g * x * p)

outerpars <- getParameters(P(trafo))
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)

## --- 6. Fit via trust() / mstrust() ------------------------------------------
myfit <- trust(obj, pouter, rinit = 1, rmax = 10)
myfit

fits <- mstrust(obj, pouter, studyname = "abFit", rinit = 1, rmax = 10,
                 fits = 200, cores = 1, iterlim = 200, sd = 0.1)

pf <- as.parframe(fits)
best <- as.parvec(pf, 1)
waterfall <- plotValues(pf)
waterfall

parscatter <- plotPars(pf)
parscatter

cat("=== mstrust() best of", nrow(pf), "starts ===\n")
cat("Best objective value:", pf$value[1], "\n\n")


## --- 7. Prediction vs. data plot ---------------------------------------------
pred_fit <- (g * x * p)(seq(0, 12, by = 0.1), best)
pl <- plot(pred_fit, data)
pl
plot_path <- file.path(out_dir, "ab_fit_plot.png")
ggplot2::ggsave(plot_path, pl, width = 7, height = 4, dpi = 120)

## --- 8. Profile likelihood ----------------------------------------------------
cat("\n=== Profile likelihood: k, q, A ===\n")
profs <- do.call(rbind, list(
  profile(obj, best, whichPar = "k", limits = c(-3, 3), method = "optimize"),
  profile(obj, best, whichPar = "q", limits = c(-3, 3), method = "optimize"),
  profile(obj, best, whichPar = "A", limits = c(-3, 3), method = "optimize")
))
print(confint(profs, val.column = "value"))

pl_prof <- plotProfile(profs)
pl_prof
profile_plot_path <- file.path(out_dir, "ab_fit_profiles.png")
ggplot2::ggsave(profile_plot_path, pl_prof, width = 7, height = 4, dpi = 120)
cat("\nProfile plot saved to", profile_plot_path, "\n")

path <- plotProfilesAndPaths(profs, c("k", "q"))
path
path_plot_path <- file.path(out_dir, "ab_fit_profile_paths.png")
ggplot2::ggsave(path_plot_path, path, width = 7, height = 4, dpi = 120)
cat("Profile-path plot (k vs q, should show the k*q=const flat direction) saved to",
    path_plot_path, "\n")

## --- 9. FIM null-space regularized re-fit -------------------------------------
## Project the Hessian at the unregularized fit `best` onto its null (k*q
## non-identifiable) direction and penalize only that -- see
## notes/group_regularization_null_space.tex and fim-portugal/functions-portugal.R.
## The identified combination k*q stays untouched; only drift along k/q's
## gauge direction gets pulled back toward `best`.
source(file.path(out_dir, "functions-portugal.R"))

H <- obj(best)$hessian[names(best), names(best)]
P <- nullSpaceProjector(H)
obj_group <- obj + groupL2(P, best, lambda = 0.1)

fits_group <- mstrust(obj_group, pouter, studyname = "abFitGroupL2", rinit = 1, rmax = 10,
                       fits = 200, cores = 1, iterlim = 200, sd = 0.1)
pf_group <- as.parframe(fits_group)
best_group <- as.parvec(pf_group, 1)

cat("\n=== groupL2-regularized fit: best of", nrow(pf_group), "starts ===\n")
cat("Raw data objective at regularized best:", obj(best_group)$value, "\n")
print(exp(best_group))

## --- 10. Profile likelihood at the regularized fit, for identifiability -------
## groupL2 belongs in the OPTIMIZATION only -- it keeps mstrust from drifting
## to an arbitrary point along the flat k/q valley, giving a reproducible
## `best_group`. For identifiability/model-reduction diagnostics we still
## need to profile the RAW data objective `obj`, not `obj_group`: profiling
## obj_group would bake the prior's artificial curvature into every
## direction, so even a genuinely non-identifiable parameter would come back
## with a finite CI, hiding exactly the signal model reduction depends on.
## Profiling obj (data only) at best_group instead gives the true structural
## picture -- flat/unbounded for k and q individually, a proper parabola for
## A -- same as Part 8, just anchored at the more reliably-found point.
cat("\n=== Profile likelihood (raw data objective, at the regularized fit): k, q, A ===\n")
profs_group <- do.call(rbind, list(
  profile(obj, best_group, whichPar = "k", limits = c(-3, 3), method = "optimize"),
  profile(obj, best_group, whichPar = "q", limits = c(-3, 3), method = "optimize"),
  profile(obj, best_group, whichPar = "A", limits = c(-3, 3), method = "optimize")
))
print(confint(profs_group, val.column = "value"))

pl_prof_group <- plotProfile(profs_group)
pl_prof_group
profile_group_plot_path <- file.path(out_dir, "ab_fit_profiles_at_groupL2_fit.png")
ggplot2::ggsave(profile_group_plot_path, pl_prof_group, width = 7, height = 4, dpi = 120)
cat("\nProfile plot (at regularized fit) saved to", profile_group_plot_path, "\n")
