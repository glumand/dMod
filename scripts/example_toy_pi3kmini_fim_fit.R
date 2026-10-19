## ---------------------------------------------------------------------------
## Toy example #6 ("PI3K-mini"): a small model deliberately built so the FIM
## null-space projector (groupL2()/nullSpaceProjector(), fim-portugal/
## functions-portugal.R) actually WINS against both no regularization and a
## plain L2 prior -- the positive counterpart to example_toy_mapk_fim_fit.R
## (where the same method fails once genuine sloppiness, not just one exact
## symmetry, is present).
##
## --- Why this model -----------------------------------------------------
##
## It is a miniature, mechanistically exact copy of the PI3K structural
## non-identifiability found in the real HGF/HepG2 model
## (HGF_HepG2/04-Results/M-22.0-comparison/05-Writeup/
## Technical-Notes_PI3K-Symmetry-and-Regularization-Diagnostics.pdf):
## an UNOBSERVED species entering the rest of the network only as a product
## with one rate constant, with its own abundance set per-condition and its
## own degradation rate untouched by the symmetry. Built here from scratch,
## on purpose, to test the method's own precondition (an isolated, exact
## symmetry with everything else genuinely well-identified) in full control,
## before meeting the same pattern embedded inside a much larger, otherwise
## sloppy model.
##
## --- Topology --------------------------------------------------------------
##
##   A -> B                (control block: k1, well time-resolved, directly
##                           observed via obsA/obsB -- has nothing to do with
##                           the symmetry below)
##   "" -> U  (rate k_gen*A)   \  U is NEVER observed. It is generated
##   U  -> "" (rate k_degU*U)  /  catalytically from A (A is not consumed)
##   Mp -> Mpp (rate k_use*Mp*U)  and degrades on its own; Mpp (observed via
##                                 obsMpp) is the only readout of U's activity.
##
##   Three conditions (cond1/cond2/cond3) differ ONLY in U's initial value
##   (U_cond1/U_cond2/U_cond3), mirroring PI3K_5/PI3K_250/PI3K_350 in the real
##   model exactly.
##
## --- The exact symmetry -----------------------------------------------------
##
##   For any c > 0: U_cond{1,2,3} -> c*U_cond{1,2,3}, k_gen -> c*k_gen,
##   k_use -> k_use/c leaves U's own ODE self-consistent (U(t) -> c*U(t)
##   exactly, since it is linear and homogeneous in U) and leaves Mp/Mpp's
##   ODEs EXACTLY unchanged (k_use*U is invariant, and that product is all
##   Mp/Mpp ever see) -- so every observable trajectory is exactly invariant.
##   k_degU does NOT load on this symmetry (same reason PI3K's own
##   k_degrade_PI3K does not in the real model: pure linear decay of U itself
##   needs no compensation) -- it stays a genuinely identified rate, shaping
##   the TIMING of Mpp's rise, not its scale.
##   A0 and Mp0 are fixed to 1 (not fitted) to avoid two OTHER, unwanted,
##   accidental scale/initial-value symmetries this design would otherwise
##   introduce (A0<->scale_A/scale_B, Mp0<->scale_Mpp) -- see conversation
##   record for the by-hand derivation of why those appear if A0/Mp0 are left
##   free. Discovered and fixed by first running this script with those two
##   free and finding TWO extra unwanted near-zero eigenvalues, not one.
##
## IMPORTANT implementation pitfall (cost real debugging time -- see
## conversation record): Y(..., attach.input = TRUE) attaches ALL raw states
## (including the hidden U!) to the prediction output. If the simulated
## "data" data.frame is built from that attached output without filtering
## down to the actual observables first, U's own (non-invariant) trajectory
## silently becomes part of what normL2() fits against, which completely
## destroys the intended symmetry without any error or warning -- the
## Hessian's eigenvalue spectrum just looks like an ordinary sloppy spread,
## no gap anywhere. Use attach.input = FALSE here for exactly this reason.
##
## --- Fair comparison: SAME regularization center for both penalties -------
##
## Both constraintL2() and groupL2() are centered at the IDENTICAL reference
## point mu_center = -1 in every log-parameter (deliberately not 0/truth, and
## deliberately not the unregularized best fit) -- so any difference between
## the two methods reflects the PENALTY GEOMETRY (diagonal/isotropic vs.
## projected/null-direction-only), not an accidental difference in anchor
## point. (An earlier version anchored groupL2() at best_unreg -- the single
## best unregularized fit -- which, along the exactly flat direction, is
## itself an arbitrary draw from the unregularized multistart's spread, and
## made the two methods' raw group-parameter estimates land at different,
## incomparable points. See conversation record.)
##
## --- What this script does, and does NOT do --------------------------------
##
## This script only builds the model, runs ONE multistart per objective
## (unregularized / L2 / FIM-group), and SAVES the resulting parframes, best
## fits, reference Hessian/eigendecomposition, and fine-grid prediction
## curves to scripts/plots/pi3kmini_fim_fit/pi3kmini_results.rds. It does not
## plot anything -- see example_toy_pi3kmini_fim_plots.R for that, which
## loads the .rds and does not need dMod's compiled model at all (predictions
## are precomputed and saved here specifically so the plotting script can run
## standalone, fast, without recompiling anything).
##
## Same repo-root/build bootstrap as scripts/example_toy_fit.R.
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
cat("Using dMod", as.character(packageVersion("dMod")), "from", find.package("dMod"), "\n\n")

set.seed(1)
out_dir <- file.path(repo_root, "scripts", "plots", "pi3kmini_fim_fit")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(tempdir())
suppressMessages(library(dplyr))

## --- 1. Model ---------------------------------------------------------------
f <- NULL
f <- addReaction(f, from = "A", to = "B", rate = "k1*A", description = "A to B")
f <- addReaction(f, from = "", to = "U", rate = "k_gen*A", description = "generate U from A (A not consumed)")
f <- addReaction(f, from = "U", to = "", rate = "k_degU*U", description = "degrade U")
f <- addReaction(f, from = "Mp", to = "Mpp", rate = "k_use*Mp*U", description = "U-catalyzed conversion")
model <- odemodel(f, modelname = "pi3kMiniModel", solver = "CppODE")
x <- Xs(model)

## Only 3 observables -- U (hidden), A, B, Mp, Mpp states are NOT attached;
## see the header note above on why attach.input = TRUE would silently break
## the intended symmetry here.
observables <- eqnvec(obsA = "scale_A*A", obsB = "scale_B*B", obsMpp = "scale_Mpp*Mpp")
g <- Y(observables, x, compile = TRUE, modelname = "pi3kMiniObs", attach.input = FALSE)

## --- 2. Parameter transformation: 3 conditions differing only in U(0) ------
innerpars <- unique(c(getParameters(g * x), getParameters(observables)))
trafo <- repar("x~x", x = innerpars)
trafo <- repar("x~1", x = "A", trafo)    # A(0) = 1 (fixed reference, avoids A0<->scale_A/B symmetry)
trafo <- repar("x~0", x = "B", trafo)    # B(0) = 0
trafo <- repar("x~1", x = "Mp", trafo)   # Mp(0) = 1 (fixed reference, avoids Mp0<->scale_Mpp symmetry)
trafo <- repar("x~0", x = "Mpp", trafo)  # Mpp(0) = 0

conditions <- c("cond1", "cond2", "cond3")
trafoList <- branch(trafo, conditions = conditions)
for (i in seq_along(conditions))
  trafoList[[i]] <- insert(trafoList[[i]], paste0("x~x_", conditions[i]), x = "U")
trafoList <- insert(trafoList, "x~exp(x)", x = .currentSymbols)  # log-parametrize everything free

p <- P(trafoList, compile = TRUE, modelname = "pi3kMiniTrafo")
outerpars <- getParameters(p)

## --- 3. "True" parameters and simulated data -------------------------------
true_pars_inner <- c(k1 = 0.3, k_gen = 0.5, k_degU = 0.2, k_use = 0.4,
                      U_cond1 = 1, U_cond2 = 3, U_cond3 = 6,
                      scale_A = 1, scale_B = 1, scale_Mpp = 1)
pouter_true <- setNames(log(true_pars_inner[outerpars]), outerpars)

times <- seq(0, 15, by = 1.13)
prediction_true <- (g * x * p)(times, pouter_true)

sigma <- 0.05
sim <- as.data.frame(prediction_true)
sim$value <- sim$value + rnorm(nrow(sim), sd = sigma)
sim$sigma <- sigma
data <- as.datalist(sim, split.by = "condition")
obj_data <- normL2(data, g * x * p)

group_pars <- c("k_gen", "k_use", "U_cond1", "U_cond2", "U_cond3")
ident_pars <- c("k1", "k_degU", "scale_A", "scale_B", "scale_Mpp")

## --- 4. Unregularized multistart (ONE run) ---------------------------------
pouter_center <- setNames(rep(-1, length(outerpars)), outerpars)
N_FITS <- 200
fits_unreg <- mstrust(obj_data, pouter_center, studyname = "pi3kmini_unreg",
                       rinit = 1, rmax = 10, fits = N_FITS, cores = 4, iterlim = 200, sd = 1)
pf_unreg <- as.parframe(fits_unreg)
best_unreg <- as.parvec(pf_unreg, 1)
cat("=== Unregularized ===\n")
cat("n within 0.1 of best:", sum(pf_unreg$value < min(pf_unreg$value) + 0.1), "/", nrow(pf_unreg), "\n")
cat("Best value:", pf_unreg$value[1], "\n")

## --- 5. Reference Hessian, null-space check --------------------------------
## Expect a clean ~12-order-of-magnitude gap isolating exactly ONE direction,
## loading equally (+-0.4472 = 1/sqrt(5)) on U_cond1/U_cond2/U_cond3/k_gen
## (same sign) and k_use (opposite sign), zero everywhere else -- the exact
## analytic symmetry direction derived in the header, recovered numerically.
H <- obj_data(pouter_center)$hessian[outerpars, outerpars]
eig <- eigen(H, symmetric = TRUE)
cat("\nEigenvalues (log10):\n"); print(round(log10(abs(eig$values) + 1e-300), 2))
v_null <- eig$vectors[, which.min(eig$values)]; names(v_null) <- outerpars
cat("\nSmallest eigenvector loadings:\n"); print(round(v_null, 4))

source(file.path(repo_root, "fim-portugal", "functions-portugal.R"))
P_null <- nullSpaceProjector(H, tol = 1e-6)

## --- 6. L2 and FIM(group) regularized objectives, SAME center -------------
lambda <- 0.5
mu_center <- setNames(rep(-1, length(outerpars)), outerpars)
obj_l2  <- obj_data + constraintL2(mu = mu_center, sigma = 1 / sqrt(2 * lambda))
obj_fim <- obj_data + groupL2(P_null, mu_center, lambda = lambda)

## --- 7. L2 multistart (ONE run) --------------------------------------------
fits_l2 <- mstrust(obj_l2, pouter_center, studyname = "pi3kmini_l2",
                    rinit = 1, rmax = 10, fits = N_FITS, cores = 4, iterlim = 200, sd = 1)
pf_l2 <- as.parframe(fits_l2)
best_l2 <- as.parvec(pf_l2, 1)
cat("\n=== L2-regularized ===\n")
cat("n within 0.1 of best (penalized objective):",
    sum(pf_l2$value < min(pf_l2$value) + 0.1), "/", nrow(pf_l2), "\n")

## --- 8. FIM(group) multistart (ONE run) ------------------------------------
fits_fim <- mstrust(obj_fim, pouter_center, studyname = "pi3kmini_fim",
                     rinit = 1, rmax = 10, fits = N_FITS, cores = 4, iterlim = 200, sd = 1)
pf_fim <- as.parframe(fits_fim)
best_fim <- as.parvec(pf_fim, 1)
cat("\n=== FIM(group)-regularized ===\n")
cat("n within 0.1 of best (penalized objective):",
    sum(pf_fim$value < min(pf_fim$value) + 0.1), "/", nrow(pf_fim), "\n")

chi2_unreg <- obj_data(best_unreg[outerpars])$value
chi2_l2    <- obj_data(best_l2[outerpars])$value
chi2_fim   <- obj_data(best_fim[outerpars])$value
cat("\n=== raw DATA chi2 (unpenalized) at each method's best fit ===\n")
cat("unreg:", chi2_unreg, " l2:", chi2_l2, " fim:", chi2_fim, "\n")

## --- 9. Precompute fine-grid prediction curves at each best fit -----------
## Saved here (not recomputed by the plotting script) so that script never
## needs to touch dMod's compiled model at all.
times_fine <- seq(0, 15, by = 0.1)
pred_unreg <- as.data.frame((g * x * p)(times_fine, best_unreg[outerpars]))
pred_l2    <- as.data.frame((g * x * p)(times_fine, best_l2[outerpars]))
pred_fim   <- as.data.frame((g * x * p)(times_fine, best_fim[outerpars]))

## --- 10. Save everything the plotting script needs -------------------------
results <- list(
  pf_unreg = pf_unreg, pf_l2 = pf_l2, pf_fim = pf_fim,
  best_unreg = best_unreg, best_l2 = best_l2, best_fim = best_fim,
  chi2_unreg = chi2_unreg, chi2_l2 = chi2_l2, chi2_fim = chi2_fim,
  pred_unreg = pred_unreg, pred_l2 = pred_l2, pred_fim = pred_fim,
  sim = sim, times = times, times_fine = times_fine,
  H = H, eig = eig, P_null = P_null,
  true_pars_inner = true_pars_inner, pouter_true = pouter_true,
  outerpars = outerpars, group_pars = group_pars, ident_pars = ident_pars,
  lambda = lambda, mu_center = mu_center, sigma = sigma, N_FITS = N_FITS
)
results_path <- file.path(out_dir, paste0("pi3kmini_results_sigma", sigma, ".rds"))
saveRDS(results, results_path)
cat("\nSaved parframes/results to", results_path, "\n")
cat("Run example_toy_pi3kmini_fim_plots.R to generate plots from this file.\n")
