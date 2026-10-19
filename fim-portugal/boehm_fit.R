## ---------------------------------------------------------------------------
## A REAL, published model as the FIM-null-space-regularization test case:
## Boehm_JProteomeRes2014, a JAK2/STAT5 phosphorylation model (PEtab
## benchmark, BenchmarkModels/Boehm_JProteomeRes2014) -- 8 dynamic states, 9
## estimated parameters (log10 scale, box [-5, 5]), 1 condition, 48 data
## points. Unlike fim-portugal/oscillator_fit.R (a synthetic, engineered
## worst case for MULTISTART itself), this script targets the OTHER failure
## mode the regularization in functions-portugal.R/notes/
## group_regularization_null_space.tex is built for: practical
## non-identifiability -- a near-flat direction in the Hessian at the fit,
## which lets many different parameter combinations reach essentially the
## same likelihood. That is a real, previously-documented property of this
## model (k_exp_hetero and k_imp_homo are known to be poorly identified from
## this data), not something engineered in for this script.
##
## Structure:
##   Part 1-2: import the PEtab model, unregularized multistart -- shows a
##     genuinely bad waterfall (most starts land near-optimal, but a
##     meaningful tail lands in distinct, much worse local optima) AND, more
##     to the point, a wide spread of PARAMETER values among the
##     near-optimal fits themselves (same likelihood, different parameters).
##   Part 3: Hessian eigenanalysis at the best raw fit identifies exactly
##     which directions are sloppy (near-zero eigenvalues) -- cross-checked
##     against which raw parameters actually showed the wide spread.
##   Part 4: groupL2-regularized multistart, same start distribution/box.
##     Compares waterfall and per-parameter spread against Part 2 -- the
##     regularization should tighten the sloppy parameters' spread by an
##     order of magnitude or more while leaving the achieved likelihood (the
##     RAW data objective, not the penalized one) essentially unchanged.
##   Part 5: summary table.
##   Part 6: profile likelihood, raw vs. groupL2. Verified empirically (not
##     assumed) that profile()'s TWO methods behave completely differently
##     here: method="optimize" (reoptimizes all other parameters at each
##     step, ignores the Hessian) already works fine on the RAW fit with no
##     regularization, and correctly shows k_exp_hetero/k_imp_homo swinging
##     over ~10 orders of magnitude for ~0 likelihood change while
##     well-identified parameters stay tightly bounded -- profiling this
##     model was never actually blocked. But the default, much cheaper
##     method="integrate" uses the FULL Hessian as its step metric, and
##     because that Hessian has 2 near-zero eigenvalues (from the same 2
##     parameters), it fails immediately for EVERY parameter, even
##     perfectly-identified ones -- an artifact of this model's Hessian
##     conditioning, not of any one parameter's identifiability. groupL2
##     fixes exactly that: the regularized Hessian has no near-zero
##     eigenvalues, method="integrate" starts working, and for
##     well-identified parameters it reproduces method="optimize"'s raw
##     profile almost exactly (much faster). For the two null-space
##     parameters themselves, the regularized method="integrate" profile
##     comes out artificially narrow -- the penalty's own curvature, not
##     real information -- so it must NOT be read as an identifiability
##     result there (same caveat as ab_fit.R Part 10: profile the RAW
##     objective, not the penalized one, whenever the question is
##     identifiability rather than optimizer convenience).
##
## This required fixing TWO real bugs in importPEtab() along the way (see
## R/PEtabInterface.R): its returned `obj` was a bare closure that (a)
## dropped `raw_obj`'s class/attributes, breaking mstrust() outright
## ("no applicable method for 'mname'"), and (b) discarded the PEtab-baked
## `fixed` values whenever a caller passed its OWN `fixed=` (as profile()
## always does to move the profiled parameter out of `pars`), breaking
## every reoptimization step with "Missing parameters: ratio, cyt, nuc,
## specC17". Both fixed; tests/testthat/test-petab.R still passes.
##
## PEtab import needs `python-libsbml` importable by SOME Python
## interpreter (see R/SBMLinterface.R's .dmod_libsbml_python() /
## vignettes/petab.Rmd) -- reticulate's own auto-provisioning needs network
## access and did not work in the sandbox this script was developed in,
## where `pip install --user --break-system-packages python-libsbml` into
## the system python3 was used instead (see below). Set
## DMOD_LIBSBML_PYTHON explicitly if your system python3 is not the right
## one.
##
## Builds/loads the dMod that lives in THIS checkout (DESCRIPTION Version,
## from repo_root) into a local, gitignored library under .Rlib/, so it
## never picks up whatever dMod happens to be installed in the default R
## library. Invoke as `Rscript fim-portugal/boehm_fit.R` from anywhere; the
## first run compiles the package, later runs reuse it. Same repo-root/
## build bootstrap as fim-portugal/oscillator_fit.R and ab_fit.R.
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

## If no libsbml-capable Python is already configured, try the system
## python3 -- covers the common case where `pip install --user
## python-libsbml` (or a system package) already satisfies it, without
## forcing reticulate's network-dependent auto-provisioning.
if (!nzchar(Sys.getenv("DMOD_LIBSBML_PYTHON"))) {
  candidate <- Sys.which("python3")
  if (nzchar(candidate)) {
    ok <- identical(suppressWarnings(system2(candidate, c("-c", shQuote("import libsbml")),
                                              stdout = FALSE, stderr = FALSE)), 0L)
    if (ok) Sys.setenv(DMOD_LIBSBML_PYTHON = candidate)
  }
}

set.seed(1)
out_dir <- file.path(repo_root, "fim-portugal")
setwd(tempdir())

source(file.path(out_dir, "functions-portugal.R"))

## --- 1. Import the PEtab model -----------------------------------------------
yaml <- file.path(repo_root, "BenchmarkModels", "Boehm_JProteomeRes2014",
                   "Boehm_JProteomeRes2014.yaml")
cat("Importing Boehm_JProteomeRes2014 (PEtab) ...\n")
t0 <- Sys.time()
pp <- importPEtab(yaml, solver = "CppODE", compile = TRUE, cores = 4)
cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t0, units = "secs")))

obj      <- pp$obj
pouter   <- pp$bestfit          # published nominal fit, used as the multistart center
parlower <- pp$parlower         # PEtab-declared bounds, log10 scale: [-5, 5] for all 9 pars
parupper <- pp$parupper

cat("\nEstimated parameters (", length(pouter), "):\n", sep = "")
print(pouter)
cat("\nn conditions:", length(pp$dataList),
    " n data points:", sum(vapply(pp$dataList, nrow, 0L)), "\n")

## --- 2. Unregularized multistart: bad waterfall + non-identifiable spread ----
## sd = 2 (log10 scale, box is +-5 wide) is wide enough that a meaningful
## fraction of starts land in distinct, clearly-worse local optima (a real
## bad waterfall, not engineered) while the near-optimal cluster remains
## large enough to characterize per-parameter spread. iterlim = 1000 is
## generous headroom -- fits either converge in ~200-300 iterations or land
## in a worse basin; verified this is not an iteration-budget artifact.
n_fits <- 200
sd_use <- 2
iterlim_use <- 1000

cat("\n=== Part 2: unregularized mstrust (", n_fits, " starts) ===\n", sep = "")
fits_raw <- mstrust(obj, pouter, studyname = "boehmRaw", rinit = 1, rmax = 10,
                     fits = n_fits, cores = min(10, parallel::detectCores()),
                     iterlim = iterlim_use, sd = sd_use,
                     parupper = parupper, parlower = parlower)
pf_raw <- as.parframe(fits_raw)
best_raw <- as.parvec(pf_raw, 1)

v_raw <- sort(pf_raw$value)
cat(sprintf("converged: %d/%d   best: %.4f   worst: %.4f\n",
            sum(pf_raw$converged), nrow(pf_raw), v_raw[1], v_raw[length(v_raw)]))
cat(sprintf("fraction within 2 of best: %.2f   within 10: %.2f\n",
            mean(v_raw < v_raw[1] + 2), mean(v_raw < v_raw[1] + 10)))

waterfall_raw <- plotValues(pf_raw) + ggplot2::ggtitle("Boehm: unregularized multistart")
waterfall_raw
waterfall_raw_path <- file.path(out_dir, "boehm_fit_waterfall_raw.png")
ggplot2::ggsave(waterfall_raw_path, waterfall_raw, width = 9, height = 5, dpi = 120)

## Per-parameter spread WITHIN the near-optimal cluster (value < best + 2):
## the whole point is that a tight likelihood band can still contain a wide
## band of parameter values along non-identifiable directions.
parnames <- names(pouter)
cluster_raw <- pf_raw[pf_raw$value < v_raw[1] + 2, ]
spread_raw <- setNames(vapply(parnames, function(p) sd(cluster_raw[[p]]), numeric(1)), parnames)
cat("\nPer-parameter sd within the near-optimal cluster (n=", nrow(cluster_raw), "), raw:\n", sep = "")
print(sort(spread_raw, decreasing = TRUE))

## plotPars()'s default view boxplots one representative fit per waterfall
## jump, which is the wrong summary here -- the point is the spread of
## parameter values WITHIN the near-optimal cluster itself (same
## likelihood, different parameters), so plot that cluster's raw values
## directly instead.
cluster_long_raw <- do.call(rbind, lapply(parnames, function(p)
  data.frame(name = p, parvalue = cluster_raw[[p]])))
parscatter_raw <- ggplot2::ggplot(cluster_long_raw, ggplot2::aes(x = name, y = parvalue)) +
  ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 1.2) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(x = NULL, y = "log10 parameter value",
                title = "Boehm: unregularized -- parameter values within the near-optimal cluster",
                subtitle = sprintf("%d fits with objective < best + 2", nrow(cluster_raw)))
parscatter_raw
parscatter_raw_path <- file.path(out_dir, "boehm_fit_parscatter_raw.png")
ggplot2::ggsave(parscatter_raw_path, parscatter_raw, width = 9, height = 6, dpi = 120)

## --- 3. Hessian eigenanalysis at the best raw fit ----------------------------
## Cross-check: do the near-zero eigenvalues' directions actually line up
## with the parameters that showed wide empirical spread above? (They
## should -- that correspondence is the entire premise of projecting the
## penalty onto the Hessian's null space instead of an ordinary ridge.)
H <- obj(best_raw)$hessian[names(best_raw), names(best_raw)]
basis <- nullSpaceBasis(H)
cat("\nHessian eigenvalues at best raw fit (decreasing):\n")
print(basis$values)
cat("\nNull/sloppy directions (", sum(basis$is_zero), " of ", length(basis$values), "):\n", sep = "")
print(round(basis$V_null, 3))

P <- nullSpaceProjector(H)

## --- 4. groupL2-regularized multistart, same start distribution/box ---------
## lambda = 1: chosen from a lambda in {0.01, 0.1, 1, 10} sweep during
## development -- 1 already tightens the sloppy parameters' cluster spread
## by roughly an order of magnitude with no visible cost to the best-found
## likelihood; 10 tightens further but starts to look like it is clamping
## the parameters to best_raw rather than letting the optimizer find them.
lambda_use <- 1
obj_group <- obj + groupL2(P, best_raw, lambda = lambda_use)

cat("\n=== Part 4: groupL2-regularized mstrust (lambda=", lambda_use, ", ",
    n_fits, " starts) ===\n", sep = "")
fits_group <- mstrust(obj_group, pouter, studyname = "boehmGroup", rinit = 1, rmax = 10,
                       fits = n_fits, cores = min(10, parallel::detectCores()),
                       iterlim = iterlim_use, sd = sd_use,
                       parupper = parupper, parlower = parlower)
pf_group <- as.parframe(fits_group)

## Score every regularized fit on the RAW data objective, not obj_group --
## groupL2 is an optimization aid, not part of the quantity we report
## (same reasoning as ab_fit.R Part 10: profiling/reporting must use the
## unpenalized likelihood or a genuinely flat direction would look
## artificially well-determined).
raw_vals_group <- vapply(seq_len(nrow(pf_group)), function(i) obj(as.parvec(pf_group, i))$value, numeric(1))
o <- order(raw_vals_group)
pf_group_byraw <- pf_group[o, ]
pf_group_byraw$value <- raw_vals_group[o]   # overwrite with the raw-objective value for plotting/reporting
v_group <- pf_group_byraw$value
best_group <- as.parvec(pf_group, o[1])

cat(sprintf("converged: %d/%d   best (raw obj): %.4f   worst (raw obj): %.4f\n",
            sum(pf_group$converged), nrow(pf_group), v_group[1], v_group[length(v_group)]))
cat(sprintf("fraction within 2 of best: %.2f   within 10: %.2f\n",
            mean(v_group < v_group[1] + 2), mean(v_group < v_group[1] + 10)))
cat("(The waterfall SHAPE is expected to look nearly identical to Part 2's --\n",
    "groupL2 targets spread WITHIN a basin, not which of several genuinely\n",
    "distinct basins a start converges to; the payoff shows up in the\n",
    "per-parameter spread below and in the parscatter plots, not here.)\n")

waterfall_group <- plotValues(pf_group_byraw) +
  ggplot2::ggtitle("Boehm: groupL2-regularized multistart (raw-objective values)")
waterfall_group
waterfall_group_path <- file.path(out_dir, "boehm_fit_waterfall_group.png")
ggplot2::ggsave(waterfall_group_path, waterfall_group, width = 9, height = 5, dpi = 120)

cluster_group <- pf_group_byraw[v_group < v_group[1] + 2, ]
spread_group <- setNames(vapply(parnames, function(p) sd(cluster_group[[p]]), numeric(1)), parnames)
cat("\nPer-parameter sd within the near-optimal cluster (n=", nrow(cluster_group), "), groupL2-regularized:\n", sep = "")
print(sort(spread_group, decreasing = TRUE))

cluster_long_group <- do.call(rbind, lapply(parnames, function(p)
  data.frame(name = p, parvalue = cluster_group[[p]])))
parscatter_group <- ggplot2::ggplot(cluster_long_group, ggplot2::aes(x = name, y = parvalue)) +
  ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 1.2, color = "darkred") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(x = NULL, y = "log10 parameter value",
                title = "Boehm: groupL2-regularized -- parameter values within the near-optimal cluster",
                subtitle = sprintf("%d fits with (raw-objective) value < best + 2", nrow(cluster_group)))
parscatter_group
parscatter_group_path <- file.path(out_dir, "boehm_fit_parscatter_group.png")
ggplot2::ggsave(parscatter_group_path, parscatter_group, width = 9, height = 6, dpi = 120)

## --- 6. Profile likelihood: does groupL2 actually help? -----------------------
## Tested empirically rather than assumed. profile()'s two methods behave
## completely differently on this model:
##   - method="optimize" reoptimizes every OTHER parameter at each step and
##     never touches the Hessian -- already works fine on the RAW,
##     unregularized fit, and correctly shows the two null-space parameters
##     swinging over many orders of magnitude for ~0 likelihood change.
##   - method="integrate" (the default, much cheaper -- no reoptimization,
##     just a Hessian-guided predictor/corrector step) uses the FULL 9x9
##     Hessian as its step metric, so its 2 near-zero eigenvalues (from
##     k_exp_hetero/k_imp_homo) break it immediately for EVERY parameter,
##     including perfectly-identified ones -- a global conditioning
##     problem, not a per-parameter one.
## groupL2 raises those 2 eigenvalues to 2*lambda (Part 3/4's P, added to
## the Hessian), which is exactly what makes method="integrate" usable
## again. Checked on 2 null-space parameters and 2 well-identified controls.
profile_pars <- c("k_exp_hetero", "k_imp_homo", "k_phos", "Epo_degradation_BaF3")
profile_limits <- c(-5, 5)

cat("\n=== Part 6: profile likelihood -- method x regularization ===\n")
profile_results <- list()
for (wp in profile_pars) {

  cat(sprintf("\n--- %s ---\n", wp))

  t0 <- Sys.time()
  prof_opt_raw <- tryCatch(
    profile(obj, best_raw, whichPar = wp, limits = profile_limits, method = "optimize", cores = 1),
    error = function(e) { cat("  optimize/raw FAILED:", conditionMessage(e), "\n"); NULL })
  t_opt_raw <- as.numeric(Sys.time() - t0, units = "secs")

  t0 <- Sys.time()
  prof_int_raw <- tryCatch(
    profile(obj, best_raw, whichPar = wp, limits = profile_limits, method = "integrate", cores = 1),
    error = function(e) { cat("  integrate/raw FAILED:", conditionMessage(e), "\n"); NULL })
  t_int_raw <- as.numeric(Sys.time() - t0, units = "secs")

  t0 <- Sys.time()
  prof_int_group <- tryCatch(
    profile(obj_group, best_raw, whichPar = wp, limits = profile_limits, method = "integrate", cores = 1),
    error = function(e) { cat("  integrate/groupL2 FAILED:", conditionMessage(e), "\n"); NULL })
  t_int_group <- as.numeric(Sys.time() - t0, units = "secs")
  ## Score the regularized-Hessian profile path on the RAW objective --
  ## obj_group's own $value column bakes in the penalty, which is not the
  ## quantity of interest (same reasoning as Part 4).
  if (!is.null(prof_int_group))
    prof_int_group$value <- vapply(seq_len(nrow(prof_int_group)),
                                    function(i) obj(unlist(prof_int_group[i, names(best_raw)]))$value,
                                    numeric(1))

  summarize <- function(p, t) if (is.null(p)) "FAILED" else
    sprintf("n=%3d  par=[%7.3f,%7.3f]  value=[%7.2f,%7.2f]  %5.1fs",
            nrow(p), min(p[[wp]]), max(p[[wp]]), min(p$value), max(p$value), t)
  cat("  optimize,  raw obj    :", summarize(prof_opt_raw,   t_opt_raw),   "\n")
  cat("  integrate, raw obj    :", summarize(prof_int_raw,   t_int_raw),   "\n")
  cat("  integrate, groupL2 obj:", summarize(prof_int_group, t_int_group), "(raw-obj value shown)\n")

  profile_results[[wp]] <- list(opt_raw = prof_opt_raw, int_raw = prof_int_raw, int_group = prof_int_group)
}

## Plot: for one null-space parameter and one well-identified control,
## overlay the ground-truth (optimize, raw) profile against the
## regularized-Hessian (integrate, groupL2) one -- expect near-perfect
## overlap for the identified parameter and a visibly too-narrow
## regularized curve for the null-space one.
plot_profile_comparison <- function(wp) {
  res <- profile_results[[wp]]
  rows <- list()
  if (!is.null(res$opt_raw))
    rows[[length(rows) + 1]] <- data.frame(par = res$opt_raw[[wp]], value = res$opt_raw$value,
                                            source = "optimize, raw obj (ground truth)")
  if (!is.null(res$int_group))
    rows[[length(rows) + 1]] <- data.frame(par = res$int_group[[wp]], value = res$int_group$value,
                                            source = "integrate, groupL2 obj (raw-obj value)")
  d <- do.call(rbind, rows)
  ggplot2::ggplot(d, ggplot2::aes(x = par, y = value, color = source)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = obj(best_raw)$value + qchisq(0.95, 1), linetype = "dashed") +
    ggplot2::coord_cartesian(ylim = c(obj(best_raw)$value - 1, obj(best_raw)$value + 20)) +
    ggplot2::labs(x = wp, y = "raw data objective", color = NULL,
                  title = paste0("Boehm profile: ", wp)) +
    ggplot2::theme(legend.position = "bottom")
}

profile_plot_identified <- plot_profile_comparison("k_phos")
profile_plot_identified
profile_plot_identified_path <- file.path(out_dir, "boehm_fit_profile_k_phos.png")
ggplot2::ggsave(profile_plot_identified_path, profile_plot_identified, width = 7, height = 5, dpi = 120)

profile_plot_sloppy <- plot_profile_comparison("k_exp_hetero")
profile_plot_sloppy
profile_plot_sloppy_path <- file.path(out_dir, "boehm_fit_profile_k_exp_hetero.png")
ggplot2::ggsave(profile_plot_sloppy_path, profile_plot_sloppy, width = 7, height = 5, dpi = 120)

cat("\nProfile plots saved:\n  ", profile_plot_identified_path, "\n  ", profile_plot_sloppy_path, "\n", sep = "")

## --- 7. Summary: raw vs. groupL2-regularized ---------------------------------
cat("\n=== Summary: raw vs. groupL2-regularized (lambda=", lambda_use, ") ===\n", sep = "")
summary_tab <- data.frame(
  quantity = c("best (raw objective)", "worst (raw objective)",
               "fraction within 2 of best", "fraction within 10 of best",
               paste0("sd(", names(sort(spread_raw, decreasing = TRUE))[1], ") in cluster"),
               paste0("sd(", names(sort(spread_raw, decreasing = TRUE))[2], ") in cluster")),
  raw = c(v_raw[1], v_raw[length(v_raw)],
          mean(v_raw < v_raw[1] + 2), mean(v_raw < v_raw[1] + 10),
          sort(spread_raw, decreasing = TRUE)[1], sort(spread_raw, decreasing = TRUE)[2]),
  groupL2 = c(v_group[1], v_group[length(v_group)],
              mean(v_group < v_group[1] + 2), mean(v_group < v_group[1] + 10),
              spread_group[names(sort(spread_raw, decreasing = TRUE))[1]],
              spread_group[names(sort(spread_raw, decreasing = TRUE))[2]])
)
print(summary_tab, row.names = FALSE)

cat("\nPlots saved:\n  ", waterfall_raw_path, "\n  ", parscatter_raw_path,
    "\n  ", waterfall_group_path, "\n  ", parscatter_group_path,
    "\n  ", profile_plot_identified_path, "\n  ", profile_plot_sloppy_path, "\n", sep = "")
