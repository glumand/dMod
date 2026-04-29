## inst/examples/runPEtabBenchmarks.R
## ---------------------------------------------------------------------------
## Run dMod's PEtab importer against the four bundled benchmark models in
## BenchmarkModels/ and report the negative log-likelihood at each problem's
## published-optimum nominalValue. Self-contained — runs interactively
## (`source(...)`) or from the shell:
##
##   Rscript inst/examples/runPEtabBenchmarks.R                # all four
##   Rscript inst/examples/runPEtabBenchmarks.R Boehm          # one (prefix-match)
##   Rscript inst/examples/runPEtabBenchmarks.R Boehm Elowitz  # two
##
## Env vars:
##   DMOD_BENCHMARKMODELS  override path to BenchmarkModels/
##                         (default: <repo-root>/BenchmarkModels, walked up
##                         from this script's location)
##   DMOD_BENCH_KEEP       set to "1" to keep generated .c/.cpp/.so artifacts
##                         (default: cleaned up after each model)
##
## Note: large models (e.g. Zheng_PNAS2012) can take several minutes to
## compile sensitivity equations. There's no in-script timeout because R's
## setTimeLimit doesn't reach gcc/g++ via system2(); use the shell's
## `timeout` if you want a hard wall-clock cap, e.g.:
##   timeout 300 Rscript inst/examples/runPEtabBenchmarks.R Zheng
##
## Each model runs in its own tempdir; failures are caught and reported per
## model so one broken model doesn't kill the suite. The script exits 0 if
## every requested model produced a finite -log L (within tolerance of the
## published optimum, when known); otherwise exits 1.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(dMod))

## --- locate BenchmarkModels/ ----------------------------------------------

.bench_dir <- function() {
  envp <- Sys.getenv("DMOD_BENCHMARKMODELS", unset = "")
  if (nzchar(envp) && dir.exists(envp)) return(normalizePath(envp))
  # When sourced via `source()`, sys.frames is populated; when run via
  # Rscript, we walk up from the script path in commandArgs(trailingOnly=FALSE).
  here <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grepl("^--file=", args)])
    if (length(f)) normalizePath(f) else
      normalizePath(sys.frames()[[1L]]$ofile)
  }, error = function(e) "")
  candidates <- c(
    if (nzchar(here)) file.path(dirname(here), "..", "..", "BenchmarkModels"),
    file.path(getwd(), "BenchmarkModels"),
    file.path(getwd(), "..", "BenchmarkModels"),
    "/home/simon/Documents/Projects/dMod/BenchmarkModels"
  )
  for (p in candidates) if (nzchar(p) && dir.exists(p)) return(normalizePath(p))
  stop("Could not locate BenchmarkModels/. Set DMOD_BENCHMARKMODELS.")
}


## --- benchmark catalogue --------------------------------------------------
##
## Reference values are -log L at the *published* parameter optimum
## (parameters.tsv `nominalValue` column). Source: Hass et al., Bioinformatics
## 35(17):3073-3082 (2019), "Benchmark problems for dynamic modeling of
## intracellular processes", Suppl. Tab. S2 — and verified against dMod once
## the importer matched the published value. Set ref_nll = NA for models
## where we have no validated reference yet; the script will then just print
## the value rather than checking it.

.benchmarks <- list(
  Boehm_JProteomeRes2014 = list(
    yaml = "Boehm_JProteomeRes2014.yaml",
    ref_nll = 138.22,
    tol = 0.5,
    notes = "JAK/STAT5; <power/> + <assignmentRule> for BaF3_Epo + log10 pars"
  ),
  Elowitz_Nature2000 = list(
    yaml = "Elowitz_Nature2000.yaml",
    ref_nll = NA_real_,
    tol = NA_real_,
    notes = "Repressilator; uses ln(2) (libsbml L3 -> rewritten to log)"
  ),
  Fujita_SciSignal2010 = list(
    yaml = "Fujita_SciSignal2010.yaml",
    ref_nll = NA_real_,
    tol = NA_real_,
    notes = "EGF stimulation; uses <piecewise> (NOT yet supported by importer)"
  ),
  Zheng_PNAS2012 = list(
    yaml = "Zheng_PNAS2012.yaml",
    ref_nll = NA_real_,
    tol = NA_real_,
    notes = "Phosphorylation grid; uses <functionDefinition> (auto-inlined)"
  )
)


## --- per-model runner -----------------------------------------------------

.run_one <- function(model_id, spec, bench_dir) {

  yaml_path <- file.path(bench_dir, model_id, spec$yaml)
  if (!file.exists(yaml_path)) {
    return(list(model = model_id, status = "missing", error =
                paste("YAML not found:", yaml_path)))
  }

  wd <- tempfile(paste0("bench_", model_id, "_"))
  dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
  setwd(wd)

  cat(sprintf("\n=== %s ===\n", model_id))
  cat(sprintf("  notes: %s\n", spec$notes))
  cat(sprintf("  workdir: %s\n", wd))

  modelname <- paste0("bench_", gsub("[^A-Za-z0-9_]", "_", model_id))
  t0 <- Sys.time()
  petab <- tryCatch(
    importPEtab(yaml_path, solver = "deSolve", modelname = modelname),
    error = function(e) e)
  if (inherits(petab, "error")) {
    return(list(model = model_id, status = "import_failed",
                error = conditionMessage(petab),
                t_total = as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  t_imp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  cat(sprintf("  imported in %.1fs:  pouter=%d  fixed=%d  obs=%d  cond=%d  meas=%d\n",
              t_imp, length(petab$pouter), length(petab$fixed),
              length(petab$observables), nrow(petab$sub_cond_map),
              sum(vapply(petab$data, nrow, 0L))))

  t1 <- Sys.time()
  out <- tryCatch(
    petab$obj(petab$pouter, fixed = petab$fixed, deriv = FALSE),
    error = function(e) e)
  if (inherits(out, "error")) {
    return(list(model = model_id, status = "obj_failed",
                error = conditionMessage(out),
                t_total = as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  t_obj <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  nll <- out$value / 2
  cat(sprintf("  -2 log L = %.6f   (-log L = %.6f)   t_obj = %.2fs\n",
              out$value, nll, t_obj))

  status <- "ok"
  if (!is.na(spec$ref_nll)) {
    diff <- abs(nll - spec$ref_nll)
    cat(sprintf("  ref     = %.4f   |Δ| = %.4f   tol = %.4f   %s\n",
                spec$ref_nll, diff, spec$tol,
                if (diff < spec$tol) "PASS" else "FAIL"))
    if (diff >= spec$tol) status <- "ref_mismatch"
  } else {
    cat("  ref     = (none recorded — pin this value if you trust it)\n")
  }

  list(model = model_id, status = status,
       value = out$value, nll = nll, ref = spec$ref_nll,
       t_import = t_imp, t_obj = t_obj,
       t_total = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}


## --- main -----------------------------------------------------------------

.cleanup_artifacts <- function() {
  for (pat in c("bench_*", "*.c", "*.cpp", "*.o", "*.so", "*_model.csv",
                "reactions_for_*"))
    unlink(pat)
}

.main <- function(selected = character(0)) {

  bench_dir <- .bench_dir()
  cat(sprintf("BenchmarkModels/ at %s\n", bench_dir))

  keep_art <- identical(Sys.getenv("DMOD_BENCH_KEEP", unset = ""), "1")

  to_run <- if (length(selected) == 0L) names(.benchmarks) else {
    matches <- unique(unlist(lapply(selected, function(q)
      grep(q, names(.benchmarks), ignore.case = TRUE, value = TRUE))))
    if (length(matches) == 0L)
      stop("No benchmark matches: ", paste(selected, collapse = ", "),
           "\nAvailable: ", paste(names(.benchmarks), collapse = ", "))
    matches
  }

  cat(sprintf("Will run %d model(s): %s\n",
              length(to_run), paste(to_run, collapse = ", ")))

  results <- lapply(to_run, function(m) {
    r <- .run_one(m, .benchmarks[[m]], bench_dir)
    if (!keep_art) .cleanup_artifacts()
    r
  })

  cat("\n", strrep("=", 72), "\n", sep = "")
  cat("SUMMARY\n")
  cat(strrep("=", 72), "\n", sep = "")
  fmt <- "  %-28s  %-14s  %12s  %10s  %s\n"
  cat(sprintf(fmt, "model", "status", "-log L", "ref", "t_total"))
  cat(sprintf(fmt, strrep("-", 28), strrep("-", 14),
              strrep("-", 12), strrep("-", 10), strrep("-", 8)))
  for (r in results) {
    cat(sprintf(fmt,
                r$model,
                r$status,
                if (is.null(r$nll)) "—" else sprintf("%.4f", r$nll),
                if (is.null(r$ref) || is.na(r$ref)) "—" else sprintf("%.4f", r$ref),
                if (is.null(r$t_total)) "—" else sprintf("%.1fs", r$t_total)))
    if (!is.null(r$error))
      cat(sprintf("      ↳ %s\n", substr(r$error, 1, 200)))
  }

  ok <- vapply(results, function(r) r$status == "ok", logical(1))
  cat(sprintf("\n%d/%d passed.\n", sum(ok), length(ok)))
  invisible(results)
}


## --- entrypoint -----------------------------------------------------------

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  res <- .main(args)
  ok <- vapply(res, function(r) r$status == "ok", logical(1))
  if (!all(ok)) quit(status = 1L)
} else {
  ## Interactive use:
  ##   results <- .main()                # all four
  ##   results <- .main("Boehm")         # one (prefix-matched)
  cat("Loaded. Call: results <- .main()  -- or  .main(\"Boehm\", \"Elowitz\")\n")
}
