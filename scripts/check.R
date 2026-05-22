## scripts/check.R
## ---------------------------------------------------------------------------
## Local equivalent of the GitHub Actions R-CMD-check job. Use this before
## pushing so the CI run isn't your first signal that something's broken.
##
## Mirrors `.github/workflows/R-CMD-check.yaml`:
##   build_args = c("--no-manual", "--no-build-vignettes",
##                  "--compact-vignettes=gs+qpdf")
##   error_on   = "warning"   (CI fails on warnings, not just errors)
##
## Invocation (from any platform):
##   Rscript scripts/check.R                # full R CMD check
##   Rscript scripts/check.R --tests-only   # just devtools::test() (faster)
##
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
tests_only <- "--tests-only" %in% args

# Locate the script regardless of how it was invoked (Rscript, source(), ...)
.script_path <- (function() {
  full <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", full, value = TRUE)
  if (length(hit)) return(sub("^--file=", "", hit[1]))
  if (!is.null(sys.frames()) && length(sys.frames()) >= 1L) {
    of <- sys.frame(1)$ofile
    if (!is.null(of)) return(of)
  }
  NA_character_
})()

repo_root <- if (is.na(.script_path)) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(.script_path), ".."),
                winslash = "/", mustWork = TRUE)
}
setwd(repo_root)
cat("[check] repo root:", repo_root, "\n")

need <- c("devtools", "rcmdcheck", "pkgbuild", "testthat")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing))
  stop("Missing required packages: ", paste(missing, collapse = ", "),
       "\nInstall with: install.packages(c(\"",
       paste(missing, collapse = "\", \""), "\"))")

if (tests_only) {
  cat("[check] tests-only mode (devtools::test)\n")
  res <- devtools::test(stop_on_failure = TRUE)
  invisible(res)
  quit(status = 0L)
}

cat("[check] running rcmdcheck (mirrors CI flags)\n")
# Note: dMod is NOT a CRAN package (it depends on CppODE via Remotes:), so
# `--as-cran` is intentionally absent — its CRAN-submission checks would
# otherwise hard-warn on the GitHub remote and the missing Authors@R, neither
# of which is actionable for an internal package.
# `--ignore-vignettes` skips the package-vignettes check; the petab vignette
# needs the libsbml Python venv that CI doesn't provision.
# _R_CHECK_FORCE_SUGGESTS_=false: don't hard-fail when packages in Suggests
# are not installed locally. JuliaCall/openxlsx are heavyweight optional
# deps. CI installs Suggests via setup-r-dependencies.
res <- rcmdcheck::rcmdcheck(
  path        = ".",
  args        = c("--no-manual", "--ignore-vignettes"),
  build_args  = c("--no-manual", "--no-build-vignettes",
                  "--compact-vignettes=gs+qpdf"),
  error_on    = "warning",
  check_dir   = "check",
  env         = c("_R_CHECK_FORCE_SUGGESTS_" = "false")
)

cat("\n[check] errors  :", length(res$errors),
    "\n[check] warnings:", length(res$warnings),
    "\n[check] notes   :", length(res$notes), "\n")

if (length(res$errors) || length(res$warnings)) quit(status = 1L)
quit(status = 0L)
