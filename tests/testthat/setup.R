# Run by testthat before any test file. Snapshot cwd and restore it at the
# end of the run so wayward setwd() inside individual tests cannot break
# testthat's relative-path lookups for the rest of the session.
.dmod_initial_wd <- getwd()
if (requireNamespace("withr", quietly = TRUE)) {
  withr::defer(setwd(.dmod_initial_wd), testthat::teardown_env())
}

# Resolve fixture directories that live OUTSIDE the installed package
# (PEtabTests/ and BenchmarkModels/ are .Rbuildignore'd because of their
# size). Tests look them up via DMOD_PETABTESTS / DMOD_BENCHMARKMODELS env
# vars; we set those here once, before any test_that block runs `setwd()`
# / `withr::local_dir(tempdir())` and breaks cwd-relative discovery.
#
# Walk up from the current wd looking for a dir of the given name. Caps at
# 8 levels so a missing fixture cannot loop to filesystem root forever.
.dmod_find_fixture <- function(name, start = .dmod_initial_wd) {
  here <- normalizePath(start, mustWork = FALSE, winslash = "/")
  for (i in seq_len(8)) {
    cand <- file.path(here, name)
    if (dir.exists(cand)) return(normalizePath(cand, winslash = "/"))
    parent <- dirname(here)
    if (parent == here) break
    here <- parent
  }
  ""
}

local({
  for (spec in list(
    c(env = "DMOD_PETABTESTS",     dir = "PEtabTests"),
    c(env = "DMOD_BENCHMARKMODELS", dir = "BenchmarkModels")
  )) {
    if (!nzchar(Sys.getenv(spec[["env"]], unset = ""))) {
      hit <- .dmod_find_fixture(spec[["dir"]])
      if (nzchar(hit)) do.call(Sys.setenv, setNames(list(hit), spec[["env"]]))
    }
  }
})
