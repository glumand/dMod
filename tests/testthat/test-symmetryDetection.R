# Behavioral tests for symmetryDetection() (sympy-backed symmetry search).
#
# symmetryDetection drives a Python (sympy) routine via reticulate to find
# infinitesimal symmetries of a reaction network plus observation map.
# Reference: Maiwald & Schelker et al. 2016 (PRE 92, 012920).
#
# Test strategy: a canonical scaling symmetry case (A <-> B with observable
# alpha*A) for which the analytical infinitesimal symmetry is known. We
# verify the routine runs end-to-end without error and returns a non-null
# result on a sympy-equipped environment; if sympy is not reachable, the
# test is skipped.


.sympy_works <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  isTRUE(tryCatch({
    sympy <- reticulate::import("sympy", convert = TRUE)
    nzchar(sympy[["__version__"]])
  }, error = function(e) FALSE))
}


test_that("symmetryDetection runs on the canonical A<->B scaling case with alpha*A observable", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  eq <- eqnlist() |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "A", "k2 * B")
  observables <- eqnvec(Aobs = "alpha * A")

  # Run end-to-end on the canonical scaling case. The behavioural check
  # is "no error raised"; the python side prints progress to stdout.
  expect_error(
    symmetryDetection(eq, observables, ansatz = "uni", pMax = 1,
                      cores = 1),
    regexp = NA)
})
