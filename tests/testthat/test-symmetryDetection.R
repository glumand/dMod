context("symmetryDetection")

# Probe whether reticulate can materialise sympy/scipy/numpy. CI provisions
# the env via reticulate::py_require() in .onLoad(), but on a fresh worker
# the first call triggers a Python+wheel download — skip the test there
# rather than fail it.
.sympy_works <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  isTRUE(tryCatch({
    sympy <- reticulate::import("sympy", convert = TRUE)
    nzchar(sympy[["__version__"]])
  }, error = function(e) FALSE))
}

test_that("symmetryDetection finds the scaling symmetry of A<->B with Aobs = alpha*A", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "A", "B", "k1*A")
  reactions <- addReaction(reactions, "B", "A", "k2*B")
  observables <- eqnvec(Aobs = "alpha * A")

  # The Python implementation prints results rather than returning them.
  # Python's sys.stdout writes directly to the OS fd, bypassing R's text
  # connections — base::capture.output() cannot see it. reticulate's
  # py_capture_output redirects sys.stdout to a StringIO for the duration
  # of `expr` and returns the captured text.
  out <- reticulate::py_capture_output({
    symmetryDetection(reactions, observables)
  })

  # For this overparameterised observation, at least one transformation
  # must be found and the type must be 'scaling'.
  expect_match(out, "transformation\\(s\\) found")
  expect_match(out, "scaling", fixed = TRUE)
})
