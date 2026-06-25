#' Package initialization
#'
#' Declares the Python dependencies needed by SBML / PEtab import/export and
#' by `symmetryDetection()`. `reticulate::py_require()` does not start Python
#' here — it just records the requirement so the first call into Python (e.g.
#' via `import_sbml()` or `symmetryDetection()`) provisions an env that has
#' the listed packages available. Users with an existing libsbml install can
#' bypass reticulate entirely by setting the `DMOD_LIBSBML_PYTHON` environment
#' variable; see `.dmod_libsbml_python()` in `R/SBMLinterface.R`.
#'
#' @keywords internal
#' @importFrom reticulate py_require
#' @importFrom stats setNames predict density median dnorm
#' @importFrom utils globalVariables
#' @noRd
.onLoad <- function(libname, pkgname) {
  reticulate::py_require(c("python-libsbml", "sympy", "scipy", "numpy", "symengine"))
}
