#' Package initialization
#'
#' Declares the Python dependencies needed by SBML / PEtab import/export.
#' `reticulate::py_require()` does not start Python here — it just records
#' the requirement so the first call into Python (e.g. via `import_sbml()`)
#' provisions an env that has `python-libsbml` available. Users with an
#' existing libsbml install can bypass reticulate entirely by setting the
#' `DMOD_LIBSBML_PYTHON` environment variable; see `.dmod_libsbml_python()`
#' in `R/SBMLinterface.R`.
#'
#' @keywords internal
#' @importFrom reticulate py_require
#' @noRd
.onLoad <- function(libname, pkgname) {
  reticulate::py_require("python-libsbml")
}
