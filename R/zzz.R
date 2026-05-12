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
#' @importFrom stats setNames predict density median
#' @importFrom utils globalVariables
#' @noRd
.onLoad <- function(libname, pkgname) {
  reticulate::py_require(c("python-libsbml", "sympy", "scipy", "numpy"))

  # C++/OpenMP objective-function defaults.
  #   dMod.objfn.cpp     TRUE  -> normL2 / constraintL2 / datapointL2 use the
  #                              C++ kernels by default. Flip to FALSE to fall
  #                              back to the R reference path.
  #   dMod.objfn.threads 1     -> per-call OpenMP thread count for the C++
  #                              kernels. Multistart / profile orchestrators
  #                              (mstrust, profile) set this to 1 before fork
  #                              to avoid oversubscription; user can opt into
  #                              inner threading by overriding it.
  if (is.null(getOption("dMod.objfn.cpp")))
    options(dMod.objfn.cpp = TRUE)
  if (is.null(getOption("dMod.objfn.threads")))
    options(dMod.objfn.threads = 1L)
}
