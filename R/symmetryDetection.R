
#' Search for symmetries in the loaded model
#'
#' @description This function follows the method published in \[1\]. The
#'   underlying Python implementation (originally written for Python 2.7 /
#'   sympy 0.7) was ported to Python 3 and is now driven through
#'   `reticulate` instead of the archived `rPython` package. The Python
#'   environment is provisioned via `reticulate::py_require()` in
#'   `.onLoad()` and pulls in `sympy`, `scipy`, and `numpy`.
#'
#' @param f object containing the ODE for which `as.eqnvec()` is defined
#' @param obsvect vector of observation functions
#' @param prediction vector containing prediction to be tested
#' @param initial vector containing initial values
#' @param ansatz type of infinitesimal ansatz used for the analysis (uni, par, multi)
#' @param pMax maximal degree of infinitesimal ansatz
#' @param inputs specify the input variables
#' @param fixed variables to concider fixed
#' @param cores maximal number of cores used for the analysis
#' @param allTrafos do not remove transformations with a common parameter factor
#' @return NULL
#'
#' @references \[1\]
#' <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.012920>
#'
#' @examples
#' \dontrun{
#' eq <- NULL
#' eq <- addReaction(eq, "A", "B", "k1*A")
#' eq <- addReaction(eq, "B", "A", "k2*B")
#'
#' observables <- eqnvec(Aobs = "alpha * A")
#'
#' symmetryDetection(eq, observables)
#'
#' }
#' @export
symmetryDetection <- function(f, obsvect = NULL, prediction = NULL,
                              initial = NULL, ansatz = 'uni', pMax = 2, inputs = NULL, fixed = NULL,
                              cores = 1, allTrafos = FALSE){

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for symmetryDetection().")
  }

  f <- as.eqnvec(f)

  f <- as.character(lapply(1:length(f), function(i)
    paste(names(f)[i],'=',f[i])))

  obsvect <- as.character(lapply(1:length(obsvect), function(i)
    paste(names(obsvect)[i],'=',obsvect[i])))

  if (!is.null(prediction)) {
    prediction <- as.character(lapply(1:length(prediction), function(i)
      paste(names(prediction)[i],'=',prediction[i])))
  }

  if (!is.null(initial)) {
    initial <- as.character(lapply(1:length(initial), function(i)
      paste(names(initial)[i],'=',initial[i])))
  }

  # Make the script directory importable, then load the entry-point module.
  # The Python sources import each other unqualified (`from functions import *`),
  # so they must live on sys.path; using reticulate::source_python() would not
  # achieve that for the transitive imports.
  code_dir <- system.file("code", package = "dMod")
  sys <- reticulate::import("sys", convert = TRUE)
  if (!(code_dir %in% sys$path)) sys$path <- c(code_dir, sys$path)

  sd <- reticulate::import("symmetryDetection", convert = TRUE)

  sd$symmetryDetectiondMod(f, obsvect, prediction, initial,
                           ansatz, as.integer(pMax),
                           inputs, fixed,
                           as.integer(cores), allTrafos)
}
