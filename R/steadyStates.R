#' Calculate analytical steady states
#'
#' Computes symbolic steady-state expressions tailored to parameter estimation
#' via the AlyssaPetit method (see references). Calls a Python script via
#' reticulate (Python 3.x). Three solver versions are selectable via the
#' \code{version} argument; see that parameter for details.
#'
#' @param model Either name of the csv-file or the eqnlist of the model.
#' @param file Name of the file to which the steady-state equations are saved.
#' @param rates Character vector, flux vector of the system
#' @param forcings Character vector with the names of the forcings
#' @param givenCQs (Unnamed) Character vector with conserved quantities. Use the
#'   format c("A + pA = totA", "B + pB = totB"). The format c("A + pA", "B +
#'   pB") works also. If NULL, conserved quantities are automatically calculated.
#' @param neglect Character vector with names of states and parameters that must
#'   not be used for solving the steady-state equations
#' @param sparsifyLevel Numeric. Upper bound for length of linear combinations
#'   used to simplify the stoichiometric matrix. Used by versions \code{"1.0"}
#'   and \code{"1.1"}; ignored by \code{"1.2"}.
#' @param outputFormat Define the output format. By default "R" generating dMod
#'   compatible output. To obtain an output appropriate for d2d \[3\] "M" must be
#'   selected.
#' @param testSteady Character, if "T" the correctness of the obtained steady
#'   states is numerically checked (this can be very time intensive). If "F" this
#'   is skipped. Ignored for version "1.0" (always tests).
#' @param walltime Integer, wall-clock budget in seconds for the solver
#'   (default `0` = unlimited). Version "1.2" only.
#' @param simplify Final-simplification mode. One of `TRUE` (default,
#'   `sympy.simplify` once per expression), `FALSE` (skip; fastest, bulkier
#'   output), or `"full"` (aggressive `cancel`/`posify`/`simplify`/`factor`
#'   pipeline; slower but more compact). Version "1.2" only.
#' @param solveQuadratic Logical. When `TRUE`, the solver attempts a closed-form
#'   quadratic state-side resolution (positive root of `a*X^2 + b*X + c = 0`)
#'   before falling back to a flux-parameter pivot for any cycle whose final
#'   ODE is quadratic-in-self after upstream linear substitutions. Keeps the
#'   pivoted flux parameter (and the chain of substituted upstream rate
#'   constants) out of `mysteadies`, at the cost of emitting `sqrt(...)`
#'   expressions. Default `FALSE`; some workflows cannot consume `sqrt(...)`
#'   in their parameter trafos. Version "1.2" only.
#' @param version Character, AlyssaPetit backend version. One of
#'   \code{"1.0"} (original), \code{"1.1"} (adds \code{testSteady}), or
#'   \code{"1.2"} (default; sink-cluster detection, \code{walltime},
#'   priority-table cycle breaking, end-of-pipeline \code{simplify} toggle,
#'   optional quadratic state-side solve).
#'
#' @return Named character vector of steady-state equations (dMod compatible).
#'
#' @references \[1\]
#' <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4863410/>
#' @references \[2\]
#' <https://github.com/marcusrosenblatt/AlyssaPetit>
#' @references \[3\]
#' <https://github.com/Data2Dynamics/d2d>
#'
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#'
#' @export
#' @importFrom utils write.table
#' @importFrom reticulate import_from_path py_require
#' @example inst/examples/steadystates.R
steadyStates <- function(model, file = NULL, rates = NULL, forcings = NULL,
                         givenCQs = NULL, neglect = NULL, sparsifyLevel = NULL,
                         outputFormat = "R", testSteady = "T",
                         walltime = 0L, simplify = TRUE, solveQuadratic = FALSE,
                         version = "1.2") {

  # Validate version
  version <- match.arg(version, choices = c("1.0", "1.1", "1.2"))

  # Default sparsifyLevel depends on version: v1.0/v1.1 still use sparsify
  # (default 2), v1.2 ignores it (default 0, avoids the info print).
  if (is.null(sparsifyLevel)) {
    sparsifyLevel <- if (version == "1.2") 0 else 2
  }

  # Check if model is an equation list
  if (inherits(model, "eqnlist")) {
    if (is.null(file)) file <- "reactions_for_Alyssa"
    write.eqnlist(model, file = paste0(file, "_model.csv"))
    model <- paste0(file, "_model.csv")
  }
  if (!is.null(givenCQs) && length(names(givenCQs)) > 0)
    stop("givenCQs must not have names. Please unname() them.")

  # Ensure Python dependencies are available
  reticulate::py_require("numpy")
  reticulate::py_require("sympy")
  # v1.2 uses scipy.optimize.linprog for structural sink-cluster detection
  # (states whose combined mass leaks monotonically and must therefore be 0).
  # Older versions don't need scipy.
  if (version == "1.2") reticulate::py_require("scipy")

  pymodule <- paste0("AlyssaPetit_ver", gsub("\\.", "_", version))
  ap <- reticulate::import_from_path(pymodule,
                                     path = system.file("code", package = "dMod"))

  # Version-specific Python signatures:
  #   v1.0: Alyssa(filename, injections, givenCQs, neglect, sparsifyLevel, outputFormat)
  #   v1.1: Alyssa(filename, injections, givenCQs, neglect, sparsifyLevel, outputFormat, testSteady)
  #   v1.2: Alyssa(filename, injections, givenCQs, neglect, sparsifyLevel, outputFormat, testSteady, walltime, simplify, solveQuadratic)
  #        — v1.2 additionally runs structural sink-cluster detection a priori,
  #          and (when `solveQuadratic=TRUE`) attempts a closed-form quadratic
  #          state-side solve before resorting to flux-parameter pivots.
  if (version == "1.0") {
    if (testSteady == "F")
      message("Note: version 1.0 does not support testSteady='F', test will always run.")
    if (isTRUE(solveQuadratic))
      message("Note: version 1.0 does not support solveQuadratic=TRUE, ignored.")
    m_ss <- ap$Alyssa(model, as.list(forcings), as.list(givenCQs),
                      as.list(neglect), sparsifyLevel, outputFormat)

  } else if (version == "1.1") {
    if (isTRUE(solveQuadratic))
      message("Note: version 1.1 does not support solveQuadratic=TRUE, ignored.")
    m_ss <- ap$Alyssa(model, as.list(forcings), as.list(givenCQs),
                      as.list(neglect), sparsifyLevel, outputFormat,
                      testSteady)

  } else {
    # v1.2
    # simplify can be TRUE / FALSE / "full" — pass through untouched so the
    # Python side sees either a Python bool or the literal string "full".
    if (is.character(simplify)) {
      simplify <- match.arg(tolower(simplify), choices = "full")
      simplify_arg <- simplify
    } else {
      simplify_arg <- as.logical(simplify)
    }
    m_ss <- ap$Alyssa(model,
                      injections     = as.list(forcings),
                      givenCQs       = as.list(givenCQs),
                      neglect        = as.list(neglect),
                      sparsifyLevel  = as.integer(sparsifyLevel),
                      outputFormat   = outputFormat,
                      testSteady     = testSteady,
                      walltime       = as.integer(walltime),
                      simplify       = simplify_arg,
                      solveQuadratic = as.logical(solveQuadratic))
  }

  if (is.null(m_ss) || identical(m_ss, 0L)) return(0)

  # All versions return a list of "lhs=rhs" strings.
  # Parse into named character vector (dMod format).
  m_ssChar <- do.call(c, lapply(strsplit(m_ss, "="), function(eq) {
    out <- trimws(eq[2])
    names(out) <- trimws(eq[1])
    return(out)
  }))

  if (length(m_ssChar) == 0) return(0)

  # Write steady states to disk
  if (!is.null(file) && is.character(file))
    saveRDS(object = m_ssChar, file = file)

  return(m_ssChar)
}
