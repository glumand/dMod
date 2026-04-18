#' Calculate analytical steady states.
#'
#' @description This function follows the method published in \[1\]. Find the
#'   latest version of the tool and some examples under \[2\]. The determined
#'   steady-state solution is tailored to parameter estimation. Please note that
#'   kinetic parameters might be fixed for solution of steady-state equations.
#'   Note that additional parameters might be introduced to ensure positivity of
#'   the solution.
#' @description The function calls a python script via the reticulate package.
#'   Use python3.x
#' @description Three versions of the AlyssaPetit backend are available:
#'   \itemize{
#'     \item \code{"1.0"}: Original version. No steady-state test toggle.
#'     \item \code{"1.1"}: Adds the \code{testSteady} option to skip
#'       verification.
#'     \item \code{"1.2"} (default): Like v1.1 with structural sink-cluster
#'       detection (states without a protected mass source are set to 0 a
#'       priori), optional \code{walltime}, direct linear solves instead of
#'       \code{sympy.solve()} in the remaining-equations phase, and an
#'       end-of-pipeline \code{simplify} toggle. Cycle breaking uses a
#'       global priority table (idea and sort criteria ported from Severin
#'       Bang's Julia reimplementation) which replaces the previous
#'       first-cycle \code{GetBestPair} heuristic AND the numpy
#'       sparsification step — the table's `type` classification
#'       implicitly penalises cycle-creating resolutions, so an a-priori
#'       compacted stoichiometry is no longer needed.
#'   }
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
#' @param sparsifyLevel numeric, Upper bound for length of linear combinations
#'   used for simplifying the stoichiometric matrix. Only used for versions
#'   \code{"1.0"} and \code{"1.1"}. Ignored (no-op) in \code{"1.2"}, which
#'   relies on the priority-table cycle-breaking heuristic instead — kept for
#'   signature compatibility. Default changed to \code{0} for v1.2.
#' @param outputFormat Define the output format. By default "R" generating dMod
#'   compatible output. To obtain an output appropriate for d2d \[3\] "M" must be
#'   selected.
#' @param testSteady Character, if "T" the correctness of the obtained steady
#'   states is numerically checked (this can be very time intensive). If "F" this
#'   is skipped. Ignored for version "1.0" (always tests).
#' @param walltime integer, total wall-clock time budget in seconds for the
#'   solver (default 0 = unlimited). Used for version "1.2". When the budget
#'   is reached during the final simplification pass, any remaining
#'   expressions are left in their cancelled-but-not-simplified form (still
#'   mathematically correct, just bulkier).
#' @param simplify Controls the final simplification pass. One of:
#'   \itemize{
#'     \item \code{TRUE} (default): each final steady-state expression is
#'       passed through \code{sympy.simplify} once, at the end of the
#'       pipeline.
#'     \item \code{FALSE}: skip the final simplification entirely — fastest,
#'       but output may be bulkier (expressions stay in \code{cancel()}
#'       normal form).
#'     \item \code{"full"}: aggressive pipeline per expression —
#'       \code{cancel} then \code{posify} then
#'       \code{simplify(ratio = oo, rational = TRUE)} then separate
#'       \code{factor} of numerator and denominator. Noticeably slower,
#'       but often yields more compact formulas (especially when
#'       positivity of rate constants enables further cancellation).
#'   }
#'   Simplification is never applied inside the solve loop (which was the
#'   main bottleneck in previous releases). Used for version "1.2" only.
#'   For very large models prefer \code{simplify = FALSE} up-front, or
#'   rely on \code{walltime} to abort the simplify pass gracefully once
#'   enough expressions have been processed.
#' @param version Character, which AlyssaPetit backend to use. One of
#'   \code{"1.0"}, \code{"1.1"}, or \code{"1.2"} (default).
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
                         walltime = 0L, simplify = TRUE, version = "1.2") {

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
  #   v1.2: Alyssa(filename, injections, givenCQs, neglect, sparsifyLevel, outputFormat, testSteady, walltime, simplify)
  #        — v1.2 additionally runs structural sink-cluster detection a priori.
  if (version == "1.0") {
    if (testSteady == "F")
      message("Note: version 1.0 does not support testSteady='F', test will always run.")
    m_ss <- ap$Alyssa(model, as.list(forcings), as.list(givenCQs),
                      as.list(neglect), sparsifyLevel, outputFormat)

  } else if (version == "1.1") {
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
                      injections    = as.list(forcings),
                      givenCQs      = as.list(givenCQs),
                      neglect       = as.list(neglect),
                      sparsifyLevel = as.integer(sparsifyLevel),
                      outputFormat  = outputFormat,
                      testSteady    = testSteady,
                      walltime      = as.integer(walltime),
                      simplify      = simplify_arg)
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
