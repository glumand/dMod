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
#' @description Since v2.0 the tool handles cycles that cannot be removed
#'   analytically by computing a minimal Feedback Vertex Set (FVS), attempting
#'   polynomial solving for any symbol in the FVS equations (preferring the
#'   lowest polynomial degree), and providing symbolic global stability
#'   conditions (1D/2D) when closed-form solutions are not available.
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
#'   used for simplifying the stoichiometric matrix
#' @param maxPolyDegree integer, Maximum polynomial degree for which closed-form
#'   solutions of FVS species are attempted. Default 2. Degrees 3-4 are
#'   supported but computationally expensive.
#' @param outputFormat Define the output format. By default "R" generating dMod
#'   compatible output. To obtain an output appropriate for d2d \[3\] "M" must be
#'   selected.
#' @param testSteady Boolean, if "T" the correctness of the obtained steady
#'   states is numerically checked (this can be very time intensive). If "F" this
#'   is skipped.
#'
#' @return List with components:
#'   \describe{
#'     \item{equations}{Named character vector of steady-state equations (dMod compatible).}
#'     \item{fvs_species}{Character vector of FVS species names (empty if no FVS needed).}
#'     \item{fvs_results}{Character vector of FVS solution details.}
#'     \item{fvs_unsolved}{Character vector of FVS species that could not be solved in closed form.}
#'     \item{convergence_conditions}{Character vector of symbolic stability conditions (only for unsolved FVS with n <= 2).}
#'   }
#'   For backward compatibility, if no FVS species are present, returns only the
#'   named character vector of equations (as in v1.x).
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
                         givenCQs = NULL, neglect = NULL, sparsifyLevel = 2,
                         maxPolyDegree = 2L, outputFormat = "R",
                         testSteady = "T") {
  
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
  
  # Calculate steady states via AlyssaPetit v2.0
  ap <- reticulate::import_from_path("AlyssaPetit_ver2_0",
                                     path = system.file("code", package = "dMod"))
  m_ss <- ap$Alyssa(model, as.list(forcings), as.list(givenCQs), as.list(neglect),
                    sparsifyLevel, as.integer(maxPolyDegree), outputFormat)
  
  if (is.null(m_ss) || identical(m_ss, 0L)) return(0)
  
  # v2.0 returns a dict; extract components
  equations     <- m_ss$equations
  fvs_species   <- m_ss$fvs_species
  fvs_results   <- m_ss$fvs_results
  fvs_unsolved  <- m_ss$fvs_unsolved
  conv_conds    <- m_ss$convergence_conditions
  
  # Parse equations into named character vector (dMod format)
  if (length(equations) > 0) {
    m_ssChar <- do.call(c, lapply(strsplit(equations, "="), function(eq) {
      out <- trimws(eq[2])
      names(out) <- trimws(eq[1])
      return(out)
    }))
  } else {
    return(0)
  }
  
  # Write steady states to disk
  if (!is.null(file) && is.character(file))
    saveRDS(object = m_ssChar, file = file)
  
  # return
  if (length(fvs_species) == 0) {
    return(m_ssChar)
  } else {
    return(list(
      equations              = m_ssChar,
      fvs_species            = fvs_species,
      fvs_results            = fvs_results,
      fvs_unsolved           = fvs_unsolved,
      convergence_conditions = conv_conds
    ))
  }
}

