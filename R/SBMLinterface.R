#' Import an SMBL model
#'
#' Requires AMICI https://github.com/ICB-DCM/AMICI/ and dependencies to be installed on your system
#' Big thanks go to Daniel Weindl!
#'
#' Kinetic laws coming out of SBML are **extensive** (amount/time) by SBML
#' convention, whereas dMod stores rates in **concentration-style**. On import,
#' each kinetic law `K` is divided by the volume of its home compartment (the
#' compartment shared by the educts, or the product compartment for pure
#' synthesis) to produce `rate_dMod = K / V`. Combined with the volume-ratio
#' factors emitted by [getFluxes()], this preserves the SBML semantics:
#' `K_SBML = rate_dMod * V`.
#'
#' @param modelpath Path to the sbml file
#' @param amicipath Path to your amici-python-installation, e.g.: AMICIPATH/python
#'
#' @return list of eqnlist, parameters and inits
#' @export
#' @importFrom rjson fromJSON
#' @importFrom stringr str_replace_all
import_sbml <- function(modelpath, amicipath = NULL) {

  importscript <- system.file("code/sbmlAmicidMod.py", package = "dMod")
  tmpfile_json <- tempfile()
  modelpath <- normalizePath(modelpath, mustWork = TRUE)

  # Call the virtualenv's python directly — that python knows its own
  # site-packages via pyvenv.cfg, no `source activate` gymnastics needed.
  # `amicipath` (when set) is prepended to PYTHONPATH so users can mix in a
  # libsbml install that lives outside the venv.
  venv_python <- path.expand("~/.virtualenvs/amici/bin/python")
  if (!file.exists(venv_python))
    stop("dMod expects a Python virtualenv at ~/.virtualenvs/amici/. ",
         "Create one with `python3 -m venv ~/.virtualenvs/amici && ",
         "~/.virtualenvs/amici/bin/pip install python-libsbml`.")

  pyenv <- if (!is.null(amicipath))
             paste0("PYTHONPATH=",
                    paste(c(amicipath,
                            Sys.getenv("PYTHONPATH", unset = "")),
                          collapse = ":"))
           else character(0)

  status <- system2(venv_python,
                    args = c(shQuote(importscript),
                             shQuote(modelpath),
                             shQuote(tmpfile_json)),
                    env = pyenv)
  if (status != 0L || !file.exists(tmpfile_json))
    stop("SBML import failed (exit ", status, "). ",
         "Check that ~/.virtualenvs/amici/ has python-libsbml installed.")
  json_content <- rjson::fromJSON(file = tmpfile_json)

  S <- do.call(cbind, json_content[["S"]])
  S[S==0] <- NA

  # libsbml L3 emits natural log as `ln(x)`, but both R's stats::D() and the
  # C math library expect `log(x)` for natural log. Apply this normalisation
  # to every formula channel — rates, initial-value expressions, and
  # AssignmentRule RHSs — at a single point. The leading boundary
  # `(^|[^A-Za-z0-9_.])` prevents matching `eln`, `arcln`, etc.
  .normalise_formula <- function(s) {
    if (length(s) == 0L) return(s)
    s <- stringr::str_replace_all(s, "\\*\\*", "^")
    s <- stringr::str_replace_all(s, "(^|[^A-Za-z0-9_.])ln\\(", "\\1log(")
    s
  }

  v <- .normalise_formula(json_content[["v"]])

  states <- json_content[["stateNames"]]

  # Build compartment records + state-to-compartment map from the libsbml payload.
  compartments <- list()
  compartmentOf <- character(0)
  comp_json <- json_content[["compartments"]]
  spc_json  <- json_content[["speciesCompartments"]]
  if (!is.null(comp_json) && length(comp_json) > 0L) {
    for (c in comp_json) {
      # Use the SBML compartment ID as the volume expression: it appears as a
      # parameter symbol, so numeric assignment happens through the standard
      # parameter transformation. Leaves cancellation opportunities with any
      # `compartment * ...` factors still embedded in the kinetic law.
      compartments[[c$id]] <- list(volume = c$id, rule = NULL)
    }
    if (!is.null(spc_json) && length(spc_json) > 0L) {
      compartmentOf <- unlist(spc_json)
      compartmentOf <- compartmentOf[intersect(names(compartmentOf), states)]
    }
  }
  if (length(compartmentOf) == 0L) {
    compartments <- NULL
    compartmentOf <- NULL
  }

  # Normalize each kinetic law: rate_dMod = K / V_home. The home compartment is
  # determined from the educt rows of S; pure-synthesis reactions fall back to
  # the product compartment. Non-factorable kinetic laws retain a `/V` term in
  # the rate string, which is mathematically correct if aesthetically ugly.
  if (!is.null(compartmentOf) && !is.null(S)) {
    for (i in seq_along(v)) {
      row_i <- S[i, ]
      educt_idx <- which(!is.na(row_i) & row_i < 0)
      product_idx <- which(!is.na(row_i) & row_i > 0)
      home_cids <- if (length(educt_idx) > 0)   unique(compartmentOf[states[educt_idx]])
                   else if (length(product_idx) > 0) unique(compartmentOf[states[product_idx]])
                   else character(0)
      if (length(home_cids) == 1L) {
        v[i] <- paste0("(", v[i], ")/(", compartments[[home_cids]]$volume, ")")
      } else if (length(home_cids) > 1L) {
        warning(sprintf("Reaction %d spans compartments (%s); kinetic law stored as-is.",
                        i, paste(home_cids, collapse = ", ")))
      }
    }
  }

  pars <- setNames(json_content[["p"]], json_content[["parameterNames"]])
  x0 <- setNames(.normalise_formula(json_content[["x0"]]),
                 json_content[["stateNames"]])

  # Inline AssignmentRules from the SBML model. Each rule `lhs := rhs` becomes
  # an algebraic substitution applied to all rates and species initials. PEtab
  # benchmark models (e.g. Boehm_JProteomeRes2014) use rules to encode
  # time-varying inputs like `BaF3_Epo := 1.25e-7 * exp(-k * time)`.
  # Iterate to a fixed point so chained rules resolve. After inlining, the
  # LHS symbols are no longer free parameters and are dropped from `pars`.
  rules <- json_content[["assignmentRules"]]
  if (length(rules)) {
    rule_lhs <- names(rules)
    rule_rhs <- .normalise_formula(unlist(rules, use.names = FALSE))
    # Wrap each RHS in parens so substitution into a sub-expression keeps
    # operator precedence intact (e.g. `1.25e-7 * exp(...)` inside `a * lhs`).
    rule_rhs <- paste0("(", rule_rhs, ")")
    max_iter <- length(rule_lhs) + 1L
    for (it in seq_len(max_iter)) {
      new_v  <- replaceSymbols(rule_lhs, rule_rhs, v)
      new_x0 <- replaceSymbols(rule_lhs, rule_rhs, x0)
      if (identical(new_v, v) && identical(new_x0, x0)) break
      v <- new_v; x0 <- new_x0
    }
    pars <- pars[setdiff(names(pars), rule_lhs)]
  }

  reactions <- eqnlist(smatrix = S, states = states, rates = v,
                       compartments = compartments, compartmentOf = compartmentOf)

  observables <- json_content[["observables"]]

  out <- list(reactions = reactions, pars = pars, inits = x0,
              observables = observables, assignmentRules = rules)
  return(out)
}


#' Export an eqnlist to an SBML Level 3 file
#'
#' Serialises an [eqnlist] plus parameter values and initial concentrations to
#' SBML Level 3 Version 2. Each reaction's kinetic law is emitted as
#' `V_home * rate_dMod` to restore SBML's extensive-flux convention
#' (`K_SBML = rate_dMod * V`). Requires a Python environment with `libsbml`
#' installed; the default location matches the one used by [import_sbml()].
#'
#' @param eqnlist Object of class [eqnlist] to export.
#' @param parameters Named numeric vector of parameter values (including any
#'   compartment-size parameters referenced in `eqnlist$compartments`). Pass
#'   `NULL` to write parameters without values.
#' @param inits Named numeric *or* character vector of initial values keyed
#'   by state name. Numeric entries (or character entries that parse as
#'   numeric) are written as `initialConcentration`; non-numeric character
#'   entries are emitted as `<initialAssignment>` formulas and let the SBML
#'   simulator resolve the expression against `parameters` at sim time.
#'   Missing states default to 0.
#' @param filepath Path to the SBML output file.
#' @param model_id SBML model identifier. Defaults to `"dMod_export"`.
#' @param amicipath Optional `PYTHONPATH` entry prepended to the python call,
#'   matching the knob on [import_sbml()].
#' @return `filepath`, invisibly.
#' @export
#' @importFrom rjson toJSON
export_sbml <- function(eqnlist, parameters = NULL, inits = NULL, filepath,
                         model_id = "dMod_export", amicipath = NULL) {

  stopifnot(is.eqnlist(eqnlist))
  if (is.null(eqnlist$compartments) || is.null(eqnlist$compartmentOf))
    stop("`eqnlist` must have populated compartments/compartmentOf. Use the updated constructor.")

  # Compartments: numeric size when the volume expression parses as numeric;
  # otherwise 1.0 and rely on the homonymous parameter for the symbolic value.
  comp_list <- lapply(names(eqnlist$compartments), function(cid) {
    entry <- eqnlist$compartments[[cid]]
    vol <- entry$volume
    size <- suppressWarnings(as.numeric(vol))
    if (is.na(size)) size <- 1.0
    list(id = cid, size = size, spatialDimensions = 3L)
  })

  species_list <- lapply(eqnlist$states, function(st) {
    raw <- if (!is.null(inits) && st %in% names(inits)) inits[[st]] else 0
    num <- suppressWarnings(as.numeric(raw))
    base <- list(id = st, compartment = unname(eqnlist$compartmentOf[[st]]))
    if (!is.na(num)) {
      c(base, list(initialConcentration = num))
    } else {
      # symbolic: emit as <initialAssignment>; the formula may reference any
      # parameter declared on the SBML side (incl. compartment-volume IDs).
      c(base, list(initialAssignment = as.character(raw)))
    }
  })

  param_list <- list()
  if (!is.null(parameters)) {
    for (nm in names(parameters))
      param_list[[length(param_list) + 1L]] <- list(id = nm, value = as.numeric(parameters[[nm]]))
  }

  # Each row of the stoichiometric matrix becomes one reaction. The kinetic
  # law is `rate * V_home` — the α-bridge identity in the export direction.
  smatrix <- eqnlist$smatrix
  rxn_list <- lapply(seq_len(nrow(smatrix)), function(i) {
    row_i <- smatrix[i, ]
    # `which()` on a named vector preserves names, which would propagate
    # through lapply() into a *named* list — rjson then serialises it as
    # a JSON object, breaking the array-of-dicts contract dmodToSbml.py
    # expects. unname() the indices.
    educt_idx <- unname(which(!is.na(row_i) & row_i < 0))
    product_idx <- unname(which(!is.na(row_i) & row_i > 0))

    educts <- lapply(educt_idx, function(j)
      list(species = eqnlist$states[j], stoich = as.numeric(abs(row_i[j]))))
    products <- lapply(product_idx, function(j)
      list(species = eqnlist$states[j], stoich = as.numeric(row_i[j])))

    home_idx <- if (length(educt_idx) > 0L) educt_idx else product_idx
    home_cid <- unique(unname(eqnlist$compartmentOf[eqnlist$states[home_idx]]))[1]
    home_vol <- eqnlist$compartments[[home_cid]]$volume
    kinetic_law <- paste0("(", home_vol, ") * (", eqnlist$rates[i], ")")

    list(id = paste0("r", i),
         reactants = educts,
         products = products,
         kineticLaw = kinetic_law)
  })

  spec <- list(modelId = model_id,
               compartments = comp_list,
               species = species_list,
               parameters = param_list,
               reactions = rxn_list,
               outfile = normalizePath(filepath, mustWork = FALSE))

  spec_json <- tempfile(fileext = ".json")
  writeLines(rjson::toJSON(spec), spec_json)

  script <- system.file("code/dmodToSbml.py", package = "dMod")
  venv_python <- path.expand("~/.virtualenvs/amici/bin/python")
  if (!file.exists(venv_python))
    stop("dMod expects a Python virtualenv at ~/.virtualenvs/amici/.")
  pyenv <- if (!is.null(amicipath))
             paste0("PYTHONPATH=",
                    paste(c(amicipath,
                            Sys.getenv("PYTHONPATH", unset = "")),
                          collapse = ":"))
           else character(0)
  status <- system2(venv_python, args = c(shQuote(script), shQuote(spec_json)),
                    env = pyenv)
  if (status != 0L) stop("SBML export failed (exit ", status, ").")

  invisible(filepath)
}
