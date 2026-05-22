#' Import an SBML model
#'
#' Reads an SBML Level 3 file via a Python helper (`inst/code/sbmlImport.py`)
#' that uses `python-libsbml`. The Python environment is provisioned
#' automatically by `reticulate` on first use — no manual venv setup is
#' required. Users who want to point dMod at an existing interpreter
#' (e.g. a hand-managed venv or conda env) can set the
#' `DMOD_LIBSBML_PYTHON` environment variable to its absolute path; that
#' bypasses reticulate entirely.
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
#'
#' @return list of eqnlist, parameters and inits
#' @export
#' @importFrom rjson fromJSON
#' @importFrom stringr str_replace_all
import_sbml <- function(modelpath) {

  importscript <- system.file("code/sbmlImport.py", package = "dMod")
  tmpfile_json <- tempfile()
  modelpath <- normalizePath(modelpath, mustWork = TRUE)

  venv_python <- .dmod_libsbml_python()

  status <- system2(venv_python,
                    args = c(shQuote(importscript),
                             shQuote(modelpath),
                             shQuote(tmpfile_json)))
  if (status != 0L || !file.exists(tmpfile_json))
    stop("SBML import failed (exit ", status, "). ",
         "Check that the libsbml virtualenv has python-libsbml installed.")
  json_content <- rjson::fromJSON(file = tmpfile_json)

  # `S` from the python side is a list-of-lists with shape
  # `[n_species][n_reactions]`. rjson::fromJSON collapses fully-nested
  # arrays whose inner length is 1 into a bare vector, which breaks
  # `cbind`'s list-of-columns contract. We reshape against the known
  # state/reaction counts so the resulting matrix is always
  # `[n_reactions rows x n_states cols]`, matching `eqnlist`'s convention.
  S_raw <- json_content[["S"]]
  n_states <- length(json_content[["stateNames"]])
  n_rxns   <- length(json_content[["v"]])
  S <- if (is.list(S_raw)) do.call(cbind, S_raw)
       else matrix(unlist(S_raw), nrow = n_rxns, ncol = n_states)
  S[S == 0] <- NA

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
      # Compartments with size = 1 (and no rule) carry no symbolic content —
      # storing them as the literal "1" keeps the compartment ID out of the
      # kinetic laws, which is what dMod's roundtrip expects when the source
      # eqnlist had volume "1". Otherwise use the SBML compartment ID as the
      # volume symbol so the trafo can override it.
      trivial <- !is.null(c$size) && is.numeric(c$size) && isTRUE(c$size == 1)
      compartments[[c$id]] <- list(volume = if (trivial) "1" else c$id,
                                   rule   = NULL)
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
        home_vol <- compartments[[home_cids]]$volume
        if (!identical(home_vol, "1"))
          v[i] <- paste0("(", v[i], ")/(", home_vol, ")")
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

  # --- rate rules ---
  # `<rateRule variable="X">` defines dX/dt = rhs. Per SBML spec, X cannot
  # also be produced/consumed by reactions, so the new column is independent
  # of existing kinetic laws. dC/dt is intensive already, so we append the
  # rate AFTER the volume-division loop above (no /V wrap). RateRules on
  # non-species (parameter / compartment) are skipped with a warning.
  rate_rules <- json_content[["rateRules"]]
  if (length(rate_rules)) {
    rr_lhs <- names(rate_rules)
    rr_rhs <- .normalise_formula(unlist(rate_rules, use.names = FALSE))
    # Inline assignment-rule LHSs into the rate-rule RHSs too, so a RateRule
    # that references a rule-defined symbol doesn't carry it as a free var.
    if (length(rules)) {
      for (it in seq_len(max_iter)) {
        new_rr <- replaceSymbols(rule_lhs, rule_rhs, rr_rhs)
        if (identical(new_rr, rr_rhs)) break
        rr_rhs <- new_rr
      }
    }
    for (k in seq_along(rr_lhs)) {
      var <- rr_lhs[k]; rhs <- rr_rhs[k]
      if (!var %in% states) {
        warning(sprintf(
          "RateRule on `%s` ignored: only species RateRules are supported (got non-species; promote to a state in the SBML or rewrite as a reaction).",
          var), call. = FALSE)
        next
      }
      # Append a virtual reaction: stoichiometry +1 on `var`, 0 elsewhere;
      # rate string = rhs. After eqnlist construction this contributes
      # +rhs to the RHS row of `var`, which is exactly the rate rule.
      new_col <- rep(NA_real_, length(states))
      new_col[match(var, states)] <- 1
      S <- if (is.null(S)) matrix(new_col, nrow = 1L)
           else rbind(S, new_col)
      v <- c(v, rhs)
    }
  }

  reactions <- eqnlist(smatrix = S, states = states, rates = v,
                       compartments = compartments, compartmentOf = compartmentOf)

  observables <- json_content[["observables"]]

  # --- events ---
  # SBML <event> -> dMod eventlist (one row per <eventAssignment>). Triggers
  # of the form `time >=/== T` (numeric or symbolic T) populate `time`; other
  # triggers fall back to a root expression of the form `lhs - (rhs)` so the
  # ODE solver can detect the zero crossing. Method is always "replace"
  # (SBML eventAssignments are assignment-style by spec).
  events_json <- json_content[["events"]]
  events_df <- NULL
  if (length(events_json)) {
    rows <- list()
    cmp_pat <- "^\\s*(.+?)\\s*(>=|>|<=|<|==|!=)\\s*(.+?)\\s*$"
    for (ev in events_json) {
      tt <- ev[["triggerTime"]]
      tf <- ev[["triggerFormula"]]
      tt_num <- suppressWarnings(as.numeric(tt))
      time_val <- if (length(tt) && !is.null(tt) && !is.na(tt_num)) tt_num
                  else if (length(tt) && nzchar(tt)) as.character(tt)
                  else NA
      root_val <- NA_character_
      if (is.na(time_val) && length(tf) && nzchar(tf)) {
        m <- regmatches(tf, regexec(cmp_pat, tf))[[1L]]
        root_val <- if (length(m) == 4L)
                      paste0("(", m[2L], ") - (", m[4L], ")")
                    else
                      tf
      }
      for (a in ev$assignments) {
        rows[[length(rows) + 1L]] <- data.frame(
          var    = a$variable,
          time   = time_val,
          value  = .normalise_formula(a$formula),
          root   = root_val,
          method = "replace",
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(rows)) {
      events_df <- as.eventlist(do.call(rbind, rows))
    }
  }

  out <- list(reactions = reactions, pars = pars, inits = x0,
              observables = observables, assignmentRules = rules,
              events = events_df)
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
#' @param modelID SBML model identifier. Defaults to `"dMod_export"`.
#' @return `filepath`, invisibly.
#' @export
#' @importFrom rjson toJSON
export_sbml <- function(eqnlist, parameters = NULL, inits = NULL, filepath,
                         modelID = "dMod_export") {

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

  spec <- list(modelId = modelID,
               compartments = comp_list,
               species = species_list,
               parameters = param_list,
               reactions = rxn_list,
               outfile = normalizePath(filepath, mustWork = FALSE))

  spec_json <- tempfile(fileext = ".json")
  writeLines(rjson::toJSON(spec), spec_json)

  script <- system.file("code/dmodToSbml.py", package = "dMod")
  venv_python <- .dmod_libsbml_python()
  status <- system2(venv_python, args = c(shQuote(script), shQuote(spec_json)))
  if (status != 0L) stop("SBML export failed (exit ", status, ").")

  invisible(filepath)
}


.dmod_libsbml_python <- function() {
  # Explicit override wins (existing user envs, conda, CI with prebuilt
  # interpreters, ...). Skip reticulate provisioning entirely.
  override <- Sys.getenv("DMOD_LIBSBML_PYTHON", unset = "")
  python <- if (nzchar(override)) {
    if (!file.exists(override))
      stop("DMOD_LIBSBML_PYTHON=", override, " does not exist.")
    override
  } else {
    # `python-libsbml` was declared via reticulate::py_require() in
    # .onLoad(). py_exe() materialises the managed env (downloads Python +
    # installs the requirement on first call) and returns the interpreter
    # path used by the system2() calls below.
    tryCatch(reticulate::py_exe(), error = function(e) {
      stop("Could not provision a Python with python-libsbml via ",
           "reticulate (", conditionMessage(e), "). Set ",
           "DMOD_LIBSBML_PYTHON to point at a Python interpreter that ",
           "has python-libsbml installed.")
    })
  }

  # Probe `import libsbml` once per session. For the override path this
  # catches a wrong interpreter early; for the reticulate path it is a
  # cheap sanity check that the requirement actually resolved. Cached via
  # an env var so the ~30 ms python spawn does not repeat across
  # import_sbml / export_sbml calls.
  if (!identical(Sys.getenv("DMOD_LIBSBML_OK", unset = ""), "1")) {
    status <- suppressWarnings(
      system2(python, args = c("-c", shQuote("import libsbml")),
              stdout = FALSE, stderr = FALSE))
    if (status != 0L)
      stop("Python at ", python, " could not `import libsbml` ",
           "(status ", status, "). ",
           if (nzchar(override))
             "Install python-libsbml into that env."
           else
             "reticulate did not provision python-libsbml as expected.")
    Sys.setenv(DMOD_LIBSBML_OK = "1")
  }
  python
}
