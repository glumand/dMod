## R/PEtabInterface.R
## ---------------------------------------------------------------------------
## PEtab v1 importer / exporter, layered on top of dMod's existing
## SBML interface (R/SBMLinterface.R) and high-level APIs (Y, P, normL2, ...).
##
## Public API: importPEtab, exportPEtab, read_petab_yaml, read_petab_tables.
## Internal helpers prefixed `.petab_*` are unexported and may change shape.
##
## See https://petab.readthedocs.io/en/latest/v1/documentation_data_format.html
## ---------------------------------------------------------------------------


## --- low-level YAML / TSV readers ------------------------------------------

#' Read a PEtab YAML manifest
#'
#' Parses a PEtab v1 YAML file and resolves the paths of the referenced
#' tables. No SBML or TSV is touched at this stage.
#'
#' @param yaml_path Path to the PEtab YAML manifest.
#' @return A list with `base_dir`, `format_version`, `parameter_file`, and a
#'   one-element `problems` list whose entry holds resolved (absolute) paths
#'   for `sbml_file`, `condition_file`, `measurement_file`, `observable_file`.
#'   PEtab allows multiple problems and multiple files per slot; the current
#'   reader supports a single problem with one file per slot and errors
#'   otherwise.
#' @export
#' @importFrom yaml read_yaml
read_petab_yaml <- function(yaml_path) {

  yaml_path <- normalizePath(yaml_path, mustWork = TRUE)
  base_dir  <- dirname(yaml_path)
  m <- yaml::read_yaml(yaml_path)

  fv <- m$format_version %||% 1L
  if (!fv %in% c(1L, 1, "1"))
    stop("PEtab format_version ", fv, " not supported (only v1).")

  if (length(m$problems) != 1L)
    stop("Only single-problem PEtab YAML is supported (got ",
         length(m$problems), ").")

  prob <- m$problems[[1]]
  pick_one <- function(x, slot) {
    if (length(x) == 0L) stop("YAML problem missing `", slot, "`.")
    if (length(x) > 1L) stop("Only one ", slot, " per problem is supported.")
    x[[1]]
  }
  resolve <- function(p) normalizePath(file.path(base_dir, p), mustWork = TRUE)

  list(
    base_dir       = base_dir,
    format_version = 1L,
    parameter_file = resolve(m$parameter_file),
    problems = list(list(
      sbml_file        = resolve(pick_one(prob$sbml_files, "sbml_files")),
      condition_file   = resolve(pick_one(prob$condition_files, "condition_files")),
      measurement_file = resolve(pick_one(prob$measurement_files, "measurement_files")),
      observable_file  = resolve(pick_one(prob$observable_files, "observable_files"))
    ))
  )
}


#' Read PEtab TSV tables (no SBML)
#'
#' Reads the four PEtab tables referenced by a YAML manifest into base R
#' data frames. Useful for inspecting a problem without invoking AMICI.
#'
#' @param yaml_path Path to the PEtab YAML manifest.
#' @return A list with `parameters`, `conditions`, `measurements`,
#'   `observables` (data frames) and `sbml_path` (character).
#' @export
read_petab_tables <- function(yaml_path) {

  m  <- read_petab_yaml(yaml_path)
  pr <- m$problems[[1]]

  list(
    parameters   = .petab_read_tsv(m$parameter_file),
    conditions   = .petab_read_tsv(pr$condition_file),
    measurements = .petab_read_tsv(pr$measurement_file),
    observables  = .petab_read_tsv(pr$observable_file),
    sbml_path    = pr$sbml_file
  )
}


.petab_read_tsv <- function(path) {
  utils::read.delim(path, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE,
                    check.names = FALSE,
                    na.strings = c("", "NA"),
                    strip.white = TRUE)
}


## --- per-table parsers (unit-testable, no SBML side effects) ---------------

# Internal: parameters.tsv -> dMod-shaped pieces.
# Returns:
#   pouter         named numeric, estimated parameters (nominalValue, scale-applied)
#   lower / upper  named numeric, bounds (scale-applied)
#   fixed          named numeric, non-estimated parameters (scale-applied)
#   scales         named character, "lin"/"log"/"log10" per parameterId
#                  (covers both pouter and fixed)
.petab_parse_parameters <- function(df) {

  required <- c("parameterId", "parameterScale", "lowerBound", "upperBound",
                "nominalValue", "estimate")
  miss <- setdiff(required, colnames(df))
  if (length(miss))
    stop("parameters.tsv missing required column(s): ",
         paste(miss, collapse = ", "))

  scales <- setNames(df$parameterScale, df$parameterId)
  if (any(!scales %in% c("lin", "log", "log10")))
    stop("Unknown parameterScale(s): ",
         paste(unique(scales[!scales %in% c("lin", "log", "log10")]),
               collapse = ", "))

  # PEtab semantics: nominalValue / lowerBound / upperBound are given on the
  # *parameter scale* the user picked. dMod's outer parameters live on that
  # same scale; the back-transform exp/10^ happens inside the trafo (see
  # .petab_build_trafo). So we keep the values as written and only translate
  # `lin` numerics to numeric — no log/log10 conversion here.
  to_num <- function(col) suppressWarnings(as.numeric(col))

  pouter_idx <- which(df$estimate == 1)
  fixed_idx  <- which(df$estimate == 0)

  pouter <- setNames(to_num(df$nominalValue[pouter_idx]),
                     df$parameterId[pouter_idx])
  fixed  <- setNames(to_num(df$nominalValue[fixed_idx]),
                     df$parameterId[fixed_idx])
  lower  <- setNames(to_num(df$lowerBound[pouter_idx]),
                     df$parameterId[pouter_idx])
  upper  <- setNames(to_num(df$upperBound[pouter_idx]),
                     df$parameterId[pouter_idx])

  list(pouter = pouter, lower = lower, upper = upper,
       fixed = fixed, scales = scales)
}


# Internal: observables.tsv -> dMod-shaped pieces.
# Returns:
#   obs       named character, observableFormula keyed by observableId.
#             *Not* yet substituted: observableParameterK_<id> / noiseParameterK_<id>
#             placeholders are left intact for per-row substitution downstream.
#   noise     named character or named numeric, noiseFormula per observableId.
#             Numeric strings are kept as numeric (constant noise) so the
#             objective dispatcher can pick the normL2 fast path.
#   obs_trafo / noise_dist  one of "lin"/"log"/"log10" / "normal"/"laplace"/"log-normal"
#                           per observableId; PEtab defaults are "lin"/"normal".
.petab_parse_observables <- function(df) {

  if (!"observableId" %in% colnames(df))
    stop("observables.tsv missing required column `observableId`.")
  if (!"observableFormula" %in% colnames(df))
    stop("observables.tsv missing required column `observableFormula`.")

  ids <- df$observableId
  if (anyDuplicated(ids))
    stop("Duplicate observableId(s) in observables.tsv: ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "))

  obs <- setNames(df$observableFormula, ids)

  noise_raw <- if ("noiseFormula" %in% colnames(df)) df$noiseFormula
               else rep("1", length(ids))
  noise_raw[is.na(noise_raw)] <- "1"
  names(noise_raw) <- ids

  trafo <- if ("observableTransformation" %in% colnames(df))
             df$observableTransformation else rep("lin", length(ids))
  trafo[is.na(trafo)] <- "lin"
  names(trafo) <- ids

  dist <- if ("noiseDistribution" %in% colnames(df))
            df$noiseDistribution else rep("normal", length(ids))
  dist[is.na(dist)] <- "normal"
  names(dist) <- ids

  if (any(!trafo %in% c("lin", "log", "log10")))
    stop("Unknown observableTransformation(s).")
  if (any(!dist %in% c("normal", "laplace", "log-normal")))
    stop("Unknown noiseDistribution(s).")

  list(obs = obs, noise = noise_raw, obs_trafo = trafo, noise_dist = dist)
}


# Internal: classify each non-id column of conditions.tsv as
#   "init"        — a species/state initial-value override (column name matches state)
#   "compartment" — a compartment volume override
#   "parameter"   — a parameter override (the catch-all, includes condition-only pars)
#
# `sbml_states`        names of species in the imported eqnlist
# `sbml_compartments`  names from reactions$compartments
# `sbml_pars`          names of SBML parameters from import_sbml()$pars
.petab_parse_conditions <- function(df, sbml_states = character(),
                                    sbml_compartments = character(),
                                    sbml_pars = character(),
                                    obs_inner = character()) {

  if (!"conditionId" %in% colnames(df))
    stop("conditions.tsv missing required column `conditionId`.")

  cond_ids <- df$conditionId
  if (anyDuplicated(cond_ids))
    stop("Duplicate conditionId(s) in conditions.tsv: ",
         paste(unique(cond_ids[duplicated(cond_ids)]), collapse = ", "))

  override_cols <- setdiff(colnames(df), c("conditionId", "conditionName"))

  col_kind <- vapply(override_cols, function(cn) {
    if (cn %in% sbml_states)              "init"
    else if (cn %in% sbml_compartments)   "compartment"
    else if (cn %in% sbml_pars)           "parameter"
    else if (cn %in% obs_inner)           "parameter"  # observable inner par
    else { warning("Condition column `", cn,
                   "` does not match any SBML or observable symbol; treating as parameter.")
           "parameter" }
  }, character(1))

  grid <- df
  rownames(grid) <- as.character(cond_ids)
  list(grid = grid, col_kind = col_kind, override_cols = override_cols)
}


# Internal: build a stable, short, name-safe key for sub-conditions whose
# observable/noise parameter strings differ within a single PEtab condition.
# Raw PEtab values can contain ";" "." etc., which dMod's modelname sanitiser
# would mangle.
.petab_subcond_hash <- function(s) {
  if (length(s) == 0L || all(is.na(s) | s == "")) return(rep("", length(s)))
  vapply(s, function(x) {
    if (is.na(x) || identical(x, "")) ""
    else substr(digest::digest(x, algo = "md5"), 1L, 8L)
  }, character(1))
}


# Internal: parse measurements.tsv into a long data.frame, splitting
# (simulationConditionId, observableParameters, noiseParameters) tuples that
# vary within a single sim condition into sub-conditions.
#
# Returns:
#   data        as.datalist-ready data.frame with columns
#                 name, time, value, sigma, condition
#               sigma defaults to NA where noiseFormula is symbolic, or to
#               the numeric noise constant from observables.tsv otherwise.
#   sub_cond_map  data.frame of sub_condition assignment, one row per
#               (orig_sim, sub_cond) pair, with the obs/noise param string.
#   peq_map     character vector cond_id -> preeq_cond_id (or "" if none).
.petab_parse_measurements <- function(df, obs_meta) {

  needed <- c("observableId", "simulationConditionId", "time", "measurement")
  miss <- setdiff(needed, colnames(df))
  if (length(miss))
    stop("measurements.tsv missing required column(s): ",
         paste(miss, collapse = ", "))

  m <- df
  m$preequilibrationConditionId <-
    if ("preequilibrationConditionId" %in% colnames(m))
      ifelse(is.na(m$preequilibrationConditionId), "",
             m$preequilibrationConditionId)
    else rep("", nrow(m))

  # read.delim infers numeric type when a column is all numeric — but PEtab
  # observable/noise parameters can be either numeric or symbol strings, and
  # downstream substitution code requires character. Cast explicitly.
  m$observableParameters <-
    if ("observableParameters" %in% colnames(m))
      ifelse(is.na(m$observableParameters), "",
             as.character(m$observableParameters))
    else rep("", nrow(m))

  m$noiseParameters <-
    if ("noiseParameters" %in% colnames(m))
      ifelse(is.na(m$noiseParameters), "",
             as.character(m$noiseParameters))
    else rep("", nrow(m))

  # Sub-condition key: (simCondId, peq, obsParStr, noiseParStr).
  # Empty strings collapse to a single key per sim cond. Hashes only kick in
  # when at least one of obs/noise param strings is non-empty AND varies.
  obs_hash <- .petab_subcond_hash(m$observableParameters)
  noi_hash <- .petab_subcond_hash(m$noiseParameters)
  peq_hash <- .petab_subcond_hash(m$preequilibrationConditionId)

  sub_cond <- m$simulationConditionId
  needs_split <- (obs_hash != "") | (noi_hash != "") | (peq_hash != "")
  if (any(needs_split)) {
    suffix <- ifelse(!needs_split, "",
                     paste0("__",
                            ifelse(peq_hash == "", "", paste0(peq_hash, "_")),
                            obs_hash,
                            ifelse(noi_hash == "", "", paste0("_", noi_hash))))
    sub_cond <- paste0(m$simulationConditionId, suffix)
  }
  m$sub_condition <- sub_cond

  # Build the sub_cond_map (one row per unique sub_cond).
  sub_cond_map <- unique(m[, c("simulationConditionId",
                               "preequilibrationConditionId",
                               "observableParameters", "noiseParameters",
                               "sub_condition")])
  rownames(sub_cond_map) <- NULL
  if (anyDuplicated(sub_cond_map$sub_condition))
    stop("Internal: sub-condition map has duplicate sub_condition rows. ",
         "This indicates inconsistent observable/noise parameter strings.")

  # Sigma per row: take the observable's noiseFormula, apply per-row
  # observable / noise parameter substitutions, then try to evaluate as a
  # constant. If all symbols are eliminated and the result is finite, use it
  # as sigma directly (fast path through normL2's data sigma column). If
  # symbols remain (case 0015's "noise"), leave sigma = NA — the err model
  # built by .petab_build_error_fn handles the per-condition formula.
  sigma <- vapply(seq_len(nrow(m)), function(i) {
    obsId <- m$observableId[i]
    f <- obs_meta$noise[obsId]
    if (is.na(f) || !nzchar(f)) return(NA_real_)
    f <- .petab_substitute_param_string(f, m$observableParameters[i],
                                        prefix = "observableParameter")
    f <- .petab_substitute_param_string(f, m$noiseParameters[i],
                                        prefix = "noiseParameter")
    .petab_eval_constant(f)
  }, numeric(1))

  # Apply observable transformations to data values. PEtab convention: the
  # measurement column is on the linear scale even when
  # observableTransformation != "lin"; the transformation is applied inside
  # the likelihood. We achieve this by wrapping the observable formula on
  # the simulation side (see .petab_build_observation_fn) and transforming
  # the data on this side, so the residual is computed on the chosen scale
  # and dMod's normL2 fast path keeps working.
  trafo_per_row <- obs_meta$obs_trafo[m$observableId]
  val <- as.numeric(m$measurement)
  log_idx   <- which(trafo_per_row == "log")
  log10_idx <- which(trafo_per_row == "log10")
  if (length(log_idx))   val[log_idx]   <- log(val[log_idx])
  if (length(log10_idx)) val[log10_idx] <- log10(val[log10_idx])

  data_df <- data.frame(
    name      = m$observableId,
    time      = as.numeric(m$time),
    value     = val,
    sigma     = as.numeric(sigma),
    condition = m$sub_condition,
    stringsAsFactors = FALSE
  )

  list(data = data_df,
       sub_cond_map = sub_cond_map,
       has_preeq = any(sub_cond_map$preequilibrationConditionId != ""))
}


# Internal: substitute PEtab parameter placeholders (observableParameterK_<id>
# or noiseParameterK_<id>) inside a *string* `formula` with the K-th value of
# `repls_str` (";"-separated).  Used for both data-side sigma evaluation and
# error-model construction.
.petab_substitute_param_string <- function(formula, repls_str, prefix) {
  if (length(formula) != 1L) {
    return(vapply(formula, .petab_substitute_param_string, character(1),
                  repls_str = repls_str, prefix = prefix))
  }
  if (is.na(formula) || !nzchar(formula)) return(formula)
  if (is.na(repls_str) || !nzchar(repls_str)) return(formula)
  parts <- trimws(strsplit(repls_str, ";", fixed = TRUE)[[1]])
  for (k in seq_along(parts)) {
    pat <- sprintf("\\b%s%d_[A-Za-z][A-Za-z0-9_]*\\b", prefix, k)
    formula <- gsub(pat, parts[k], formula, perl = TRUE)
  }
  formula
}


# Internal: try to evaluate a formula string as a constant numeric. Supports
# pure literals ("0.5"), arithmetic ("0.5 + 2"), and basic transcendentals
# ("log(2)"). Returns NA_real_ if the result is non-numeric, non-finite, or
# if any free symbol remains (eval would error in baseenv()).
.petab_eval_constant <- function(formula) {
  if (is.na(formula) || !nzchar(formula)) return(NA_real_)
  num <- suppressWarnings(as.numeric(formula))
  if (!is.na(num)) return(num)
  v <- tryCatch(eval(parse(text = formula), envir = baseenv()),
                error = function(e) NULL)
  if (is.numeric(v) && length(v) == 1L && is.finite(v)) v else NA_real_
}


## --- core trafo / observation / objective builders --------------------------

# Build the per-condition parameter trafo (parfn).
#
# Inputs:
#   sub_cond_map   data.frame from .petab_parse_measurements
#   conditions     data.frame from .petab_parse_conditions$grid
#   col_kind       named character from .petab_parse_conditions
#   override_cols  character vector of override column names
#   inits          named character; symbolic species initial expressions
#                  from import_sbml()$inits
#   sbml_pars      named numeric; SBML parameter defaults from
#                  import_sbml()$pars
#   states         character vector of state names
#   inner_pars     character vector of inner ODE parameters (kinetic rates,
#                  compartment volumes, anything in the rates' getSymbols
#                  minus states/time)
#   pouter_names   character vector of estimated parameter names
#   fixed          named numeric of fixed parameters
#   scales         named character "lin"/"log"/"log10" per parameterId
#   obs_inner      character vector of inner observable parameters that
#                  appear in observable / noise formulas (for substitution)
#   reactions      eqnlist (used for steadyStates fallback if pre-eq present)
#
# Returns a parfn that maps outer pars (estimated + fixed) to the inner-side
# parameter set (states' initial values, inner_pars, obs_inner).
.petab_build_trafo <- function(sub_cond_map, conditions, col_kind, override_cols,
                               inits, sbml_pars, states, inner_pars,
                               pouter_names, fixed, scales,
                               obs_inner = character(),
                               reactions = NULL,
                               compile = TRUE,
                               modelname = "petab_trafo",
                               preeqMethod = c("analytic", "numeric")) {

  preeqMethod <- match.arg(preeqMethod)

  # Inner side that the trafo must produce per condition.
  inner_targets <- unique(c(states, inner_pars, obs_inner))

  # Default trafo: every inner target is mapped to itself, except states which
  # default to their SBML initial expression. After this baseline we apply
  # condition overrides and observable/noise parameter substitutions, then
  # the parameter-scale chain rule.
  build_default <- function() {
    base <- setNames(inner_targets, inner_targets)
    for (st in intersect(states, names(inits))) {
      v <- inits[[st]]
      base[st] <- as.character(v)
    }
    base
  }

  apply_scale_chain_rule <- function(tr) {
    log_pars   <- names(scales)[scales == "log"   & names(scales) %in% pouter_names]
    log10_pars <- names(scales)[scales == "log10" & names(scales) %in% pouter_names]
    if (length(log_pars))   tr <- repar("x ~ exp(x)", tr, x = log_pars)
    if (length(log10_pars)) tr <- repar("x ~ 10^(x)", tr, x = log10_pars)
    tr
  }

  apply_row_overrides <- function(tr, cond_id, scope) {
    # scope: "all" (apply every override column),
    #        "state-only" (only init-kind columns; parameter overrides handled
    #         elsewhere — for the Pequil pre/post stages where parameter
    #         overrides are baked into the rates).
    if (!cond_id %in% rownames(conditions)) return(tr)
    for (cn in override_cols) {
      v <- conditions[cond_id, cn]
      if (is.na(v) || is.null(v)) next
      v <- trimws(as.character(v))
      if (v == "") next
      kind <- col_kind[[cn]]
      if (kind == "init") {
        tr[cn] <- v
      } else if (identical(scope, "all")) {
        tr <- repar(paste0(cn, " ~ ", v), tr)
      }
    }
    tr
  }

  apply_petab_param_subs <- function(tr, obs_str, noi_str) {
    if (length(obs_str) && nzchar(obs_str))
      tr <- .petab_apply_param_substitution(tr, obs_str, prefix = "observableParameter")
    if (length(noi_str) && nzchar(noi_str))
      tr <- .petab_apply_param_substitution(tr, noi_str, prefix = "noiseParameter")
    tr
  }

  # Pre-compute analytic steady states once if requested. Reused across all
  # sub-conditions whose preeqId is non-empty.
  mysteadies <- NULL
  analytic_ok <- FALSE
  if (preeqMethod == "analytic" &&
      any(nzchar(sub_cond_map$preequilibrationConditionId))) {
    mysteadies <- tryCatch(
      suppressMessages(steadyStates(reactions)),
      error = function(e) {
        warning("steadyStates() failed: ", conditionMessage(e),
                ". Falling back to numeric pre-equilibration via Pequil().",
                call. = FALSE)
        NULL
      })
    if (!is.null(mysteadies) && is.character(mysteadies)) {
      resolved <- names(mysteadies)[unname(mysteadies) != names(mysteadies)]
      unresolved <- setdiff(states, resolved)
      analytic_ok <- length(unresolved) == 0L
      if (!analytic_ok)
        warning("steadyStates() left state(s) ",
                paste(unresolved, collapse = ", "),
                " unresolved (RHS = LHS). Falling back to numeric ",
                "pre-equilibration via Pequil() for sub-conditions with a ",
                "preequilibrationConditionId.", call. = FALSE)
    }
  }

  parfns <- list()

  for (sub in sub_cond_map$sub_condition) {

    row_idx <- which(sub_cond_map$sub_condition == sub)
    sim     <- sub_cond_map$simulationConditionId[row_idx]
    peq     <- sub_cond_map$preequilibrationConditionId[row_idx]
    obs_str <- sub_cond_map$observableParameters[row_idx]
    noi_str <- sub_cond_map$noiseParameters[row_idx]

    has_peq <- length(peq) && nzchar(peq)

    if (!has_peq) {
      # ----- single-stage Pexpl path (Stage-1 path) ---------------------
      tr <- build_default()
      tr <- apply_row_overrides(tr, sim, "all")
      tr <- apply_petab_param_subs(tr, obs_str, noi_str)
      tr <- apply_scale_chain_rule(tr)
      parfns[[sub]] <- P(tr, condition = sub, compile = compile,
                         modelname = paste(modelname,
                                           sanitizeConditions(sub),
                                           sep = "_"))
      next
    }

    # ----- preeq path -------------------------------------------------------
    if (analytic_ok) {
      # Inject steadyStates expressions as state initials, with peq-row's
      # parameter overrides substituted into each SS RHS, then apply sim-row
      # overrides on top.
      ss <- mysteadies
      if (peq %in% rownames(conditions)) {
        for (cn in override_cols) {
          v <- conditions[peq, cn]
          if (is.na(v) || is.null(v)) next
          v <- trimws(as.character(v))
          if (v == "") next
          kind <- col_kind[[cn]]
          if (kind != "init")
            ss <- replaceSymbols(cn, v, ss)
          else
            ss[cn] <- v
        }
      }
      tr <- build_default()
      for (st in intersect(states, names(ss))) tr[st] <- as.character(ss[[st]])
      tr <- apply_row_overrides(tr, sim, "all")
      tr <- apply_petab_param_subs(tr, obs_str, noi_str)
      tr <- apply_scale_chain_rule(tr)
      parfns[[sub]] <- P(tr, condition = sub, compile = compile,
                         modelname = paste(modelname,
                                           sanitizeConditions(sub),
                                           sep = "_"))
      next
    }

    # ----- Pequil composition (numeric preeq) -----------------------------
    # 1. Substitute peq-row's parameter overrides directly into the reaction
    #    rates → peq_reactions. Parameter-scale chain rule must be applied
    #    *before* peq's parameter overrides, since peq overrides come from
    #    conditions.tsv (linear-scale literals) and would otherwise also be
    #    wrapped in exp/10^.
    peq_reactions <- reactions
    peq_param_subs <- list()
    peq_state_subs <- list()
    if (peq %in% rownames(conditions)) {
      for (cn in override_cols) {
        v <- conditions[peq, cn]
        if (is.na(v) || is.null(v)) next
        v <- trimws(as.character(v))
        if (v == "") next
        kind <- col_kind[[cn]]
        if (kind == "init") peq_state_subs[[cn]] <- v
        else                 peq_param_subs[[cn]] <- v
      }
    }
    if (length(peq_param_subs)) {
      peq_reactions$rates <- replaceSymbols(
        names(peq_param_subs), unlist(peq_param_subs), peq_reactions$rates)
    }

    # 2. p_pre: outer → c(state inits with peq state overrides applied,
    #    inner pars identity, observable/noise placeholders pre-substituted).
    tr_pre <- build_default()
    for (st in names(peq_state_subs)) tr_pre[st] <- peq_state_subs[[st]]
    tr_pre <- apply_petab_param_subs(tr_pre, obs_str, noi_str)
    tr_pre <- apply_scale_chain_rule(tr_pre)
    p_pre <- P(tr_pre, condition = sub, compile = compile,
               modelname = paste(modelname, sanitizeConditions(sub),
                                 "pre", sep = "_"))

    # 3. p_eq: integrate peq_reactions to steady state. attach.input ensures
    #    parameters and observable/noise placeholders flow through unchanged.
    p_eq <- Pequil(peq_reactions, condition = sub, attach.input = TRUE,
                   compile = compile,
                   modelname = paste(modelname, sanitizeConditions(sub),
                                     "eq", sep = "_"))

    # 4. p_post: identity for states (using SS values from p_eq), apply
    #    sim-row overrides (which include parameter overrides like a
    #    different k1 in case 0009, and may re-override a state like B=0
    #    in case 0010).
    tr_post <- setNames(inner_targets, inner_targets)
    tr_post <- apply_row_overrides(tr_post, sim, "all")
    p_post <- P(tr_post, condition = sub, compile = compile,
                modelname = paste(modelname, sanitizeConditions(sub),
                                  "post", sep = "_"))

    parfns[[sub]] <- p_post * p_eq * p_pre
  }

  Reduce(`+`, parfns)
}


# Internal: apply a ";"-separated PEtab parameter replacement string to a
# trafo by overriding every inner-target whose name matches
# <prefix>1_*, <prefix>2_*, ... with the corresponding entry in `repls`.
.petab_apply_param_substitution <- function(tr, repls, prefix) {
  parts <- trimws(strsplit(repls, ";", fixed = TRUE)[[1]])
  inner <- names(tr)
  pat   <- paste0("^", prefix, "([0-9]+)_(.+)$")
  for (sym in inner) {
    m <- regmatches(sym, regexec(pat, sym))[[1]]
    if (length(m) == 3L) {
      k <- as.integer(m[2])
      if (k >= 1L && k <= length(parts)) tr[sym] <- parts[k]
    }
  }
  tr
}


# Build the observation function `g`. Observables with `observableTransformation`
# of `log` or `log10` get wrapped at construction time (`obs_b -> log10(B)`),
# which keeps dMod's normL2 fast path on PEtab `{log,log10} × normal` cases:
# the residual is computed on the transformed scale on both sides of the
# subtraction (see .petab_parse_measurements for the data side).
.petab_build_observation_fn <- function(obs, obs_trafo, reactions,
                                        compile = TRUE,
                                        modelname = "petab_obs") {
  obs_eqn <- mapply(function(formula, trafo) {
    if (identical(trafo, "log"))   sprintf("log(%s)",   formula)
    else if (identical(trafo, "log10")) sprintf("log10(%s)", formula)
    else formula
  }, obs, obs_trafo[names(obs)], SIMPLIFY = TRUE, USE.NAMES = TRUE)
  Y(g          = as.eqnvec(obs_eqn),
    f          = as.eqnvec(reactions),
    attach.input = FALSE,
    compile    = compile,
    modelname  = modelname)
}


# Build the error model when at least one (sub-condition × observable) noise
# formula carries free symbols after PEtab parameter substitution.
#
# Strategy: build a single, condition-unspecific Y over the *original*
# noiseFormula in PEtab placeholder form (e.g. `noiseParameter1_obs_a`) plus
# the observable formulas as f-states. The placeholders are inner_targets
# of the trafo (.petab_build_trafo adds them to obs_inner) and the trafo's
# per-sub-condition `.petab_apply_param_substitution` rewrites them to the
# concrete row value (numeric literal or outer-parameter symbol). At runtime
# the post-trafo inner parameter vector pinner — which normL2 hands to the
# err model — therefore already carries the substituted value under the
# placeholder name, so a single Y suffices.
#
# Returns NULL when every (sub_cond, observable) noise formula evaluates to
# a constant after substitution — sigma is then carried by the data column
# and normL2's fast path handles the likelihood without an error model.
.petab_build_error_fn <- function(obs_meta, sub_cond_map, reactions,
                                  compile = TRUE, modelname = "petab_err") {

  any_symbolic <- any(vapply(seq_len(nrow(sub_cond_map)), function(ix) {
    obs_str <- sub_cond_map$observableParameters[ix]
    noi_str <- sub_cond_map$noiseParameters[ix]
    any(vapply(names(obs_meta$noise), function(obsId) {
      f <- obs_meta$noise[[obsId]]
      f <- .petab_substitute_param_string(f, obs_str, "observableParameter")
      f <- .petab_substitute_param_string(f, noi_str, "noiseParameter")
      is.na(.petab_eval_constant(f))
    }, logical(1)))
  }, logical(1)))
  if (!any_symbolic) return(NULL)

  # f for err Y: ODE states + observable formulas (so a noise formula could
  # reference an observable, e.g. relative noise `sigma * obs_a`). We pass the
  # *untransformed* observable formulas — the obs trafo (log/log10) doesn't
  # propagate to the noise model; PEtab sigma already lives on the
  # transformed scale by convention.
  obs_eqnvec <- as.eqnvec(setNames(unname(unlist(obs_meta$obs)),
                                   names(obs_meta$obs)))
  reactions_eqnvec <- as.eqnvec(reactions)

  Y(g            = as.eqnvec(obs_meta$noise),
    f            = c(reactions_eqnvec, obs_eqnvec),
    states       = names(obs_eqnvec),
    attach.input = FALSE,
    compile      = compile,
    modelname    = modelname)
}


# Build the prediction function `x` (Xs) for the chosen solver.
# `odemodel()` always compiles native code; the `compile` knob is only
# honoured by Y / P / compile().
.petab_build_odemodel <- function(reactions, solver,
                                  modelname = "petab_model") {
  m <- odemodel(reactions, modelname = modelname, solver = solver)
  Xs(m)
}


# Build the objective.
#
# Stage 2 supports `{lin, log, log10} × normal` noise: log/log10 enter via the
# pre-wrapped observation function (data values are matched-side transformed
# in .petab_parse_measurements), so dMod's normL2 fast path still gives the
# correct residual. Symbolic sigmas flow through the error model produced by
# .petab_build_error_fn.
#
# Non-normal distributions (laplace, log-normal) are not yet implemented;
# they would need a petabL2 objective with per-cell residual logic.
.petab_build_objective <- function(data, prd, errmodel, obs_meta) {

  if (!all(obs_meta$noise_dist == "normal"))
    stop("noiseDistribution(s) ",
         paste(unique(obs_meta$noise_dist[obs_meta$noise_dist != "normal"]),
               collapse = ", "),
         " are not yet supported. Only `normal` is implemented in Stage 2.")

  if (!all(obs_meta$obs_trafo %in% c("lin", "log", "log10")))
    stop("observableTransformation(s) ",
         paste(unique(obs_meta$obs_trafo[!obs_meta$obs_trafo %in%
                                          c("lin", "log", "log10")]),
               collapse = ", "),
         " are not recognised.")

  # PEtab's likelihood definition does not include the small-sample Bessel
  # correction that normL2 applies by default when an errmodel is present.
  # The correction also breaks numerically on the small test cases (n - p
  # can be <= 0 → sqrt of negative → NaN).
  base_obj <- normL2(data = data, x = prd, errmodel = errmodel,
                     use.bessel = FALSE, cores = 1L)

  # Data-coordinate Jacobian for log / log10 observable transformations.
  # PEtab's likelihood is on the linear y_obs:
  #   log L = log f_lin(y_obs) = log f_trafo(trafo(y_obs)) − log|d trafo / dy_obs|
  # so for log:    -2 log L includes  2 * Σ log(y_obs)
  # for log10:    -2 log L includes  2 * Σ log(y_obs · ln 10)
  # The offset depends only on the data, not on parameters, so gradients
  # and Hessians are unchanged — only the absolute objective value shifts
  # to match PEtab's -2*llh.
  jac_offset <- .petab_likelihood_offset(data, obs_meta)
  if (jac_offset == 0) return(base_obj)

  myfn <- function(..., fixed = NULL, deriv = TRUE, env = NULL) {
    out <- base_obj(..., fixed = fixed, deriv = deriv, env = env)
    out$value <- out$value + jac_offset
    attr_nm <- "data"
    if (!is.null(attr(out, attr_nm)))
      attr(out, attr_nm) <- attr(out, attr_nm) + jac_offset
    out
  }
  for (a in setdiff(names(attributes(base_obj)), "class"))
    attr(myfn, a) <- attr(base_obj, a)
  class(myfn) <- class(base_obj)
  myfn
}


# Internal: data-side log-likelihood offset for {log, log10} observables.
# `data` is a datalist whose `value` columns have already been transformed
# in .petab_parse_measurements (log or log10 applied). We invert that here
# to recover y_obs in linear units before computing the Jacobian factor.
.petab_likelihood_offset <- function(data, obs_meta) {
  off <- 0
  for (cn in names(data)) {
    df <- data[[cn]]
    for (obsId in unique(df$name)) {
      trafo <- obs_meta$obs_trafo[[obsId]]
      if (is.null(trafo) || identical(trafo, "lin")) next
      val <- df$value[df$name == obsId]
      y_obs <- if (identical(trafo, "log"))   exp(val)
               else if (identical(trafo, "log10")) 10^val
               else val
      if (identical(trafo, "log"))
        off <- off + 2 * sum(log(y_obs))
      else if (identical(trafo, "log10"))
        off <- off + 2 * sum(log(y_obs * log(10)))
    }
  }
  off
}


## --- public top-level entry point ------------------------------------------

#' Import a PEtab v1 problem into dMod
#'
#' Reads a PEtab YAML manifest plus the associated SBML model and four TSV
#' tables, and assembles a fully-composed dMod problem: prediction function,
#' observation function, parameter transformation, datalist, and objective.
#' The SBML side is delegated to [import_sbml()] (AMICI-based).
#'
#' Stage 1 supports test cases 0001–0006 of the PEtab test suite: single or
#' multi-condition models, condition-table parameter overrides (numeric or
#' parameter-symbol valued), per-row `observableParameters`, and Gaussian
#' likelihood on linear observables. Pre-equilibration, non-Gaussian noise,
#' observable transformations (log/log10) and symbolic noise formulas
#' currently raise a clear "Stage 2" error.
#'
#' @param yaml_path Path to the PEtab YAML manifest.
#' @param solver Required: one of `"deSolve"` or `"CppODE"`. Forwarded to
#'   [odemodel()].
#' @param amicipath Optional `PYTHONPATH` for AMICI; forwarded to
#'   [import_sbml()].
#' @param compile Logical. If `TRUE` (default) the generated trafo,
#'   observation function, and ODE model are compiled to native code. Set to
#'   `FALSE` for inspection-only use.
#' @param modelname Optional base modelname for the generated native files.
#'   Defaults to the YAML basename.
#' @param preeqMethod How to compute pre-equilibration steady states for
#'   measurement rows that carry a `preequilibrationConditionId`. One of
#'   `"analytic"` (default) — solve symbolically via [steadyStates()] and
#'   inject the closed-form expressions as state initial values, falling back
#'   to `"numeric"` automatically when the symbolic solver leaves any state
#'   unresolved — or `"numeric"` — always integrate the pre-equilibration
#'   condition to steady state via [Pequil()] and feed the result into the
#'   simulation phase.
#' @return A list with class `"petab_problem"` holding `prd`, `p`, `g`,
#'   `err`, `data`, `obj`, `pouter`, `lower`, `upper`, `fixed`,
#'   `condition.grid`, `observables`, `reactions`, `inits`, `model_id`,
#'   `source_yaml`, `sub_cond_map`, `obs_meta`, `param_meta`. The `inits`
#'   slot is a named character vector of (possibly symbolic) species initial
#'   expressions sourced from SBML; [exportPEtab()] re-emits non-numeric
#'   entries as `<initialAssignment>` elements so the roundtrip preserves
#'   parameter-bound initials. The `pouter` vector carries
#'   `attr(., "petab_scales")` recording each parameter's PEtab scale,
#'   which the exporter reads back.
#' @export
#' @example inst/examples/PEtabInterface.R
importPEtab <- function(yaml_path, solver, amicipath = NULL,
                        compile = TRUE, modelname = NULL,
                        preeqMethod = c("analytic", "numeric")) {

  preeqMethod <- match.arg(preeqMethod)

  if (missing(solver))
    stop("Argument `solver` is required (one of \"deSolve\", \"CppODE\").")
  solver <- match.arg(solver, c("deSolve", "CppODE"))

  yaml_path <- normalizePath(yaml_path, mustWork = TRUE)
  if (is.null(modelname))
    modelname <- sub("\\.ya?ml$", "", basename(yaml_path), ignore.case = TRUE)
  modelname <- gsub("[^A-Za-z0-9_]", "_", modelname)

  # 1) read tables and SBML
  tables <- read_petab_tables(yaml_path)
  sbml   <- import_sbml(tables$sbml_path, amicipath = amicipath)

  # 2) parse PEtab parameter / observable tables (independent of SBML)
  param_meta <- .petab_parse_parameters(tables$parameters)
  obs_meta   <- .petab_parse_observables(tables$observables)

  # 3) determine inner-parameter set of the model + observables
  states    <- sbml$reactions$states
  rate_syms <- unique(unlist(lapply(sbml$reactions$rates, function(r)
                getSymbols(r, exclude = c(states, "time")))))
  inner_pars <- setdiff(rate_syms, states)

  obs_syms   <- unique(unlist(lapply(obs_meta$obs, function(f)
                  getSymbols(f, exclude = c(states, "time")))))
  noise_syms <- unique(unlist(lapply(obs_meta$noise, function(f) {
                  if (is.na(suppressWarnings(as.numeric(f)))) getSymbols(f)
                  else character(0)
                })))
  obs_inner <- unique(c(obs_syms, noise_syms))
  obs_inner <- setdiff(obs_inner, c(states, inner_pars, "time"))

  # 4) parse condition / measurement tables (need symbol pool to type-check)
  cond_info <- .petab_parse_conditions(
                 tables$conditions,
                 sbml_states       = states,
                 sbml_compartments = names(sbml$reactions$compartments %||% list()),
                 sbml_pars         = names(sbml$pars),
                 obs_inner         = obs_inner)
  meas_info <- .petab_parse_measurements(tables$measurements, obs_meta)

  pouter_names <- names(param_meta$pouter)
  fixed        <- param_meta$fixed

  # SBML-default-as-fixed gap: any inner symbol referenced by rates / inits /
  # observables that is neither in parameters.tsv nor a state must come from
  # SBML's default parameter table (e.g. compartment volumes that PEtab leaves
  # unspecified). We pull those from sbml$pars and append to `fixed`. Symbols
  # without a default raise a clear error.
  init_syms <- unique(unlist(lapply(sbml$inits, function(e) {
    if (is.character(e)) getSymbols(e) else character(0)
  })))
  required <- unique(c(inner_pars, obs_inner, init_syms))
  required <- setdiff(required, c(states, "time", pouter_names, names(fixed)))
  # observableParameterK_<obsId> / noiseParameterK_<obsId> placeholders are
  # PEtab spec sentinels that get bound per-sub-condition by the trafo
  # (see .petab_apply_param_substitution). They never appear as outer pars.
  is_petab_placeholder <- grepl("^(observable|noise)Parameter[0-9]+_", required)
  required <- required[!is_petab_placeholder]
  if (length(required)) {
    miss <- setdiff(required, names(sbml$pars))
    if (length(miss))
      stop("Symbol(s) referenced by the model but not in parameters.tsv ",
           "and not in SBML defaults: ",
           paste(miss, collapse = ", "))
    extra_fixed <- as.numeric(sbml$pars[required])
    names(extra_fixed) <- required
    fixed <- c(fixed, extra_fixed)
  }

  # 5) build prediction / observation / trafo / objective
  outdir <- getwd()  # native files land here, per dMod convention

  ode <- .petab_build_odemodel(sbml$reactions, solver = solver,
                               modelname = paste0(modelname, "_ode"))
  g <- .petab_build_observation_fn(obs_meta$obs, obs_meta$obs_trafo,
                                   sbml$reactions,
                                   compile = compile,
                                   modelname = paste0(modelname, "_obs"))
  err <- .petab_build_error_fn(obs_meta, meas_info$sub_cond_map,
                               sbml$reactions, compile = compile,
                               modelname = paste0(modelname, "_err"))
  p  <- .petab_build_trafo(
          sub_cond_map  = meas_info$sub_cond_map,
          conditions    = cond_info$grid,
          col_kind      = cond_info$col_kind,
          override_cols = cond_info$override_cols,
          inits         = sbml$inits,
          sbml_pars     = sbml$pars,
          states        = states,
          inner_pars    = inner_pars,
          pouter_names  = pouter_names,
          fixed         = fixed,
          scales        = param_meta$scales,
          obs_inner     = obs_inner,
          reactions     = sbml$reactions,
          compile       = compile,
          modelname     = paste0(modelname, "_trafo"),
          preeqMethod   = preeqMethod)

  # Datalist: split data into per-condition data.frames keyed by sub_condition
  dat <- as.datalist(meas_info$data, split.by = "condition")

  prd <- g * ode * p

  obj <- .petab_build_objective(data = dat, prd = prd, errmodel = err,
                                obs_meta = obs_meta)

  # 6) assemble petab_problem
  pouter <- param_meta$pouter
  attr(pouter, "petab_scales") <- param_meta$scales[names(pouter)]

  out <- list(
    prd            = prd,
    p              = p,
    g              = g,
    err            = err,
    data           = dat,
    obj            = obj,
    pouter         = pouter,
    lower          = param_meta$lower,
    upper          = param_meta$upper,
    fixed          = fixed,
    condition.grid = cond_info$grid,
    observables    = obs_meta$obs,
    reactions      = sbml$reactions,
    inits          = sbml$inits,
    model_id       = modelname,
    source_yaml    = yaml_path,
    sub_cond_map   = meas_info$sub_cond_map,
    obs_meta       = obs_meta,
    param_meta     = param_meta
  )
  class(out) <- "petab_problem"
  out
}


#' Print method for `petab_problem`
#' @param x A `petab_problem`.
#' @param ... Unused.
#' @export
print.petab_problem <- function(x, ...) {
  cat("<petab_problem ", x$model_id, ">\n", sep = "")
  cat("  source:        ", x$source_yaml, "\n", sep = "")
  cat("  conditions:    ", nrow(x$sub_cond_map),
      " (", length(unique(x$sub_cond_map$simulationConditionId)),
      " sim, ", length(unique(x$sub_cond_map$preequilibrationConditionId[
        x$sub_cond_map$preequilibrationConditionId != ""])),
      " preeq)\n", sep = "")
  cat("  observables:   ", paste(names(x$observables), collapse = ", "), "\n",
      sep = "")
  cat("  measurements:  ", sum(vapply(x$data, nrow, 0L)), "\n", sep = "")
  cat("  pouter (n=", length(x$pouter), "): ",
      paste(names(x$pouter), collapse = ", "), "\n", sep = "")
  if (length(x$fixed) > 0L)
    cat("  fixed  (n=", length(x$fixed), "): ",
        paste(names(x$fixed), collapse = ", "), "\n", sep = "")
  invisible(x)
}


## --- exporter --------------------------------------------------------------

#' Export a dMod `petab_problem` back to PEtab v1
#'
#' Writes the five PEtab tables (parameters, observables, conditions,
#' measurements) plus the SBML model and a YAML manifest, recovering the
#' problem produced by [importPEtab()] up to algebraic equivalence on the
#' SBML side. Sub-conditions synthesised by the importer are collapsed back
#' to their PEtab condition + per-row `observableParameters` /
#' `noiseParameters` representation.
#'
#' Symbolic species initials (SBML `<initialAssignment>` elements) survive
#' the roundtrip: they are carried as-is on `petab$inits` and re-emitted by
#' [export_sbml()]. The reimported objective therefore matches the original
#' value at the same `pouter`.
#'
#' Lossy steps documented:
#' \itemize{
#'   \item AssignmentRules / RateRules are not yet read by [import_sbml()];
#'         models that rely on them are out of scope.
#'   \item Parameter scales survive only via `attr(petab$pouter,
#'         "petab_scales")` set by the importer; hand-built problems lacking
#'         this attribute default to `"lin"`.
#' }
#'
#' @param petab A `petab_problem` produced by [importPEtab()] (or hand-built
#'   list with the same slot names).
#' @param dir Output directory; created if missing.
#' @param model_id SBML model identifier (defaults to `petab$model_id` or
#'   `"dMod_export"`).
#' @param amicipath Forwarded to [export_sbml()].
#' @param overwrite Logical. If `FALSE` (default) errors when files already
#'   exist in `dir`.
#' @return Path to the written YAML manifest, invisibly.
#' @export
#' @importFrom yaml write_yaml
exportPEtab <- function(petab, dir, model_id = NULL, amicipath = NULL,
                        overwrite = FALSE) {

  stopifnot(inherits(petab, "petab_problem") || is.list(petab))
  if (is.null(model_id)) model_id <- petab$model_id %||% "dMod_export"

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  paths <- list(
    parameters   = file.path(dir, paste0("parameters_",   model_id, ".tsv")),
    observables  = file.path(dir, paste0("observables_",  model_id, ".tsv")),
    conditions   = file.path(dir, paste0("conditions_",   model_id, ".tsv")),
    measurements = file.path(dir, paste0("measurements_", model_id, ".tsv")),
    sbml         = file.path(dir, paste0(model_id, ".xml")),
    yaml         = file.path(dir, paste0(model_id, ".yaml"))
  )
  if (!overwrite) {
    existing <- vapply(paths, file.exists, logical(1))
    if (any(existing))
      stop("Output file(s) already exist; pass overwrite = TRUE to replace: ",
           paste(unlist(paths)[existing], collapse = ", "))
  }

  # parameters.tsv: pouter (estimate=1) + fixed (estimate=0)
  scales <- attr(petab$pouter, "petab_scales")
  if (is.null(scales))
    scales <- setNames(rep("lin", length(petab$pouter)), names(petab$pouter))

  pouter_df <- data.frame(
    parameterId    = names(petab$pouter),
    parameterScale = unname(scales[names(petab$pouter)]),
    lowerBound     = unname(petab$lower[names(petab$pouter)]),
    upperBound     = unname(petab$upper[names(petab$pouter)]),
    nominalValue   = unname(petab$pouter),
    estimate       = 1L,
    stringsAsFactors = FALSE
  )
  fixed_df <- if (length(petab$fixed)) data.frame(
    parameterId    = names(petab$fixed),
    parameterScale = "lin",
    lowerBound     = -Inf,
    upperBound     = Inf,
    nominalValue   = unname(petab$fixed),
    estimate       = 0L,
    stringsAsFactors = FALSE
  ) else NULL

  utils::write.table(rbind(pouter_df, fixed_df), paths$parameters,
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  # observables.tsv
  om <- petab$obs_meta
  obs_df <- data.frame(
    observableId             = names(om$obs),
    observableFormula        = unname(om$obs),
    observableTransformation = unname(om$obs_trafo[names(om$obs)]),
    noiseFormula             = unname(om$noise[names(om$obs)]),
    noiseDistribution        = unname(om$noise_dist[names(om$obs)]),
    stringsAsFactors = FALSE
  )
  utils::write.table(obs_df, paths$observables, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")

  # conditions.tsv: collapse sub_cond_map back to original conditions
  scm <- petab$sub_cond_map
  uniq_sims <- unique(scm$simulationConditionId)
  cond_df <- petab$condition.grid[uniq_sims, , drop = FALSE]
  cond_df$conditionId <- uniq_sims
  cond_df <- cond_df[, c("conditionId",
                         setdiff(colnames(cond_df), "conditionId")),
                     drop = FALSE]
  utils::write.table(cond_df, paths$conditions, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")

  # measurements.tsv: long-format with re-attached obs/noise param strings
  meas_rows <- do.call(rbind, lapply(names(petab$data), function(sub) {
    df <- petab$data[[sub]]
    ix <- match(sub, scm$sub_condition)
    data.frame(
      observableId               = df$name,
      simulationConditionId      = scm$simulationConditionId[ix],
      time                       = df$time,
      measurement                = df$value,
      preequilibrationConditionId = scm$preequilibrationConditionId[ix],
      observableParameters       = scm$observableParameters[ix],
      noiseParameters            = scm$noiseParameters[ix],
      stringsAsFactors = FALSE
    )
  }))
  # drop empty optional cols
  for (col in c("preequilibrationConditionId", "observableParameters",
                "noiseParameters")) {
    if (all(meas_rows[[col]] == "")) meas_rows[[col]] <- NULL
  }
  utils::write.table(meas_rows, paths$measurements, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")

  # SBML. Symbolic initials (whose strings reference parameters from
  # parameters.tsv or compartment volumes) get re-emitted as
  # <initialAssignment> elements; numeric initials use initialConcentration.
  # See export_sbml() for the dispatch on character vs. numeric `inits`.
  all_pars <- c(petab$pouter, petab$fixed)
  inits <- petab$inits %||%
           setNames(rep(0, length(petab$reactions$states)),
                    petab$reactions$states)
  export_sbml(petab$reactions, parameters = all_pars, inits = inits,
              filepath = paths$sbml, model_id = model_id,
              amicipath = amicipath)

  # YAML manifest
  manifest <- list(
    format_version = 1L,
    parameter_file = basename(paths$parameters),
    problems = list(list(
      sbml_files        = list(basename(paths$sbml)),
      condition_files   = list(basename(paths$conditions)),
      measurement_files = list(basename(paths$measurements)),
      observable_files  = list(basename(paths$observables))
    ))
  )
  yaml::write_yaml(manifest, paths$yaml)

  invisible(paths$yaml)
}


# `%||%` is base R since 4.4; relied on per dMod's existing usage.
