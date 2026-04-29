## R/PEtabInterface.R
## ---------------------------------------------------------------------------
## PEtab v1 importer / exporter, layered on top of dMod's existing
## SBML interface (R/SBMLinterface.R) and high-level APIs (Y, P, normL2, ...).
##
## Public API: importPEtab, exportPEtab, exportPEtabObject,
##             read_petab_yaml, read_petab_tables.
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

  # PEtab v1 spec: nominalValue / lowerBound / upperBound are written on the
  # *linear* scale, regardless of parameterScale. dMod's outer parameters
  # live on the chosen parameter scale (the trafo applies `10^x` / `exp(x)`
  # in apply_scale_chain_rule), so we pre-transform here:
  #   log10 → log10(.),   log → log(.),   lin → identity.
  to_num <- function(col) suppressWarnings(as.numeric(col))
  apply_fwd_scale <- function(values, ids) {
    sc <- scales[ids]
    out <- values
    log_idx   <- which(sc == "log")
    log10_idx <- which(sc == "log10")
    if (length(log_idx))   out[log_idx]   <- log(values[log_idx])
    if (length(log10_idx)) out[log10_idx] <- log10(values[log10_idx])
    out
  }

  pouter_idx <- which(df$estimate == 1)
  fixed_idx  <- which(df$estimate == 0)

  pouter_ids <- df$parameterId[pouter_idx]
  fixed_ids  <- df$parameterId[fixed_idx]

  pouter <- setNames(apply_fwd_scale(to_num(df$nominalValue[pouter_idx]),
                                     pouter_ids),
                     pouter_ids)
  # `fixed` parameters are passed straight through to the trafo as numeric
  # constants on the *inner* (linear) scale — no scale wrapping in the
  # trafo, so we keep linear values here regardless of parameterScale.
  fixed  <- setNames(to_num(df$nominalValue[fixed_idx]), fixed_ids)
  lower  <- setNames(apply_fwd_scale(to_num(df$lowerBound[pouter_idx]),
                                     pouter_ids),
                     pouter_ids)
  upper  <- setNames(apply_fwd_scale(to_num(df$upperBound[pouter_idx]),
                                     pouter_ids),
                     pouter_ids)

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

  # When the per-row noise resolves to a finite numeric, the value is carried
  # by the `sigma` column and consumed via normL2's fast path; the trafo's
  # placeholder substitution is dead code (errmodel will be NULL for such
  # observables). Collapse those rows' `noiseParameters` strings to a
  # uniform literal so they share one sub-condition rather than spawning one
  # compiled parameter trafo per unique sigma value. The literal preserves
  # the original ";"-separated arity so any K-th placeholder reference still
  # resolves cleanly.
  numeric_noise <- !is.na(sigma) & nzchar(m$noiseParameters)
  if (any(numeric_noise)) {
    n_parts <- lengths(strsplit(m$noiseParameters[numeric_noise],
                                ";", fixed = TRUE))
    m$noiseParameters[numeric_noise] <- vapply(
      pmax(n_parts, 1L),
      function(n) paste(rep("1", n), collapse = ";"),
      character(1))
  }

  # Sub-condition key: (simCondId, peq, obsParStr, noiseParStr).
  # A suffix is applied only when a single simCondId carries more than one
  # distinct tuple — otherwise the original simCondId is kept as the
  # sub_condition name (no synthetic ``cond__hash`` artifacts).
  obs_hash <- .petab_subcond_hash(m$observableParameters)
  noi_hash <- .petab_subcond_hash(m$noiseParameters)
  peq_hash <- .petab_subcond_hash(m$preequilibrationConditionId)

  tuple_key <- paste(peq_hash, obs_hash, noi_hash, sep = "")
  sub_cond  <- m$simulationConditionId
  for (sc in unique(m$simulationConditionId)) {
    ix <- which(m$simulationConditionId == sc)
    if (length(unique(tuple_key[ix])) <= 1L) next
    for (k in unique(tuple_key[ix])) {
      ix2 <- ix[tuple_key[ix] == k]
      ph  <- peq_hash[ix2[1L]]
      oh  <- obs_hash[ix2[1L]]
      nh  <- noi_hash[ix2[1L]]
      parts <- c(if (nzchar(ph)) ph,
                 if (nzchar(oh)) oh,
                 if (nzchar(nh)) nh)
      sub_cond[ix2] <- paste0(sc, "__", paste(parts, collapse = "_"))
    }
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


# Build the underlying odemodel and the Xs() prediction function for the
# chosen solver. `compile` is forwarded to cOde::funC / CppODE::CppODE so
# the importer can defer linking until a single batched compile().
.petab_build_odemodel <- function(reactions, solver,
                                  modelname = "petab_model",
                                  compile = TRUE) {
  m <- odemodel(reactions, modelname = modelname, solver = solver,
                compile = compile)
  list(odemodel = m, x = Xs(m))
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
#' @param cores Number of parallel compilation jobs (Unix only) forwarded to
#'   [compile()]. The importer batches all generated source files into a
#'   single `compile()` call so this directly controls native-build
#'   concurrency.
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
#' @return A list with class `"PEtabProblem"` holding `dataList`,
#'   `reactions`, `odemodel`, `g`, `x`, `p`, `e`, `prd` (the composite
#'   `g * x * p`), `obj`, `bestfit`, `parlower`, `parupper`. The `obj`
#'   closure has the PEtab fixed parameters baked in, so calling
#'   `obj(bestfit)` evaluates the likelihood at the current estimate
#'   without the user having to pass `fixed = ...`. The `bestfit` vector
#'   carries `attr(., "petab_scales")` recording each parameter's PEtab
#'   scale (the exporter reads it back). Internal metadata needed for
#'   the round-trip exporter (`fixed`, `inits`, `model_id`,
#'   `source_yaml`, `sub_cond_map`, `obs_meta`, `param_meta`) lives on
#'   `attr(., "petab_meta")`.
#' @export
#' @example inst/examples/PEtabInterface.R
importPEtab <- function(yaml_path, solver, amicipath = NULL,
                        compile = TRUE, cores = 1L, modelname = NULL,
                        preeqMethod = c("analytic", "numeric")) {

  preeqMethod <- match.arg(preeqMethod)
  cores <- as.integer(cores)
  if (length(cores) != 1L || is.na(cores) || cores < 1L)
    stop("`cores` must be a single positive integer.")

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
  #
  # Defer per-call compilation: each helper writes generated source files but
  # leaves linking to a single batched compile() so we can honour `cores` and
  # avoid one g++ invocation per sub-condition.
  outdir <- getwd()  # native files land here, per dMod convention

  ode_pair <- .petab_build_odemodel(sbml$reactions, solver = solver,
                                    modelname = paste0(modelname, "_ode"),
                                    compile = FALSE)
  odeobj <- ode_pair$odemodel
  x      <- ode_pair$x
  g <- .petab_build_observation_fn(obs_meta$obs, obs_meta$obs_trafo,
                                   sbml$reactions,
                                   compile = FALSE,
                                   modelname = paste0(modelname, "_obs"))
  e <- .petab_build_error_fn(obs_meta, meas_info$sub_cond_map,
                             sbml$reactions, compile = FALSE,
                             modelname = paste0(modelname, "_err"))
  p <- .petab_build_trafo(
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
          compile       = FALSE,
          modelname     = paste0(modelname, "_trafo"),
          preeqMethod   = preeqMethod)

  # Datalist: split data into per-condition data.frames keyed by sub_condition
  dataList <- as.datalist(meas_info$data, split.by = "condition")

  prd <- g * x * p

  if (isTRUE(compile)) {
    if (is.null(e)) dMod::compile(prd, cores = cores)
    else            dMod::compile(prd, e, cores = cores)
  }

  raw_obj <- .petab_build_objective(data = dataList, prd = prd, errmodel = e,
                                    obs_meta = obs_meta)

  # Bake `fixed` into the objective so users only pass the estimated bestfit
  # vector. `fixed = NULL` (the default) lets the closure inject the PEtab
  # fixed values; an explicit `fixed = ...` overrides per call.
  baked_fixed <- fixed
  obj <- function(pars, fixed = NULL, ...) {
    if (is.null(fixed)) fixed <- baked_fixed
    raw_obj(pars, fixed = fixed, ...)
  }

  # 6) assemble PEtabProblem
  bestfit <- param_meta$pouter
  attr(bestfit, "petab_scales") <- param_meta$scales[names(bestfit)]

  out <- list(
    dataList  = dataList,
    reactions = sbml$reactions,
    odemodel  = odeobj,
    g         = g,
    x         = x,
    p         = p,
    e         = e,
    prd       = prd,
    obj       = obj,
    bestfit   = bestfit,
    parlower  = param_meta$lower,
    parupper  = param_meta$upper
  )
  attr(out, "petab_meta") <- list(
    fixed        = fixed,
    inits        = sbml$inits,
    model_id     = modelname,
    source_yaml  = yaml_path,
    sub_cond_map = meas_info$sub_cond_map,
    obs_meta     = obs_meta,
    param_meta   = param_meta
  )
  class(out) <- "PEtabProblem"
  out
}


#' Print method for `PEtabProblem`
#' @param x A `PEtabProblem`.
#' @param ... Unused.
#' @export
print.PEtabProblem <- function(x, ...) {
  meta  <- attr(x, "petab_meta") %||% list()
  scm   <- meta$sub_cond_map
  obs   <- meta$obs_meta$obs
  fixed <- meta$fixed
  cat("<PEtabProblem ", meta$model_id %||% "", ">\n", sep = "")
  if (!is.null(meta$source_yaml))
    cat("  source:        ", meta$source_yaml, "\n", sep = "")
  if (!is.null(scm))
    cat("  conditions:    ", nrow(scm),
        " (", length(unique(scm$simulationConditionId)),
        " sim, ", length(unique(scm$preequilibrationConditionId[
          scm$preequilibrationConditionId != ""])),
        " preeq)\n", sep = "")
  if (!is.null(obs))
    cat("  observables:   ", paste(names(obs), collapse = ", "), "\n",
        sep = "")
  cat("  measurements:  ", sum(vapply(x$dataList, nrow, 0L)), "\n", sep = "")
  cat("  bestfit (n=", length(x$bestfit), "): ",
      paste(names(x$bestfit), collapse = ", "), "\n", sep = "")
  if (length(fixed) > 0L)
    cat("  fixed   (n=", length(fixed), "): ",
        paste(names(fixed), collapse = ", "), "\n", sep = "")
  invisible(x)
}


## --- exporter --------------------------------------------------------------
##
## The forward direction is the symbolic inverse of `.petab_build_trafo`
## (line 441). Given a dMod parfn `p`, we read its per-condition LHS=RHS
## eqnvec via `getEquations(p)`, strip the parameter-scale chain rule
## (replacing `10^(op)` / `exp(op)` subexpressions with bare `op`), and
## classify each LHS — the resulting (a) constants land as SBML defaults
## or initialConcentrations, (b) constant symbolic mappings land as
## conditions.tsv columns or initialAssignments, (c) per-condition-varying
## mappings land as conditions.tsv override columns, (d) identity mappings
## are the importer's default and need no emission. The single source of
## truth is `getEquations(p)`; `attr(data, "condition.grid")` is ignored.

# Internal: does `e` syntactically match the bare symbol `op_sym`, possibly
# wrapped in redundant parentheses (`(X)` parses as a call to `(`)? dMod's
# `repar` emits `10^(K)` which parses as `^(10, (K))` — i.e. the exponent
# slot is the parenthesis-call, not the bare K symbol. We treat both forms
# as equivalent.
.petab_is_bare <- function(e, op_sym) {
  if (is.symbol(e) && identical(e, op_sym)) return(TRUE)
  if (is.call(e) && length(e) == 2L &&
      identical(e[[1L]], as.symbol("(")) &&
      .petab_is_bare(e[[2L]], op_sym)) return(TRUE)
  FALSE
}

# Compensate every occurrence of a single log/log10-scaled outer parameter
# `op` inside `expr` for the importer's chain rule wrap (which substitutes
# `op → 10^(op)` for log10 / `op → exp(op)` for log on every RHS). Two
# cases per occurrence:
#
#   - Direct exponent of `10^(.)` for log10 (or argument of `exp(.)` for
#     log): leave the inner symbol bare. After import, chain rule re-wraps
#     it so the original `10^(op)` form is reproduced verbatim.
#
#   - Anywhere else (compound expression, bare symbol elsewhere): wrap
#     with `log10(.)` (or `log(.)`). After chain rule, `log10(10^(op))`
#     simplifies (numerically) to `op`, leaving the surrounding compound
#     expression untouched. This makes round-trip work for any
#     algebraic combination, e.g. `10^(KM + 5)` or `K1 * K2 + offset`.
#
# Mixed wrap (e.g. `10^(K) + K`) is no longer ambiguous — the clean wrap
# stays clean, the bare K gets compensated independently.
.petab_compensate_chain_rule <- function(expr, op, scale) {
  op_sym  <- as.symbol(op)
  if (scale == "log10") {
    inv_head     <- as.symbol("log10")
    is_clean_wrap <- function(e) {
      is.call(e) && length(e) == 3L && identical(e[[1L]], as.symbol("^")) &&
      is.numeric(e[[2L]]) && length(e[[2L]]) == 1L && e[[2L]] == 10 &&
      .petab_is_bare(e[[3L]], op_sym)
    }
  } else {  # "log"
    inv_head     <- as.symbol("log")
    is_clean_wrap <- function(e) {
      is.call(e) && length(e) == 2L && identical(e[[1L]], as.symbol("exp")) &&
      .petab_is_bare(e[[2L]], op_sym)
    }
  }
  walk <- function(e) {
    if (is_clean_wrap(e)) return(op_sym)
    if (is.call(e)) {
      for (i in seq_along(e)[-1L]) e[[i]] <- walk(e[[i]])
      return(e)
    }
    if (is.symbol(e) && identical(e, op_sym))
      return(call(as.character(inv_head), op_sym))
    e
  }
  walk(expr)
}

# Compensate a single RHS string for the importer's chain rule, for every
# log/log10-scaled outer parameter at once. Pure numeric literals pass
# through unchanged.
.petab_strip_param_scale <- function(rhs_str, scales) {
  s <- as.character(rhs_str)
  if (!is.na(suppressWarnings(as.numeric(s)))) return(s)
  expr <- tryCatch(parse(text = s, keep.source = FALSE)[[1L]],
                   error = function(e)
                     stop("Cannot parse trafo RHS `", s, "`: ",
                          conditionMessage(e), call. = FALSE))
  log10_pars <- names(scales)[scales == "log10"]
  log_pars   <- names(scales)[scales == "log"]
  for (op in log10_pars)
    expr <- .petab_compensate_chain_rule(expr, op, "log10")
  for (op in log_pars)
    expr <- .petab_compensate_chain_rule(expr, op, "log")
  paste(deparse(expr, width.cutoff = 500L), collapse = "")
}

# Apply the strip to every per-condition eqnvec.
.petab_strip_trafo <- function(eqs, scales) {
  lapply(eqs, function(eqv) {
    nm  <- names(eqv)
    raw <- as.character(unclass(eqv))
    out <- vapply(raw,
                  function(rhs) .petab_strip_param_scale(rhs, scales),
                  character(1))
    names(out) <- nm
    out
  })
}

# Classify a single LHS across conditions:
#   "missing"        — at least one condition has no entry for this LHS
#   "identity"       — every condition's RHS == LHS (importer's build_default)
#   "const_numeric"  — every condition's RHS is the same numeric literal
#   "const_symbolic" — every condition's RHS is the same symbolic formula
#   "varying"        — RHS differs across conditions
.petab_classify_lhs <- function(stripped_eqs, lhs, conds) {
  rhs <- vapply(conds, function(c) {
    e <- stripped_eqs[[c]]
    v <- if (lhs %in% names(e)) e[[lhs]] else NA_character_
    if (length(v) == 0L || is.na(v)) NA_character_ else as.character(v)
  }, character(1))
  if (any(is.na(rhs))) return(list(kind = "missing"))
  if (all(rhs == lhs)) return(list(kind = "identity"))
  # Try to collapse closed-form numeric expressions (e.g. "10^0" → 1) so they
  # land as SBML defaults / initialConcentrations rather than conditions.tsv
  # columns referencing literal arithmetic.
  vals <- vapply(rhs, .petab_eval_constant, numeric(1))
  if (length(unique(rhs)) == 1L) {
    if (!is.na(vals[[1L]]))
      return(list(kind = "const_numeric", value = unname(vals[[1L]])))
    return(list(kind = "const_symbolic", formula = unname(rhs[[1L]])))
  }
  # Per-condition values: if every cell evaluates to a constant, store the
  # numeric values; otherwise pass the raw formulas through.
  if (all(!is.na(vals)))
    return(list(kind = "varying",
                per_cond = setNames(as.character(unname(vals)), conds)))
  list(kind = "varying", per_cond = rhs)
}

# Decompose the per-condition trafo into PEtab v1 building blocks.
# Returns:
#   conditions_df   data.frame with conditionId column + override columns
#                   (NA where no override). Rownames = condition names.
#   inits           named list keyed by state — numeric (initialConcentration)
#                   or character (initialAssignment formula).
#   sbml_extra_pars named numeric — additional fixed-default parameters that
#                   must appear in SBML so condition.tsv columns and
#                   collapsed inner_pars resolve to a declared SId.
.petab_decompose_trafo <- function(eqs, states, inner_pars, obs_inner,
                                   pouter_names, fixed, scales) {

  conds        <- names(eqs)
  stripped     <- .petab_strip_trafo(eqs, scales)
  inner_targets <- unique(c(states, inner_pars, obs_inner))

  inits          <- list()
  cond_overrides <- list()
  sbml_extra_pars <- numeric(0)

  declare_extra <- function(nm, val) {
    if (!nm %in% c(pouter_names, names(fixed)))
      sbml_extra_pars[[nm]] <<- val
  }

  for (lhs in inner_targets) {
    cls    <- .petab_classify_lhs(stripped, lhs, conds)
    state  <- lhs %in% states

    switch(cls$kind,
      missing = stop(sprintf(
        "Trafo does not cover %s `%s` for at least one condition. Every state and inner parameter must appear on the trafo's LHS.",
        if (state) "state" else "inner parameter", lhs), call. = FALSE),

      identity = {
        if (state) {
          # State LHS = bare RHS of the same name. The trafo will resolve
          # this at runtime via an outer parameter (or fixed) of the same
          # name. Emit as an SBML <initialAssignment> referencing that
          # parameter; the importer's chain rule will then reproduce
          # `state = 10^(state)` (when the parameter is on log10 scale).
          if (lhs %in% pouter_names || lhs %in% names(fixed))
            inits[[lhs]] <- lhs
          # otherwise leave to the post-loop default (0) — undeclared
          # identity on a state means "init = 0 unless the user pouter has
          # a name match", and the trafo's getSymbols would already have
          # raised it as an undeclared symbol.
        } else {
          declare_extra(lhs, 1)
        }
      },

      const_numeric = {
        if (state) inits[[lhs]] <- cls$value
        else       declare_extra(lhs, cls$value)
      },

      const_symbolic = {
        if (state) {
          inits[[lhs]] <- cls$formula
        } else {
          # Constant-across-conditions inner_par mapping. Emit as a
          # conditions.tsv column with the same value in every row, and
          # register the inner par as an SBML SId (placeholder default).
          cond_overrides[[lhs]] <- setNames(rep(cls$formula, length(conds)),
                                            conds)
          declare_extra(lhs, 1)
        }
      },

      varying = {
        cond_overrides[[lhs]] <- cls$per_cond[conds]
        if (state) {
          # SBML default for a per-condition state init: a placeholder; the
          # importer reads the per-row override anyway.
          v0 <- suppressWarnings(as.numeric(cls$per_cond[[1L]]))
          inits[[lhs]] <- if (!is.na(v0)) v0 else 0
        } else {
          declare_extra(lhs, 1)
        }
      }
    )

    # Free-symbol sanity check on the post-strip RHS.
    if (cls$kind %in% c("const_symbolic", "varying")) {
      rhs_strs <- if (cls$kind == "varying") cls$per_cond else cls$formula
      free <- unique(unlist(lapply(rhs_strs, function(s) {
        if (!is.na(suppressWarnings(as.numeric(s)))) character(0)
        else getSymbols(s, exclude = c(states, "time"))
      })))
      declared <- c(pouter_names, names(fixed), names(sbml_extra_pars))
      miss <- setdiff(free, declared)
      if (length(miss))
        stop(sprintf(
          "Trafo RHS for `%s` references undeclared symbol(s): %s. Add to `pouter` or `fixed`.",
          lhs, paste(miss, collapse = ", ")), call. = FALSE)
    }
  }

  # Every state must have an init recorded (default 0).
  for (st in setdiff(states, names(inits))) inits[[st]] <- 0

  conditions_df <- data.frame(conditionId      = conds,
                              row.names        = conds,
                              stringsAsFactors = FALSE)
  for (cn in names(cond_overrides))
    conditions_df[[cn]] <- unname(cond_overrides[[cn]])

  list(conditions_df   = conditions_df,
       inits           = inits,
       sbml_extra_pars = sbml_extra_pars)
}


#' Export a dMod problem to PEtab v1
#'
#' Symbolically decomposes a dMod parameter transformation `p` and writes
#' the corresponding PEtab v1 problem (parameters / observables / conditions
#' / measurements TSVs, an SBML model, and a YAML manifest). The decomposer
#' inverts the importer's [importPEtab()] trafo construction so the
#' roundtrip preserves outer parameter names, scales, and per-condition
#' overrides; the imported PEtab problem is fitted on the same outer
#' parameter vector as the dMod-native problem.
#'
#' @section Algorithm:
#' Given `p`, the exporter reads `getEquations(p)` (one eqnvec per
#' condition keyed by inner-side LHS), strips the parameter-scale chain
#' rule (replacing `10^(op)` with `op` for each `op` declared on log10
#' scale; symmetric for `log`), and classifies each LHS. The
#' classification drives where the mapping is emitted:
#' \itemize{
#'   \item Constant-numeric state inits → SBML `<initialConcentration>`.
#'   \item Constant-symbolic state inits → SBML `<initialAssignment>`.
#'   \item Constant-symbolic inner-parameter mappings → conditions.tsv
#'         column with the same value in every row.
#'   \item Per-condition-varying mappings → conditions.tsv column with
#'         per-row values.
#'   \item Identity mappings (RHS == LHS) → emitted nowhere; the importer's
#'         `build_default` reproduces them.
#' }
#' Inner parameters that survive only as conditions.tsv overrides (or
#' collapsed-numeric SBML defaults) are appended to `fixed` with placeholder
#' values so the SBML model declares them as SIds and the importer's
#' required-symbol check finds them.
#'
#' @section Limitations:
#' \itemize{
#'   \item Pre-equilibration is not yet roundtrippable through `p`. Use
#'         [exportPEtabObject()] on the original `petabProblem` if you
#'         need it.
#'   \item Per-row data sigmas are not preserved; supply an explicit
#'         `errors` formula (e.g. parameterised noise) if you need
#'         non-constant noise.
#'   \item Fixed parameters are always written with `parameterScale = "lin"`
#'         (PEtab+dMod convention; cf. `exportPEtabObject` line 1387).
#'         Wrapping a fixed-name in `10^(.)` inside the trafo will trigger
#'         the mixed-wrap error.
#' }
#'
#' @param data A [datalist] (or a list of data.frames keyed by condition,
#'   each with `name`, `time`, `value` columns; or a long-format data.frame
#'   that [as.datalist()] accepts). Only used as the measurement source —
#'   `attr(data, "condition.grid")` is ignored.
#' @param reactions An [eqnlist] describing the ODE network.
#' @param observables Observable formulas keyed by observableId. Accepts a
#'   named character vector, an [eqnvec], or an observation function produced
#'   by [Y()] (formulas read from `attr(observables, "equations")`).
#' @param p A `parfn` produced by [P()]. Required: the symbolic trafo is
#'   the source of truth for parameters.tsv, conditions.tsv, and
#'   `<initialAssignment>` formulas.
#' @param pouter Named numeric vector of estimated outer parameters, on the
#'   chosen `parameterScale`. Names become `parameterId`s in parameters.tsv.
#' @param errors Noise formulas keyed by observableId. Accepts a named
#'   character vector, [eqnvec], or a Y-built error function. If `NULL`,
#'   defaults to `"1"` per observable (constant unit noise).
#' @param lower,upper Named numeric vectors of bounds, on the same scale as
#'   `pouter`. If `NULL`, written as `-Inf`/`Inf`.
#' @param fixed Optional named numeric vector of non-estimated parameters
#'   on the linear scale.
#' @param parameterScale Scalar `"lin"`/`"log"`/`"log10"` (broadcast to all
#'   names in `pouter`) or a named character vector keyed by parameterId.
#'   Defaults to `"log10"` (matches dMod's `insert("x ~ 10^X", ...)` idiom).
#' @param observableTransformation Scalar or named character —
#'   `"lin"`/`"log"`/`"log10"` per observableId.
#' @param noiseDistribution Scalar or named character — `"normal"`/
#'   `"laplace"`/`"log-normal"` per observableId.
#' @param model_id SBML model identifier; defaults to `"dMod_export"`.
#' @param amicipath Forwarded to [export_sbml()].
#' @param dir Output directory; created if missing.
#' @param overwrite Whether to overwrite existing PEtab files in `dir`.
#' @return Path to the written YAML manifest, invisibly.
#' @seealso [exportPEtabObject()] for the lower-level entry that takes a
#'   pre-assembled `petabProblem` list. [importPEtab()] for the inverse.
#' @export
exportPEtab <- function(data, reactions, observables, p, pouter,
                        errors = NULL,
                        lower = NULL, upper = NULL, fixed = NULL,
                        parameterScale = "log10",
                        observableTransformation = "lin",
                        noiseDistribution = "normal",
                        model_id = "dMod_export",
                        amicipath = NULL,
                        dir, overwrite = FALSE) {

  ## --- 1. datalist normalisation ------------------------------------------
  if (is.data.frame(data)) data <- as.datalist(data)
  if (is.list(data) && !inherits(data, "datalist")) data <- as.datalist(data)
  if (!inherits(data, "datalist"))
    stop("`data` must be a datalist, list of data.frames, or long-format data.frame.")

  ## --- 2. extract per-condition trafo from p ------------------------------
  if (!is.function(p) || is.null(attr(p, "mappings")))
    stop("`p` must be a parfn produced by P() — needed to decompose the parameter trafo.")
  eqs <- getEquations(p)
  if (!is.list(eqs)) eqs <- list(eqs)
  conds <- names(eqs)
  if (is.null(conds) || any(!nzchar(conds)))
    stop("`getEquations(p)` returned an unnamed list — every condition must have a name.")

  ## --- 3. reactions / observables / errors --------------------------------
  if (!inherits(reactions, "eqnlist"))
    stop("`reactions` must be an eqnlist (the ODE network with stoichiometry).")
  states <- reactions$states

  pull_eqns <- function(obj, what) {
    if (is.function(obj)) {
      eqns <- attr(obj, "equations")
      if (is.null(eqns))
        stop(sprintf(
          "`%s` is a function but carries no `equations` attribute.", what))
      if (is.list(eqns)) eqns <- eqns[[1L]]
      return(setNames(as.character(eqns), names(eqns)))
    }
    if (inherits(obj, "eqnvec"))
      return(setNames(as.character(unclass(obj)), names(obj)))
    if (is.character(obj) && !is.null(names(obj))) return(obj)
    stop(sprintf("`%s` must be a named character vector, eqnvec, or Y-built fn.", what))
  }
  obs_eqns <- pull_eqns(observables, "observables")
  if (is.null(names(obs_eqns)) || any(!nzchar(names(obs_eqns))))
    stop("`observables` entries must be named (observableId -> formula).")
  obs_ids <- names(obs_eqns)

  # Default error model: if the datalist carries a non-NA `sigma` column,
  # encode per-row sigmas via PEtab's `noiseParameter1_<obsId>` placeholder
  # so per-row noise survives the round-trip (the importer reads
  # `noiseParameters` column from measurements.tsv and substitutes it back
  # in). Falls back to constant `"1"` only if the datalist has no sigma.
  has_per_row_sigma <- any(vapply(data, function(d) {
    "sigma" %in% colnames(d) && any(!is.na(d$sigma))
  }, logical(1)))
  if (is.null(errors)) {
    err_eqns <- if (has_per_row_sigma)
      setNames(paste0("noiseParameter1_", obs_ids), obs_ids)
    else
      setNames(rep("1", length(obs_ids)), obs_ids)
  } else {
    err_eqns <- pull_eqns(errors, "errors")
  }
  miss <- setdiff(obs_ids, names(err_eqns))
  if (length(miss))
    stop("`errors` is missing entries for: ", paste(miss, collapse = ", "))
  err_eqns <- err_eqns[obs_ids]
  use_per_row_noise <- is.null(errors) && has_per_row_sigma

  ## --- 4. inner_pars / obs_inner (mirror importer logic) ------------------
  rate_syms <- unique(unlist(lapply(reactions$rates, function(r)
                getSymbols(r, exclude = c(states, "time")))))
  inner_pars <- setdiff(rate_syms, states)
  obs_syms <- unique(unlist(lapply(obs_eqns, function(f)
                getSymbols(f, exclude = c(states, "time")))))
  noise_syms <- unique(unlist(lapply(err_eqns, function(f) {
                if (is.na(suppressWarnings(as.numeric(f)))) getSymbols(f)
                else character(0)
              })))
  obs_inner <- setdiff(unique(c(obs_syms, noise_syms)),
                       c(states, inner_pars, "time"))
  # PEtab placeholders (`observableParameter<k>_<obsId>` /
  # `noiseParameter<k>_<obsId>`) are spec sentinels bound per-row in
  # measurements.tsv, not real inner parameters — strip them so the
  # trafo decomposer doesn't expect a mapping for them. Mirror of the
  # importer's filter at the required-symbol gap-fill (line 935-938).
  obs_inner <- obs_inner[!grepl(
    "^(observable|noise)Parameter[0-9]+_", obs_inner)]

  ## --- 5. validate pouter / fixed / scales --------------------------------
  pouter_ids <- names(pouter)
  if (is.null(pouter_ids) || any(!nzchar(pouter_ids)))
    stop("`pouter` must be a named numeric vector.")
  if (is.null(fixed)) fixed <- numeric(0)
  if (length(fixed) > 0L && (is.null(names(fixed)) || any(!nzchar(names(fixed)))))
    stop("`fixed` must be a named numeric vector.")
  if (anyDuplicated(c(pouter_ids, names(fixed))))
    stop("Parameter ids overlap between `pouter` and `fixed`.")

  broadcast_scale <- function(s, ids) {
    if (length(s) == 1L && (is.null(names(s)) || !nzchar(names(s))))
      return(setNames(rep(unname(s), length(ids)), ids))
    miss <- setdiff(ids, names(s))
    if (length(miss))
      stop("`parameterScale` is missing entries for: ", paste(miss, collapse = ", "))
    s[ids]
  }
  scales_pouter <- broadcast_scale(parameterScale, pouter_ids)
  bad <- setdiff(unique(scales_pouter), c("lin", "log", "log10"))
  if (length(bad))
    stop("Unknown parameterScale(s): ", paste(bad, collapse = ", "))

  # Fixed parameters MUST be linear (PEtab+dMod convention; see
  # exportPEtabObject:1387 — fixed always written as parameterScale="lin").
  scales_fixed <- setNames(rep("lin", length(fixed)), names(fixed))
  scales_all   <- c(scales_pouter, scales_fixed)

  ## --- 6. bounds ----------------------------------------------------------
  if (is.null(lower)) lower <- setNames(rep(-Inf, length(pouter_ids)), pouter_ids)
  if (is.null(upper)) upper <- setNames(rep( Inf, length(pouter_ids)), pouter_ids)
  if (!setequal(names(lower), pouter_ids))
    stop("`lower` must be named like `pouter`.")
  if (!setequal(names(upper), pouter_ids))
    stop("`upper` must be named like `pouter`.")
  lower <- lower[pouter_ids]; upper <- upper[pouter_ids]

  ## --- 7. decompose trafo -------------------------------------------------
  decomp <- .petab_decompose_trafo(
    eqs          = eqs,
    states       = states,
    inner_pars   = inner_pars,
    obs_inner    = obs_inner,
    pouter_names = pouter_ids,
    fixed        = fixed,
    scales       = scales_all)

  fixed_full <- c(fixed, decomp$sbml_extra_pars)
  if (anyDuplicated(names(fixed_full)))
    stop("Internal error: duplicate fixed parameter id after decomposition: ",
         paste(names(fixed_full)[duplicated(names(fixed_full))], collapse = ", "))

  ## --- 8. observable metadata --------------------------------------------
  broadcast_obs <- function(x, ids, what, allowed) {
    if (length(x) == 1L && (is.null(names(x)) || !nzchar(names(x))))
      x <- setNames(rep(unname(x), length(ids)), ids)
    else {
      miss <- setdiff(ids, names(x))
      if (length(miss))
        stop(sprintf("`%s` is missing entries for: %s", what,
                     paste(miss, collapse = ", ")))
      x <- x[ids]
    }
    bad <- setdiff(unique(x), allowed)
    if (length(bad))
      stop(sprintf("Unknown %s value(s): %s", what, paste(bad, collapse = ", ")))
    x
  }
  obs_trafo  <- broadcast_obs(observableTransformation, obs_ids,
                              "observableTransformation",
                              c("lin", "log", "log10"))
  noise_dist <- broadcast_obs(noiseDistribution, obs_ids,
                              "noiseDistribution",
                              c("normal", "laplace", "log-normal"))
  obs_meta <- list(obs        = obs_eqns,
                   noise      = err_eqns,
                   obs_trafo  = obs_trafo,
                   noise_dist = noise_dist)

  ## --- 9. drive condition set from p (not from data); split sub-conds -----
  ## When per-row sigmas vary inside a single condition, encode them by
  ## splitting the data into one sub-condition per unique sigma value (the
  ## noiseParameters string). Sub-condition naming matches the importer's
  ## `<sim_cond>__<noi_hash>` convention so a roundtrip yields the same
  ## sub_cond_map shape on both sides.
  data_conds <- names(data)
  miss_in_data <- setdiff(conds, data_conds)
  miss_in_p    <- setdiff(data_conds, conds)
  if (length(miss_in_p))
    warning("data has condition(s) not covered by p (dropped): ",
            paste(miss_in_p, collapse = ", "))

  data_filtered <- list()
  scm_rows      <- list()

  for (c in intersect(conds, data_conds)) {
    d <- data[[c]]
    if (use_per_row_noise && "sigma" %in% colnames(d) && nrow(d) > 0L) {
      sig_str <- ifelse(is.na(d$sigma), "1", formatC(d$sigma, digits = 17,
                                                     format = "g"))
      noi_h   <- .petab_subcond_hash(sig_str)
      uniq_h  <- unique(noi_h)
      single  <- length(uniq_h) <= 1L
      for (h in uniq_h) {
        sel  <- which(noi_h == h)
        sub_name <- if (single || h == "") c else paste0(c, "__", h)
        d_sub <- d[sel, , drop = FALSE]
        d_sub$sigma <- NA_real_  # importer rebuilds from noiseParameters
        data_filtered[[sub_name]] <- d_sub
        scm_rows[[length(scm_rows) + 1L]] <- data.frame(
          simulationConditionId       = c,
          preequilibrationConditionId = "",
          observableParameters        = "",
          noiseParameters             = sig_str[sel[[1L]]],
          sub_condition               = sub_name,
          stringsAsFactors            = FALSE)
      }
    } else {
      data_filtered[[c]] <- d
      scm_rows[[length(scm_rows) + 1L]] <- data.frame(
        simulationConditionId       = c,
        preequilibrationConditionId = "",
        observableParameters        = "",
        noiseParameters             = "",
        sub_condition               = c,
        stringsAsFactors            = FALSE)
    }
  }
  for (c in miss_in_data) {
    data_filtered[[c]] <- data.frame(
      name = character(0), time = numeric(0),
      value = numeric(0), sigma = numeric(0),
      stringsAsFactors = FALSE)
    scm_rows[[length(scm_rows) + 1L]] <- data.frame(
      simulationConditionId       = c,
      preequilibrationConditionId = "",
      observableParameters        = "",
      noiseParameters             = "",
      sub_condition               = c,
      stringsAsFactors            = FALSE)
  }
  attr(data_filtered, "class") <- attr(data, "class")  # preserve datalist
  sub_cond_map <- do.call(rbind, scm_rows)
  rownames(sub_cond_map) <- NULL

  ## --- 10. tag pouter with scales attr (consumed by exportPEtabObject) ---
  attr(pouter, "petab_scales") <- scales_pouter

  petab <- list(
    pouter         = pouter,
    lower          = lower,
    upper          = upper,
    fixed          = fixed_full,
    obs_meta       = obs_meta,
    sub_cond_map   = sub_cond_map,
    condition.grid = decomp$conditions_df,
    data           = data_filtered,
    reactions      = reactions,
    inits          = decomp$inits,
    model_id       = model_id)

  exportPEtabObject(petab, dir, model_id = model_id,
                    amicipath = amicipath, overwrite = overwrite)
}


#' Export a dMod `petabProblem` back to PEtab v1
#'
#' Low-level exporter that writes the four PEtab tables (parameters,
#' observables, conditions, measurements) plus the SBML model and a YAML
#' manifest, given a fully-populated `petabProblem`-shaped list (i.e. an
#' object as produced by [importPEtab()]). Sub-conditions synthesised by the
#' importer are collapsed back to their PEtab condition + per-row
#' `observableParameters` / `noiseParameters` representation.
#'
#' If you have only the dMod-native pieces (datalist, odemodel, observation
#' / parameter / prediction functions, and a numeric `pouter`) and never went
#' through `importPEtab()`, use the higher-level [exportPEtab()] adapter,
#' which synthesises the missing PEtab metadata and dispatches here.
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
#' @param petab A `petabProblem` produced by [importPEtab()] (or hand-built
#'   list with the same slot names).
#' @param dir Output directory; created if missing.
#' @param model_id SBML model identifier (defaults to `petab$model_id` or
#'   `"dMod_export"`).
#' @param amicipath Forwarded to [export_sbml()].
#' @param overwrite Logical. If `FALSE` (default) errors when files already
#'   exist in `dir`.
#' @return Path to the written YAML manifest, invisibly.
#' @seealso [exportPEtab()] for a dMod-native entry point.
#' @export
#' @importFrom yaml write_yaml
exportPEtabObject <- function(petab, dir, model_id = NULL, amicipath = NULL,
                              overwrite = FALSE) {

  stopifnot(inherits(petab, "PEtabProblem") || is.list(petab))

  meta <- attr(petab, "petab_meta") %||% list()
  fixed        <- meta$fixed        %||% petab$fixed        %||% numeric(0)
  inits_meta   <- meta$inits        %||% petab$inits
  obs_meta     <- meta$obs_meta     %||% petab$obs_meta
  sub_cond_map <- meta$sub_cond_map %||% petab$sub_cond_map
  param_meta   <- meta$param_meta   %||% petab$param_meta
  cond_grid    <- attr(petab$dataList %||% petab$data, "condition.grid")
  if (is.null(cond_grid)) cond_grid <- petab$condition.grid
  bestfit      <- petab$bestfit %||% petab$pouter
  parlower     <- petab$parlower %||% petab$lower
  parupper     <- petab$parupper %||% petab$upper
  data_list    <- petab$dataList %||% petab$data
  reactions    <- petab$reactions
  if (is.null(model_id))
    model_id <- meta$model_id %||% petab$model_id %||% "dMod_export"

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

  # parameters.tsv: bestfit (estimate=1) + fixed (estimate=0)
  scales <- attr(bestfit, "petab_scales")
  if (is.null(scales))
    scales <- setNames(rep("lin", length(bestfit)), names(bestfit))

  # PEtab v1 spec: nominalValue / lowerBound / upperBound are written on the
  # linear scale. Internally dMod stores parameters on the parameter scale,
  # so we invert the importer's forward transform: log10 → 10^x, log → exp(x).
  apply_inv_scale <- function(values, ids) {
    sc <- scales[ids]
    out <- values
    log_idx   <- which(sc == "log")
    log10_idx <- which(sc == "log10")
    if (length(log_idx))   out[log_idx]   <- exp(values[log_idx])
    if (length(log10_idx)) out[log10_idx] <- 10 ^ values[log10_idx]
    out
  }

  est_ids <- names(bestfit)
  est_df <- data.frame(
    parameterId    = est_ids,
    parameterScale = unname(scales[est_ids]),
    lowerBound     = unname(apply_inv_scale(parlower[est_ids], est_ids)),
    upperBound     = unname(apply_inv_scale(parupper[est_ids], est_ids)),
    nominalValue   = unname(apply_inv_scale(bestfit, est_ids)),
    estimate       = 1L,
    stringsAsFactors = FALSE
  )
  fixed_df <- if (length(fixed)) data.frame(
    parameterId    = names(fixed),
    parameterScale = "lin",
    lowerBound     = -Inf,
    upperBound     = Inf,
    nominalValue   = unname(fixed),
    estimate       = 0L,
    stringsAsFactors = FALSE
  ) else NULL

  utils::write.table(rbind(est_df, fixed_df), paths$parameters,
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  # observables.tsv
  om <- obs_meta
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
  scm <- sub_cond_map
  uniq_sims <- unique(scm$simulationConditionId)
  cond_df <- cond_grid[uniq_sims, , drop = FALSE]
  cond_df$conditionId <- uniq_sims
  cond_df <- cond_df[, c("conditionId",
                         setdiff(colnames(cond_df), "conditionId")),
                     drop = FALSE]
  utils::write.table(cond_df, paths$conditions, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")

  # measurements.tsv: long-format with re-attached obs/noise param strings
  meas_rows <- do.call(rbind, lapply(names(data_list), function(sub) {
    df <- data_list[[sub]]
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
  all_pars <- c(bestfit, fixed)
  inits <- inits_meta %||%
           setNames(rep(0, length(reactions$states)), reactions$states)
  export_sbml(reactions, parameters = all_pars, inits = inits,
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
