## R/PEtabInterface.R
## ---------------------------------------------------------------------------
## PEtab v1/v2 importer / exporter, layered on top of dMod's existing
## SBML interface (R/SBMLinterface.R) and high-level APIs (Y, P, normL2, ...).
##
## Public API: importPEtab, exportPEtab, exportPEtabObject,
##             read_petab_yaml, read_petab_tables.
## Internal helpers prefixed `.petab_*` are unexported and may change shape.
##
## v2 strategy: the YAML reader dispatches on `format_version`. The v2 path
## reads the new schema (long conditions, experiments table, combined
## noiseDistribution, explicit observable/noise placeholders, no
## parameterScale) and translates it into the internal shapes that the
## v1-shaped `.petab_parse_*` helpers consume -- so the trafo / observation /
## error / objective builders are unchanged. The exporter's `format_version`
## argument symmetrically chooses the output shape (default "2.0.0").
##
## v1 spec: https://petab.readthedocs.io/en/latest/v1/documentation_data_format.html
## v2 spec: https://petab.readthedocs.io/en/latest/v2/documentation_data_format.html
## ---------------------------------------------------------------------------

# Internal: classify a YAML format_version string/number as the major version
# integer dMod cares about. v1 accepts 1 / "1" / "1.0.0"; v2 accepts 2 /
# "2.0.0" / "2.x.y"; everything else is an error.
.petab_major_version <- function(fv) {
  if (length(fv) == 0L) return(1L)  # v1 default for unversioned YAMLs
  s <- as.character(fv)[1L]
  m <- regmatches(s, regexec("^([1-9][0-9]*)", s))[[1L]]
  if (length(m) < 2L)
    stop("Unrecognised PEtab format_version: ", s)
  as.integer(m[2L])
}


## --- low-level YAML / TSV readers ------------------------------------------

#' Read a PEtab YAML manifest
#'
#' Parses a PEtab v1 or v2 YAML file and resolves the paths of the referenced
#' tables. No SBML or TSV is touched at this stage.
#'
#' Dispatch is on the YAML `format_version` key: values starting with "1"
#' use the v1 schema (`problems` list, `sbml_files` per problem); values
#' starting with "2" use the v2 schema (flat top-level `model_files`,
#' optional `experiment_files` and `mapping_files`).
#'
#' @param yamlPath Path to the PEtab YAML manifest.
#' @return A list with `baseDir`, `formatVersion` (integer major version),
#'   `parameterFile`, and a one-element `problems` list whose entry holds
#'   resolved (absolute) paths for `sbmlFile`, `conditionFile`,
#'   `measurementFile`, `observableFile`. v2 manifests additionally fill
#'   `experimentFile`, `mappingFile` (both possibly `NULL`), `modelID`
#'   (the first / canonical model id), and `models` (a named character
#'   vector mapping every model id from `model_files` to its resolved SBML
#'   path; length 1 in single-model problems). PEtab allows multiple
#'   problems per file; the reader supports a single problem entry but any
#'   number of v2 model_files inside it.
#' @export
#' @importFrom yaml read_yaml
read_petab_yaml <- function(yamlPath) {

  yamlPath <- normalizePath(yamlPath, mustWork = TRUE)
  baseDir  <- dirname(yamlPath)
  m <- yaml::read_yaml(yamlPath)

  major <- .petab_major_version(m$format_version)
  if (!major %in% c(1L, 2L))
    stop("PEtab format_version ", m$format_version,
         " not supported (only v1 and v2).")

  resolve <- function(p) normalizePath(file.path(baseDir, p), mustWork = TRUE)

  if (major == 1L) {
    if (length(m$problems) != 1L)
      stop("Only single-problem PEtab YAML is supported (got ",
           length(m$problems), ").")
    prob <- m$problems[[1]]
    pick_one <- function(x, slot) {
      if (length(x) == 0L) stop("YAML problem missing `", slot, "`.")
      if (length(x) > 1L) stop("Only one ", slot, " per problem is supported.")
      x[[1]]
    }
    sbml_path <- resolve(pick_one(prob$sbml_files, "sbml_files"))
    return(list(
      baseDir         = baseDir,
      formatVersion   = 1L,
      parameterFile   = resolve(m$parameter_file),
      problems = list(list(
        sbmlFile        = sbml_path,
        conditionFile   = resolve(pick_one(prob$condition_files, "condition_files")),
        measurementFile = resolve(pick_one(prob$measurement_files, "measurement_files")),
        observableFile  = resolve(pick_one(prob$observable_files, "observable_files")),
        experimentFile  = NULL,
        mappingFile     = NULL,
        modelID         = NA_character_,
        # `model` is a synthetic key for the v1 single-model case; v1 has no
        # modelId column on measurements, so this only ever appears as the
        # default fall-through in `importPEtab()`.
        models          = setNames(sbml_path, "model")
      ))
    ))
  }

  ## --- v2 -----------------------------------------------------------------
  pick_one_top <- function(x, slot, optional = FALSE) {
    if (length(x) == 0L) {
      if (optional) return(NULL)
      stop("YAML missing required `", slot, "`.")
    }
    if (length(x) > 1L)
      stop("Only one ", slot, " per problem is supported (got ",
           length(x), ").")
    if (is.list(x)) x[[1L]] else as.character(x)[1L]
  }

  param_file <- pick_one_top(m$parameter_files, "parameter_files")

  models <- m$model_files
  if (length(models) == 0L)
    stop("YAML missing required `model_files`.")
  if (is.null(names(models)) || any(!nzchar(names(models))))
    stop("`model_files` entries must be a named mapping `<modelId>: { location: ... }`.")

  resolved_models <- character(length(models))
  names(resolved_models) <- names(models)
  for (i in seq_along(models)) {
    entry <- models[[i]]
    mid   <- names(models)[i]
    if (is.null(entry$location))
      stop("model_files entry `", mid, "` missing `location`.")
    language <- tolower(as.character(entry$language %||% "sbml"))
    if (!identical(language, "sbml"))
      stop("Only SBML models are supported by dMod's PEtab importer (got `",
           language, "` for model id `", mid, "`).")
    resolved_models[i] <- resolve(entry$location)
  }

  obs_file  <- pick_one_top(m$observable_files,   "observable_files")
  meas_file <- pick_one_top(m$measurement_files,  "measurement_files")
  cond_file <- pick_one_top(m$condition_files,    "condition_files",  optional = TRUE)
  exp_file  <- pick_one_top(m$experiment_files,   "experiment_files", optional = TRUE)
  map_file  <- pick_one_top(m$mapping_files,      "mapping_files",    optional = TRUE)

  list(
    baseDir         = baseDir,
    formatVersion   = 2L,
    parameterFile   = resolve(param_file),
    problems = list(list(
      sbmlFile        = unname(resolved_models[1L]),
      conditionFile   = if (!is.null(cond_file)) resolve(cond_file) else NULL,
      measurementFile = resolve(meas_file),
      observableFile  = resolve(obs_file),
      experimentFile  = if (!is.null(exp_file)) resolve(exp_file) else NULL,
      mappingFile     = if (!is.null(map_file)) resolve(map_file) else NULL,
      modelID         = names(resolved_models)[1L],
      models          = resolved_models
    ))
  )
}


#' Read PEtab TSV tables (no SBML)
#'
#' Reads the PEtab tables referenced by a YAML manifest into base R data
#' frames. Useful for inspecting a problem without touching the SBML model.
#'
#' For v2 manifests, `experiments` and `mapping` slots are populated when
#' the corresponding files are present, otherwise `NULL`. The `conditions`
#' slot is also `NULL` when v2 declares no condition file.
#'
#' @param yamlPath Path to the PEtab YAML manifest.
#' @return A list with `parameters`, `conditions`, `measurements`,
#'   `observables`, optional `experiments` and `mapping` (data frames or
#'   `NULL`), `sbmlPath` (the first / canonical SBML path), `sbmlPaths`
#'   (named character vector keyed by modelId; length 1 for single-model
#'   problems), and `formatVersion` (integer).
#' @export
read_petab_tables <- function(yamlPath) {

  m  <- read_petab_yaml(yamlPath)
  pr <- m$problems[[1]]

  list(
    parameters    = .petab_read_tsv(m$parameterFile),
    conditions    = if (!is.null(pr$conditionFile))
                      .petab_read_tsv(pr$conditionFile) else NULL,
    measurements  = .petab_read_tsv(pr$measurementFile),
    observables   = .petab_read_tsv(pr$observableFile),
    experiments   = if (!is.null(pr$experimentFile))
                      .petab_read_tsv(pr$experimentFile) else NULL,
    mapping       = if (!is.null(pr$mappingFile))
                      .petab_read_tsv(pr$mappingFile) else NULL,
    sbmlPath      = pr$sbmlFile,
    sbmlPaths     = pr$models %||% setNames(pr$sbmlFile, pr$modelID),
    formatVersion = m$formatVersion
  )
}


.petab_read_tsv <- function(path) {
  utils::read.delim(path, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE,
                    check.names = FALSE,
                    na.strings = c("", "NA"),
                    strip.white = TRUE)
}


## --- v2 -> v1 normalizer ---------------------------------------------------
##
## Translates the four/six v2 table shapes into the v1-shapes that the
## `.petab_parse_*` helpers consume:
##   parameters: synthesise `parameterScale = "lin"`, coerce
##               `estimate` from logical/text to {1,0}.
##   observables: split combined `noiseDistribution` into v1's
##               `observableTransformation` + `noiseDistribution`; rewrite
##               named placeholders into v1 `<prefix>Parameter<k>_<obsId>`
##               sentinels.
##   conditions: pivot long (`conditionId`,`targetId`,`targetValue`) to
##               wide.
##   measurements: rewrite `experimentId` -> {`simulationConditionId`,
##               `preequilibrationConditionId`} via the experiments table.
## The mapping table (when present) is applied as a textual rewrite of
## `petabEntityId` -> `modelEntityId` across every string column we pass to
## the parsers, so the imported SBML's symbol names line up.
.petab_v2_normalize_tables <- function(tables) {

  ## --- 1. mapping table: build petab->model substitution -----------------
  mapping <- tables$mapping
  apply_mapping <- function(s) s
  if (!is.null(mapping) && nrow(mapping) > 0L) {
    if (!all(c("petabEntityId", "modelEntityId") %in% colnames(mapping)))
      stop("mapping.tsv must have columns `petabEntityId` and `modelEntityId`.")
    aliases <- mapping[!is.na(mapping$modelEntityId) &
                       nzchar(mapping$modelEntityId), , drop = FALSE]
    if (nrow(aliases) > 0L) {
      pat <- paste0("\\b", aliases$petabEntityId, "\\b")
      repl <- aliases$modelEntityId
      apply_mapping <- function(s) {
        if (length(s) == 0L) return(s)
        out <- s
        for (i in seq_along(pat))
          out <- gsub(pat[i], repl[i], out, perl = TRUE)
        out
      }
    }
  }

  ## --- 2. parameters: add parameterScale, coerce estimate ----------------
  par_df <- tables$parameters
  if (!"parameterId" %in% colnames(par_df))
    stop("parameters.tsv missing required column `parameterId`.")
  if (!"parameterScale" %in% colnames(par_df))
    par_df$parameterScale <- "lin"
  if ("estimate" %in% colnames(par_df)) {
    e <- par_df$estimate
    if (is.logical(e)) {
      par_df$estimate <- as.integer(e)
    } else {
      etxt <- tolower(trimws(as.character(e)))
      par_df$estimate <- ifelse(etxt %in% c("true", "1"), 1L,
                         ifelse(etxt %in% c("false", "0"), 0L, NA_integer_))
      if (any(is.na(par_df$estimate)))
        stop("parameters.tsv `estimate` must be true/false (or 1/0); got: ",
             paste(unique(e[is.na(par_df$estimate)]), collapse = ", "))
    }
  }

  ## --- 3. observables: split noiseDistribution, rewrite placeholders -----
  obs_df <- tables$observables
  if (!"observableId" %in% colnames(obs_df))
    stop("observables.tsv missing required column `observableId`.")
  if (!"observableFormula" %in% colnames(obs_df))
    stop("observables.tsv missing required column `observableFormula`.")
  # read.delim infers numeric type when a column parses entirely as numbers;
  # we substitute strings into these formulas, so coerce to character first.
  for (col in c("observableFormula", "noiseFormula",
                "observablePlaceholders", "noisePlaceholders")) {
    if (col %in% colnames(obs_df))
      obs_df[[col]] <- as.character(obs_df[[col]])
  }

  v2_dist <- if ("noiseDistribution" %in% colnames(obs_df))
               obs_df$noiseDistribution else rep("normal", nrow(obs_df))
  v2_dist[is.na(v2_dist) | !nzchar(v2_dist)] <- "normal"
  bad <- setdiff(unique(v2_dist),
                 c("normal", "log-normal", "laplace", "log-laplace"))
  if (length(bad))
    stop("Unsupported v2 noiseDistribution(s): ",
         paste(bad, collapse = ", "),
         ". dMod uses L2 likelihoods only; supported v2 values are ",
         "normal, log-normal, laplace, log-laplace.")
  obs_df$observableTransformation <- ifelse(
    v2_dist %in% c("log-normal", "log-laplace"), "log", "lin")
  obs_df$noiseDistribution <- ifelse(
    v2_dist %in% c("laplace", "log-laplace"), "laplace", "normal")

  rewrite_placeholders <- function(formula, ph_str, prefix, obs_id) {
    if (is.na(formula) || !nzchar(formula)) return(formula)
    if (is.na(ph_str) || !nzchar(ph_str))   return(formula)
    parts <- trimws(strsplit(ph_str, ";", fixed = TRUE)[[1L]])
    parts <- parts[nzchar(parts)]
    for (k in seq_along(parts)) {
      pat <- sprintf("\\b%s\\b", parts[k])
      formula <- gsub(pat, sprintf("%s%d_%s", prefix, k, obs_id),
                      formula, perl = TRUE)
    }
    formula
  }
  if ("observablePlaceholders" %in% colnames(obs_df)) {
    obs_df$observableFormula <- vapply(seq_len(nrow(obs_df)), function(i)
      rewrite_placeholders(obs_df$observableFormula[i],
                           obs_df$observablePlaceholders[i],
                           "observableParameter",
                           obs_df$observableId[i]), character(1))
  }
  if ("noiseFormula" %in% colnames(obs_df) &&
      "noisePlaceholders" %in% colnames(obs_df)) {
    obs_df$noiseFormula <- vapply(seq_len(nrow(obs_df)), function(i)
      rewrite_placeholders(obs_df$noiseFormula[i],
                           obs_df$noisePlaceholders[i],
                           "noiseParameter",
                           obs_df$observableId[i]), character(1))
  }
  obs_df$observableFormula <- apply_mapping(obs_df$observableFormula)
  if ("noiseFormula" %in% colnames(obs_df))
    obs_df$noiseFormula <- apply_mapping(obs_df$noiseFormula)

  ## --- 4. conditions: long -> wide --------------------------------------
  cond_v2 <- tables$conditions
  if (is.null(cond_v2) || nrow(cond_v2) == 0L) {
    cond_wide <- data.frame(conditionId = character(0),
                            stringsAsFactors = FALSE)
  } else {
    needed <- c("conditionId", "targetId", "targetValue")
    miss <- setdiff(needed, colnames(cond_v2))
    if (length(miss))
      stop("v2 conditions.tsv missing column(s): ",
           paste(miss, collapse = ", "))
    cond_v2$targetId    <- apply_mapping(as.character(cond_v2$targetId))
    cond_v2$targetValue <- apply_mapping(as.character(cond_v2$targetValue))
    uniq_cids <- unique(cond_v2$conditionId)
    uniq_tids <- unique(cond_v2$targetId)
    cond_wide <- data.frame(conditionId = uniq_cids,
                            stringsAsFactors = FALSE)
    for (tid in uniq_tids)
      cond_wide[[tid]] <- NA_character_
    for (i in seq_len(nrow(cond_v2))) {
      r <- match(cond_v2$conditionId[i], cond_wide$conditionId)
      tid <- cond_v2$targetId[i]
      if (!is.na(cond_wide[r, tid]) &&
          !identical(cond_wide[r, tid], cond_v2$targetValue[i]))
        stop("Duplicate (conditionId,targetId) in v2 conditions.tsv: ",
             cond_v2$conditionId[i], " / ", tid)
      cond_wide[r, tid] <- cond_v2$targetValue[i]
    }
  }

  ## --- 5. experiments -> (simCondId, preeqCondId) per measurement -------
  meas <- tables$measurements
  required_meas <- c("observableId", "time", "measurement")
  miss <- setdiff(required_meas, colnames(meas))
  if (length(miss))
    stop("measurements.tsv missing required column(s): ",
         paste(miss, collapse = ", "))

  # `modelId` is preserved for the multi-model dispatch in importPEtab().
  # Empty/NA cells are filled later (the importer knows the canonical default).

  exp_df <- tables$experiments
  exp_map <- list()  # experimentId -> list(sim=..., preeq=...)

  exp_ids_in_meas <- if ("experimentId" %in% colnames(meas))
                       unique(meas$experimentId[!is.na(meas$experimentId) &
                                                nzchar(meas$experimentId)])
                     else character(0)

  if (!is.null(exp_df) && nrow(exp_df) > 0L) {
    needed <- c("experimentId", "time", "conditionId")
    miss <- setdiff(needed, colnames(exp_df))
    if (length(miss))
      stop("experiments.tsv missing column(s): ",
           paste(miss, collapse = ", "))
    for (eid in unique(exp_df$experimentId)) {
      sub <- exp_df[exp_df$experimentId == eid, , drop = FALSE]
      times <- suppressWarnings(as.numeric(sub$time))
      # Accept "-inf"/"-Inf" parsing as -Inf
      txt <- tolower(trimws(as.character(sub$time)))
      times[txt %in% c("-inf", "-infinity")] <- -Inf
      times[txt %in% c("inf", "+inf", "infinity")] <- Inf
      ord <- order(times)
      times <- times[ord]
      cids  <- as.character(sub$conditionId)[ord]
      if (length(times) == 1L) {
        if (!(times[1L] == 0 || is.infinite(times[1L]) && times[1L] > 0))
          stop("Experiment `", eid, "` single-period time must be 0 or +inf, got ",
               times[1L], ".")
        exp_map[[eid]] <- list(sim = cids[1L], preeq = "")
      } else if (length(times) == 2L) {
        if (!is.infinite(times[1L]) || times[1L] > 0 || times[2L] != 0)
          stop("Experiment `", eid, "` two-period pattern must be ",
               "(time=-inf, condA) + (time=0, condB); got times (",
               paste(times, collapse = ", "), ").")
        exp_map[[eid]] <- list(sim = cids[2L], preeq = cids[1L])
      } else {
        stop("Experiment `", eid, "` has ", length(times),
             " periods; dMod supports at most 2 (preequilibration + ",
             "simulation).")
      }
    }
  }

  if ("experimentId" %in% colnames(meas)) {
    eid <- ifelse(is.na(meas$experimentId), "", meas$experimentId)
    sim   <- character(nrow(meas))
    preeq <- character(nrow(meas))
    for (i in seq_along(eid)) {
      e <- eid[i]
      if (!nzchar(e)) {
        # v2 spec: empty experimentId means "use model as-is" -- synthesise a
        # sentinel sim condition so the trafo stage has something to key on.
        sim[i]   <- "_dmod_default_condition_"
        preeq[i] <- ""
        next
      }
      if (is.null(exp_map[[e]]))
        stop("measurements.tsv references unknown experimentId `", e, "`.")
      sim[i]   <- exp_map[[e]]$sim
      preeq[i] <- exp_map[[e]]$preeq
    }
    meas$simulationConditionId        <- sim
    meas$preequilibrationConditionId  <- preeq
    meas$experimentId <- NULL

    # If we synthesised a default condition, append a no-override row to
    # cond_wide so .petab_parse_conditions sees it.
    if (any(sim == "_dmod_default_condition_") &&
        !"_dmod_default_condition_" %in% cond_wide$conditionId) {
      new_row <- cond_wide[NA_integer_, , drop = FALSE][1L, , drop = FALSE]
      new_row$conditionId <- "_dmod_default_condition_"
      cond_wide <- rbind(cond_wide, new_row)
    }
  } else {
    # No experimentId column at all; spec says empty means model-as-is.
    meas$simulationConditionId <- "_dmod_default_condition_"
    meas$preequilibrationConditionId <- ""
    if (!"_dmod_default_condition_" %in% cond_wide$conditionId) {
      new_row <- cond_wide[NA_integer_, , drop = FALSE][1L, , drop = FALSE]
      new_row$conditionId <- "_dmod_default_condition_"
      cond_wide <- rbind(cond_wide, new_row)
    }
  }

  # observableParameters / noiseParameters columns survive unchanged; their
  # placeholder substitution machinery already keys on the v1-style sentinels
  # we wrote into the observable/noise formulas above.
  for (col in c("observableParameters", "noiseParameters")) {
    if (col %in% colnames(meas))
      meas[[col]] <- ifelse(is.na(meas[[col]]), "", as.character(meas[[col]]))
  }

  list(parameters   = par_df,
       observables  = obs_df,
       conditions   = cond_wide,
       measurements = meas,
       sbmlPath      = tables$sbmlPath,
       sbmlPaths     = tables$sbmlPaths,
       formatVersion = 1L)
}


## --- per-table parsers (unit-testable, no SBML side effects) ---------------

# Internal: parameters.tsv -> dMod-shaped pieces.
# Returns:
#   pouter         named numeric, estimated parameters (nominalValue, scale-applied)
#   lower / upper  named numeric, bounds (scale-applied)
#   fixed          named numeric, non-estimated parameters (scale-applied)
#   scales         named character, "lin"/"log"/"log10" per parameterId
#                  (covers both pouter and fixed)
#   priors         NULL, or list(mu, sigma) with named numerics on the
#                  optimizer's view (i.e. on parameterScale). Built from
#                  `priorDistribution` / `priorParameters` (v2) or
#                  `objectivePriorType` / `objectivePriorParameters` (v1).
#                  dMod's only objective is normL2; only `parameterScaleNormal`
#                  (and `normal` on `lin` parameters, where the two coincide)
#                  is mapped to a `constraintL2` term. Other distributions raise
#                  a clear error.
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
  #   log10 -> log10(.),   log -> log(.),   lin -> identity.
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
  # constants on the *inner* (linear) scale -- no scale wrapping in the
  # trafo, so we keep linear values here regardless of parameterScale.
  fixed  <- setNames(to_num(df$nominalValue[fixed_idx]), fixed_ids)
  lower  <- setNames(apply_fwd_scale(to_num(df$lowerBound[pouter_idx]),
                                     pouter_ids),
                     pouter_ids)
  upper  <- setNames(apply_fwd_scale(to_num(df$upperBound[pouter_idx]),
                                     pouter_ids),
                     pouter_ids)

  priors <- .petab_parse_priors(df, scales)

  list(pouter = pouter, lower = lower, upper = upper,
       fixed = fixed, scales = scales, priors = priors)
}


# Internal: read PEtab prior columns and produce a named (mu, sigma) pair
# usable by `constraintL2`. Returns NULL when no prior is declared.
#
# Accepted column names (auto-detected, v2 takes priority over v1):
#   priorDistribution / priorParameters         (v2)
#   objectivePriorType / objectivePriorParameters (v1)
#
# Accepted distributions: `parameterScaleNormal` (always), and `normal`
# *only* on parameters whose `parameterScale = lin`. On non-lin parameters
# `normal` would impose the prior on the linear-scale value, not on the
# parameter as the optimizer sees it, and would silently misbehave.
# `uniform` is treated as "no prior". Anything else is rejected.
.petab_parse_priors <- function(df, scales) {

  pick <- function(a, b)
    if (a %in% colnames(df)) a else if (b %in% colnames(df)) b else NA_character_

  dist_col <- pick("priorDistribution",  "objectivePriorType")
  pars_col <- pick("priorParameters",    "objectivePriorParameters")
  if (is.na(dist_col) || is.na(pars_col)) return(NULL)

  pd <- as.character(df[[dist_col]])
  pp <- as.character(df[[pars_col]])
  est <- if ("estimate" %in% colnames(df)) df$estimate == 1L else rep(TRUE, nrow(df))
  active <- !is.na(pd) & nzchar(pd) & tolower(pd) != "uniform" & est
  if (!any(active)) return(NULL)

  mu <- numeric(0)
  sg <- numeric(0)
  for (i in which(active)) {
    pid <- df$parameterId[i]
    d   <- pd[i]
    sc  <- scales[[pid]]
    if (!d %in% c("normal", "parameterScaleNormal"))
      stop("Unsupported priorDistribution `", d, "` for parameter `", pid,
           "`. dMod only supports L2 priors (`parameterScaleNormal`, ",
           "or `normal` when parameterScale = lin).")
    if (d == "normal" && !identical(sc, "lin"))
      stop("priorDistribution `normal` for parameter `", pid,
           "` (parameterScale = `", sc, "`) is ambiguous: dMod's L2 prior ",
           "acts on the optimizer's parameter (the parameterScale-",
           "transformed value). Use `parameterScaleNormal` to declare a ",
           "prior on that view, or set parameterScale = `lin`.")
    parts <- trimws(strsplit(pp[i], ";", fixed = TRUE)[[1L]])
    if (length(parts) != 2L)
      stop("priorParameters for `", pid, "` must be two `;`-separated ",
           "numbers (`mu;sigma`); got `", pp[i], "`.")
    mu_i <- suppressWarnings(as.numeric(parts[1L]))
    sg_i <- suppressWarnings(as.numeric(parts[2L]))
    if (is.na(mu_i) || is.na(sg_i) || sg_i <= 0)
      stop("priorParameters for `", pid, "` must parse as `mu;sigma` with ",
           "positive sigma; got `", pp[i], "`.")
    mu[pid] <- mu_i
    sg[pid] <- sg_i
  }
  if (length(mu) == 0L) return(NULL)
  list(mu = mu, sigma = sg)
}


# Internal: observables.tsv -> dMod-shaped pieces.
# Returns:
#   obs       named character, observableFormula keyed by observableId.
#             *Not* yet substituted: observableParameterK_<id> / noiseParameterK_<id>
#             placeholders are left intact for per-row substitution.
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
#   "init"        -- a species/state initial-value override (column name matches state)
#   "compartment" -- a compartment volume override
#   "parameter"   -- a parameter override (the catch-all, includes condition-only pars)
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


# Internal: build a readable, name-safe label for sub-conditions whose
# observable/noise parameter strings differ within a single PEtab condition.
# The raw parameter string (e.g. "sd_pSTAT5A_rel") is preferred over an
# opaque hash so users see what's actually different between sub-conditions.
# Punctuation collapses to "_"; over-long labels get a short md5 tail to
# stay file-system-friendly.
.petab_subcond_label <- function(s) {
  if (length(s) == 0L || all(is.na(s) | s == "")) return(rep("", length(s)))
  vapply(s, function(x) {
    if (is.na(x) || identical(x, "")) return("")
    lab <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
    lab <- gsub("^_+|_+$", "", lab)
    if (!nzchar(lab))
      lab <- substr(digest::digest(x, algo = "md5"), 1L, 8L)
    if (nchar(lab) > 32L) {
      tag <- substr(digest::digest(x, algo = "md5"), 1L, 6L)
      lab <- paste0(substr(lab, 1L, 24L), "_", tag)
    }
    lab
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

  # read.delim infers numeric type when a column is all numeric, but PEtab
  # observable/noise parameters can be either numeric or symbol strings.
  # Cast to character so the substitution step has a uniform input type.
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
  # symbols remain (case 0015's "noise"), leave sigma = NA -- the err model
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

  # Sub-condition assignment.
  #
  # PEtab `observableParameters` / `noiseParameters` columns substitute
  # placeholders `<prefix>K_<observableId>` -- i.e. they are *observable-
  # specific*. Two rows with different observableIds substitute disjoint
  # placeholder sets and can therefore coexist in a single trafo without
  # interference (e.g. Boehm: `noiseParameter1_pSTAT5A_rel = sd_pSTAT5A_rel`,
  # `noiseParameter1_pSTAT5B_rel = sd_pSTAT5B_rel`, ... in one sub-condition).
  #
  # We merge whenever every (simCondId, preeq) group is observable-consistent:
  # for each observableId X in the group, all rows with `observableId == X`
  # share the same (obsParStr, noiseParStr) tuple. Result: one sub-condition
  # per (simCondId, preeq) pair, carrying a per-observable substitution map.
  #
  # If consistency fails (same obsId carries different tuples within the
  # same group -- e.g. replicate-specific scaling, rare), we fall back to
  # the per-tuple split: one sub-condition per unique (peq, obsParStr,
  # noiseParStr) tuple, with a wildcard `"*"` substitution applied to all
  # placeholders.
  obs_lab <- .petab_subcond_label(m$observableParameters)
  noi_lab <- .petab_subcond_label(m$noiseParameters)
  peq_lab <- .petab_subcond_label(m$preequilibrationConditionId)

  m$sub_condition  <- m$simulationConditionId  # default; overwritten below
  obs_subs_by_sub  <- list()
  noi_subs_by_sub  <- list()

  # Group by (simCondId, preeq); decide merge-vs-split per group.
  group_key <- paste(m$simulationConditionId, m$preequilibrationConditionId,
                     sep = "")
  for (gk in unique(group_key)) {
    ix <- which(group_key == gk)
    sc   <- m$simulationConditionId[ix[1L]]
    peq  <- m$preequilibrationConditionId[ix[1L]]
    plab <- peq_lab[ix[1L]]

    # Per observableId in this group: distinct (obsPar, noisePar) tuples.
    by_obs <- split(ix, m$observableId[ix])
    consistent <- all(vapply(by_obs, function(rs) {
      length(unique(paste(m$observableParameters[rs],
                          m$noiseParameters[rs],
                          sep = ""))) <= 1L
    }, logical(1)))

    if (consistent) {
      # ----- merge: one sub-condition for the whole (sc, peq) group -------
      sub_name <- if (nzchar(plab)) paste0(sc, "__", plab) else sc
      m$sub_condition[ix] <- sub_name
      obs_map <- vapply(by_obs, function(rs) m$observableParameters[rs[1L]],
                        character(1))
      noi_map <- vapply(by_obs, function(rs) m$noiseParameters[rs[1L]],
                        character(1))
      obs_map <- obs_map[nzchar(obs_map)]
      noi_map <- noi_map[nzchar(noi_map)]
      if (length(obs_map)) obs_subs_by_sub[[sub_name]] <- obs_map
      if (length(noi_map)) noi_subs_by_sub[[sub_name]] <- noi_map
      next
    }

    # ----- split: per-tuple fallback (rare; same obsId carries multiple
    # distinct (obsPar, noisePar) within one (sc, peq) group) ------------
    tuple_key <- paste(peq_lab[ix], obs_lab[ix], noi_lab[ix], sep = "")
    keys_here <- unique(tuple_key)
    suffixes <- vapply(keys_here, function(k) {
      jx <- ix[tuple_key == k][1L]
      parts <- c(if (nzchar(peq_lab[jx])) peq_lab[jx],
                 if (nzchar(obs_lab[jx])) obs_lab[jx],
                 if (nzchar(noi_lab[jx])) noi_lab[jx])
      paste(parts, collapse = "_")
    }, character(1))
    if (anyDuplicated(suffixes)) {
      tags <- vapply(keys_here, function(k)
        substr(digest::digest(k, algo = "md5"), 1L, 6L), character(1))
      suffixes <- paste(suffixes, tags, sep = "_")
    }
    for (i in seq_along(keys_here)) {
      jx <- ix[tuple_key == keys_here[i]]
      sub_name <- paste0(sc, "__", suffixes[i])
      m$sub_condition[jx] <- sub_name
      o <- m$observableParameters[jx[1L]]
      n <- m$noiseParameters[jx[1L]]
      if (nzchar(o)) obs_subs_by_sub[[sub_name]] <- c("*" = o)
      if (nzchar(n)) noi_subs_by_sub[[sub_name]] <- c("*" = n)
    }
  }

  # Build the sub_cond_map (one row per unique sub_cond). The list-valued
  # substitution maps live in attributes so the data.frame stays trivially
  # serialisable; consumers (build_trafo, build_error_fn, exporter) read
  # them via attr(, "obs_subs") / attr(, "noi_subs").
  sub_cond_map <- unique(m[, c("simulationConditionId",
                               "preequilibrationConditionId",
                               "sub_condition")])
  rownames(sub_cond_map) <- NULL
  if (anyDuplicated(sub_cond_map$sub_condition))
    stop("Internal: sub-condition map has duplicate sub_condition rows.")
  attr(sub_cond_map, "obs_subs") <- obs_subs_by_sub
  attr(sub_cond_map, "noi_subs") <- noi_subs_by_sub

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
#   reactions      eqnlist (used for Pequil composition if pre-eq present)
#
# Returns a parfn that maps outer pars (estimated + fixed) to the inner-side
# parameter set (states' initial values, inner_pars, obs_inner).
.petab_build_trafo <- function(sub_cond_map, conditions, col_kind, override_cols,
                               inits, sbml_pars, states, inner_pars,
                               pouter_names, fixed, scales,
                               obs_inner = character(),
                               reactions = NULL,
                               compile = TRUE,
                               modelname = "petab_trafo") {

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
    #         elsewhere -- for the Pequil pre/post stages where parameter
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

  apply_petab_param_subs <- function(tr, obs_subs, noi_subs) {
    tr <- .petab_apply_subs_map(tr, obs_subs, prefix = "observableParameter")
    tr <- .petab_apply_subs_map(tr, noi_subs, prefix = "noiseParameter")
    tr
  }

  obs_subs_attr <- attr(sub_cond_map, "obs_subs") %||% list()
  noi_subs_attr <- attr(sub_cond_map, "noi_subs") %||% list()

  parfns <- list()

  for (sub in sub_cond_map$sub_condition) {

    row_idx <- which(sub_cond_map$sub_condition == sub)
    sim     <- sub_cond_map$simulationConditionId[row_idx]
    peq     <- sub_cond_map$preequilibrationConditionId[row_idx]
    obs_subs <- obs_subs_attr[[sub]]
    noi_subs <- noi_subs_attr[[sub]]

    has_peq <- length(peq) && nzchar(peq)

    if (!has_peq) {
      # ----- single-stage Pexpl path (Stage-1 path) ---------------------
      tr <- build_default()
      tr <- apply_row_overrides(tr, sim, "all")
      tr <- apply_petab_param_subs(tr, obs_subs, noi_subs)
      tr <- apply_scale_chain_rule(tr)
      parfns[[sub]] <- P(tr, condition = sub, compile = compile,
                         modelname = paste(modelname,
                                           sanitizeConditions(sub),
                                           sep = "_"))
      next
    }

    # ----- Pequil composition --------------------------------------------
    # 1. Substitute peq-row's parameter overrides directly into the reaction
    #    rates -> peq_reactions. Parameter-scale chain rule must be applied
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

    # 2. p_pre: outer -> c(state inits with peq state overrides applied,
    #    inner pars identity, observable/noise placeholders pre-substituted).
    tr_pre <- build_default()
    for (st in names(peq_state_subs)) tr_pre[st] <- peq_state_subs[[st]]
    tr_pre <- apply_petab_param_subs(tr_pre, obs_subs, noi_subs)
    tr_pre <- apply_scale_chain_rule(tr_pre)

    # Pre-equilibration must integrate the FULL network to steady state, with
    # every conserved moiety fixed by its initial total. When the network has
    # conserved moieties we therefore run Pequil with expressInTotals = TRUE
    # (which keeps all moiety species dynamical and pins each `total_i`) and
    # supply each total from the scale-corrected species initials. The
    # default expressInTotals = FALSE would instead hold a pivot species at
    # its initial value, violating mass conservation.
    peq_totals <- getTotals(peq_reactions)
    eq_in_totals <- length(peq_totals) > 0L &&
      all(vapply(peq_totals, function(e)
        all(getSymbols(e) %in% names(tr_pre)), logical(1)))
    if (eq_in_totals) {
      for (tn in names(peq_totals)) {
        sp <- getSymbols(peq_totals[[tn]])
        tr_pre[tn] <- replaceSymbols(sp, paste0("(", tr_pre[sp], ")"),
                                     peq_totals[[tn]])
      }
    }

    p_pre <- P(tr_pre, condition = sub, compile = compile,
               modelname = paste(modelname, sanitizeConditions(sub),
                                 "pre", sep = "_"))

    # 3. p_eq: integrate peq_reactions to steady state. attach.input ensures
    #    parameters and observable/noise placeholders flow through unchanged.
    p_eq <- Pequil(peq_reactions, condition = sub, attach.input = TRUE,
                   expressInTotals = eq_in_totals,
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
# `obs_id`: when non-NULL, only placeholders whose observableId suffix equals
# `obs_id` are substituted -- this lets a single sub-condition trafo carry
# per-observable PEtab placeholder substitutions (e.g. Boehm's three
# `noiseParameter1_<obsId>` slots). NULL substitutes every matching placeholder
# regardless of obsId.
.petab_apply_param_substitution <- function(tr, repls, prefix, obs_id = NULL) {
  parts <- trimws(strsplit(repls, ";", fixed = TRUE)[[1]])
  inner <- names(tr)
  pat <- if (is.null(obs_id))
           paste0("^", prefix, "([0-9]+)_(.+)$")
         else
           paste0("^", prefix, "([0-9]+)_",
                  gsub("([][\\\\.|()^$*+?{}])", "\\\\\\1", obs_id), "$")
  for (sym in inner) {
    m <- regmatches(sym, regexec(pat, sym))[[1]]
    if (length(m) >= 2L) {
      k <- as.integer(m[2])
      if (k >= 1L && k <= length(parts)) tr[sym] <- parts[k]
    }
  }
  tr
}


# Internal: apply a per-observable substitution map to a trafo. `subs` is a
# named character vector keyed by observableId (or by "*" meaning "apply to
# all matching placeholders regardless of observableId" -- the per-tuple
# fallback semantic). Each value is a ";"-separated PEtab replacement string.
.petab_apply_subs_map <- function(tr, subs, prefix) {
  if (length(subs) == 0L) return(tr)
  for (key in names(subs)) {
    repl <- subs[[key]]
    if (!nzchar(repl)) next
    tr <- .petab_apply_param_substitution(
      tr, repl, prefix,
      obs_id = if (identical(key, "*")) NULL else key)
  }
  tr
}


# Build the observation function `g`. Observables with `observableTransformation`
# of `log` or `log10` get wrapped at construction time (`obs_b -> log10(B)`),
# which keeps dMod's normL2 fast path on PEtab `{log,log10} * normal` cases:
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


# Build the error model when at least one (sub-condition * observable) noise
# formula carries free symbols after PEtab parameter substitution.
#
# Strategy: build a single, condition-unspecific Y over the *original*
# noiseFormula in PEtab placeholder form (e.g. `noiseParameter1_obs_a`) plus
# the observable formulas as f-states. The placeholders are inner_targets
# of the trafo (.petab_build_trafo adds them to obs_inner) and the trafo's
# per-sub-condition `.petab_apply_param_substitution` rewrites them to the
# concrete row value (numeric literal or outer-parameter symbol). At runtime
# the post-trafo inner parameter vector pinner -- which normL2 hands to the
# err model -- therefore already carries the substituted value under the
# placeholder name, so a single Y suffices.
#
# Returns NULL when every (sub_cond, observable) noise formula evaluates to
# a constant after substitution -- sigma is then carried by the data column
# and normL2's fast path handles the likelihood without an error model.
.petab_build_error_fn <- function(obs_meta, sub_cond_map, reactions,
                                  compile = TRUE, modelname = "petab_err") {

  obs_subs <- attr(sub_cond_map, "obs_subs") %||% list()
  noi_subs <- attr(sub_cond_map, "noi_subs") %||% list()

  # For each (sub-condition, observableId), evaluate the noise formula after
  # the per-row substitutions. If anything stays symbolic, we need a Y-based
  # error model; otherwise sigma is carried by the data column and normL2's
  # fast path handles it.
  pick_str <- function(map, sub, obsId) {
    m <- map[[sub]]
    if (is.null(m)) return("")
    if (obsId %in% names(m)) return(unname(m[obsId]))
    if ("*" %in% names(m))   return(unname(m["*"]))
    ""
  }
  any_symbolic <- any(vapply(sub_cond_map$sub_condition, function(sub) {
    any(vapply(names(obs_meta$noise), function(obsId) {
      f <- obs_meta$noise[[obsId]]
      f <- .petab_substitute_param_string(
        f, pick_str(obs_subs, sub, obsId), "observableParameter")
      f <- .petab_substitute_param_string(
        f, pick_str(noi_subs, sub, obsId), "noiseParameter")
      is.na(.petab_eval_constant(f))
    }, logical(1)))
  }, logical(1)))
  if (!any_symbolic) return(NULL)

  # f for err Y: ODE states + observable formulas (so a noise formula could
  # reference an observable, e.g. relative noise `sigma * obs_a`). We pass the
  # *untransformed* observable formulas -- the obs trafo (log/log10) doesn't
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
# the importer can defer linking until a single batched compile(). `events`
# (an eventlist or NULL) is what `import_sbml()` reads from <event> blocks.
.petab_build_odemodel <- function(reactions, solver,
                                  modelname = "petab_model",
                                  compile = TRUE,
                                  events = NULL) {
  m <- odemodel(reactions, modelname = modelname, solver = solver,
                events = events, compile = compile)
  list(odemodel = m, x = Xs(m))
}


# Internal: take a NULL-condition obsfn (e.g. a Y() result) and rewire its
# `mappings` / `conditions` attributes so it carries explicit `conds`, all
# routing to the same X2Y kernel. Used by the multi-model importer to
# combine per-model error functions via `+.fn`, which requires named
# mappings on both sides.
.obsfn_with_conditions <- function(fn, conds) {
  if (is.null(fn) || !length(conds)) return(NULL)
  if (!inherits(fn, "obsfn")) stop("not an obsfn")
  base <- attr(fn, "mappings")[[1L]]
  attr(fn, "mappings")   <- setNames(rep(list(base), length(conds)), conds)
  attr(fn, "conditions") <- conds
  fn
}


# Internal: shared "build per-model dMod pieces" routine. Used by the
# single- and multi-model branches of importPEtab. `meas_m` is the slice
# of the (already-normalised) measurements table belonging to one model;
# `obs_meta_full` is the global observables parse, which we filter down
# to the observables that actually appear in `meas_m` so the model's `g`
# / `e` only compile what they need.
.petab_build_model_pieces <- function(sbml, meas_m, conditions_df,
                                      obs_meta_full, param_meta,
                                      modelname, solver, compile,
                                      sub_cond_prefix = "") {

  obs_used <- unique(meas_m$observableId)
  obs_meta <- list(
    obs        = obs_meta_full$obs[obs_used],
    obs_trafo  = obs_meta_full$obs_trafo[obs_used],
    noise      = obs_meta_full$noise[obs_used],
    noise_dist = obs_meta_full$noise_dist[obs_used]
  )

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

  cond_info <- .petab_parse_conditions(
                 conditions_df,
                 sbml_states       = states,
                 sbml_compartments = names(sbml$reactions$compartments %||% list()),
                 sbml_pars         = names(sbml$pars),
                 obs_inner         = obs_inner)
  meas_info <- .petab_parse_measurements(meas_m, obs_meta)

  # In multi-model setups, sub-condition keys built from PEtab columns
  # (simulationConditionId, observableParameters, noiseParameters) are not
  # guaranteed to be disjoint across models; prefixing with the modelId
  # prevents collisions when per-model trafos and datalists are summed.
  if (nzchar(sub_cond_prefix)) {
    rn <- function(x) paste0(sub_cond_prefix, x)
    meas_info$sub_cond_map$sub_condition <- rn(meas_info$sub_cond_map$sub_condition)
    meas_info$data$condition             <- rn(meas_info$data$condition)
    rename_attr_keys <- function(scm, key) {
      a <- attr(scm, key)
      if (length(a)) names(a) <- rn(names(a))
      attr(scm, key) <- a
      scm
    }
    meas_info$sub_cond_map <- rename_attr_keys(meas_info$sub_cond_map, "obs_subs")
    meas_info$sub_cond_map <- rename_attr_keys(meas_info$sub_cond_map, "noi_subs")
  }

  pouter_names <- names(param_meta$pouter)
  fixed        <- param_meta$fixed

  # SBML-default-as-fixed gap (see importPEtab() comment).
  init_syms <- unique(unlist(lapply(sbml$inits, function(e) {
    if (is.character(e)) getSymbols(e) else character(0)
  })))
  required <- unique(c(inner_pars, obs_inner, init_syms))
  required <- setdiff(required, c(states, "time", pouter_names, names(fixed)))
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

  ode_pair <- .petab_build_odemodel(sbml$reactions, solver = solver,
                                    modelname = paste0(modelname, "_ode"),
                                    compile = compile,
                                    events = sbml$events)
  g <- .petab_build_observation_fn(obs_meta$obs, obs_meta$obs_trafo,
                                   sbml$reactions,
                                   compile = compile,
                                   modelname = paste0(modelname, "_obs"))
  e <- .petab_build_error_fn(obs_meta, meas_info$sub_cond_map,
                             sbml$reactions, compile = compile,
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
          compile       = compile,
          modelname     = paste0(modelname, "_trafo"))

  dataList <- as.datalist(meas_info$data, split.by = "condition")

  list(
    odemodel     = ode_pair$odemodel,
    x            = ode_pair$x,
    g            = g, e = e, p = p,
    dataList     = dataList,
    fixed        = fixed,
    obs_meta     = obs_meta,
    sub_cond_map = meas_info$sub_cond_map,
    cond_grid    = cond_info$grid,
    col_kind     = cond_info$col_kind,
    sbml         = sbml,
    sub_conds    = unique(meas_info$sub_cond_map$sub_condition)
  )
}


# Build the objective.
#
# Supports `{lin, log, log10} * normal` noise: log/log10 enter via the
# pre-wrapped observation function (data values are matched-side
# transformed), so normL2 still gives the correct residual. Symbolic sigmas
# flow through the error model.
#
# Non-normal distributions (laplace, log-normal) are not implemented; they
# would need a per-cell residual objective.
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
  # can be <= 0 -> sqrt of negative -> NaN).
  base_obj <- normL2(data = data, x = prd, errmodel = errmodel,
                     use.bessel = FALSE, cores = 1L)

  # Data-coordinate Jacobian for log / log10 observable transformations.
  # PEtab's likelihood is on the linear y_obs:
  #   log L = log f_lin(y_obs) = log f_trafo(trafo(y_obs)) - log|d trafo / dy_obs|
  # so for log:    -2 log L includes  2 * Sigma log(y_obs)
  # for log10:    -2 log L includes  2 * Sigma log(y_obs * ln 10)
  # The offset depends only on the data, not on parameters, so gradients
  # and Hessians are unchanged -- only the absolute objective value shifts
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

#' Import a PEtab v1 or v2 problem into dMod
#'
#' Reads a PEtab YAML manifest plus the associated SBML model and TSV tables,
#' and assembles a fully-composed dMod problem: prediction function,
#' observation function, parameter transformation, datalist, and objective.
#' The SBML side is delegated to [import_sbml()] (libsbml-based).
#'
#' Dispatch is on `format_version`. v2 manifests are translated to the
#' internal v1 shapes so the trafo / observation / objective pipeline is
#' shared.
#'
#' v2 specifics handled here: long-format conditions are pivoted to wide,
#' the `experiments.tsv` `(time, conditionId)` sequences are mapped to
#' dMod's `(simulationConditionId, preequilibrationConditionId)` pair (max
#' two periods), the combined `noiseDistribution` (`normal`/`log-normal`/
#' `laplace`/`log-laplace`) is split into observable transformation +
#' distribution, explicit `observablePlaceholders` / `noisePlaceholders`
#' are rewritten to v1's `observableParameter${k}_${id}` sentinels, and an
#' optional `mapping.tsv` is applied as a textual rewrite of PEtab entity
#' IDs into model entity IDs.
#'
#' L2 priors on parameters are read from `priorDistribution` /
#' `priorParameters` (v2) or `objectivePriorType` /
#' `objectivePriorParameters` (v1) and added to `obj` as a `constraintL2`
#' term. `parameterScaleNormal` is mapped directly; `normal` is mapped on
#' parameters with `parameterScale = lin` (where the two coincide). All
#' other prior distributions raise an error: dMod's only objective is
#' `normL2`.
#'
#' Out of scope: non-SBML model languages (BNGL / CellML / PySB), and
#' multi-period experiments (>2 periods).
#'
#' @param yamlPath Path to the PEtab YAML manifest.
#' @param solver Required: one of `"deSolve"` or `"CppODE"`. Forwarded to
#'   [odemodel()].
#' @param compile Logical. If `TRUE` (default) the generated trafo,
#'   observation function, and ODE model are compiled to native code. Set to
#'   `FALSE` for inspection-only use.
#' @param cores Number of parallel compilation jobs (Unix only) forwarded to
#'   [compile()]. The importer batches all generated source files into a
#'   single `compile()` call so this directly controls native-build
#'   concurrency.
#' @param modelname Optional base modelname for the generated native files.
#'   Defaults to the YAML basename.
#' @return A list with class `"PEtabProblem"` holding `dataList`,
#'   `reactions`, `odemodel`, `g`, `x`, `p`, `e`, `prd` (the composite
#'   `g * x * p`), `obj`, `bestfit`, `parlower`, `parupper`. The `obj`
#'   closure has the PEtab fixed parameters baked in, so calling
#'   `obj(bestfit)` evaluates the likelihood at the published MLE
#'   without the user having to pass `fixed = ...`. `bestfit` is the
#'   PEtab `nominalValue` for every estimated parameter -- for benchmark
#'   problems this is the published maximum-likelihood estimate
#'   transported through the manifest; for user-authored problems it is
#'   the current best estimate (and serves as the optimizer's starting
#'   point). The vector carries `attr(., "petab_scales")` recording each
#'   parameter's PEtab scale (the exporter reads it back). Internal
#'   metadata used by the exporter (`fixed`, `inits`, `modelID`,
#'   `source_yaml`, `sub_cond_map`, `obs_meta`, `param_meta`) lives on
#'   `attr(., "petab_meta")`.
#' @export
#' @example inst/examples/PEtabInterface.R
importPEtab <- function(yamlPath, solver,
                        compile = TRUE, cores = 1L, modelname = NULL) {

  cores <- as.integer(cores)
  if (length(cores) != 1L || is.na(cores) || cores < 1L)
    stop("`cores` must be a single positive integer.")

  if (missing(solver))
    stop("Argument `solver` is required (one of \"deSolve\", \"CppODE\").")
  solver <- match.arg(solver, c("deSolve", "CppODE"))

  yamlPath <- normalizePath(yamlPath, mustWork = TRUE)
  if (is.null(modelname))
    modelname <- sub("\\.ya?ml$", "", basename(yamlPath), ignore.case = TRUE)
  modelname <- gsub("[^A-Za-z0-9_]", "_", modelname)

  tables <- read_petab_tables(yamlPath)
  if (identical(tables$formatVersion, 2L))
    tables <- .petab_v2_normalize_tables(tables)

  sbml_paths <- tables$sbmlPaths
  if (is.null(sbml_paths) || !length(sbml_paths))
    sbml_paths <- setNames(tables$sbmlPath, "model")
  default_mid <- names(sbml_paths)[1L]

  meas <- tables$measurements
  if (!"modelId" %in% colnames(meas)) {
    if (length(sbml_paths) > 1L)
      stop("measurements.tsv has no `modelId` column but the YAML declares ",
           length(sbml_paths), " models. Add a `modelId` column referencing ",
           "the appropriate model_files key per row.")
    meas$modelId <- default_mid
  } else {
    nas <- is.na(meas$modelId) | !nzchar(meas$modelId)
    if (any(nas)) {
      if (length(sbml_paths) > 1L)
        stop("measurements.tsv has missing `modelId` values; in a multi-model ",
             "problem every measurement row must specify its modelId.")
      meas$modelId[nas] <- default_mid
    }
    bad <- !meas$modelId %in% names(sbml_paths)
    if (any(bad))
      stop("measurements.tsv references unknown modelId(s) not declared in ",
           "the YAML's model_files: ",
           paste(unique(meas$modelId[bad]), collapse = ", "))
  }
  tables$measurements <- meas

  param_meta <- .petab_parse_parameters(tables$parameters)
  obs_meta   <- .petab_parse_observables(tables$observables)

  # Compilation of the generated source files is deferred to a single batched
  # compile() at the end so `cores` controls native-build concurrency and we
  # avoid one g++ invocation per sub-cond.
  per_model <- list()
  for (mid in names(sbml_paths)) {
    meas_m <- meas[meas$modelId == mid, , drop = FALSE]
    if (nrow(meas_m) == 0L) {
      warning("Model `", mid, "` has no measurements; skipping.",
              call. = FALSE)
      next
    }
    meas_m$modelId <- NULL
    sbml_m <- import_sbml(unname(sbml_paths[mid]))
    suffix <- if (length(sbml_paths) == 1L) ""
              else paste0("__", gsub("[^A-Za-z0-9_]", "_", mid))
    sc_prefix <- if (length(sbml_paths) == 1L) ""
                 else paste0(gsub("[^A-Za-z0-9_]", "_", mid), "__")
    pieces <- .petab_build_model_pieces(
      sbml             = sbml_m,
      meas_m           = meas_m,
      conditions_df    = tables$conditions,
      obs_meta_full    = obs_meta,
      param_meta       = param_meta,
      modelname        = paste0(modelname, suffix),
      solver           = solver,
      compile          = FALSE,
      sub_cond_prefix  = sc_prefix)
    pieces$modelID <- mid
    per_model[[mid]] <- pieces
  }
  if (length(per_model) == 0L)
    stop("No measurements matched any declared model.")

  # For multi-model problems each `prd_M = g_M * x_M * p_M` carries `p_M`'s
  # sub-conds as conditions, so `Reduce(\`+\`, prd_list)` is disjoint.
  multi_model <- length(per_model) > 1L

  if (!multi_model) {
    pm  <- per_model[[1L]]
    odeobj   <- pm$odemodel
    x        <- pm$x
    g        <- pm$g
    e        <- pm$e
    p        <- pm$p
    dataList <- pm$dataList
    fixed    <- pm$fixed
    sbml     <- pm$sbml
    sub_cond_map <- pm$sub_cond_map
    cond_grid    <- pm$cond_grid
    col_kind     <- pm$col_kind
    obs_meta_used <- pm$obs_meta
    prd <- g * x * p
  } else {
    # Compose per-model prd's, then sum via +.fn (mappings disjoint by
    # construction). Each prd_M's conditions come from p_M (multi-cond),
    # so the result has all sub-conds across all models.
    prd_list <- lapply(per_model, function(pm) pm$g * pm$x * pm$p)
    prd <- Reduce(`+`, prd_list)

    # Error model: each per-model `e_M` is a NULL-condition Y (one Y per
    # model). To combine via `+.fn` we attach explicit conditions = M's
    # sub-conds; +.fn then dispatches each sub-cond to the right kernel.
    e_pieces <- Filter(Negate(is.null),
      lapply(per_model, function(pm)
        if (is.null(pm$e)) NULL
        else .obsfn_with_conditions(pm$e, pm$sub_conds)))
    e <- if (length(e_pieces)) Reduce(`+`, e_pieces) else NULL

    # Combined dataList (sub-cond names are disjoint by construction).
    dataList <- do.call(c, lapply(per_model, `[[`, "dataList"))
    class(dataList) <- "datalist"

    # Combined fixed: union across models. Conflicting defaults (the same
    # symbol declared in multiple SBMLs with different values) are an error;
    # add them under different names if they really mean different things.
    fixed <- Reduce(function(a, b) {
      ov <- intersect(names(a), names(b))
      if (length(ov)) {
        diff <- ov[abs(a[ov] - b[ov]) > 1e-12]
        if (length(diff))
          stop("Conflicting SBML default(s) for symbol(s) `",
               paste(diff, collapse = "`, `"), "` across models.")
      }
      c(a, b[setdiff(names(b), names(a))])
    }, lapply(per_model, `[[`, "fixed"))

    # Use a representative odemodel / sbml / sub_cond_map / cond_grid /
    # col_kind for the returned PEtabProblem slots. Multi-model fits surface
    # all per-model pieces under `attr(., "petab_meta")$models` for callers
    # that need the disaggregated view.
    odeobj <- per_model[[1L]]$odemodel
    sbml   <- per_model[[1L]]$sbml
    x <- prd  # x is conceptually the underlying ODE predictor; `prd`
              # already absorbs the per-model dispatch, so expose that.
    g <- NULL
    p <- NULL
    sub_cond_map <- do.call(rbind, lapply(per_model, `[[`, "sub_cond_map"))
    cond_grid    <- per_model[[1L]]$cond_grid
    col_kind     <- per_model[[1L]]$col_kind
    obs_meta_used <- obs_meta
  }

  if (isTRUE(compile)) {
    if (is.null(e)) dMod::compile(prd, cores = cores)
    else            dMod::compile(prd, e, cores = cores)
  }

  raw_obj <- .petab_build_objective(data = dataList, prd = prd, errmodel = e,
                                    obs_meta = obs_meta_used)

  # PEtab L2 priors: parameterScaleNormal (and `normal` on lin) become a
  # constraintL2 term added to the objective. constraintL2 reads parameter
  # values from c(pars, fixed), so this is independent of which parameters
  # end up in `baked_fixed` below.
  if (!is.null(param_meta$priors)) {
    prior_obj <- constraintL2(mu       = param_meta$priors$mu,
                              sigma    = param_meta$priors$sigma,
                              attr.name = "prior")
    raw_obj <- raw_obj + prior_obj
  }

  # Symbols pulled from sbml$pars (the "default-as-fixed" gap) that turn
  # out to be inner targets of `p` rather than outer parameters are
  # determined by the trafo (typically via conditions.tsv overrides). They
  # need an SBML <parameter> entry so targetIds resolve, but they are NOT
  # genuine fixed parameters of the PEtab problem; passing them as
  # `fixed` to the trafo would just be noise. Move them to a separate
  # `sbml_only_pars` slot so a re-export preserves the SBML declaration
  # without polluting parameters.tsv. In multi-model setups we union the
  # outer parameters from every per-model trafo.
  outer_p <- if (multi_model)
               unique(unlist(lapply(per_model, function(pm) getParameters(pm$p))))
             else
               getParameters(p)
  sbml_only_pars <- fixed[setdiff(names(fixed), outer_p)]
  fixed          <- fixed[intersect(names(fixed), outer_p)]

  # Bake `fixed` into the objective so users only pass the estimated bestfit
  # vector. `fixed = NULL` (the default) lets the closure inject the PEtab
  # fixed values; an explicit `fixed = ...` overrides per call.
  baked_fixed <- fixed
  obj <- function(pars, fixed = NULL, ...) {
    if (is.null(fixed)) fixed <- baked_fixed
    raw_obj(pars, fixed = fixed, ...)
  }

  bestfit <- param_meta$pouter
  attr(bestfit, "petab_scales") <- param_meta$scales[names(bestfit)]

  reactions_top <- if (multi_model)
                     lapply(per_model, function(pm) pm$sbml$reactions)
                   else sbml$reactions

  out <- list(
    dataList  = dataList,
    reactions = reactions_top,
    odemodel  = if (multi_model) lapply(per_model, `[[`, "odemodel") else odeobj,
    g         = if (multi_model) lapply(per_model, `[[`, "g") else g,
    x         = if (multi_model) lapply(per_model, `[[`, "x") else x,
    p         = if (multi_model) lapply(per_model, `[[`, "p") else p,
    e         = e,
    prd       = prd,
    obj       = obj,
    bestfit   = bestfit,
    parlower  = param_meta$lower,
    parupper  = param_meta$upper
  )
  attr(out, "petab_meta") <- list(
    fixed          = fixed,
    sbml_only_pars = sbml_only_pars,
    inits          = if (multi_model) lapply(per_model, function(pm) pm$sbml$inits)
                     else sbml$inits,
    modelID        = if (multi_model) names(per_model) else modelname,
    sourceYaml     = yamlPath,
    sub_cond_map   = sub_cond_map,
    obs_meta       = obs_meta,
    param_meta     = param_meta,
    # Preserve the wide-format condition grid (one row per conditionId,
    # one column per overridden target) so `exportPEtabObject` can
    # reconstruct conditions.tsv. Without this, condition-side state-init
    # / parameter overrides survive only inside the per-condition trafo
    # functions and are invisible to the low-level exporter.
    cond_grid      = cond_grid,
    col_kind       = col_kind,
    # Multi-model dispatch table (NULL for single-model problems): a list
    # keyed by modelId carrying per-model {odemodel, x, g, e, p, dataList,
    # sub_conds, sbml, fixed, obs_meta} for callers that need the
    # disaggregated view.
    models         = if (multi_model) per_model else NULL
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
  models <- meta$models
  mid_label <- if (length(meta$modelID) > 1L)
                 paste(meta$modelID, collapse = ", ")
               else (meta$modelID %||% "")
  cat("<PEtabProblem ", mid_label, ">\n", sep = "")
  if (!is.null(meta$sourceYaml))
    cat("  source:        ", meta$sourceYaml, "\n", sep = "")
  if (!is.null(models))
    cat("  models:        ", length(models), " (",
        paste(names(models), collapse = ", "), ")\n", sep = "")
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
  if (!is.null(meta$param_meta$priors))
    cat("  priors  (n=", length(meta$param_meta$priors$mu), "): ",
        paste(names(meta$param_meta$priors$mu), collapse = ", "), "\n",
        sep = "")
  invisible(x)
}


## --- exporter --------------------------------------------------------------
##
## The forward direction is the symbolic inverse of `.petab_build_trafo`
## (line 441). Given a dMod parfn `p`, we read its per-condition LHS=RHS
## eqnvec via `getEquations(p)`, strip the parameter-scale chain rule
## (replacing `10^(op)` / `exp(op)` subexpressions with bare `op`), and
## classify each LHS -- the resulting (a) constants land as SBML defaults
## or initialConcentrations, (b) constant symbolic mappings land as
## conditions.tsv columns or initialAssignments, (c) per-condition-varying
## mappings land as conditions.tsv override columns, (d) identity mappings
## are the importer's default and need no emission. The single source of
## truth is `getEquations(p)`; `attr(data, "condition.grid")` is ignored.

# Internal: does `e` syntactically match the bare symbol `op_sym`, possibly
# wrapped in redundant parentheses (`(X)` parses as a call to `(`)? dMod's
# `repar` emits `10^(K)` which parses as `^(10, (K))` -- i.e. the exponent
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
# `op -> 10^(op)` for log10 / `op -> exp(op)` for log on every RHS). Two
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
# Mixed wrap (e.g. `10^(K) + K`) is no longer ambiguous -- the clean wrap
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
#   "missing"        -- at least one condition has no entry for this LHS
#   "identity"       -- every condition's RHS == LHS (importer's build_default)
#   "const_numeric"  -- every condition's RHS is the same numeric literal
#   "const_symbolic" -- every condition's RHS is the same symbolic formula
#   "varying"        -- RHS differs across conditions
.petab_classify_lhs <- function(stripped_eqs, lhs, conds) {
  rhs <- vapply(conds, function(c) {
    e <- stripped_eqs[[c]]
    v <- if (lhs %in% names(e)) e[[lhs]] else NA_character_
    if (length(v) == 0L || is.na(v)) NA_character_ else as.character(v)
  }, character(1))
  if (any(is.na(rhs))) return(list(kind = "missing"))
  if (all(rhs == lhs)) return(list(kind = "identity"))
  # Try to collapse closed-form numeric expressions (e.g. "10^0" -> 1) so they
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
#   inits           named list keyed by state -- numeric (initialConcentration)
#                   or character (initialAssignment formula).
#   sbml_extra_pars named numeric -- additional fixed-default parameters that
#                   must appear in SBML so condition.tsv columns and
#                   collapsed inner_pars resolve to a declared SId.
.petab_decompose_trafo <- function(eqs, states, inner_pars, obs_inner,
                                   pouter_names, fixed, scales,
                                   formatVersion = 1L) {

  conds        <- names(eqs)
  # v1: strip 10^(...) / exp(...) wraps so they land in parameters.tsv as
  # `parameterScale = log10/log` + linear nominalValue. v2 has no
  # parameterScale column -- the trafo's wraps stay intact and end up in
  # conditions.tsv:targetValue / SBML <initialAssignment>, where v2 expects
  # parameter transformations to live.
  stripped     <- if (identical(formatVersion, 1L))
                    .petab_strip_trafo(eqs, scales)
                  else
                    lapply(eqs, function(eqv) {
                      out <- as.character(unclass(eqv))
                      names(out) <- names(eqv); out
                    })
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
          # otherwise leave to the post-loop default (0) -- undeclared
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


#' Export a dMod problem to PEtab v1 or v2
#'
#' Symbolically decomposes a dMod parameter transformation `p` and writes
#' the corresponding PEtab problem (parameters / observables / conditions
#' / measurements TSVs, an SBML model, and a YAML manifest; v2 additionally
#' writes an experiments TSV). Output format is selected via
#' `formatVersion` (default `"2.0.0"`). Parameter scales, fixed values,
#' and per-condition overrides are read off `p` and written into the
#' corresponding PEtab tables.
#'
#' @section Limitations:
#' Pre-equilibration cannot be expressed through the trafo `p` alone — use
#' [exportPEtabObject()] for those problems. Fixed parameters are always
#' written with `parameterScale = "lin"`.
#'
#' @param data A [datalist] (or a list of data.frames keyed by condition,
#'   each with `name`, `time`, `value` columns; or a long-format data.frame
#'   that [as.datalist()] accepts). Only used as the measurement source --
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
#' @param observableTransformation Scalar or named character --
#'   `"lin"`/`"log"`/`"log10"` per observableId.
#' @param noiseDistribution Scalar or named character -- `"normal"`/
#'   `"laplace"`/`"log-normal"` per observableId.
#' @param modelID SBML model identifier; defaults to `"dMod_export"`.
#' @param dir Output directory; created if missing.
#' @param formatVersion Output format. One of `"2.0.0"` (default) or
#'   `"1"`. Forwarded to [exportPEtabObject()].
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
                        modelID = "dMod_export",
                        formatVersion = "2.0.0",
                        dir, overwrite = FALSE) {

  ## --- 1. datalist normalisation ------------------------------------------
  if (is.data.frame(data)) data <- as.datalist(data)
  if (is.list(data) && !inherits(data, "datalist")) data <- as.datalist(data)
  if (!inherits(data, "datalist"))
    stop("`data` must be a datalist, list of data.frames, or long-format data.frame.")

  ## --- 2. extract per-condition trafo from p ------------------------------
  if (!is.function(p) || is.null(attr(p, "mappings")))
    stop("`p` must be a parfn produced by P() -- needed to decompose the parameter trafo.")
  eqs <- getEquations(p)
  if (!is.list(eqs)) eqs <- list(eqs)
  conds <- names(eqs)
  if (is.null(conds) || any(!nzchar(conds)))
    stop("`getEquations(p)` returned an unnamed list -- every condition must have a name.")

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
  # measurements.tsv, not real inner parameters -- strip them so the
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

  # v2 has no parameterScale column -- the trafo `p` already encodes the
  # scale via `10^(...)` / `exp(...)` wraps in its equations. Force lin so
  # later steps do not strip those wraps.
  fv_major <- .petab_major_version(formatVersion)
  if (fv_major == 2L) {
    if (any(scales_pouter != "lin"))
      warning("`parameterScale` is ignored for PEtab v2 (the parameter ",
              "transformation lives in `p`, which keeps its `10^(...)` / ",
              "`exp(...)` wraps in conditions.tsv:targetValue or SBML ",
              "<initialAssignment>). Forcing `lin` on disk.",
              call. = FALSE)
    scales_pouter <- setNames(rep("lin", length(pouter_ids)), pouter_ids)
  }

  # Fixed parameters MUST be linear (PEtab+dMod convention; see
  # exportPEtabObject:1387 -- fixed always written as parameterScale="lin").
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
    eqs           = eqs,
    states        = states,
    inner_pars    = inner_pars,
    obs_inner     = obs_inner,
    pouter_names  = pouter_ids,
    fixed         = fixed,
    scales        = scales_all,
    formatVersion = fv_major)

  # Split decomposer-emitted SBML defaults into two buckets:
  #   - bound: targets of a conditions.tsv override (the trafo determines
  #     them on import); they need an SBML <parameter> entry so targetIds
  #     resolve, but NOT a parameters.tsv row (their nominal "1" is never
  #     used).
  #   - unbound: identity / collapsed-numeric inner symbols without an
  #     override (e.g. `s -> 10^0 = 1`); these become genuine fixed outer
  #     parameters of the imported p, so they belong in parameters.tsv.
  override_targets <- setdiff(colnames(decomp$conditions_df), "conditionId")
  bound_idx        <- names(decomp$sbml_extra_pars) %in% override_targets
  sbml_only_pars   <- decomp$sbml_extra_pars[bound_idx]
  sbml_unbound     <- decomp$sbml_extra_pars[!bound_idx]
  fixed_full <- c(fixed, sbml_unbound)
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

  data_filtered    <- list()
  scm_rows         <- list()
  obs_subs_by_sub  <- list()
  noi_subs_by_sub  <- list()

  for (c in intersect(conds, data_conds)) {
    d <- data[[c]]
    if (use_per_row_noise && "sigma" %in% colnames(d) && nrow(d) > 0L) {
      sig_str <- ifelse(is.na(d$sigma), "1", formatC(d$sigma, digits = 17,
                                                     format = "g"))
      noi_h   <- .petab_subcond_label(sig_str)
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
          sub_condition               = sub_name,
          stringsAsFactors            = FALSE)
        noi_subs_by_sub[[sub_name]] <- c("*" = sig_str[sel[[1L]]])
      }
    } else {
      data_filtered[[c]] <- d
      scm_rows[[length(scm_rows) + 1L]] <- data.frame(
        simulationConditionId       = c,
        preequilibrationConditionId = "",
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
      sub_condition               = c,
      stringsAsFactors            = FALSE)
  }
  attr(data_filtered, "class") <- attr(data, "class")  # preserve datalist
  sub_cond_map <- do.call(rbind, scm_rows)
  rownames(sub_cond_map) <- NULL
  attr(sub_cond_map, "obs_subs") <- obs_subs_by_sub
  attr(sub_cond_map, "noi_subs") <- noi_subs_by_sub

  ## --- 10. tag pouter with scales attr (consumed by exportPEtabObject) ---
  attr(pouter, "petab_scales") <- scales_pouter

  petab <- list(
    pouter         = pouter,
    lower          = lower,
    upper          = upper,
    fixed          = fixed_full,
    sbml_only_pars = sbml_only_pars,
    obs_meta       = obs_meta,
    sub_cond_map   = sub_cond_map,
    condition.grid = decomp$conditions_df,
    data           = data_filtered,
    reactions      = reactions,
    inits          = decomp$inits,
    modelID        = modelID)

  exportPEtabObject(petab, dir, modelID = modelID,
                    formatVersion = formatVersion,
                    overwrite = overwrite)
}


#' Export a dMod `petabProblem` to a PEtab problem on disk
#'
#' Low-level exporter that writes the PEtab tables (parameters, observables,
#' conditions, measurements; v2 additionally writes experiments) plus the
#' SBML model and a YAML manifest, given a fully-populated `petabProblem`-
#' shaped list (i.e. an object as produced by [importPEtab()]).
#' Sub-conditions synthesised by the importer are collapsed back to their
#' PEtab condition + per-row `observableParameters` / `noiseParameters`
#' representation.
#'
#' If you have only the dMod-native pieces (datalist, odemodel, observation
#' / parameter / prediction functions, and a numeric `pouter`) and never went
#' through `importPEtab()`, use the higher-level [exportPEtab()] adapter,
#' which synthesises the missing PEtab metadata and dispatches here.
#'
#' Symbolic species initials (SBML `<initialAssignment>` elements) on
#' `petab$inits` are written verbatim through [export_sbml()] into the
#' new SBML model.
#'
#' Lossy steps documented:
#' \itemize{
#'   \item Parameter scales survive only via `attr(petab$pouter,
#'         "petab_scales")` set by the importer; hand-built problems lacking
#'         this attribute default to `"lin"`.
#'   \item v2 has no `parameterScale` column. When exporting to v2, scales
#'         other than `"lin"` are linearised on disk (values *and* bounds)
#'         with a one-time warning.
#'   \item v2 has no `log10` observable transformation. When exporting an
#'         observable with `obs_trafo == "log10"`, dMod warns and emits
#'         `noiseDistribution = "log-normal"` (mathematically `log`).
#'   \item State-init overrides that live only on the trafo `p` (i.e.
#'         `cond_grid` columns whose target is a state initial concentration,
#'         baked in by the v1 importer) are not written to the v2 long-form
#'         conditions table. Exporting such a problem to v2 therefore drops
#'         those overrides; export to `formatVersion = "1"` keeps them.
#'         The native [exportPEtab()] is unaffected because state inits
#'         there flow through the trafo + SBML `<initialAssignment>` path.
#'   \item SBML `<algebraicRule>` elements are silently skipped on import
#'         (dMod has no DAE solver), so any SBML written back from such a
#'         problem lacks them.
#' }
#'
#' @param petab A `petabProblem` produced by [importPEtab()] (or hand-built
#'   list with the same slot names).
#' @param dir Output directory; created if missing.
#' @param modelID SBML model identifier (defaults to `petab$modelID` or
#'   `"dMod_export"`).
#' @param formatVersion Output format. One of `"2.0.0"` (default) or
#'   `"1"`. v2 writes the new YAML schema with `model_files` /
#'   `experiment_files` and a long-format conditions table.
#' @param overwrite Logical. If `FALSE` (default) errors when files already
#'   exist in `dir`.
#' @return Path to the written YAML manifest, invisibly.
#' @seealso [exportPEtab()] for a dMod-native entry point.
#' @export
#' @importFrom yaml write_yaml
exportPEtabObject <- function(petab, dir, modelID = NULL,
                              formatVersion = "2.0.0",
                              overwrite = FALSE) {

  stopifnot(inherits(petab, "PEtabProblem") || is.list(petab))
  major <- .petab_major_version(formatVersion)
  if (!major %in% c(1L, 2L))
    stop("`formatVersion` must be '1' or '2.0.0' (got ", formatVersion, ").")

  meta <- attr(petab, "petab_meta") %||% list()
  fixed        <- meta$fixed        %||% petab$fixed        %||% numeric(0)
  # SBML-only parameters: declared in the SBML model so conditions.tsv
  # targetIds resolve, but excluded from parameters.tsv (their values are
  # determined by per-condition overrides, not the nominal default).
  sbml_only    <- meta$sbml_only_pars %||% petab$sbml_only_pars %||% numeric(0)
  inits_meta   <- meta$inits        %||% petab$inits
  obs_meta     <- meta$obs_meta     %||% petab$obs_meta
  sub_cond_map <- meta$sub_cond_map %||% petab$sub_cond_map
  param_meta   <- meta$param_meta   %||% petab$param_meta
  # Prefer the original PEtab conditions table (saved by importPEtab) over
  # the dMod-internal condition.grid, which only echoes conditionIds and
  # loses the override columns (a0, k1, ...) that live on the per-condition
  # trafo. exportPEtab synthesises its own cond_grid via the trafo
  # decomposer and stashes it in petab$condition.grid -- that path is also
  # honoured.
  cond_grid <- meta$cond_grid
  if (is.null(cond_grid))
    cond_grid <- petab$condition.grid
  if (is.null(cond_grid))
    cond_grid <- attr(petab$dataList %||% petab$data, "condition.grid")
  bestfit      <- petab$bestfit %||% petab$pouter
  parlower     <- petab$parlower %||% petab$lower
  parupper     <- petab$parupper %||% petab$upper
  data_list    <- petab$dataList %||% petab$data
  reactions    <- petab$reactions
  if (is.null(modelID))
    modelID <- meta$modelID %||% petab$modelID %||% "dMod_export"

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  paths <- list(
    parameters   = file.path(dir, paste0("parameters_",   modelID, ".tsv")),
    observables  = file.path(dir, paste0("observables_",  modelID, ".tsv")),
    conditions   = file.path(dir, paste0("conditions_",   modelID, ".tsv")),
    experiments  = file.path(dir, paste0("experiments_",  modelID, ".tsv")),
    measurements = file.path(dir, paste0("measurements_", modelID, ".tsv")),
    sbml         = file.path(dir, paste0(modelID, ".xml")),
    yaml         = file.path(dir, paste0(modelID, ".yaml"))
  )
  must_exist <- c("parameters", "observables", "conditions",
                  "measurements", "sbml", "yaml")
  if (major == 2L) must_exist <- c(must_exist, "experiments")
  if (!overwrite) {
    existing <- vapply(paths[must_exist], file.exists, logical(1))
    if (any(existing))
      stop("Output file(s) already exist; pass overwrite = TRUE to replace: ",
           paste(unlist(paths[must_exist])[existing], collapse = ", "))
  }

  scales <- attr(bestfit, "petab_scales")
  if (is.null(scales))
    scales <- setNames(rep("lin", length(bestfit)), names(bestfit))

  est_ids <- names(bestfit)

  if (major == 1L) {
    # v1 nominalValue / lowerBound / upperBound live on the linear scale
    # regardless of parameterScale. Internally pouter is on the parameter
    # scale, so we invert the importer's forward transform here.
    apply_inv_scale <- function(values, ids) {
      sc <- scales[ids]
      out <- values
      log_idx   <- which(sc == "log")
      log10_idx <- which(sc == "log10")
      if (length(log_idx))   out[log_idx]   <- exp(values[log_idx])
      if (length(log10_idx)) out[log10_idx] <- 10 ^ values[log10_idx]
      out
    }
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
  } else {
    # v2 has no parameterScale column. The parameter transformation lives
    # on the trafo `p` (and via export -- in conditions.tsv:targetValue or
    # SBML <initialAssignment>), so the literal value of the outer
    # parameter is what goes into nominalValue. No inverse-scaling and no
    # warning: the optimisation scale is preserved through the wraps in
    # the model description, exactly as v2 expects.
    #
    # v2 also drops `estimate = false` rows: the canonical PEtab v2 test
    # suite (e.g. case 0010) keeps non-estimated parameters in the SBML
    # <parameter> table only and lists just estimated parameters here.
    # Their nominal values are still written to SBML below via
    # `c(bestfit, fixed, sbml_only)`, so the importer's
    # SBML-default-as-fixed gap recovers them on round-trip.
    est_df <- data.frame(
      parameterId  = est_ids,
      lowerBound   = unname(parlower[est_ids]),
      upperBound   = unname(parupper[est_ids]),
      nominalValue = unname(bestfit),
      estimate     = "true",
      stringsAsFactors = FALSE
    )
    fixed_df <- NULL
  }

  # Round-trip priors (constraintL2 -> priorDistribution / priorParameters).
  # We emit `parameterScaleNormal` on every prior'd row regardless of the
  # original PEtab distribution: dMod's L2 prior acts on the parameter as
  # the optimizer sees it (i.e. on parameterScale), so this is the
  # unambiguous v2 spelling. For `lin` parameters this coincides with
  # `normal`, so it round-trips a `normal` lin prior without a semantic
  # shift.
  priors <- param_meta$priors
  if (!is.null(priors)) {
    pd <- rep(NA_character_, nrow(est_df))
    pp <- rep(NA_character_, nrow(est_df))
    hits <- match(names(priors$mu), est_df$parameterId)
    keep <- !is.na(hits)
    pd[hits[keep]] <- "parameterScaleNormal"
    pp[hits[keep]] <- sprintf("%g;%g",
                              priors$mu[keep], priors$sigma[keep])
    if (major == 1L) {
      est_df$objectivePriorType       <- pd
      est_df$objectivePriorParameters <- pp
      if (!is.null(fixed_df)) {
        fixed_df$objectivePriorType       <- NA_character_
        fixed_df$objectivePriorParameters <- NA_character_
      }
    } else {
      est_df$priorDistribution <- pd
      est_df$priorParameters   <- pp
    }
  }

  utils::write.table(rbind(est_df, fixed_df), paths$parameters,
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  ## --- observables.tsv ---------------------------------------------------
  om <- obs_meta
  if (major == 1L) {
    obs_df <- data.frame(
      observableId             = names(om$obs),
      observableFormula        = unname(om$obs),
      observableTransformation = unname(om$obs_trafo[names(om$obs)]),
      noiseFormula             = unname(om$noise[names(om$obs)]),
      noiseDistribution        = unname(om$noise_dist[names(om$obs)]),
      stringsAsFactors = FALSE
    )
  } else {
    obs_ids <- names(om$obs)
    obs_trafo  <- om$obs_trafo[obs_ids]
    noise_dist <- om$noise_dist[obs_ids]
    if (any(obs_trafo == "log10"))
      warning("PEtab v2 dropped the `log10` observable transformation; ",
              "emitting `log-normal` for: ",
              paste(obs_ids[obs_trafo == "log10"], collapse = ", "),
              ". Likelihood values shift by the log10/log Jacobian.",
              call. = FALSE)
    v2_dist <- ifelse(obs_trafo %in% c("log", "log10") &
                      noise_dist == "normal", "log-normal",
              ifelse(obs_trafo == "lin" & noise_dist == "normal", "normal",
              ifelse(obs_trafo %in% c("log", "log10") &
                     noise_dist == "laplace", "log-laplace",
              ifelse(obs_trafo == "lin" & noise_dist == "laplace", "laplace",
                     NA_character_))))
    if (any(is.na(v2_dist)))
      stop("Cannot encode (obs_trafo, noise_dist) into v2 noiseDistribution: ",
           paste(sprintf("%s=(%s,%s)", obs_ids[is.na(v2_dist)],
                         obs_trafo[is.na(v2_dist)],
                         noise_dist[is.na(v2_dist)]), collapse = "; "))
    obs_form <- unname(om$obs[obs_ids])
    noi_form <- unname(om$noise[obs_ids])
    scan_placeholders <- function(formula, prefix, obs_id) {
      if (is.na(formula) || !nzchar(formula)) return("")
      pat <- sprintf("%s[0-9]+_%s", prefix, obs_id)
      m <- regmatches(formula, gregexpr(pat, formula, perl = TRUE))[[1L]]
      if (!length(m)) return("")
      idx <- as.integer(sub(sprintf("^%s([0-9]+)_.*", prefix), "\\1", m))
      uniq <- m[!duplicated(idx)]
      uniq <- uniq[order(idx[!duplicated(idx)])]
      paste(uniq, collapse = ";")
    }
    obs_ph   <- vapply(seq_along(obs_ids), function(i)
      scan_placeholders(obs_form[i], "observableParameter", obs_ids[i]),
      character(1))
    noise_ph <- vapply(seq_along(obs_ids), function(i)
      scan_placeholders(noi_form[i], "noiseParameter", obs_ids[i]),
      character(1))
    obs_df <- data.frame(
      observableId           = obs_ids,
      observableFormula      = obs_form,
      observablePlaceholders = obs_ph,
      noiseFormula           = noi_form,
      noiseDistribution      = v2_dist,
      noisePlaceholders      = noise_ph,
      stringsAsFactors = FALSE
    )
  }
  utils::write.table(obs_df, paths$observables, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")

  ## --- conditions.tsv ---------------------------------------------------
  scm <- sub_cond_map
  uniq_sims  <- unique(scm$simulationConditionId)
  uniq_preeq <- setdiff(unique(scm$preequilibrationConditionId), "")
  uniq_conds <- unique(c(uniq_sims, uniq_preeq))

  if (major == 1L) {
    cond_df <- cond_grid[uniq_sims, , drop = FALSE]
    cond_df$conditionId <- uniq_sims
    cond_df <- cond_df[, c("conditionId",
                           setdiff(colnames(cond_df), "conditionId")),
                       drop = FALSE]
  } else {
    # long-format: one row per (conditionId, targetId) override.
    # Skip the dMod-internal `condition` column whose values just echo
    # conditionId -- it's a `as.datalist` artifact, not a real override.
    rows <- list()
    cg_cols <- setdiff(colnames(cond_grid), c("conditionId", "condition"))
    for (cid in uniq_conds) {
      if (!cid %in% rownames(cond_grid)) next
      for (tid in cg_cols) {
        v <- cond_grid[cid, tid]
        if (is.na(v) || !nzchar(trimws(as.character(v)))) next
        # Self-mapping (value == conditionId) is also a dMod artifact.
        if (identical(as.character(v), cid)) next
        rows[[length(rows) + 1L]] <- data.frame(
          conditionId  = cid,
          targetId     = tid,
          targetValue  = as.character(v),
          stringsAsFactors = FALSE)
      }
    }
    cond_df <- if (length(rows)) do.call(rbind, rows) else
      data.frame(conditionId = character(0), targetId = character(0),
                 targetValue = character(0), stringsAsFactors = FALSE)
  }
  utils::write.table(cond_df, paths$conditions, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")

  ## --- experiments.tsv (v2 only) ----------------------------------------
  # Mapping from (sim, preeq) sub-condition tuple to a synthesised
  # experimentId. Used to rewrite measurements.tsv and to populate
  # experiments.tsv.
  uniq_pairs <- unique(scm[, c("simulationConditionId",
                               "preequilibrationConditionId"),
                           drop = FALSE])
  exp_id_for <- function(sim, preeq) {
    if (identical(sim, "_dmod_default_condition_")) return("")
    if (!nzchar(preeq)) sim else paste0(preeq, "__", sim)
  }
  uniq_pairs$experimentId <- mapply(
    exp_id_for,
    uniq_pairs$simulationConditionId,
    uniq_pairs$preequilibrationConditionId)

  if (major == 2L) {
    exp_rows <- list()
    for (i in seq_len(nrow(uniq_pairs))) {
      eid <- uniq_pairs$experimentId[i]
      if (!nzchar(eid)) next  # default condition: no experiment row
      sim <- uniq_pairs$simulationConditionId[i]
      pre <- uniq_pairs$preequilibrationConditionId[i]
      if (nzchar(pre))
        exp_rows[[length(exp_rows) + 1L]] <- data.frame(
          experimentId = eid, time = "-inf", conditionId = pre,
          stringsAsFactors = FALSE)
      exp_rows[[length(exp_rows) + 1L]] <- data.frame(
        experimentId = eid, time = "0", conditionId = sim,
        stringsAsFactors = FALSE)
    }
    exp_df <- if (length(exp_rows)) do.call(rbind, exp_rows) else
      data.frame(experimentId = character(0), time = character(0),
                 conditionId = character(0), stringsAsFactors = FALSE)
    utils::write.table(exp_df, paths$experiments, sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "")
  }

  ## --- measurements.tsv -------------------------------------------------
  scm_obs_subs <- attr(scm, "obs_subs") %||% list()
  scm_noi_subs <- attr(scm, "noi_subs") %||% list()
  pick_subs <- function(sub, obsId, map) {
    m <- map[[sub]]
    if (is.null(m)) return("")
    if (obsId %in% names(m)) return(unname(m[obsId]))
    if ("*" %in% names(m))   return(unname(m["*"]))
    ""
  }
  meas_rows <- do.call(rbind, lapply(names(data_list), function(sub) {
    df <- data_list[[sub]]
    ix <- match(sub, scm$sub_condition)
    sim <- scm$simulationConditionId[ix]
    pre <- scm$preequilibrationConditionId[ix]
    obs_par <- vapply(df$name, function(o) pick_subs(sub, o, scm_obs_subs),
                      character(1))
    noi_par <- vapply(df$name, function(o) pick_subs(sub, o, scm_noi_subs),
                      character(1))
    base <- data.frame(
      observableId         = df$name,
      time                 = df$time,
      measurement          = df$value,
      observableParameters = obs_par,
      noiseParameters      = noi_par,
      stringsAsFactors = FALSE)
    if (major == 1L) {
      base$simulationConditionId       <- sim
      base$preequilibrationConditionId <- pre
    } else {
      pair_ix <- which(uniq_pairs$simulationConditionId == sim &
                       uniq_pairs$preequilibrationConditionId == pre)
      base$experimentId <- uniq_pairs$experimentId[pair_ix]
    }
    base
  }))
  if (major == 1L) {
    base_cols <- c("observableId", "simulationConditionId", "time",
                   "measurement", "preequilibrationConditionId",
                   "observableParameters", "noiseParameters")
    meas_rows <- meas_rows[, base_cols, drop = FALSE]
  } else {
    base_cols <- c("observableId", "experimentId", "time", "measurement",
                   "observableParameters", "noiseParameters")
    meas_rows <- meas_rows[, base_cols, drop = FALSE]
  }
  drop_optional <- c("preequilibrationConditionId", "observableParameters",
                     "noiseParameters", "experimentId")
  for (col in intersect(drop_optional, colnames(meas_rows))) {
    if (all(meas_rows[[col]] == "")) meas_rows[[col]] <- NULL
  }
  utils::write.table(meas_rows, paths$measurements, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")

  ## --- SBML --------------------------------------------------------------
  all_pars <- c(bestfit, fixed, sbml_only)
  inits <- inits_meta %||%
           setNames(rep(0, length(reactions$states)), reactions$states)
  export_sbml(reactions, parameters = all_pars, inits = inits,
              filepath = paths$sbml, modelID = modelID)

  ## --- YAML manifest -----------------------------------------------------
  if (major == 1L) {
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
  } else {
    manifest <- list(
      format_version    = "2.0.0",
      parameter_files   = list(basename(paths$parameters)),
      model_files       = setNames(
        list(list(location = basename(paths$sbml), language = "sbml")),
        modelID),
      observable_files  = list(basename(paths$observables)),
      measurement_files = list(basename(paths$measurements)),
      condition_files   = list(basename(paths$conditions)),
      experiment_files  = list(basename(paths$experiments))
    )
  }
  yaml::write_yaml(manifest, paths$yaml)

  invisible(paths$yaml)
}


# `%||%` is base R since 4.4; relied on per dMod's existing usage.
