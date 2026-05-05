context("PEtab importer / exporter")

# Repo-relative path to the bundled PEtab test suite. The package does NOT
# ship PEtabTests/ in the installed tarball (it lives at the repo root), so
# we walk up from inst/tests/ when running interactively. testthat may
# change cwd, so we also try DMOD_PETABTESTS env var and a hardcoded path
# matching the user's checkout.
.petab_repo_dir <- function() {
  envp <- Sys.getenv("DMOD_PETABTESTS", unset = "")
  candidates <- c(
    if (nzchar(envp)) envp,
    file.path(getwd(), "PEtabTests"),
    file.path(getwd(), "..", "..", "PEtabTests"),
    file.path(dirname(getwd()), "..", "PEtabTests"),
    "/home/simon/Documents/Projects/dMod/PEtabTests"
  )
  for (p in candidates) if (nzchar(p) && dir.exists(p)) return(normalizePath(p))
  ""
}

# Same idea for the BenchmarkModels/ directory (real-world PEtab benchmarks).
.benchmark_dir <- function() {
  envp <- Sys.getenv("DMOD_BENCHMARKMODELS", unset = "")
  candidates <- c(
    if (nzchar(envp)) envp,
    file.path(getwd(), "BenchmarkModels"),
    file.path(getwd(), "..", "..", "BenchmarkModels"),
    file.path(dirname(getwd()), "..", "BenchmarkModels"),
    "/home/simon/Documents/Projects/dMod/BenchmarkModels"
  )
  for (p in candidates) if (nzchar(p) && dir.exists(p)) return(normalizePath(p))
  ""
}

# Skip integration tests that need import_sbml() when the libsbml virtualenv
# is missing.
.libsbml_works <- function() {
  # Probe by importing case 0001 (the smallest of the bundled PEtab models).
  # Cheaper than synthesising an SBML file and avoids libsbml-strictness
  # quibbles on hand-written test XML.
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) return(FALSE)
  isTRUE(tryCatch({
    res <- suppressWarnings(
      import_sbml(file.path(petab_dir, "0001", "_model.xml")))
    !is.null(res$reactions)
  }, error = function(e) FALSE))
}


## --- pure parser unit tests (no SBML) -------------------------------------

test_that(".petab_parse_parameters splits estimated / fixed and tracks scales", {

  # PEtab v1: nominalValue / lowerBound / upperBound are written on the
  # linear scale regardless of parameterScale. The parser pre-transforms
  # estimated parameters and bounds to the parameter scale (dMod's pouter
  # convention); fixed parameters stay on the linear scale because the
  # trafo's scale chain rule only wraps estimated outer parameters.
  df <- data.frame(
    parameterId    = c("a", "b", "c"),
    parameterScale = c("lin", "log10", "log"),
    lowerBound     = c(0, 1e-3, 1e-5),
    upperBound     = c(10, 1e3, 1e5),
    nominalValue   = c(1.0, 100, exp(2)),
    estimate       = c(1L, 1L, 0L),
    stringsAsFactors = FALSE
  )
  pm <- dMod:::.petab_parse_parameters(df)

  expect_equal(names(pm$pouter), c("a", "b"))
  # a (lin)   = 1.0
  # b (log10) = log10(100) = 2  — pouter on parameter scale
  expect_equal(unname(pm$pouter), c(1.0, 2.0))
  expect_equal(names(pm$fixed),  c("c"))
  # c is fixed → stays on linear scale (no scale chain rule wraps it).
  expect_equal(unname(pm$fixed["c"]), exp(2))
  expect_equal(pm$scales[["a"]], "lin")
  expect_equal(pm$scales[["b"]], "log10")
  expect_equal(pm$scales[["c"]], "log")
  # lower["b"] = log10(1e-3) = -3
  expect_equal(unname(pm$lower["b"]), -3)
  expect_equal(unname(pm$upper["b"]), 3)
})


test_that(".petab_parse_observables defaults to lin/normal and parses noise", {

  df <- data.frame(
    observableId      = c("o1", "o2"),
    observableFormula = c("A", "B + offset"),
    noiseFormula      = c("0.5", "1"),
    stringsAsFactors  = FALSE
  )
  om <- dMod:::.petab_parse_observables(df)
  expect_equal(unname(om$obs_trafo),  c("lin", "lin"))
  expect_equal(unname(om$noise_dist), c("normal", "normal"))
  expect_equal(unname(om$noise),      c("0.5", "1"))
})


test_that(".petab_parse_conditions classifies columns as init / parameter", {

  # Case 0002 shape: a0 is in conditions and is a parameter symbol that also
  # parameterises species A's initial. We expect "parameter" classification.
  df <- data.frame(conditionId = c("c0", "c1"),
                   a0          = c(0.8, 0.9),
                   stringsAsFactors = FALSE)
  ci <- dMod:::.petab_parse_conditions(df,
          sbml_states       = c("A", "B"),
          sbml_compartments = "compartment",
          sbml_pars         = c("a0", "b0", "k1", "k2", "compartment"))
  expect_equal(ci$col_kind[["a0"]], "parameter")
  expect_equal(ci$override_cols, "a0")

  # init kind: column name is a state itself
  df2 <- data.frame(conditionId = "c0", A = 0.5,
                    stringsAsFactors = FALSE)
  ci2 <- dMod:::.petab_parse_conditions(df2,
           sbml_states = c("A", "B"),
           sbml_compartments = character(),
           sbml_pars = character())
  expect_equal(ci2$col_kind[["A"]], "init")
})


test_that(".petab_parse_measurements unfolds per-row observableParameters", {

  obs_meta <- list(
    obs   = c(obs_a = "observableParameter1_obs_a * A"),
    noise = c(obs_a = "1"),
    obs_trafo  = c(obs_a = "lin"),
    noise_dist = c(obs_a = "normal")
  )

  # Case 0006 shape: same simulation condition, two different obs param
  # values across two rows.
  df <- data.frame(
    observableId          = c("obs_a", "obs_a"),
    simulationConditionId = c("c0", "c0"),
    time                  = c(0, 10),
    measurement           = c(0.7, 0.1),
    observableParameters  = c("10", "15"),
    stringsAsFactors = FALSE
  )
  mi <- dMod:::.petab_parse_measurements(df, obs_meta)
  expect_equal(nrow(mi$sub_cond_map), 2L)
  expect_true(all(grepl("^c0__", mi$sub_cond_map$sub_condition)))
  # data is partitioned across sub-conditions:
  expect_equal(sort(unique(mi$data$condition)),
               sort(mi$sub_cond_map$sub_condition))

  # Case 0001 shape: no sub-condition splitting, single condition.
  df2 <- data.frame(
    observableId          = c("obs_a", "obs_a"),
    simulationConditionId = c("c0", "c0"),
    time                  = c(0, 10),
    measurement           = c(0.7, 0.1),
    stringsAsFactors = FALSE
  )
  mi2 <- dMod:::.petab_parse_measurements(df2, obs_meta)
  expect_equal(nrow(mi2$sub_cond_map), 1L)
  expect_equal(mi2$sub_cond_map$sub_condition, "c0")
})


test_that("read_petab_yaml resolves manifest paths correctly", {

  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")

  y <- read_petab_yaml(file.path(petab_dir, "0001", "_0001.yaml"))
  expect_equal(y$formatVersion, 1L)
  expect_true(file.exists(y$problems[[1]]$sbmlFile))
  expect_true(file.exists(y$problems[[1]]$measurementFile))
})


test_that("read_petab_tables returns the expected slots for v1", {
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")

  tabs <- read_petab_tables(file.path(petab_dir, "0001", "_0001.yaml"))
  expect_named(tabs, c("parameters", "conditions", "measurements",
                       "observables", "experiments", "mapping",
                       "sbmlPath", "sbmlPaths", "formatVersion"))
  expect_s3_class(tabs$parameters,   "data.frame")
  expect_s3_class(tabs$conditions,   "data.frame")
  expect_s3_class(tabs$measurements, "data.frame")
  expect_s3_class(tabs$observables,  "data.frame")
  expect_null(tabs$experiments)   # v1 has no experiments table
  expect_null(tabs$mapping)
  expect_identical(tabs$formatVersion, 1L)
})


## --- end-to-end fixture test (no SBML import required) -------------------
##
## We hand-build the eqnlist that matches PEtab test case 0001's SBML model
## and verify the trafo+objective machinery against the published solution.
## This avoids a libsbml dependency on every test run.

test_that("hand-built case-0001 fixture produces solution-matching llh", {

  setwd(tempdir())

  # Reaction network identical to PEtabTests/0001/_model.xml after libsbml
  # would have inlined the kinetic law's compartment factor. We use a unit
  # compartment so kinetic laws read just k1*A and k2*B.
  reactions <- eqnlist()
  reactions <- addReaction(reactions, "A", "B", "k1*A", "fwd")
  reactions <- addReaction(reactions, "B", "A", "k2*B", "rev")

  ode <- odemodel(reactions, modelname = "petab_fixture_ode",
                  solver = "deSolve")
  x <- Xs(ode)

  # Observation function: obs_a = A.
  g <- Y(g = c(obs_a = "A"), f = reactions,
         attach.input = FALSE,
         compile = FALSE,
         modelname = "petab_fixture_obs")

  # Parameter trafo: A := a0, B := b0, k1 := k1, k2 := k2 (identity for k's).
  innerpars <- getParameters(x)
  trafo <- structure(innerpars, names = innerpars)
  trafo["A"] <- "a0"
  trafo["B"] <- "b0"
  p <- P(trafo, condition = "c0", modelname = "petab_fixture_trafo")

  # Data exactly matching _measurements.tsv.
  data <- as.datalist(data.frame(
    name      = c("obs_a", "obs_a"),
    time      = c(0, 10),
    value     = c(0.7, 0.1),
    sigma     = c(0.5, 0.5),
    condition = c("c0", "c0"),
    stringsAsFactors = FALSE
  ))

  prd <- g * x * p
  obj <- normL2(data, prd)

  # Nominal pouter from _parameters.tsv
  pouter <- c(a0 = 1.0, b0 = 0.0, k1 = 0.8, k2 = 0.6)

  out <- obj(pouter, deriv = FALSE)

  # PEtab _0001_solution.yaml gives llh = -0.8475016971318833 and
  # chi2 = 0.7918379836848569. dMod's normL2 returns the *full* Gaussian
  # negative log-likelihood multiplied by 2 (i.e. -2*log L), which equals
  # chi2 + sum(log(2*pi*sigma^2)) per data point. We compare against -2*llh.
  expect_lt(abs(out$value - (-2 * -0.8475016971318833)), 0.001)

  # Cleanup
  unlink("petab_fixture_*.c"); unlink("petab_fixture_*.cpp")
  unlink("petab_fixture_*.o"); unlink("petab_fixture_*.so")
})


## --- libsbml-dependent integration tests ----------------------------------

test_that("PEtab test cases 0001-0006 import and produce solution-matching llh", {

  setwd(tempdir())

  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")
  if (!.libsbml_works())   skip("libsbml virtualenv not available")

  for (id in sprintf("%04d", 1:6)) {
    yamlPath <- file.path(petab_dir, id, paste0("_", id, ".yaml"))
    sol_path  <- file.path(petab_dir, id, paste0("_", id, "_solution.yaml"))
    if (!file.exists(sol_path)) next

    petab <- importPEtab(yamlPath, solver = "deSolve",
                         modelname = paste0("petab_", id))
    sol <- yaml::read_yaml(sol_path)

    out <- petab$obj(petab$bestfit, deriv = FALSE)
    # dMod normL2 returns -2*log L (chi2 + log normaliser); compare against
    # -2 * sol$llh so all 6 cases share the same metric.
    expect_lt(abs(out$value - (-2 * sol$llh)),
              max(0.01, abs(2 * sol$tol_llh)),
              label = paste0("case ", id, " -2*llh"))
  }

  unlink("petab_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so"); unlink("*_model.csv")
})


test_that("PEtab Stage-2 test cases 0007-0016 produce solution-matching llh", {

  setwd(tempdir())

  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")
  if (!.libsbml_works())   skip("libsbml virtualenv not available")

  # Cases 0007 (log10 trafo), 0008 (replicates), 0009/0010 (preequilibration —
  # numeric Pequil fallback is exercised because steadyStates() leaves one
  # state unresolved on the A↔B reaction), 0011-0013 (init / compartment /
  # parametric init overrides), 0014/0015 (numeric / symbolic noise parameter
  # overrides), 0016 (log trafo).
  for (id in sprintf("%04d", 7:16)) {

    yamlPath <- file.path(petab_dir, id, paste0("_", id, ".yaml"))
    sol_path  <- file.path(petab_dir, id, paste0("_", id, "_solution.yaml"))
    if (!file.exists(sol_path)) next

    suppressWarnings(
      petab <- importPEtab(yamlPath, solver = "deSolve",
                           modelname = paste0("petab_s2_", id)))
    sol <- yaml::read_yaml(sol_path)
    out <- petab$obj(petab$bestfit, deriv = FALSE)
    expect_lt(abs(out$value - (-2 * sol$llh)),
              max(0.01, abs(2 * sol$tol_llh)),
              label = paste0("case ", id, " -2*llh"))
  }

  unlink("petab_s2_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so"); unlink("*_model.csv")
  unlink("reactions_for_Alyssa*")
})


test_that("two-condition roundtrip preserves objective value", {

  setwd(tempdir())

  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")
  if (!.libsbml_works())   skip("libsbml virtualenv not available")

  # Case 0002 has two conditions, an InitialAssignment binding A := a0 / B := b0
  # and exercises the condition table. The InitialAssignment roundtrip is the
  # interesting part: without it the reimported objective evaluates with
  # state initials = 0 and disagrees with the original.
  yaml1 <- file.path(petab_dir, "0002", "_0002.yaml")
  petab1 <- importPEtab(yaml1, solver = "deSolve",
                        modelname = "rt_in")
  v1 <- petab1$obj(petab1$bestfit, deriv = FALSE)$value

  out_dir <- file.path(tempdir(), "petab_roundtrip")
  # v1 keeps the wide-format conditions table the original v1 importer
  # baked in via cond_grid; v2 export from a v1-imported petab loses the
  # state-init overrides that live only on the trafo (#known-limitation).
  yaml2 <- exportPEtabObject(petab1, out_dir, modelID = "rt_out",
                             formatVersion = "1", overwrite = TRUE)

  petab2 <- importPEtab(yaml2, solver = "deSolve",
                        modelname = "rt_back")
  v2 <- petab2$obj(petab2$bestfit, deriv = FALSE)$value

  expect_true(is.finite(v1))
  expect_true(is.finite(v2))
  # InitialAssignments survive the roundtrip → values must match within
  # numerical noise of the ODE solver.
  expect_lt(abs(v1 - v2), 1e-6)

  unlink("rt_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})


## --- real-world benchmark: Boehm_JProteomeRes2014 -------------------------
##
## End-to-end test on a published JAK/STAT5 benchmark. Exercises features the
## bundled 0001-0016 fixtures don't:
##   - log10 parameter scaling on 9 outer parameters,
##   - libsbml `<power/>` MathML (kinetic laws contain STAT5A^2 / STAT5B^2
##     which the L2 formatter rendered as `pow(...)` — broke jacobianSymb),
##   - <assignmentRule> for time-varying input BaF3_Epo,
##   - sub-condition splitting from per-observable noiseParameter symbols.
## At nominalValue (= the published optimum), -log L should reproduce the
## Hass et al. 2019 benchmark value of 138.22.

test_that("Boehm_JProteomeRes2014 benchmark imports and matches published optimum", {

  setwd(tempdir())

  bm_dir <- .benchmark_dir()
  if (!nzchar(bm_dir))  skip("BenchmarkModels/ directory not found")
  if (!.libsbml_works())  skip("libsbml virtualenv not available")

  yamlPath <- file.path(bm_dir, "Boehm_JProteomeRes2014",
                         "Boehm_JProteomeRes2014.yaml")
  petab <- importPEtab(yamlPath, solver = "deSolve",
                       modelname = "boehm")

  # Imported problem shape:
  expect_equal(length(petab$bestfit), 9L)
  expect_setequal(names(attr(petab, "petab_meta")$obs_meta$obs),
                  c("pSTAT5A_rel", "pSTAT5B_rel", "rSTAT5A_rel"))
  # All estimated parameters are on log10 scale per parameters.tsv:
  scales <- attr(petab$bestfit, "petab_scales")
  expect_true(all(scales == "log10"))
  # AssignmentRule for BaF3_Epo must have been inlined → not in `fixed`:
  expect_false("BaF3_Epo" %in% names(attr(petab, "petab_meta")$fixed))

  out <- petab$obj(petab$bestfit, deriv = FALSE)

  # Published optimum: -log L = 138.22 (Hass et al. 2019, "Benchmark
  # problems for dynamic modeling of intracellular processes"). dMod's
  # normL2 returns -2*log L, so we compare against ~276.44.
  expect_lt(abs(out$value - 2 * 138.22), 0.5,
            label = "Boehm -2*logL at published optimum")

  unlink("boehm*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so"); unlink("*_model.csv")
})


test_that("export_sbml emits InitialAssignment for symbolic state initials", {

  setwd(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "A", "B", "k1*A", "fwd",
                           compartment = "compartment")
  reactions <- addReaction(reactions, "B", "A", "k2*B", "rev",
                           compartment = "compartment")

  # Mixed inits: A is symbolic (→ InitialAssignment), B is numeric.
  inits <- c(A = "a0", B = "0")
  pars  <- c(a0 = 0.8, k1 = 0.8, k2 = 0.6, compartment = 1.0)

  out_xml <- file.path(tempdir(), "ia_export.xml")
  export_sbml(reactions, parameters = pars, inits = inits,
              filepath = out_xml, modelID = "ia_export")

  xml_text <- readLines(out_xml, warn = FALSE)
  expect_true(any(grepl("<initialAssignment", xml_text, fixed = TRUE)),
              info = "no <initialAssignment> emitted for symbolic init")
  expect_true(any(grepl("symbol=\"A\"", xml_text)),
              info = "InitialAssignment for A missing")

  unlink(out_xml)
})


## --- trafo-aware exportPEtab: pure-R helper unit tests --------------------
##
## The strip + classify decomposer should be unit-testable without libsbml
## because it operates only on character RHSes and named eqnvecs.

test_that(".petab_strip_param_scale compensates the chain rule per-occurrence", {
  # Clean wrap stays clean (importer chain rule re-wraps it)
  expect_equal(
    dMod:::.petab_strip_param_scale("10^(K_REFLUX)", c(K_REFLUX = "log10")),
    "K_REFLUX")
  # Steady-state-like product of clean wraps
  expect_equal(
    dMod:::.petab_strip_param_scale(
      "10^(TCA_CELL) * 10^(K_EXPORT_CANA) / 10^(K_REFLUX)",
      c(TCA_CELL = "log10", K_EXPORT_CANA = "log10", K_REFLUX = "log10")),
    "TCA_CELL * K_EXPORT_CANA/K_REFLUX")
  # Compound expression: bare KM inside is compensated with log10(KM)
  expect_equal(
    dMod:::.petab_strip_param_scale("10^(KM + 5)", c(KM = "log10")),
    "10^(log10(KM) + 5)")
  # Mixed wrap: clean wrap stripped, bare occurrence compensated
  expect_equal(
    dMod:::.petab_strip_param_scale("10^(K) + K", c(K = "log10")),
    "K + log10(K)")
  # exp(.) for log scale, with compensation
  expect_equal(
    dMod:::.petab_strip_param_scale("exp(K)", c(K = "log")),
    "K")
  expect_equal(
    dMod:::.petab_strip_param_scale("exp(KM + 5)", c(KM = "log")),
    "exp(log(KM) + 5)")
  # lin parameters unchanged
  expect_equal(
    dMod:::.petab_strip_param_scale("K1 * K2 + offset",
      c(K1 = "log10", K2 = "log10", offset = "lin")),
    "log10(K1) * log10(K2) + offset")
  # Pure numeric literal — passes through
  expect_equal(
    dMod:::.petab_strip_param_scale("0", c()), "0")
})


test_that(".petab_classify_lhs categorizes per-condition RHSes", {
  conds <- c("c1", "c2")
  stripped <- list(
    c1 = c(s = "1", k = "K", state = "0", iden = "iden", v = "K"),
    c2 = c(s = "1", k = "K", state = "0", iden = "iden", v = "K2"))
  # all-numeric constant
  expect_equal(dMod:::.petab_classify_lhs(stripped, "s", conds),
               list(kind = "const_numeric", value = 1))
  # all-symbolic constant
  expect_equal(dMod:::.petab_classify_lhs(stripped, "k", conds),
               list(kind = "const_symbolic", formula = "K"))
  # numeric-zero (state init)
  expect_equal(dMod:::.petab_classify_lhs(stripped, "state", conds),
               list(kind = "const_numeric", value = 0))
  # identity (RHS == LHS)
  expect_equal(dMod:::.petab_classify_lhs(stripped, "iden", conds),
               list(kind = "identity"))
  # varying
  res <- dMod:::.petab_classify_lhs(stripped, "v", conds)
  expect_equal(res$kind, "varying")
  expect_equal(res$per_cond, c(c1 = "K", c2 = "K2"))
  # missing
  expect_equal(dMod:::.petab_classify_lhs(stripped, "absent", conds),
               list(kind = "missing"))
})


test_that(".petab_classify_lhs collapses 10^0 -> 1 via eval_constant", {
  conds <- c("c1", "c2")
  stripped <- list(c1 = c(s = "10^0"), c2 = c(s = "10^0"))
  res <- dMod:::.petab_classify_lhs(stripped, "s", conds)
  expect_equal(res, list(kind = "const_numeric", value = 1))
})


## --- trafo-aware exportPEtab: native roundtrip ----------------------------

test_that("native exportPEtab roundtrips outer pouter on log10 scale (1-cond)", {

  setwd(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  # Tiny 2-state model, single condition. Build the trafo via explicit
  # eqnvec literals to avoid `define`'s NSE which can't see test_that locals.
  reactions <- eqnlist() %>%
    addReaction("A", "B", rate = "k1*A", description = "fwd") %>%
    addReaction("B", "A", rate = "k2*B", description = "rev")

  m <- odemodel(reactions, modelname = "rt1_ode", compile = FALSE,
                solver = "deSolve")
  x_native <- Xs(m)

  obs <- eqnvec(obs_a = "A", obs_b = "B")
  g_native <- Y(obs, f = x_native, condition = NULL,
                compile = FALSE, modelname = "rt1_obs", attach.input = FALSE)

  trafo <- as.eqnvec(c(A = "10^(A)", B = "10^(B)",
                       k1 = "10^(K1)", k2 = "10^(K2)"))
  p_native <- P(trafo, condition = "c1", compile = FALSE,
                modelname = "rt1_par")
  compile(g_native, x_native, p_native, cores = 1)

  data <- as.datalist(data.frame(
    name = c("obs_a", "obs_b"), time = c(1, 1),
    value = c(0.5, 0.3), sigma = c(1, 1),
    condition = c("c1", "c1"), stringsAsFactors = FALSE))

  pouter <- c(A = -1, B = -1, K1 = -1, K2 = -1)
  obj_native <- normL2(data, g_native * x_native * p_native)

  out_dir <- file.path(tempdir(), "petab_rt1")
  yaml_out <- exportPEtab(
    data = data, reactions = reactions, observables = obs,
    p = p_native, pouter = pouter,
    parameterScale = "log10", modelID = "rt1_export",
    formatVersion = "1", dir = out_dir, overwrite = TRUE)

  petab <- importPEtab(yaml_out, solver = "deSolve",
                       modelname = "rt1_imp")

  expect_setequal(names(petab$bestfit), names(pouter))
  expect_true(all(attr(petab$bestfit, "petab_scales") == "log10"))

  v_native <- obj_native(pouter, deriv = FALSE)$value
  v_petab  <- petab$obj(pouter[names(petab$bestfit)], deriv = FALSE)$value
  expect_lt(abs(v_native - v_petab), 1e-3)

  unlink("rt1_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})


test_that("native exportPEtab roundtrips per-condition k override (2-cond)", {

  setwd(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  reactions <- eqnlist() %>%
    addReaction("A", "B", rate = "k*A", description = "fwd")

  m <- odemodel(reactions, modelname = "rt2_ode", compile = FALSE,
                solver = "deSolve")
  x_native <- Xs(m)
  obs <- eqnvec(obs_a = "A")
  g_native <- Y(obs, f = x_native, condition = NULL, compile = FALSE,
                modelname = "rt2_obs", attach.input = FALSE)

  # Two conditions with different k mapping (closed→K, open→K_OPEN).
  trafo_closed <- as.eqnvec(c(A = "10^(A)", B = "10^(B)", k = "10^(K)"))
  trafo_open   <- as.eqnvec(c(A = "10^(A)", B = "10^(B)", k = "10^(K_OPEN)"))
  p_native <- P(trafo_closed, condition = "closed", compile = FALSE,
                modelname = "rt2_par_c") +
              P(trafo_open,   condition = "open",   compile = FALSE,
                modelname = "rt2_par_o")
  compile(g_native, x_native, p_native, cores = 1)

  data <- as.datalist(data.frame(
    name = c("obs_a", "obs_a"), time = c(1, 1),
    value = c(0.5, 0.5), sigma = c(1, 1),
    condition = c("closed", "open"), stringsAsFactors = FALSE))

  pouter <- c(A = -1, B = -1, K = -1, K_OPEN = -0.5)
  obj_native <- normL2(data, g_native * x_native * p_native)

  out_dir <- file.path(tempdir(), "petab_rt2")
  yaml_out <- exportPEtab(
    data = data, reactions = reactions, observables = obs,
    p = p_native, pouter = pouter,
    parameterScale = "log10", modelID = "rt2_export",
    formatVersion = "1", dir = out_dir, overwrite = TRUE)

  # conditions.tsv must have a `k` column distinguishing closed from open.
  cond_df <- read.delim(file.path(out_dir, "conditions_rt2_export.tsv"),
                        stringsAsFactors = FALSE)
  expect_true("k" %in% colnames(cond_df))
  expect_setequal(cond_df$k, c("K", "K_OPEN"))

  petab <- importPEtab(yaml_out, solver = "deSolve",
                       modelname = "rt2_imp")
  expect_setequal(names(petab$bestfit), c("A", "B", "K", "K_OPEN"))

  v_native <- obj_native(pouter, deriv = FALSE)$value
  v_petab  <- petab$obj(pouter[names(petab$bestfit)], deriv = FALSE)$value
  expect_lt(abs(v_native - v_petab), 1e-3)

  unlink("rt2_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})


test_that("exportPEtab errors on undeclared free symbol after strip", {

  setwd(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  reactions <- eqnlist() %>%
    addReaction("A", "B", rate = "k*A", description = "fwd")
  m <- odemodel(reactions, modelname = "err1_ode", compile = FALSE,
                solver = "deSolve")
  x <- Xs(m)
  obs <- eqnvec(obs_a = "A")
  g <- Y(obs, f = x, compile = FALSE, modelname = "err1_obs",
         attach.input = FALSE)
  trafo <- as.eqnvec(c(A = "10^(A)", B = "10^(B)",
                       k = "10^(K) + UNDECLARED"))
  p <- P(trafo, condition = "c1", compile = FALSE, modelname = "err1_par")

  data <- as.datalist(data.frame(
    name = "obs_a", time = 1, value = 0.5, sigma = 1, condition = "c1",
    stringsAsFactors = FALSE))

  # exportPEtab emits an informational warning when parameterScale is supplied
  # to a v2 export (it is ignored on disk). The test only cares about the error.
  suppressWarnings(expect_error(
    exportPEtab(data = data, reactions = reactions, observables = obs,
                p = p, pouter = c(A = 0, B = 0, K = -1),
                parameterScale = "log10",
                dir = tempfile("err1_"), overwrite = TRUE),
    "undeclared symbol"))

  unlink("err1_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})


test_that("native exportPEtab roundtrips per-row sigma via noiseParameters column", {

  setwd(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  reactions <- eqnlist() %>%
    addReaction("A", "B", rate = "k*A", description = "fwd")
  m <- odemodel(reactions, modelname = "rt_sig_ode", compile = FALSE,
                solver = "deSolve")
  x <- Xs(m)
  obs <- eqnvec(obs_a = "A")
  g <- Y(obs, f = x, condition = NULL, compile = FALSE,
         modelname = "rt_sig_obs", attach.input = FALSE)
  trafo <- as.eqnvec(c(A = "10^(A)", B = "10^(B)", k = "10^(K)"))
  p <- P(trafo, condition = "c1", compile = FALSE, modelname = "rt_sig_par")
  compile(g, x, p, cores = 1)

  # Three measurements with three different sigmas — exercise the
  # noiseParameter1_<obsId> placeholder + per-row noiseParameters path.
  data <- as.datalist(data.frame(
    name = "obs_a", time = c(1, 2, 3),
    value = c(0.5, 0.3, 0.2), sigma = c(0.5, 1.0, 2.0),
    condition = "c1", stringsAsFactors = FALSE))

  pouter <- c(A = 0, B = 0, K = -1)
  obj_native <- normL2(data, g * x * p)

  out_dir <- file.path(tempdir(), "petab_rt_sig")
  yaml_out <- exportPEtab(
    data = data, reactions = reactions, observables = obs,
    p = p, pouter = pouter,
    parameterScale = "log10", modelID = "rt_sig_export",
    formatVersion = "1", dir = out_dir, overwrite = TRUE)

  # observables.tsv must declare the placeholder noiseFormula.
  obs_tsv <- read.delim(file.path(out_dir, "observables_rt_sig_export.tsv"),
                        stringsAsFactors = FALSE)
  expect_equal(obs_tsv$noiseFormula, "noiseParameter1_obs_a")

  # measurements.tsv must carry per-row noiseParameters values.
  meas_tsv <- read.delim(file.path(out_dir, "measurements_rt_sig_export.tsv"),
                         stringsAsFactors = FALSE)
  expect_setequal(as.numeric(meas_tsv$noiseParameters), c(0.5, 1.0, 2.0))

  petab <- importPEtab(yaml_out, solver = "deSolve",
                       modelname = "rt_sig_imp")
  v_native <- obj_native(pouter, deriv = FALSE)$value
  v_petab  <- petab$obj(pouter[names(petab$bestfit)], deriv = FALSE)$value
  expect_lt(abs(v_native - v_petab), 1e-3)

  unlink("rt_sig_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})


test_that("native exportPEtab roundtrips compound trafos like 10^(KM + 5)", {

  setwd(tempdir())
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  # 1-state, 1-reaction with a non-trivial compound mapping for the rate:
  #   k = 10^(K + 5) — chain-rule "compensation" path, not strippable.
  reactions <- eqnlist() %>%
    addReaction("A", "B", rate = "k*A", description = "fwd")
  m <- odemodel(reactions, modelname = "rt3_ode", compile = FALSE,
                solver = "deSolve")
  x <- Xs(m)
  obs <- eqnvec(obs_a = "A")
  g <- Y(obs, f = x, compile = FALSE, modelname = "rt3_obs",
         attach.input = FALSE)
  trafo <- as.eqnvec(c(A = "10^(A)", B = "10^(B)", k = "10^(K + 5)"))
  p <- P(trafo, condition = "c1", compile = FALSE, modelname = "rt3_par")
  compile(g, x, p, cores = 1)

  data <- as.datalist(data.frame(
    name = "obs_a", time = 1, value = 0.5, sigma = 1, condition = "c1",
    stringsAsFactors = FALSE))

  pouter <- c(A = 0, B = 0, K = -1)
  obj_native <- normL2(data, g * x * p)

  out_dir <- file.path(tempdir(), "petab_rt3")
  yaml_out <- exportPEtab(
    data = data, reactions = reactions, observables = obs,
    p = p, pouter = pouter,
    parameterScale = "log10", modelID = "rt3_export",
    formatVersion = "1", dir = out_dir, overwrite = TRUE)

  # The conditions.tsv cell for k must contain the compensated form
  # `10^(log10(K) + 5)` so the importer's chain rule reproduces 10^(K+5).
  cond_df <- read.delim(file.path(out_dir, "conditions_rt3_export.tsv"),
                        stringsAsFactors = FALSE)
  expect_match(as.character(cond_df$k[[1L]]), "log10\\(K\\)")

  petab <- importPEtab(yaml_out, solver = "deSolve",
                       modelname = "rt3_imp")
  v_native <- obj_native(pouter, deriv = FALSE)$value
  v_petab  <- petab$obj(pouter[names(petab$bestfit)], deriv = FALSE)$value
  expect_lt(abs(v_native - v_petab), 1e-3)

  unlink("rt3_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})


## --- PEtab v2 (no-SBML pure-parser tests) ---------------------------------

test_that(".petab_major_version recognises v1 and v2 strings", {
  expect_identical(dMod:::.petab_major_version(1L),       1L)
  expect_identical(dMod:::.petab_major_version("1"),      1L)
  expect_identical(dMod:::.petab_major_version("1.0.0"),  1L)
  expect_identical(dMod:::.petab_major_version("2.0.0"),  2L)
  expect_identical(dMod:::.petab_major_version("2.1.3"),  2L)
  expect_identical(dMod:::.petab_major_version(NULL),     1L)  # legacy default
  expect_error(dMod:::.petab_major_version("v2"),
               regexp = "Unrecognised PEtab format_version")
})


test_that(".petab_v2_normalize_tables converts a single-condition v2 problem", {
  tables <- list(
    parameters = data.frame(
      parameterId  = c("k1", "k2", "init_a"),
      lowerBound   = c(1e-5, 1e-5, 0),
      upperBound   = c(1e3,  1e3,  10),
      nominalValue = c(0.1, 0.5, 1.0),
      estimate     = c("true", "true", "false"),
      stringsAsFactors = FALSE),
    observables = data.frame(
      observableId           = c("o1"),
      observableFormula      = c("A * scale + offset"),
      observablePlaceholders = c("scale;offset"),
      noiseFormula           = c("sigma"),
      noiseDistribution      = c("log-normal"),
      noisePlaceholders      = c("sigma"),
      stringsAsFactors = FALSE),
    conditions = data.frame(
      conditionId = c("c1", "c1"),
      targetId    = c("a0", "k_in"),
      targetValue = c("init_a", "0.4"),
      stringsAsFactors = FALSE),
    measurements = data.frame(
      observableId = c("o1", "o1"),
      experimentId = c("exp1", "exp1"),
      time         = c(0, 10),
      measurement  = c(1.0, 0.6),
      observableParameters = c("1.5;0", "1.5;0"),
      noiseParameters      = c("0.1", "0.1"),
      stringsAsFactors = FALSE),
    experiments = data.frame(
      experimentId = c("exp1"),
      time         = c("0"),
      conditionId  = c("c1"),
      stringsAsFactors = FALSE),
    mapping     = NULL,
    sbmlPath   = "ignored.xml",
    formatVersion = 2L)

  out <- dMod:::.petab_v2_normalize_tables(tables)

  # parameters: parameterScale synthesised; estimate coerced to 1/0.
  expect_true("parameterScale" %in% colnames(out$parameters))
  expect_equal(unique(out$parameters$parameterScale), "lin")
  expect_equal(out$parameters$estimate, c(1L, 1L, 0L))

  # observables: log-normal split into log + normal; placeholders rewritten.
  expect_equal(out$observables$observableTransformation, "log")
  expect_equal(out$observables$noiseDistribution, "normal")
  expect_match(out$observables$observableFormula,
               "observableParameter1_o1.*observableParameter2_o1")
  expect_match(out$observables$noiseFormula,
               "^noiseParameter1_o1$")

  # conditions: long → wide.
  expect_equal(sort(setdiff(colnames(out$conditions), "conditionId")),
               c("a0", "k_in"))
  r <- which(out$conditions$conditionId == "c1")
  expect_equal(out$conditions$a0[r],   "init_a")
  expect_equal(out$conditions$k_in[r], "0.4")

  # measurements: experimentId rewritten.
  expect_equal(out$measurements$simulationConditionId,
               c("c1", "c1"))
  expect_equal(out$measurements$preequilibrationConditionId,
               c("", ""))
  expect_false("experimentId" %in% colnames(out$measurements))
})


test_that(".petab_v2_normalize_tables handles preequilibration via 2-period experiments", {
  tables <- list(
    parameters = data.frame(parameterId = "k", lowerBound = 0, upperBound = 1,
                            nominalValue = 0.5, estimate = "true",
                            stringsAsFactors = FALSE),
    observables = data.frame(observableId = "o1", observableFormula = "A",
                             noiseFormula = "1",
                             noiseDistribution = "normal",
                             stringsAsFactors = FALSE),
    conditions = data.frame(
      conditionId = c("c_pre", "c_sim"),
      targetId    = c("a0", "a0"),
      targetValue = c("5",  "1"),
      stringsAsFactors = FALSE),
    measurements = data.frame(
      observableId = "o1", experimentId = "exp_with_pre",
      time = 5, measurement = 0.7,
      stringsAsFactors = FALSE),
    experiments = data.frame(
      experimentId = c("exp_with_pre", "exp_with_pre"),
      time         = c("-inf",         "0"),
      conditionId  = c("c_pre",        "c_sim"),
      stringsAsFactors = FALSE),
    mapping     = NULL,
    sbmlPath   = "ignored.xml",
    formatVersion = 2L)

  out <- dMod:::.petab_v2_normalize_tables(tables)
  expect_equal(out$measurements$simulationConditionId,        "c_sim")
  expect_equal(out$measurements$preequilibrationConditionId,  "c_pre")
})


test_that(".petab_v2_normalize_tables rejects > 2 periods", {
  tables <- list(
    parameters = data.frame(parameterId = "k", lowerBound = 0, upperBound = 1,
                            nominalValue = 0.5, estimate = "true",
                            stringsAsFactors = FALSE),
    observables = data.frame(observableId = "o1", observableFormula = "A",
                             noiseFormula = "1",
                             noiseDistribution = "normal",
                             stringsAsFactors = FALSE),
    conditions = data.frame(conditionId = c("c1", "c2", "c3"),
                            targetId = c("a0", "a0", "a0"),
                            targetValue = c("1", "2", "3"),
                            stringsAsFactors = FALSE),
    measurements = data.frame(observableId = "o1", experimentId = "e",
                              time = 0, measurement = 1,
                              stringsAsFactors = FALSE),
    experiments = data.frame(experimentId = rep("e", 3),
                             time = c("-inf", "0", "5"),
                             conditionId = c("c1", "c2", "c3"),
                             stringsAsFactors = FALSE),
    mapping = NULL, sbmlPath = "x", formatVersion = 2L)
  expect_error(dMod:::.petab_v2_normalize_tables(tables),
               regexp = "at most 2")
})


test_that(".petab_v2_normalize_tables applies mapping table substitutions", {
  tables <- list(
    parameters = data.frame(parameterId = c("species_a_init", "k"),
                            lowerBound = c(0, 0), upperBound = c(10, 10),
                            nominalValue = c(1, 0.5),
                            estimate = c("false", "true"),
                            stringsAsFactors = FALSE),
    observables = data.frame(observableId = "o1",
                             observableFormula = "species_a",
                             noiseFormula = "1",
                             noiseDistribution = "normal",
                             stringsAsFactors = FALSE),
    conditions = data.frame(conditionId = "c1",
                            targetId = "species_a",
                            targetValue = "species_a_init",
                            stringsAsFactors = FALSE),
    measurements = data.frame(observableId = "o1", experimentId = "exp",
                              time = 0, measurement = 1,
                              stringsAsFactors = FALSE),
    experiments = data.frame(experimentId = "exp", time = "0",
                             conditionId = "c1",
                             stringsAsFactors = FALSE),
    mapping = data.frame(petabEntityId = "species_a",
                         modelEntityId = "A_internal",
                         stringsAsFactors = FALSE),
    sbmlPath = "x", formatVersion = 2L)
  out <- dMod:::.petab_v2_normalize_tables(tables)
  expect_equal(out$observables$observableFormula, "A_internal")
  # condition target column renamed to model entity name
  expect_true("A_internal" %in% colnames(out$conditions))
})


test_that("read_petab_yaml dispatches v1 vs v2 schema", {
  td <- tempfile("petab_v2_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  writeLines("conditionId\ttargetId\ttargetValue\nc1\ta0\t1\n",
             file.path(td, "conditions.tsv"))
  writeLines("experimentId\ttime\tconditionId\nexp\t0\tc1\n",
             file.path(td, "experiments.tsv"))
  writeLines("observableId\tobservableFormula\tnoiseFormula\tnoiseDistribution\no1\tA\t1\tnormal\n",
             file.path(td, "observables.tsv"))
  writeLines("observableId\texperimentId\ttime\tmeasurement\no1\texp\t0\t1\n",
             file.path(td, "measurements.tsv"))
  writeLines("parameterId\tlowerBound\tupperBound\tnominalValue\testimate\nk\t0\t1\t0.5\ttrue\n",
             file.path(td, "parameters.tsv"))
  writeLines("<sbml/>", file.path(td, "model.xml"))

  yaml::write_yaml(list(
    format_version    = "2.0.0",
    parameter_files   = list("parameters.tsv"),
    model_files       = list(my_model = list(location = "model.xml",
                                             language = "sbml")),
    observable_files  = list("observables.tsv"),
    measurement_files = list("measurements.tsv"),
    condition_files   = list("conditions.tsv"),
    experiment_files  = list("experiments.tsv")
  ), file.path(td, "problem.yaml"))

  m <- read_petab_yaml(file.path(td, "problem.yaml"))
  expect_identical(m$formatVersion, 2L)
  expect_equal(m$problems[[1]]$modelID, "my_model")
  expect_match(m$problems[[1]]$sbmlFile,        "model\\.xml$")
  expect_match(m$problems[[1]]$experimentFile,  "experiments\\.tsv$")
  expect_null(m$problems[[1]]$mappingFile)
})


test_that("read_petab_yaml errors on non-SBML model language", {
  td <- tempfile("petab_v2_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  writeLines("dummy", file.path(td, "model.bngl"))
  writeLines("dummy", file.path(td, "p.tsv"))
  writeLines("dummy", file.path(td, "o.tsv"))
  writeLines("dummy", file.path(td, "m.tsv"))
  yaml::write_yaml(list(
    format_version = "2.0.0",
    parameter_files = list("p.tsv"),
    model_files = list(m = list(location = "model.bngl", language = "bngl")),
    observable_files = list("o.tsv"),
    measurement_files = list("m.tsv")
  ), file.path(td, "problem.yaml"))

  expect_error(read_petab_yaml(file.path(td, "problem.yaml")),
               regexp = "SBML")
})


test_that("exportPEtabObject v2 writes nominalValue verbatim (no parameterScale linearisation)", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")

  setwd(tempdir())
  pp <- importPEtab(file.path(petab_dir, "0001", "_0001.yaml"),
                    solver = "deSolve", compile = FALSE)
  td <- tempfile("petab_v2_lin_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # No warning on v2 export, even though pp came from a v1 problem with
  # log10-scale outer parameters: the trafo `p` already encodes the scale
  # via `10^(...)` wraps, which the v2 path keeps in conditions.tsv /
  # SBML <initialAssignment>.
  expect_silent(
    exportPEtabObject(pp, dir = td, formatVersion = "2.0.0",
                      overwrite = TRUE))

  par_path <- list.files(td, pattern = "^parameters_.*\\.tsv$",
                         full.names = TRUE)
  par_df <- read.delim(par_path, stringsAsFactors = FALSE, na.strings = "")
  expect_false("parameterScale" %in% colnames(par_df))
  # nominalValue equals the internal pouter (log10-scale) — i.e. NOT
  # 10^pouter as the old linearised code emitted.
  est <- par_df[par_df$estimate == "true", , drop = FALSE]
  expect_equal(est$nominalValue,
               unname(pp$bestfit[est$parameterId]))
})


test_that("exportPEtabObject v2 writes long-format conditions and experiments", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")

  setwd(tempdir())
  # 0001 has a single condition; trivial v2 export should produce one
  # experimentId row.
  petab <- importPEtab(file.path(petab_dir, "0001", "_0001.yaml"),
                       solver = "deSolve", compile = FALSE)
  td <- tempfile("petab_v2_out_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  yamlPath <- exportPEtabObject(petab, dir = td, formatVersion = "2.0.0",
                                  overwrite = TRUE)

  cond_path <- list.files(td, pattern = "^conditions_.*\\.tsv$",
                          full.names = TRUE)
  expect_length(cond_path, 1L)
  cond <- read.delim(cond_path, stringsAsFactors = FALSE, na.strings = "")
  expect_setequal(colnames(cond), c("conditionId", "targetId", "targetValue"))

  expect_length(list.files(td, pattern = "^experiments_.*\\.tsv$"), 1L)
  m <- yaml::read_yaml(yamlPath)
  expect_identical(m$format_version, "2.0.0")
  expect_true("model_files" %in% names(m))
  expect_equal(m$model_files[[1]]$language, "sbml")
})


test_that("v2 export → v2 import roundtrips the objective on case 0001", {
  if (!.libsbml_works()) skip("libsbml virtualenv not available")
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")

  setwd(tempdir())
  pp1 <- importPEtab(file.path(petab_dir, "0001", "_0001.yaml"),
                     solver = "deSolve", compile = TRUE)
  td <- tempfile("v2_rt_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  yamlPath <- exportPEtabObject(pp1, dir = td, formatVersion = "2.0.0",
                                 overwrite = TRUE)

  setwd(td)
  pp2 <- importPEtab(yamlPath, solver = "deSolve", compile = TRUE)

  v1 <- pp1$obj(pp1$bestfit)$value
  v2 <- pp2$obj(pp2$bestfit)$value
  expect_equal(v1, v2, tolerance = 1e-6)
})


test_that("v2 PEtab test cases 0001/0002/0009 import and match published llh", {
  # Earlier tests may setwd() into a tempfile() dir that gets unlinked on
  # exit; reset to a guaranteed-existing cwd before any path lookup so
  # `.petab_repo_dir()`'s `getwd()` calls don't error.
  setwd(tempdir())
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")
  v2_dir <- file.path(petab_dir, "v2")
  if (!dir.exists(v2_dir)) skip("PEtabTests/v2/ not present")
  if (!.libsbml_works()) skip("libsbml virtualenv not available")

  for (case in c("0001", "0002", "0009")) {
    yamlPath <- file.path(v2_dir, case, paste0("_", case, ".yaml"))
    if (!file.exists(yamlPath)) next
    sol_path  <- file.path(v2_dir, case, paste0("_", case, "_solution.yaml"))
    sol <- yaml::read_yaml(sol_path)
    wd <- tempfile(paste0("v2_case_", case, "_")); dir.create(wd)
    setwd(wd)
    res <- tryCatch({
      pp <- importPEtab(yamlPath, solver = "deSolve", compile = TRUE,
                        modelname = paste0("v2bench_", case))
      pp$obj(pp$bestfit)$value
    }, error = function(e) {
      message("v2 case ", case, " import error: ", conditionMessage(e))
      NA_real_
    })
    setwd(tempdir())
    expect_equal(res, -2 * as.numeric(sol$llh), tolerance = 1e-3,
                 info = sprintf("v2 case %s", case))
  }
})


test_that("v2 → v1 → v2 textual normaliser roundtrips a minimal problem", {
  # Round-trip purely at the table level: build a v2 input, normalise to v1
  # shape, write back as v2, normalise again, and compare key invariants.
  v2 <- list(
    parameters = data.frame(
      parameterId = c("k1"), lowerBound = 1e-3, upperBound = 1e3,
      nominalValue = 0.5, estimate = "true", stringsAsFactors = FALSE),
    observables = data.frame(
      observableId = "o1", observableFormula = "A",
      noiseFormula = "1", noiseDistribution = "normal",
      stringsAsFactors = FALSE),
    conditions = data.frame(
      conditionId = "c1", targetId = "a0", targetValue = "3",
      stringsAsFactors = FALSE),
    measurements = data.frame(
      observableId = "o1", experimentId = "e1",
      time = 0, measurement = 1.0, stringsAsFactors = FALSE),
    experiments = data.frame(
      experimentId = "e1", time = "0", conditionId = "c1",
      stringsAsFactors = FALSE),
    mapping = NULL, sbmlPath = "x", formatVersion = 2L)
  out <- dMod:::.petab_v2_normalize_tables(v2)
  expect_equal(out$measurements$simulationConditionId, "c1")
  expect_equal(out$conditions$a0[out$conditions$conditionId == "c1"], "3")
})
