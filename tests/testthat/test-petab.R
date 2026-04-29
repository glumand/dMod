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

# AMICI is a heavy GitHub-only dependency. Skip integration tests that need
# it if we cannot run import_sbml() at all.
.amici_works <- function() {
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
  expect_equal(y$format_version, 1L)
  expect_true(file.exists(y$problems[[1]]$sbml_file))
  expect_true(file.exists(y$problems[[1]]$measurement_file))
})


test_that("read_petab_tables returns 4 data frames", {
  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")

  tabs <- read_petab_tables(file.path(petab_dir, "0001", "_0001.yaml"))
  expect_named(tabs, c("parameters", "conditions", "measurements",
                       "observables", "sbml_path"))
  expect_s3_class(tabs$parameters,   "data.frame")
  expect_s3_class(tabs$conditions,   "data.frame")
  expect_s3_class(tabs$measurements, "data.frame")
  expect_s3_class(tabs$observables,  "data.frame")
})


## --- end-to-end fixture test (no AMICI required) --------------------------
##
## We hand-build the eqnlist that matches PEtab test case 0001's SBML model
## and verify the trafo+objective machinery against the published solution.
## This avoids an AMICI dependency on every test run.

test_that("hand-built case-0001 fixture produces solution-matching llh", {

  setwd(tempdir())

  # Reaction network identical to PEtabTests/0001/_model.xml after AMICI
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


## --- AMICI-dependent integration tests ------------------------------------

test_that("PEtab test cases 0001-0006 import and produce solution-matching llh", {

  setwd(tempdir())

  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")
  if (!.amici_works())   skip("AMICI / Python virtualenv not available")

  for (id in sprintf("%04d", 1:6)) {
    yaml_path <- file.path(petab_dir, id, paste0("_", id, ".yaml"))
    sol_path  <- file.path(petab_dir, id, paste0("_", id, "_solution.yaml"))
    if (!file.exists(sol_path)) next

    petab <- importPEtab(yaml_path, solver = "deSolve",
                         modelname = paste0("petab_", id))
    sol <- yaml::read_yaml(sol_path)

    out <- petab$obj(petab$pouter, fixed = petab$fixed, deriv = FALSE)
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
  if (!.amici_works())   skip("AMICI / Python virtualenv not available")

  # Cases 0007 (log10 trafo), 0008 (replicates), 0009/0010 (preequilibration —
  # numeric Pequil fallback is exercised because steadyStates() leaves one
  # state unresolved on the A↔B reaction), 0011-0013 (init / compartment /
  # parametric init overrides), 0014/0015 (numeric / symbolic noise parameter
  # overrides), 0016 (log trafo).
  for (id in sprintf("%04d", 7:16)) {

    yaml_path <- file.path(petab_dir, id, paste0("_", id, ".yaml"))
    sol_path  <- file.path(petab_dir, id, paste0("_", id, "_solution.yaml"))
    if (!file.exists(sol_path)) next

    suppressWarnings(
      petab <- importPEtab(yaml_path, solver = "deSolve",
                           modelname = paste0("petab_s2_", id)))
    sol <- yaml::read_yaml(sol_path)
    out <- petab$obj(petab$pouter, fixed = petab$fixed, deriv = FALSE)
    expect_lt(abs(out$value - (-2 * sol$llh)),
              max(0.01, abs(2 * sol$tol_llh)),
              label = paste0("case ", id, " -2*llh"))
  }

  unlink("petab_s2_*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so"); unlink("*_model.csv")
  unlink("reactions_for_Alyssa*")
})


test_that("preeqMethod = 'numeric' forces Pequil even when analytic would work", {

  setwd(tempdir())

  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")
  if (!.amici_works())   skip("AMICI / Python virtualenv not available")

  # Case 0009 has no symbolic conserved-quantity hint, so analytic falls back
  # automatically. The explicit "numeric" override should still match.
  yaml_path <- file.path(petab_dir, "0009", "_0009.yaml")
  sol_path  <- file.path(petab_dir, "0009", "_0009_solution.yaml")
  petab <- importPEtab(yaml_path, solver = "deSolve",
                       modelname = "petab_0009_num",
                       preeqMethod = "numeric")
  sol <- yaml::read_yaml(sol_path)
  out <- petab$obj(petab$pouter, fixed = petab$fixed, deriv = FALSE)
  expect_lt(abs(out$value - (-2 * sol$llh)),
            max(0.01, abs(2 * sol$tol_llh)))

  unlink("petab_0009_num*"); unlink("*.c"); unlink("*.cpp")
  unlink("*.o"); unlink("*.so")
})


test_that("two-condition roundtrip preserves objective value", {

  setwd(tempdir())

  petab_dir <- .petab_repo_dir()
  if (!nzchar(petab_dir)) skip("PEtabTests/ directory not found")
  if (!.amici_works())   skip("AMICI / Python virtualenv not available")

  # Case 0002 has two conditions, an InitialAssignment binding A := a0 / B := b0
  # and exercises the condition table. The InitialAssignment roundtrip is the
  # interesting part: without it the reimported objective evaluates with
  # state initials = 0 and disagrees with the original.
  yaml1 <- file.path(petab_dir, "0002", "_0002.yaml")
  petab1 <- importPEtab(yaml1, solver = "deSolve",
                        modelname = "rt_in")
  v1 <- petab1$obj(petab1$pouter, fixed = petab1$fixed, deriv = FALSE)$value

  out_dir <- file.path(tempdir(), "petab_roundtrip")
  yaml2 <- exportPEtabObject(petab1, out_dir, model_id = "rt_out", overwrite = TRUE)

  petab2 <- importPEtab(yaml2, solver = "deSolve",
                        modelname = "rt_back")
  v2 <- petab2$obj(petab2$pouter, fixed = petab2$fixed, deriv = FALSE)$value

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
  if (!.amici_works())  skip("AMICI / Python virtualenv not available")

  yaml_path <- file.path(bm_dir, "Boehm_JProteomeRes2014",
                         "Boehm_JProteomeRes2014.yaml")
  petab <- importPEtab(yaml_path, solver = "deSolve",
                       modelname = "boehm")

  # Imported problem shape:
  expect_equal(length(petab$pouter), 9L)
  expect_setequal(names(petab$observables),
                  c("pSTAT5A_rel", "pSTAT5B_rel", "rSTAT5A_rel"))
  # All estimated parameters are on log10 scale per parameters.tsv:
  scales <- attr(petab$pouter, "petab_scales")
  expect_true(all(scales == "log10"))
  # AssignmentRule for BaF3_Epo must have been inlined → not in `fixed`:
  expect_false("BaF3_Epo" %in% names(petab$fixed))

  out <- petab$obj(petab$pouter, fixed = petab$fixed, deriv = FALSE)

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
  if (!.amici_works()) skip("Python / libsbml virtualenv not available")

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
              filepath = out_xml, model_id = "ia_export")

  xml_text <- readLines(out_xml, warn = FALSE)
  expect_true(any(grepl("<initialAssignment", xml_text, fixed = TRUE)),
              info = "no <initialAssignment> emitted for symbolic init")
  expect_true(any(grepl("symbol=\"A\"", xml_text)),
              info = "InitialAssignment for A missing")

  unlink(out_xml)
})
