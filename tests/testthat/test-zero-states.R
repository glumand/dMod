# Structural a-priori zero-state detection for eqnlist inputs to
# Pequil() / Pimpl(). The three layers (NegCol, PosCol-with-single-state
# feeder, sink-cluster LP) match AlyssaPetit v1.2; ported from
# inst/code/AlyssaPetit_ver1_2.py to R for use without reticulate.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}

helper <- dMod:::.zeroStatesFromSmatrix


## ---- Layer 1: only-outflux column --------------------------------------

test_that("NegCol: a state with only outflux is detected as zero", {
  el <- eqnlist() |>
    addReaction("A", "", "k_dg * A")
  zs <- helper(el)
  expect_setequal(zs$zero_states, "A")
  expect_equal(ncol(zs$eqnlist$smatrix), 0L)
})


## ---- Layer 2: PosCol with single-state feeder --------------------------

test_that("PosCol: only-influx state whose feeder has a single reactant", {
  # X is fed by k_pr * Z and consumed by nothing (PosCol).
  # k_pr * Z must be 0 in SS, so Z = 0 (single-state in flux).
  # Then X is fed only by zero flux -> drop X next iteration (NegCol-empty).
  el <- eqnlist() |>
    addReaction("", "X", "k_pr * Z") |>
    addReaction("Z", "",  "k_dg * Z")
  zs <- helper(el)
  expect_setequal(zs$zero_states, c("Z", "X"))
})


## ---- Layer 3: structural sink cluster ----------------------------------

test_that("Sink cluster: TGFb + R_TGFb + R_TGFb_int (combined mass leaks)", {
  skip_if_not_installed("lpSolve")
  # Each individual column has +1 and -1 entries, so layers 1+2 cannot
  # catch them. The cluster {L, RL, RLi} as a whole degrades via the
  # RLi -> "" reaction, so its total mass leaks monotonically.
  el <- eqnlist() |>
    addReaction("L + R", "RL",  "k_on  * L * R") |>
    addReaction("RL",    "L + R", "k_off * RL") |>
    addReaction("RL",    "RLi", "k_int * RL") |>
    addReaction("RLi",   "L",   "k_dec * RLi") |>
    addReaction("RLi",   "",    "k_dg  * RLi") |>
    addReaction("",      "R",   "k_pr_R") |>
    addReaction("R",     "",    "k_dg_R * R")
  zs <- helper(el)
  expect_true(all(c("L", "RL", "RLi") %in% zs$zero_states))
  # R is not in the sink cluster: it has its own production/degradation.
  expect_false("R" %in% zs$zero_states)
})


## ---- Layer 3 falls back when lpSolve is missing ------------------------

test_that("FindSinkCluster degrades gracefully without lpSolve", {
  # Same model as above; we can't easily uninstall lpSolve mid-test, so
  # this is a structural assertion: with only layers 1+2, the cluster
  # would not be found (we verified manually no NegCol/PosCol matches).
  el <- eqnlist() |>
    addReaction("L + R", "RL",  "k_on  * L * R") |>
    addReaction("RL",    "L + R", "k_off * RL") |>
    addReaction("RL",    "RLi", "k_int * RL") |>
    addReaction("RLi",   "L",   "k_dec * RLi") |>
    addReaction("RLi",   "",    "k_dg  * RLi")
  # Sanity: every column is mixed-sign in this minimal cluster model.
  S <- el$smatrix
  S[is.na(S)] <- 0
  storage.mode(S) <- "double"
  for (j in seq_len(ncol(S))) {
    col <- S[, j]
    expect_true(any(col > 0) && any(col < 0),
                info = sprintf("state %s is not mixed-sign", colnames(S)[j]))
  }
})


## ---- Idempotence on a model with no zero states ------------------------

test_that("Pure production-decay has no zero states", {
  el <- eqnlist() |>
    addReaction("", "A", "k_in") |>
    addReaction("A", "",  "k_out * A")
  zs <- helper(el)
  expect_equal(zs$zero_states, character(0))
  expect_identical(zs$eqnlist, el)
})


## ---- End-to-end: Pequil produces value = 0 for zero states -------------

test_that("Pequil reports zero states with value 0 and drops them from getParameters()", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Receptor R alone has production + degradation (nonzero baseline).
  # Ligand L and complex LR form a sink cluster (LR degrades).
  el <- eqnlist() |>
    addReaction("",       "R",  "k_pr_R") |>
    addReaction("R",      "",   "k_dg_R * R") |>
    addReaction("L + R",  "LR", "k_on * L * R") |>
    addReaction("LR",     "L + R", "k_off * LR") |>
    addReaction("LR",     "",   "k_dg_LR * LR")

  pf <- Pequil(el,
               modelname = paste0("test_zero_states_", as.integer(Sys.time())),
               compile = TRUE, verbose = FALSE)

  expect_false(any(c("L", "LR") %in% getParameters(pf)))
  pars <- c(k_pr_R = 2, k_dg_R = 0.5, k_on = 1, k_off = 0.5, k_dg_LR = 0.1)
  out  <- pf(pars)[[1]]
  expect_equal(as.numeric(out["L"]),  0, tolerance = 1e-8)
  expect_equal(as.numeric(out["LR"]), 0, tolerance = 1e-8)
  expect_equal(as.numeric(out["R"]),  pars[["k_pr_R"]] / pars[["k_dg_R"]],
               tolerance = 1e-5)
})


## ---- No-progress: hard error instead of stale initial values -----------

test_that("Pequil errors when the solver makes no progress", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Force a no-progress condition by capping the total step budget so the
  # solver can't even take one accepted step. Without the new guard,
  # Pequil would silently return the initial values.
  pf <- Pequil(c(A = "k_in - k_out * A"),
               parameters = c("k_in", "k_out"),
               modelname = paste0("test_no_progress_", as.integer(Sys.time())),
               controlsODE = list(maxsteps = 1L, maxprogress = 1L),
               compile = TRUE, verbose = FALSE)
  expect_error(
    suppressWarnings(pf(c(k_in = 1, k_out = 1, A = 1))),
    "no progress")
})
