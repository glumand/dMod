# Benchmark + regression harness for symmetryDetection().
#
# Times the modular observability engine on two fixtures drawn from
# inst/examples/symmetryDetection.R and fingerprints the result (rank, dim,
# identifiable, and the rendered non-identifiability directions) so an
# optimization can be shown to preserve the analysis exactly while getting
# faster. `bench/` is .Rbuildignore'd.
#
# Usage (from the repo root, env per memory local-check-env):
#   Rscript bench/bench_symmetryDetection.R            # EGF only, cores 1 & 4
#   DMOD_BENCH_SMAD=1 Rscript bench/bench_symmetryDetection.R   # + the slow SMAD model
#   DMOD_BENCH_CORES=1,4,8 Rscript bench/bench_symmetryDetection.R
#
# Or interactively:
#   source("bench/bench_symmetryDetection.R"); res <- bench_run()
#   saveRDS(res$fp, "before.rds")   # ... make changes ...  bench_compare("before.rds")

suppressMessages(devtools::load_all(".", quiet = TRUE))

# ---- fixtures -------------------------------------------------------------

egf_model <- function() {
  reactions <- eqnlist() |>
    addReaction("EGF + EGFR", "EGF_EGFR", "k_bind * EGF * EGFR") |>
    addReaction("EGF_EGFR", "EGF + EGFR", "k_unbind * EGF_EGFR") |>
    addReaction("MEK", "pMEK", "k_phos_MEK * EGF_EGFR * MEK") |>
    addReaction("pMEK", "MEK", "k_dephos_MEK * pMEK") |>
    addReaction("ERK", "pERK", "k_phos_ERK * pMEK * ERK") |>
    addReaction("pERK", "ERK", "k_dephos_ERK * pERK")
  reactions <- customTotals(reactions, list(
    totalEGF  = "EGF + EGF_EGFR", totalEGFR = "EGFR + EGF_EGFR",
    totalMEK  = "MEK + pMEK",     totalERK  = "ERK + pERK"))
  observables <- eqnvec(pMEK_obs = "scale_pMEK * pMEK",
                        pERK_obs = "scale_pERK * pERK")
  list(f = reactions, g = observables,
       args = list(method = "observability", reduceCQ = FALSE, closedForm = TRUE))
}

smad_model <- function() {
  addRC <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "Cell", ...)
  addRE <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "extraCell", ...)
  reactions <- eqnlist() |>
    addRC("", "bool_ActD",  "0") |> addRC("", "bool_CHX", "0") |> addRC("", "bool_MG132", "0") |>
    addRE("", "TGFb", "0") |>
    addRC("", "R1mRNA",  "k_pr_R1mRNA * (1 + k_inh_R1mRNA_FB3 * FB3^nhill_R1) * (1 - bool_ActD)") |>
    addRC("", "R2mRNA",  "k_pr_R2mRNA / (1 + k_inh_R2mRNA_FB4 * FB4^nhill_R2) * (1 - bool_ActD)") |>
    addRC("", "FB2mRNA", "k_pr_FB2mRNA * C3^nhill_FB2mRNA / (Km_FB2mRNA^nhill_FB2mRNA + C3^nhill_FB2mRNA) * (1 - bool_ActD)") |>
    addRC("", "FB3mRNA", "k_pr_FB3mRNA * C3^nhill_FB3mRNA / (Km_FB3mRNA^nhill_FB3mRNA + C3^nhill_FB3mRNA) * (1 - bool_ActD)") |>
    addRC("", "FB4mRNA", "k_pr_FB4mRNA * C3^nhill_FB4mRNA / (Km_FB4mRNA^nhill_FB4mRNA + C3^nhill_FB4mRNA) * (1 - bool_ActD)") |>
    addRC("R1mRNA", "", "k_dg_R1mRNA * R1mRNA") |> addRC("R2mRNA", "", "k_dg_R2mRNA * R2mRNA") |>
    addRC("FB2mRNA", "", "k_dg_FB2 * FB2mRNA") |> addRC("FB3mRNA", "", "k_dg_FB3 * FB3mRNA") |>
    addRC("FB4mRNA", "", "k_dg_FB4 * FB4mRNA") |>
    addRC("", "R1",  "k_pr_R1 * R1mRNA * (1 - bool_CHX)") |>
    addRC("", "R2",  "k_pr_R2 * R2mRNA * (1 - bool_CHX)") |>
    addRC("", "FB2", "k_pr_FB2 * FB2mRNA * (1 - bool_CHX)") |>
    addRC("", "FB3", "k_pr_FB3 * FB3mRNA * (1 - bool_CHX)") |>
    addRC("", "FB4", "k_pr_FB4 * FB4mRNA * (1 - bool_CHX)") |>
    addRC("R1", "", "k_dg_R1 * R1 * (1 - bool_MG132)") |> addRC("R2", "", "k_dg_R2 * R2 * (1 - bool_MG132)") |>
    addRC("FB2", "", "k_dg_FB2 * FB2 * (1 - bool_MG132)") |> addRC("FB3", "", "k_dg_FB3 * FB3 * (1 - bool_MG132)") |>
    addRC("FB4", "", "k_dg_FB4 * FB4 * (1 - bool_MG132)") |>
    addRC("R1 + R2", "R1_R2", "k_act_R1_R2 * R1 * R2") |>
    addRC("R1_R2", "R1 + R2", "k_deact_R1_R2 * R1_R2") |>
    addRC("R1_R2", "", "k_dg_R1_R2 * R1_R2 * (1 - bool_MG132)") |>
    addRC("R2 + TGFb", "R2_TGFb", "k_act_R2_TGFb * TGFb * R2 / (km_R2 + R2 + TGFb)", rateCompartment = "Cell") |>
    addRC("R2_TGFb", "R2 + TGFb", "k_deact_R2_TGFb * R2_TGFb") |>
    addRC("R2_TGFb", "R2_TGFb_int", "k_int_R2_TGFb * R2_TGFb") |>
    addRC("R2_TGFb_int", "TGFb", "k_decay_R2_TGFb_int * R2_TGFb_int") |>
    addRC("R2_TGFb_int", "", "k_dg_R2_TGFb_int * R2_TGFb_int * (1 - bool_MG132)") |>
    addRC("R1 + TGFb", "R1_TGFb", "k_act_R1_TGFb * TGFb * R1 / (km_R1 + R1 + TGFb)", rateCompartment = "Cell") |>
    addRC("R1_TGFb", "R1 + TGFb", "k_deact_R1_TGFb * R1_TGFb") |>
    addRC("R1_TGFb", "R1_TGFb_int", "k_int_R1_TGFb * R1_TGFb") |>
    addRC("R1_TGFb_int", "TGFb", "k_decay_R1_TGFb_int * R1_TGFb_int") |>
    addRC("R1_TGFb_int", "", "k_dg_R1_TGFb_int * R1_TGFb_int * (1 - bool_MG132)") |>
    addRC("R1 + R2_TGFb", "R1_R2_TGFb", "k_act_R1_R2_TGFb * R1 * R2_TGFb") |>
    addRC("R2 + R1_TGFb", "R1_R2_TGFb", "k_act_R2_R1_TGFb * R2 * R1_TGFb") |>
    addRC("R1_R2 + TGFb", "R1_R2_TGFb", "k_act_R1_R2_TGFb_direct * TGFb * R1_R2 / (km_R1_R2 + R1_R2 + TGFb)", rateCompartment = "Cell") |>
    addRC("R1_R2_TGFb", "", "(k_dg_R1_R2_TGFb + k_dg_R1_R2_TGFb_FB1 * FB2) * R1_R2_TGFb * (1 - bool_MG132)") |>
    addRC("Smad2", "pSmad2", "(k_phospho_pS2 / (1 + k_inh_pSmad2_FB2 * FB2)) * Smad2 * R1_R2") |>
    addRC("Smad2", "pSmad2", "(k_phospho_pS2 / (1 + k_inh_pSmad2_FB1 * FB2)) * Smad2 * R1_R2_TGFb") |>
    addRC("pSmad2", "Smad2", "k_dephos_S2 * pSmad2") |>
    addRC("Smad3", "pSmad3", "(k_phospho_pS3 / (1 + k_inh_pSmad3_FB2 * FB2)) * Smad3 * R1_R2") |>
    addRC("Smad3", "pSmad3", "(k_phospho_pS3 / (1 + k_inh_pSmad3_FB2 * FB2)) * Smad3 * R1_R2_TGFb") |>
    addRC("pSmad3", "Smad3", "k_dephos_S3 * pSmad3") |>
    addRC("pSmad2 + pSmad3 + Smad4", "C3", "k_form_S4Coip * pSmad2 * pSmad3 * Smad4") |>
    addRC("C3", "Smad2 + pSmad3 + Smad4", "k_dissolve_C3_dp2 * C3") |>
    addRC("C3", "pSmad2 + Smad3 + Smad4", "k_dissolve_C3_dp3 * C3")
  reactions$compartments$Cell$volume      <- "1"
  reactions$compartments$extraCell$volume <- "volumeEC"
  reactions <- customTotals(reactions, list(totalSMAD2 = "Smad2 + pSmad2 + C3",
                                            totalSMAD3 = "Smad3 + pSmad3 + C3",
                                            totalSMAD4 = "Smad4 + C3"))
  observables <- eqnvec(
    R1_obs = "scale_R1 * R1", R2_obs = "scale_R2 * R2",
    pSmad2_obs = "scale_pSmad2 * (pSmad2 + C3)", pSmad3_obs = "scale_pSmad3 * (pSmad3 + C3)",
    TSmad2_obs = "scale_TSmad2 * (Smad2 + pSmad2 + C3)", TSmad3_obs = "scale_TSmad3 * (Smad3 + pSmad3 + C3)",
    Smad4_CoIP_obs = "scale_CoIP * C3", TGFBR1_mRNA_obs = "scale_R1mRNA * R1mRNA",
    TGFBR2_mRNA_obs = "scale_R2mRNA * R2mRNA", TGFb_obs = "scale_TGFb * TGFb")
  events <- eventlist() |>
    addEvent(var = "TGFb",       time = 0,   value = "init_TGFb",      method = "replace") |>
    addEvent(var = "bool_CHX",   time = -30, value = "var_bool_CHX",   method = "replace") |>
    addEvent(var = "bool_MG132", time = -30, value = "var_bool_MG132", method = "replace") |>
    addEvent(var = "bool_ActD",  time = -30, value = "var_bool_ActD",  method = "replace")
  cond.grid <- data.frame(Pertubation = c("Ctrl", "ActD", "CHX", "MG132", "R1Knd", "R2Knd"),
                          init_TGFb = 1, stringsAsFactors = FALSE)
  cond.grid$var_bool_ActD  <- ifelse(cond.grid$Pertubation == "ActD",  1, 0)
  cond.grid$var_bool_CHX   <- ifelse(cond.grid$Pertubation == "CHX",   1, 0)
  cond.grid$var_bool_MG132 <- ifelse(cond.grid$Pertubation == "MG132", 1, 0)
  cond.grid$k_pr_R1mRNA <- ifelse(cond.grid$Pertubation == "R1Knd", "k_pr_R1mRNA_R1Knd", "k_pr_R1mRNA")
  cond.grid$k_pr_R2mRNA <- ifelse(cond.grid$Pertubation == "R2Knd", "k_pr_R2mRNA_R2Knd", "k_pr_R2mRNA")
  rownames(cond.grid) <- cond.grid$Pertubation
  cond.grid$Pertubation <- NULL
  # define()/branch() re-evaluate their `x=` argument by name up to the global
  # environment (they work at top level in the example); expose the parameter set
  # there for the duration of the trafo construction, then remove it.
  assign(".smad_pars", getParameters(reactions, events), envir = globalenv())
  on.exit(suppressWarnings(rm(".smad_pars", envir = globalenv())), add = TRUE)
  cond.trafo <- eqnvec() |>
    define("x~x", x = .smad_pars) |>
    branch(table = cond.grid, apply = "insert")
  list(f = reactions, g = observables,
       args = list(method = "observability", events = events, trafo = cond.trafo,
                   forcings = c("bool_ActD","bool_CHX","bool_MG132","TGFb"),
                   equilibrate = TRUE, reduceCQ = TRUE, closedForm = TRUE))
}

# ---- fingerprint (result identity, for regression) ------------------------

fingerprint <- function(o) {
  lines <- tryCatch(
    vapply(o$nonIdentifiable, function(d) dMod:::.sym_direction_line(d), character(1)),
    error = function(e) character(0))
  list(rank = o$rank, dim = o$dim, identifiable = isTRUE(o$identifiable),
       nDir = length(o$nonIdentifiable), lines = sort(lines))
}

# ---- runner ---------------------------------------------------------------

bench_one <- function(name, model, cores) {
  args <- c(list(f = model$f, g = model$g), model$args,
            list(cores.conditions = cores, cores.GLp = cores))
  gc()
  t <- system.time(o <- do.call(symmetryDetection, args))["elapsed"]
  cat(sprintf("  %-6s cores=%d  %8.2fs  rank %d/%d  ident=%s  #dir=%d\n",
              name, cores, t, o$rank, o$dim, isTRUE(o$identifiable),
              length(o$nonIdentifiable)))
  list(name = name, cores = cores, elapsed = as.numeric(t), fp = fingerprint(o))
}

bench_run <- function() {
  coresVec <- as.integer(strsplit(Sys.getenv("DMOD_BENCH_CORES", "1,4"), ",")[[1]])
  runSmad  <- nzchar(Sys.getenv("DMOD_BENCH_SMAD"))
  results <- list()
  cat("== EGF/EGFR -> MEK/ERK (Section 8, closedForm) ==\n")
  egf <- egf_model()
  for (cc in coresVec) results[[length(results) + 1L]] <- bench_one("egf", egf, cc)
  if (runSmad) {
    cat("== TGF-beta / SMAD (Section 10, equilibrate, closedForm) ==\n")
    smad <- smad_model()
    for (cc in coresVec) results[[length(results) + 1L]] <- bench_one("smad", smad, cc)
  } else {
    cat("(SMAD skipped; set DMOD_BENCH_SMAD=1 to include the large fixture)\n")
  }
  fp <- lapply(results, function(r) r$fp)
  names(fp) <- vapply(results, function(r) paste0(r$name, "_c", r$cores), character(1))
  invisible(list(results = results, fp = fp))
}

# Compare current fingerprints against a saved baseline; stops on any mismatch.
bench_compare <- function(baselinePath) {
  before <- readRDS(baselinePath)
  now <- bench_run()$fp
  keys <- union(names(before), names(now))
  ok <- TRUE
  for (k in keys) {
    b <- before[[k]]; n <- now[[k]]
    same <- !is.null(b) && !is.null(n) &&
      identical(b$rank, n$rank) && identical(b$dim, n$dim) &&
      identical(b$identifiable, n$identifiable) && identical(b$lines, n$lines)
    cat(sprintf("  %-10s %s\n", k, if (same) "OK (identical result)" else "*** MISMATCH ***"))
    if (!same) { ok <- FALSE; cat("    before: rank ", b$rank, " lines:\n",
                                   paste("     ", b$lines, collapse = "\n"), "\n",
                                   "    now:    rank ", n$rank, " lines:\n",
                                   paste("     ", n$lines, collapse = "\n"), "\n", sep = "") }
  }
  if (!ok) stop("fingerprint mismatch: an optimization changed the analysis result")
  cat("all fingerprints identical to baseline\n")
  invisible(TRUE)
}

# auto-run only when invoked directly (Rscript bench/...), not when source()'d
if (sys.nframe() == 0L && !interactive()) {
  res <- bench_run()
  cat("\nelapsed summary:\n")
  for (r in res$results)
    cat(sprintf("  %-10s %8.2fs\n", paste0(r$name, "_c", r$cores), r$elapsed))
}
