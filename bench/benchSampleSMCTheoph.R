## ============================================================================
## Bayesian NLME (Pfad A, marginal Laplace) SMC benchmark on Theophylline.
##
## Profiles mcmc(bayesNLMEMarginal(...), sequenceType = "sequential") with the
## "good defaults" recipe used in inst/examples/theophyllineBayes.R:
##
##   1. Fit FOCEI MAP first (one nlmeFit(method = "focei") call), reuse its
##      Hessian both as the fixed Langevin preconditioner and as the
##      covariance of the Laplace cloud (priorSample).
##   2. SMC moves: moveType = "langevin", metric = "fixed",
##      langevinControl(GFixed = fit_focei$hessian).
##   3. SMC schedule: sequenceSchedule = "adaptiveEss",
##      smcControl(malaSteps = 5, essThreshold = 0.7).
##
## Each SMC evaluation still runs one FOCEI inner trust per subject, so total
## wallclock scales as
##   (nParticles * malaSteps * nLevels * N_subjects * inner_trust_iters).
## The Laplace cloud + fixed preconditioner cuts nLevels roughly in half vs
## starting from a wide prior, and the per-MALA cost is one inner trust per
## subject (no Fisher refactorisation per step).
##
## Default smoke configuration: 4 subjects, 200 particles, 5 MALA steps per
## level (~2-4 min wallclock). Set DMOD_BENCH_FULL=1 for the canonical
## 12-subject Theoph run (~15-30 min). Per-axis overrides:
##   DMOD_BENCH_PARTICLES, DMOD_BENCH_MALASTEPS, DMOD_BENCH_CORES,
##   DMOD_BENCH_ESS (essThreshold).
##
## Roundtrip check: a fixed seed (set.seed(1)) makes the run reproducible.
## The first invocation persists a baseline record
## (smc_theoph_baseline_<mode>.rds in bench/baselines/); subsequent runs
## compare posterior mean and log-evidence against the baseline and print
## PASS / FAIL. Overwrite the baseline intentionally with DMOD_BENCH_RECORD=1.
##
## Usage:   Rscript bench/benchSampleSMCTheoph.R
##          DMOD_BENCH_FULL=1 Rscript bench/benchSampleSMCTheoph.R
##          DMOD_BENCH_RECORD=1 Rscript bench/benchSampleSMCTheoph.R
##
## Outputs (in tempdir() by default; override via DMOD_BENCH_OUTDIR env):
##   smc_theoph_LAST_<mode>.rds       most recent run
##   smc_theoph_Rprof_<mode>.out/.txt raw + pretty Rprof
##
## Outputs (only on DMOD_BENCH_RECORD=1 or if missing):
##   bench/baselines/smc_theoph_baseline_<mode>.rds   immutable baseline
## ============================================================================

.dmod_root <- "/home/simon/Documents/Projects/dMod"
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(.dmod_root, quiet = TRUE)
} else {
  library(dMod)
}

baseline_dir <- file.path(.dmod_root, "bench", "baselines")
dir.create(baseline_dir, recursive = TRUE, showWarnings = FALSE)
scratch_dir  <- Sys.getenv("DMOD_BENCH_OUTDIR",
                           unset = file.path(tempdir(), "dMod_bench_smc"))
dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)
bench_mode <- if (identical(Sys.getenv("DMOD_BENCH_FULL"), "1"))
  "full" else "smoke"
rprof_out  <- file.path(scratch_dir,  paste0("smc_theoph_Rprof_",
                                              bench_mode, ".out"))
rprof_txt  <- file.path(scratch_dir,  paste0("smc_theoph_Rprof_",
                                              bench_mode, ".txt"))
last_rds   <- file.path(scratch_dir,  paste0("smc_theoph_LAST_",
                                              bench_mode, ".rds"))
base_rds   <- file.path(baseline_dir, paste0("smc_theoph_baseline_",
                                              bench_mode, ".rds"))

wd <- file.path(tempdir(), "smc_bench")
unlink(wd, recursive = TRUE, force = TRUE)
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
setwd(wd)

set.seed(1)

## Model identical to bench/bench_focei_theoph.R for cross-comparison.
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
all_subjects <- sort(unique(Theoph$Subject))

full_bench <- identical(Sys.getenv("DMOD_BENCH_FULL"), "1")
n_subjects_use <- if (full_bench) length(all_subjects) else 4L
subjects <- all_subjects[seq_len(n_subjects_use)]
Theoph <- Theoph[Theoph$Subject %in% subjects, , drop = FALSE]
N <- length(subjects)

n_particles_use <- as.integer(Sys.getenv("DMOD_BENCH_PARTICLES",
                                          unset = if (full_bench) "500" else "200"))
mala_steps_use  <- as.integer(Sys.getenv("DMOD_BENCH_MALASTEPS", unset = "5"))
ess_threshold   <- as.numeric(Sys.getenv("DMOD_BENCH_ESS", unset = "0.7"))
cores_use <- {
  user_cores <- Sys.getenv("DMOD_BENCH_CORES", unset = "")
  if (nzchar(user_cores)) as.integer(user_cores)
  else max(1L, min(n_particles_use, parallel::detectCores() - 1L))
}

doses <- vapply(subjects, function(s) {
  rec <- Theoph[Theoph$Subject == s, ][1, ]
  rec$Dose * rec$Wt
}, 0.0)

data_df <- data.frame(name = "y", time = Theoph$Time, value = Theoph$conc,
                      sigma = NA_real_, condition = Theoph$Subject,
                      stringsAsFactors = FALSE)
dlist <- as.datalist(data_df)

reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")
m <- odemodel(reactions, modelname = "theoph_smc_ode", compile = TRUE,
              solver = "CppODE", deriv2 = TRUE)
x <- Xs(m)
g <- Y(c(y = "Cc"), x, modelname = "theoph_smc_obs", compile = TRUE, deriv2 = TRUE)
err <- Y(eqnvec(y = "sigma_add"), g, attach.input = FALSE,
         compile = TRUE, modelname = "theoph_smc_err")

trafo <- eqnvec(Ka        = "exp(tka + eta_Ka)",
                V         = "exp(tv  + eta_V)",
                Cl        = "exp(tcl + eta_Cl)",
                Ag        = "Ag_init",
                Cc        = "0",
                sigma_add = "exp(log_sigma_add)")
subj_table <- data.frame(
  eta_Ka  = paste0("eta_Ka_", subjects),
  eta_V   = paste0("eta_V_",  subjects),
  eta_Cl  = paste0("eta_Cl_", subjects),
  Ag_init = doses, row.names = subjects, stringsAsFactors = FALSE)
trafos <- branch(trafo, table = subj_table, apply = "insert")
p <- P(trafos, method = "explicit", compile = TRUE, modelname = "theoph_smc_p",
       deriv2 = TRUE)
prdfn <- g * x * p

om <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
joint <- normL2(dlist, prdfn, errmodel = err, use.bessel = FALSE) +
         constraintL2(mu = 0, Omega = om)

init_bayes <- c(tka = 0.0, tv = 3.0, tcl = 0.7, log_sigma_add = log(0.5))
init_bayes <- c(init_bayes,
                setNames(rep(log(0.3), length(om$cholPars)), om$cholPars))

## ----------------------------------------------------------------------------
## 1. FOCEI MAP -- warm-start point + Hessian for the SMC preconditioner
## ----------------------------------------------------------------------------
cat(sprintf("[%s] fitting FOCEI MAP (warm-start for SMC)...\n",
            format(Sys.time(), "%H:%M:%S")))
t_focei <- system.time(
  focei_fit <- nlmeFit(joint, om, init_bayes,
                        prdfn = prdfn, data = dlist, errfn = err,
                        method = "focei",
                        control = list(focei = list(
                          innerControl = list(rtol = 1e-7, maxit = 30),
                          trustControl = list(iterlim = 80))),
                        verbose = FALSE))
cat(sprintf("  FOCEI converged in %.2fs, OFV = %.4f\n",
            t_focei[["elapsed"]], focei_fit$value))

## ----------------------------------------------------------------------------
## 2. Laplace cloud + priors centered at the MAP
## ----------------------------------------------------------------------------
mapArg    <- focei_fit$argument
mapCov    <- 2 * solve(focei_fit$hessian)         # -2 log L -> log p
mapCovSym <- (mapCov + t(mapCov)) / 2
mapL      <- chol(mapCovSym)
parNames  <- names(mapArg)

priorSample <- function(n) {
  Z <- matrix(rnorm(n * length(mapArg)), nrow = n)
  out <- sweep(Z %*% mapL, 2L, mapArg, FUN = "+")
  colnames(out) <- parNames
  out
}
priorTheta    <- constraintL2(mapArg[setdiff(parNames, om$cholPars)],
                              sigma = 5)
priorOmegaObj <- priorOmega(om, kind = "LKJHalfNormal", scaleSD = 1.0)

tgt <- bayesNLMEMarginal(
  obj           = joint, omegaSpec = om,
  prdfn         = prdfn, data = dlist, errfn = err,
  init          = init_bayes,
  priorTheta    = priorTheta,
  priorOmegaObj = priorOmegaObj,
  priorSample   = priorSample)

cat(sprintf("[%s] starting SMC (Theoph%s, N=%d subjects, n_particles=%d, malaSteps=%d, ess=%.2f, cores=%d)\n",
            format(Sys.time(), "%H:%M:%S"),
            if (full_bench) " FULL" else " smoke",
            N, n_particles_use, mala_steps_use, ess_threshold, cores_use))

## ----------------------------------------------------------------------------
## 3. SMC with fixed (MAP-curvature) Langevin preconditioner
## ----------------------------------------------------------------------------
Rprof(rprof_out, interval = 0.02)
t_smc <- system.time(
  chain <- mcmc(
    target           = tgt,
    sequenceType     = "sequential",
    sequenceSchedule = "adaptiveEss",
    populationSize   = n_particles_use,
    moveType         = "langevin",
    metric           = "fixed",
    moveControl      = langevinControl(GFixed = focei_fit$hessian),
    sequenceControl  = smcControl(malaSteps = mala_steps_use,
                                   essThreshold = ess_threshold,
                                   verbose = TRUE),
    cores            = cores_use)
)
Rprof(NULL)

prof_summary <- summaryRprof(rprof_out)

cat(sprintf("\n[done] elapsed = %.2fs   logZ = %.3f   levels = %d\n",
            t_smc[["elapsed"]], chain$logEvidence, chain$nLevels))
cat("\n--- top 10 self time (R-side) ---\n")
print(head(prof_summary$by.self, 10))
cat("\n--- top 10 total time ---\n")
print(head(prof_summary$by.total, 10))

## ----------------------------------------------------------------------------
## 4. Roundtrip cross-check vs the FOCEI MAP we already have
## ----------------------------------------------------------------------------
posterior_mean <- colMeans(chain$samples)
posterior_sd   <- apply(chain$samples, 2L, sd)
shared <- intersect(names(focei_fit$argument), names(posterior_mean))
struct_axis <- intersect(shared, names(init_bayes)[
  !names(init_bayes) %in% om$cholPars])
map_match_struct <- max(abs(posterior_mean[struct_axis] -
                            focei_fit$argument[struct_axis]) /
                        posterior_sd[struct_axis])

cat(sprintf("\nFOCEI MAP cross-check (structural axis):\n"))
for (nm in struct_axis) {
  cat(sprintf("  %-15s MAP = %+8.4f  Bayes mean = %+8.4f  diff/SD = %+5.2f\n",
              nm, focei_fit$argument[nm], posterior_mean[nm],
              (posterior_mean[nm] - focei_fit$argument[nm]) /
                posterior_sd[nm]))
}
cat(sprintf("  max |Bayes mean - MAP| / SD = %.2f\n", map_match_struct))

record <- list(
  timestamp     = Sys.time(),
  dmod_sha      = tryCatch(
    system(paste("git -C", shQuote(.dmod_root), "rev-parse HEAD"),
           intern = TRUE),
    error = function(e) NA_character_),
  full_bench    = full_bench,
  seed          = 1L,
  elapsed       = unname(t_smc[["elapsed"]]),
  foceiElapsed  = unname(t_focei[["elapsed"]]),
  nSubjects     = N,
  nParticles    = n_particles_use,
  malaSteps     = mala_steps_use,
  essThreshold  = ess_threshold,
  nLevels       = chain$nLevels,
  logEvidence   = chain$logEvidence,
  acceptRates   = chain$acceptRates,
  posteriorMean = posterior_mean,
  posteriorSD   = posterior_sd,
  foceiArgument = focei_fit$argument,
  foceiOFV      = focei_fit$value,
  mapMatchStruct = map_match_struct,
  prof_by_self  = head(prof_summary$by.self,  20),
  prof_by_total = head(prof_summary$by.total, 20))
saveRDS(record, last_rds)
capture.output(prof_summary, file = rprof_txt)
cat("\nWrote: ", last_rds, "\n", sep = "")

## ----------------------------------------------------------------------------
## 5. Baseline write / compare
## ----------------------------------------------------------------------------
write_baseline <- !file.exists(base_rds) ||
                  identical(Sys.getenv("DMOD_BENCH_RECORD"), "1")
if (write_baseline) {
  saveRDS(record, base_rds)
  cat("Wrote: ", base_rds, "  (baseline)\n", sep = "")
  cat("  (this is the first run or DMOD_BENCH_RECORD=1 was set; ",
      "subsequent runs will compare against it)\n", sep = "")
} else {
  ref <- readRDS(base_rds)
  cat("\nRoundtrip vs baseline (", base_rds, "):\n", sep = "")
  shared_pars <- intersect(names(record$posteriorMean), names(ref$posteriorMean))
  diff_mean_sd <- abs(record$posteriorMean[shared_pars] -
                      ref$posteriorMean[shared_pars]) /
                  pmax(record$posteriorSD[shared_pars], 1e-8)
  diff_logZ <- abs(record$logEvidence - ref$logEvidence)
  cat(sprintf("  max |mean(now) - mean(ref)| / SD(now) = %.2f\n",
              max(diff_mean_sd)))
  cat(sprintf("  |logZ(now) - logZ(ref)|              = %.3f\n", diff_logZ))
  # PASS criterion: posterior mean within 1 SD across all parameters,
  # log-evidence within 5. Loose because SMC has substantial Monte-Carlo
  # variance even at large particle counts.
  pass <- max(diff_mean_sd) < 1.0 && diff_logZ < 5.0
  cat(sprintf("  %s  (overwrite intentionally with DMOD_BENCH_RECORD=1)\n",
              if (pass) "PASS" else "FAIL"))
}

cat("\nOutputs written under: ", scratch_dir, "\n", sep = "")
