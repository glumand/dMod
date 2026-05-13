## ============================================================================
## FOCEI benchmark fixture (Theophylline, 12 subjects).
##
## Runs the canonical Theophylline FOCEI fit, captures wallclock + Rprof
## breakdown, persists a baseline RDS for cross-commit comparison.
##
## Usage:
##   Rscript bench/bench_focei_theoph.R
##
##   # or, to overwrite the baseline (only do this when introducing a new
##   # reference point intentionally):
##   DMOD_BENCH_RECORD=1 Rscript bench/bench_focei_theoph.R
##
## Outputs (always overwritten, in a transient scratch dir):
##   <scratch>/focei_theoph_LAST.rds       most recent run
##   <scratch>/focei_theoph_Rprof.out      raw Rprof trace
##   <scratch>/focei_theoph_Rprof.txt      summaryRprof() dump
##
##   <scratch> defaults to tempdir(); override with env DMOD_BENCH_OUTDIR.
##
## Outputs (only on DMOD_BENCH_RECORD=1, or if file missing):
##   bench/baselines/focei_theoph_pre_rewrite.rds   immutable baseline
## ============================================================================

.dmod_root <- "/home/simon/Documents/Projects/dMod"
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(.dmod_root, quiet = TRUE)
} else {
  library(dMod)
}

# Transient outputs go to a scratch dir (tempdir() by default). The immutable
# pre-rewrite baseline is the only thing that lives under bench/baselines/.
baseline_dir <- file.path(.dmod_root, "bench", "baselines")
dir.create(baseline_dir, recursive = TRUE, showWarnings = FALSE)
scratch_dir  <- Sys.getenv("DMOD_BENCH_OUTDIR",
                           unset = file.path(tempdir(), "dMod_bench_focei"))
dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)
rprof_out  <- file.path(scratch_dir,  "focei_theoph_Rprof.out")
rprof_txt  <- file.path(scratch_dir,  "focei_theoph_Rprof.txt")
last_rds   <- file.path(scratch_dir,  "focei_theoph_LAST.rds")
base_rds   <- file.path(baseline_dir, "focei_theoph_pre_rewrite.rds")

wd <- file.path(tempdir(), "focei_bench")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
oldwd <- setwd(wd); on.exit(setwd(oldwd), add = TRUE)
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
       force = TRUE)

set.seed(1)

## ----------------------------------------------------------------------------
## Model (identical to inst/examples/theophylline.R)
## ----------------------------------------------------------------------------
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
subjects <- sort(unique(Theoph$Subject))
N <- length(subjects)

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
m <- odemodel(reactions, modelname = "theoph_ode", compile = TRUE,
              solver = "CppODE", deriv2 = TRUE)
x <- Xs(m)
g <- Y(c(y = "Cc"), x, modelname = "theoph_obs", compile = TRUE, deriv2 = TRUE)
err <- Y(eqnvec(y = "sigma_add"), g, attach.input = FALSE,
         compile = TRUE, modelname = "theoph_err")

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
p <- P(trafos, method = "explicit", compile = TRUE, modelname = "theoph_p",
       deriv2 = TRUE)
model <- g * x * p

om <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
joint <- normL2(dlist, model, errmodel = err, use.bessel = FALSE) +
         constraintL2(mu = 0, Omega = om)

init <- c(tka = 0.0, tv = 3.0, tcl = 0.7, log_sigma_add = log(0.5))
init[om$cholPars] <- rep(log(0.3), length(om$cholPars))

## ----------------------------------------------------------------------------
## Run with Rprof
## ----------------------------------------------------------------------------
cat(sprintf("[%s] starting FOCEI fit (Theoph, N=%d, init = log-scale)\n",
            format(Sys.time(), "%H:%M:%S"), N))

Rprof(rprof_out, interval = 0.02, line.profiling = FALSE, memory.profiling = FALSE)
t_focei <- system.time(
  fit <- nlmeFit(joint, om, init,
                 model    = model,
                 data     = dlist,
                 errmodel = err,
                 method   = "focei",
                 control  = list(focei = list(
                   innerControl = list(rtol = 1e-7, maxit = 30),
                   correction   = "eager",
                   trustControl = list(iterlim = 100))),
                 verbose = FALSE)
)
Rprof(NULL)

prof_summary <- summaryRprof(rprof_out)

cat(sprintf("\n[done] elapsed = %.2fs   OFV = %.6f   iter = %s   converged = %s\n",
            t_focei[["elapsed"]], fit$value,
            format(fit$iterations %||% NA), format(fit$converged)))

## ----------------------------------------------------------------------------
## Persist
## ----------------------------------------------------------------------------
record <- list(
  timestamp     = Sys.time(),
  dmod_sha      = tryCatch(
    system(paste("git -C", shQuote(.dmod_root), "rev-parse HEAD"),
           intern = TRUE),
    error = function(e) NA_character_),
  elapsed       = unname(t_focei[["elapsed"]]),
  user          = unname(t_focei[["user.self"]]),
  sys           = unname(t_focei[["sys.self"]]),
  value         = fit$value,
  argument      = fit$argument,
  iterations    = fit$iterations,
  converged     = fit$converged,
  Omega         = fit$Omega,
  etaModes      = fit$etaModes,
  init          = init,
  prof_by_total = head(prof_summary$by.total, 30),
  prof_by_self  = head(prof_summary$by.self,  30),
  prof_sampling = prof_summary$sampling.time
)

saveRDS(record, last_rds)
cat("Wrote: ", last_rds, "\n")

# Persist the canonical "pre-rewrite" baseline only if missing or env asks.
write_baseline <- !file.exists(base_rds) ||
                  identical(Sys.getenv("DMOD_BENCH_RECORD"), "1")
if (write_baseline) {
  saveRDS(record, base_rds)
  cat("Wrote: ", base_rds, "  (baseline)\n")
} else {
  cat("Existing baseline kept: ", base_rds, "\n")
  cat("  (set DMOD_BENCH_RECORD=1 to overwrite)\n")
}

## ----------------------------------------------------------------------------
## Pretty Rprof summary
## ----------------------------------------------------------------------------
sink(rprof_txt)
cat("dMod FOCEI benchmark Rprof summary\n")
cat("==================================\n\n")
cat(sprintf("elapsed       : %.2f s\n", record$elapsed))
cat(sprintf("OFV           : %.6f\n",   record$value))
cat(sprintf("outer iter    : %s\n",     format(record$iterations)))
cat(sprintf("converged     : %s\n",     format(record$converged)))
cat(sprintf("sampling time : %.2f s\n", record$prof_sampling))
cat("\nTop 20 by total time:\n")
print(head(record$prof_by_total, 20))
cat("\nTop 20 by self time:\n")
print(head(record$prof_by_self, 20))
sink()
cat("Wrote: ", rprof_txt, "\n")

cat("\nTop 10 by self time:\n")
print(head(record$prof_by_self, 10))
