## End-to-end FOCEI benchmark on the Theophylline NLME dataset
##
## Intended use: source() this script interactively (e.g. from RStudio).
## Generates compiled C/C++ artefacts in tempdir().
##
## Reference values (Pinheiro & Bates 2000, NONMEM/nlmixr2 FOCEI):
##   Ka_pop  ~ 1.50  (1/h)
##   V_pop   ~ 32    (L)
##   Cl_pop  ~ 3.0   (L/h)
##   sd(eta_Ka) ~ 0.5
##   sd(eta_V)  ~ 0.1
##   sd(eta_Cl) ~ 0.3

rm(list = ls(all.names = TRUE))
.dmod_root <- "/home/simon/Documents/Projects/dMod"
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(.dmod_root, quiet = TRUE)
} else {
  library(dMod)
}
setwd(tempdir())
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
       force = TRUE)

set.seed(1)

## ---------------------------------------------------------------------------
## 1. Data: nlme::Theoph (12 subjects, oral dosing, 11 timepoints each)
## ---------------------------------------------------------------------------
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
subjects <- sort(unique(Theoph$Subject))
N <- length(subjects)
cat(sprintf("Theoph: %d subjects, %d total observations\n", N, nrow(Theoph)))

# Total dose (mg) per subject = Dose (mg/kg) x Wt (kg)
doses <- vapply(subjects, function(s) {
  rec <- Theoph[Theoph$Subject == s, ][1, ]
  rec$Dose * rec$Wt
}, 0.0)

# Build datalist: sigma column needed by normL2; use 0.5 mg/L as starting noise
# (matches typical NONMEM additive error). For exact match swap in proportional.
data_df <- data.frame(
  name      = "y",
  time      = Theoph$Time,
  value     = Theoph$conc,
  sigma     = 0.5,
  condition = Theoph$Subject,
  stringsAsFactors = FALSE
)
dlist <- as.datalist(data_df)

## ---------------------------------------------------------------------------
## 2. ODE model: 1-compartment oral PK
##      dAg/dt = -Ka * Ag                  (gut depot)
##      dCc/dt =  Ka * Ag / V - Cl/V * Cc  (central concentration)
## ---------------------------------------------------------------------------
reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")

m <- odemodel(reactions, modelname = "theoph_ode", compile = TRUE,
              solver = "CppODE", deriv2 = TRUE)
x <- Xs(m)

g <- Y(c(y = "Cc"), x, modelname = "theoph_obs", compile = TRUE,
       deriv2 = TRUE)

## ---------------------------------------------------------------------------
## 3. Per-subject parameter trafos via branch + insert
## ---------------------------------------------------------------------------
trafo <- eqnvec(
  Ka = "Ka_pop * exp(eta_Ka)",
  V  = "V_pop  * exp(eta_V)",
  Cl = "Cl_pop * exp(eta_Cl)",
  Ag = "Ag_init",
  Cc = "0"
)

subj_table <- data.frame(
  eta_Ka  = paste0("eta_Ka_",  subjects),
  eta_V   = paste0("eta_V_",   subjects),
  eta_Cl  = paste0("eta_Cl_",  subjects),
  Ag_init = doses,
  row.names = subjects,
  stringsAsFactors = FALSE
)
trafos <- branch(trafo, table = subj_table, apply = "insert")
p <- P(trafos, method = "explicit", compile = TRUE, modelname = "theoph_p",
       deriv2 = TRUE)

## ---------------------------------------------------------------------------
## 4. omegaSpec (diagonal) and joint objfn
## ---------------------------------------------------------------------------
om <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)

joint <- normL2(dlist, g * x * p) + constraintL2(mu = 0, Omega = om)

# Stage 1 (envelope-only) and Stage 2 (analytical log|H|-correction, lagged).
focei_none <- focei(joint, om, innerControl = list(rtol = 1e-6, maxit = 50))
focei_lag  <- focei(joint, om, model = g * x * p, data = dlist,
                    correction = "lagged",
                    innerControl = list(rtol = 1e-6, maxit = 50))

## ---------------------------------------------------------------------------
## 5. Initial population parameters and outer optimisation
## ---------------------------------------------------------------------------
init_outer <- c(
  Ka_pop      = 1.5,
  V_pop       = 30.0,
  Cl_pop      = 3.0,
  omega_Ka_Ka = log(0.5),
  omega_V_V   = log(0.2),
  omega_Cl_Cl = log(0.3)
)

cat(sprintf("Initial OFV (Stage 1): %.4f\n", focei_none(init_outer)$value))
cat(sprintf("Initial OFV (Stage 2): %.4f\n", focei_lag (init_outer)$value))

run_fit <- function(focei_obj, label) {
  fit <- suppressMessages(trust(
    focei_obj, init_outer,
    rinit = 1, rmax = 10, iterlim = 200,
    fterm = 1e-6, mterm = 1e-6,
    parlower = c(Ka_pop = 0.01, V_pop = 1, Cl_pop = 0.01,
                 omega_Ka_Ka = -10, omega_V_V = -10, omega_Cl_Cl = -10),
    parupper = Inf, printIter = FALSE,
    on_step = attr(focei_obj, "on_step")
  ))
  cat(sprintf("\n--- %s --------------------------------------------\n", label))
  cat(sprintf("Converged: %s after %d outer iterations, OFV: %.4f\n",
              fit$converged, fit$iterations, fit$value))
  cat("Population estimates:\n")
  print(round(fit$argument[c("Ka_pop", "V_pop", "Cl_pop")], 4))
  sd_eta <- exp(fit$argument[c("omega_Ka_Ka", "omega_V_V", "omega_Cl_Cl")])
  names(sd_eta) <- c("sd_eta_Ka", "sd_eta_V", "sd_eta_Cl")
  cat("Random-effect SDs:\n"); print(round(sd_eta, 4))
  diag_info <- attr(focei_obj(fit$argument, deriv = FALSE), "emDiag")
  cat(sprintf("Mean inner iter/subject: %.1f, refreshCount: %s\n",
              mean(diag_info$iter),
              if (is.null(diag_info$refreshCount)) "n/a"
              else as.character(diag_info$refreshCount)))
  fit
}

cat("\n========== FOCEI Theophylline Benchmark ==========\n")
fit_none <- run_fit(focei_none, "Stage 1 (envelope only)")
fit_lag  <- run_fit(focei_lag,  "Stage 2 (analytical log|H|-correction, lagged)")

cat("\nReference (Pinheiro & Bates, NONMEM-FOCEI region):\n")
cat("  Ka_pop ~ 1.50,  V_pop ~ 32,  Cl_pop ~ 3.0\n")
cat("  sd(eta_Ka) ~ 0.5, sd(eta_V) ~ 0.1, sd(eta_Cl) ~ 0.3\n\n")

invisible(list(none = fit_none, lagged = fit_lag))
