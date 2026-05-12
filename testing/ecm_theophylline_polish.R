## ECM polish benchmark on the Theophylline NLME dataset
##
## Builds the same model as testing/focei_theophylline.R, runs focei() via the
## new unified API (emObjfn + focei orchestrator), then polishes with
## ecmFit(mode="focei_warmstart"). Tabulates per-stage OFV and parameter
## shifts vs FOCEI and vs the published NONMEM/nlmixr2 reference values.
##
## Intended use: source() interactively (RStudio / shell). Generates compiled
## C/C++ artefacts in tempdir(). ECM polish at K_eta = 3 with epsQuadLevels
## = c(4, 5, 6) -> (7, 31, 105) nodes per subject; on a 12-subject problem
## this is ~3x the FOCEI cost end-to-end.

rm(list = ls(all.names = TRUE))
.dmod_root <- "/home/simon/Documents/Projects/dMod"
devtools::load_all(.dmod_root, quiet = TRUE)
setwd(tempdir())
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
       force = TRUE)

set.seed(1)

## ---- 1. Data + model (matches testing/focei_theophylline.R) -----------------
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
subjects <- sort(unique(Theoph$Subject))
N <- length(subjects)
cat(sprintf("Theoph: %d subjects, %d total observations\n", N, nrow(Theoph)))

doses <- vapply(subjects, function(s) {
  rec <- Theoph[Theoph$Subject == s, ][1, ]
  rec$Dose * rec$Wt
}, 0.0)

data_df <- data.frame(name = "y", time = Theoph$Time, value = Theoph$conc,
                      sigma = 0.5, condition = Theoph$Subject,
                      stringsAsFactors = FALSE)
dlist <- as.datalist(data_df)

reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")
m <- odemodel(reactions, modelname = "ecm_theoph", compile = TRUE,
              solver = "CppODE", deriv2 = TRUE)
x <- Xs(m)
g <- Y(c(y = "Cc"), x, modelname = "ecm_theoph_obs", compile = TRUE,
       deriv2 = TRUE)

trafo <- eqnvec(Ka = "Ka_pop * exp(eta_Ka)", V = "V_pop * exp(eta_V)",
                Cl = "Cl_pop * exp(eta_Cl)", Ag = "Ag_init", Cc = "0")
subj_table <- data.frame(
  eta_Ka  = paste0("eta_Ka_", subjects),
  eta_V   = paste0("eta_V_",  subjects),
  eta_Cl  = paste0("eta_Cl_", subjects),
  Ag_init = doses, row.names = subjects, stringsAsFactors = FALSE)
trafos <- branch(trafo, table = subj_table, apply = "insert")
p <- P(trafos, method = "explicit", compile = TRUE, modelname = "ecm_theoph_p",
       deriv2 = TRUE)
om <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
joint <- normL2(dlist, g * x * p) + constraintL2(mu = 0, Omega = om)

init <- c(Ka_pop = 1.5, V_pop = 32, Cl_pop = 3.0)
init[om$cholPars] <- c(log(0.5), log(0.1), log(0.3))

## ---- 2. FOCEI via the new unified API ---------------------------------------
cat("\n== FOCEI baseline ==\n")
em_lap <- emObjfn(joint, om, method = "laplace",
                  control = list(innerControl = list(rtol = 1e-7, maxit = 30)))
t_focei <- system.time(
  fit_focei <- focei(em_lap, init, trustControl = list(iterlim = 100))
)
print(t_focei)
print(fit_focei)

## ---- 3. ECM polish via ecmFit(mode = "focei_warmstart") --------------------
cat("\n== ECM polish ==\n")
em_qd <- emObjfn(joint, om, model = g * x * p, data = dlist,
                 method = "quadrature", control = list(level = 4L))
t_ecm <- system.time(
  fit_ecm <- ecmFit(em_qd, init,
                    mode              = "focei_warmstart",
                    epsQuadLevels   = c(4L, 5L),
                    maxEcmPerStage = 3L,
                    epsEcm           = 1e-4,
                    verbose           = TRUE)
)
print(t_ecm)
print(fit_ecm)

## ---- 4. Compare -------------------------------------------------------------
cat("\n== Parameter shifts vs FOCEI ==\n")
shift <- data.frame(
  par   = names(fit_focei$argument),
  focei = round(as.numeric(fit_focei$argument), 4),
  ecm   = round(as.numeric(fit_ecm$argument[names(fit_focei$argument)]), 4)
)
shift$delta    <- shift$ecm - shift$focei
shift$rel_pct  <- round(100 * shift$delta / pmax(abs(shift$focei), 1e-12), 2)
print(shift)

cat("\n== OFV trajectory ==\n")
cat(sprintf("  FOCEI         : %.6f\n", fit_focei$value))
cat(sprintf("  ECM stages    :\n"))
print(fit_ecm$stageTrace)
cat(sprintf("  Final ECM OFV : %.6f\n", fit_ecm$value))
cat(sprintf("  OFV delta     : %.6f\n", fit_ecm$value - fit_focei$value))

cat("\n== NONMEM/nlmixr2 reference (Pinheiro & Bates) ==\n")
cat("  Ka_pop ~ 1.50, V_pop ~ 32, Cl_pop ~ 3.0\n")
cat("  sd(eta_Ka) ~ 0.5, sd(eta_V) ~ 0.1, sd(eta_Cl) ~ 0.3\n")
