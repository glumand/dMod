## ============================================================================
## Theophylline NLME example: FOCEI in dMod (introductory walkthrough)
## ============================================================================
##
## End-to-end NLME walkthrough on the classic Theophylline dataset
## (12 subjects, 1-compartment oral PK). Demonstrates:
##
##   1. Building the dMod prediction chain (eqnlist -> odemodel -> Xs -> Y).
##   2. Per-subject parameter trafos with branch() + P(method = "explicit").
##   3. Assembling the marginal likelihood: normL2(...) + constraintL2(...).
##   4. Fitting via the unified API: nlmeFit(method = "focei").
##   5. Diagnostic plots from R/plots.R (camelCase S3 generics):
##        plotIndivs(fit)      per-subject IPRED/PRED curves + observed dots
##        plot(fit)            DV vs IPRED + DV vs PRED scatter (plot.nlmeFit)
##        plotResiduals(fit)   IWRES vs IPRED + IWRES vs TIME
##        plotHistIndivs(fit)  eta histogram + QQ vs N(0, Omega_kk)
##
## Observation model: concentrations are log-transformed and modelled with a
## constant additive Gaussian residual on the log scale,
##
##     log(DV_ij) = log(Cc_ij) + eps_ij,   eps_ij ~ N(0, sigma_add^2).
##
## On the linear scale this is the standard lognormal / proportional residual
## (constant CV) used in NONMEM / Monolix / nlmixr2 PK fits.
##
## Parameterisation matches NONMEM / Monolix / nlmixr2: structural parameters
## live on the natural-log scale,
##
##     Ka_i = exp(tka + eta_Ka_i),   V_i = exp(tv  + eta_V_i),
##     Cl_i = exp(tcl + eta_Cl_i),   sigma_add = exp(log_sigma_add).
##
## Intended use: source() interactively. The script returns ggplot objects in
## the global environment (p_indivs, p_obs_pred, p_resid, p_etas); print them
## to view. No files are written - generated C/C++ artefacts go to tempdir().
##
## Reference values (Pinheiro & Bates 2000):
##   Ka_pop  ~ 1.50 (1/h)   sd(eta_Ka) ~ 0.5   (% CV ~ 53)
##   V_pop   ~ 32   (L)     sd(eta_V)  ~ 0.1   (% CV ~ 10)
##   Cl_pop  ~ 3.0  (L/h)   sd(eta_Cl) ~ 0.3   (% CV ~ 31)
## ============================================================================

rm(list = ls(all.names = T))

if (requireNamespace("devtools", quietly = T)) {
  devtools::load_all(quiet = T)
} else {
  library(dMod)
}

oldwd <- setwd(tempdir())
on.exit(setwd(oldwd), add = T)
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = T),
       force = T)

set.seed(1)

## ----------------------------------------------------------------------------
## 1. Data
##
## Drop pre-dose rows (time == 0). The 1-cmt oral PK model predicts Cc(0) = 0
## by construction, so any t=0 measurement would force log(0). The remaining
## 120 rows are post-dose absorption / elimination samples (all DV > 0).
## ----------------------------------------------------------------------------
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
Theoph_pos <- Theoph[Theoph$Time > 0, , drop = FALSE]
subjects <- sort(unique(Theoph_pos$Subject))
N <- length(subjects)
cat(sprintf("Theoph: N=%d subjects, %d observations (after dropping pre-dose t==0)\n",
            N, nrow(Theoph_pos)))

doses <- vapply(subjects, function(s) {
  rec <- Theoph[Theoph$Subject == s, ][1, ]
  rec$Dose * rec$Wt
}, 0.0)

# DV is stored on the log scale; the additive error model below is the residual
# SD of log(DV). sigma = NA so the errfn estimates the residual SD jointly with
# the structural parameters and omegas (matches NONMEM / nlmixr2 convention).
data_df <- data.frame(
  name      = "y",         time      = Theoph_pos$Time,
  value     = log(Theoph_pos$conc), sigma     = NA_real_,
  condition = Theoph_pos$Subject,
  stringsAsFactors = FALSE)
dlist <- as.datalist(data_df)

## ----------------------------------------------------------------------------
## 2. ODE prd: 1-compartment oral PK
##      dAg/dt = -Ka * Ag                  (gut depot)
##      dCc/dt =  Ka * Ag / V - Cl/V * Cc  (central concentration)
##
## Observable y = log(Cc + 1e-9). The tiny offset only matters at t == 0 where
## CppODE always emits a row at Cc = 0; for all retained observations y is
## identical to log(Cc) at relative precision ~ 1e-9, far below the residual
## SD.
## ----------------------------------------------------------------------------
reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")
m <- odemodel(reactions, modelname = "theoph_ode", compile = F,
              solver = "CppODE", deriv2 = T)
x <- Xs(m)
g <- Y(c(y = "log(Cc + 1e-9)"), x, modelname = "theoph_obs",
       compile = F, deriv2 = T)

# Additive Gaussian residual on log(Cc). sigma_add is the residual SD on log
# scale; estimated jointly via log_sigma_add.
e <- Y(eqnvec(y = "sigma_add"), g, attach.input = FALSE, deriv2 = T,
       compile = F, modelname = "theoph_err")

## ----------------------------------------------------------------------------
## 3. Per-subject parameter trafos (log-scale, positivity guaranteed)
##      theta_i = exp(t_theta + eta_i)
##      sigma_add = exp(log_sigma_add)
## ----------------------------------------------------------------------------
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
  Ag_init = doses, row.names = subjects, stringsAsFactors = F)
trafos <- branch(trafo, table = subj_table, apply = "insert")
p <- P(trafos, method = "explicit", compile = F, modelname = "theoph_p", deriv2 = T)
prd <- g * x * p

compile(prd, e, cores = 6)

om <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
# use.bessel = FALSE keeps the marginal-likelihood OFV directly comparable to
# NONMEM / Monolix / nlmixr2 (plain ML, no n/(n-p) inflation).
obj <- normL2(dlist, prd, errmodel = e, use.bessel = FALSE) +
         constraintL2(mu = 0, Omega = om)

## Starting values close to (but not at) the Pinheiro-Bates reference.
init <- c(tka           = 0.0,
          tv            = 3.0,
          tcl           = 0.7,
          log_sigma_add = log(0.2))
init[om$cholPars] <- rep(log(0.3), length(om$cholPars))

## ----------------------------------------------------------------------------
## 4. FOCEI fit
## ----------------------------------------------------------------------------
cat("\n== nlmeFit(method = 'focei') ==\n")
fit_focei <- nlmeFit(obj, om, init,
                     prdfn    = prd,
                     data     = dlist,
                     errfn    = e,
                     method   = "focei",
                     control  = list(focei = list(
                       innerControl = list(rtol = 1e-7, maxit = 30),
                       trustControl = list(iterlim = 100))))
print(fit_focei)

## ----------------------------------------------------------------------------
## 5. Population summary (linear pop pars + %CV + OFV)
## ----------------------------------------------------------------------------
cv_pct <- function(omega2) 100 * sqrt(exp(omega2) - 1)
Omega <- fit_focei$Omega
arg   <- fit_focei$argument
pop <- c(
  Ka_pop    = exp(as.numeric(arg["tka"])),
  V_pop     = exp(as.numeric(arg["tv"])),
  Cl_pop    = exp(as.numeric(arg["tcl"])),
  sigma_add = exp(as.numeric(arg["log_sigma_add"])),
  `CV%_Ka`  = cv_pct(Omega["eta_Ka", "eta_Ka"]),
  `CV%_V`   = cv_pct(Omega["eta_V",  "eta_V" ]),
  `CV%_Cl`  = cv_pct(Omega["eta_Cl", "eta_Cl"]),
  OFV       = fit_focei$value)
cat("\n== Population summary ==\n")
print(round(pop, 3))

## ----------------------------------------------------------------------------
## 6. Diagnostic plots (returned as ggplot objects)
##
## Print any of these to view interactively, or wrap with ggplot2::ggsave()
## if you want to write them out. cowplot is optional; without it,
## plotHistIndivs() returns list(hist, qq) instead of one stacked panel.
## ----------------------------------------------------------------------------
p_indivs   <- plotIndivs(fit_focei)
p_obs_pred <- plot(fit_focei)
p_resid    <- plotResiduals(fit_focei)
p_etas     <- plotHistIndivs(fit_focei)

print(p_indivs)
print(p_obs_pred)
print(p_resid)
print(p_etas)
