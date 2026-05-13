## ============================================================================
## Warfarin PK/PD NLME example (introductory walkthrough)
## ============================================================================
##
## End-to-end PK-PD walkthrough on the classic Warfarin dataset from
## `nlmixr2data::warfarin` (32 healthy volunteers, single oral dose, plasma
## warfarin concentration cp + prothrombin complex activity pca over 5 days).
## Demonstrates:
##
##   1. Building a coupled 1-cmt-PK + indirect-response-PD model with TWO
##      observables in one dMod prdfn chain.
##   2. Running FOCEI via the unified API: nlmeFit(method = "focei").
##   3. Multi-start FOCEI via msnlmeFit() and the classic dMod multistart
##      diagnostic: parlist -> as.parframe() -> plotValues() waterfall.
##   4. Multi-observable diagnostic plots from R/plots.R that switch to a
##      facet_grid(name ~ condition) layout automatically:
##      - plotIndivs(fit)      observable x subject grid with IPRED + sigma
##        ribbon, PRED overlay, and observed dots.
##      - plot(fit)            DV vs IPRED + DV vs PRED, faceted by
##        observable in rows.
##      - plotResiduals(fit)   IWRES vs IPRED + IWRES vs TIME, observable
##        in rows.
##      - plotHistIndivs(fit)  per-eta histogram with the panel-specific
##        N(0, Omega_kk) overlaid, plus QQ-plot.
##
## Pharmacological model:
##
##   PK (1-cmt oral first-order absorption + linear elimination):
##     dGut/dt    = -Ka * Gut
##     dCenter/dt =  Ka * Gut - (Cl / V) * Center
##     cp         =  Center / V
##
##   PD (turnover / indirect response with synthesis inhibition by cp,
##       Imax fixed to 1, baseline pca0 = kin / kout):
##     dPCA/dt = kin * IC50 / (IC50 + cp) - kout * PCA
##     PCA(0)  = pca0
##
## Parameterisation: structural fixed effects on the natural-log scale, with
## diagonal random effects on V, Cl, pca0, kout. Ka and IC50 are population
## only, which keeps the FOCEI fit well-identified on N = 32.
##
##     V_i    = exp(tv    + eta_V_i)      Ka  = exp(tka)
##     Cl_i   = exp(tcl   + eta_Cl_i)     IC50 = exp(tic50)
##     pca0_i = exp(tpca0 + eta_pca0_i)
##     kout_i = exp(tkout + eta_kout_i)
##
## Observation model:
##
##   y_cp  = log(cp + 1e-9)   (cp log-transformed; constant additive Gaussian
##                              residual on log scale = constant CV on the
##                              linear scale).
##   y_pca = PCA              (linear-scale PCA; constant additive Gaussian
##                              residual). PCA is on a bounded 0-100 scale
##                              so a linear additive SD is appropriate.
##
## Intended use: source() interactively. The script returns ggplot objects in
## the global environment (p_waterfall, p_indivs, p_obs_pred, p_resid, p_etas);
## print them to view. No files are written - generated C/C++ artefacts go to
## tempdir().
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
## Drop dosing rows (evid == 1) and the four cp == 0 observations (log is
## undefined; the absorption phase has no information at strict zero). pca
## stays unchanged.
## ----------------------------------------------------------------------------
if (!requireNamespace("nlmixr2data", quietly = T))
  stop("Install the `nlmixr2data` package to access the warfarin dataset.")
data(warfarin, package = "nlmixr2data")
warfarin$id <- as.character(warfarin$id)

obs <- warfarin[warfarin$evid == 0 &
                  !(warfarin$dvid == "cp" & warfarin$dv == 0), ,
                drop = FALSE]
subjects <- sort(unique(obs$id))
N <- length(subjects)
cat(sprintf("Warfarin: N=%d subjects, %d observations (%d cp + %d pca)\n",
            N, nrow(obs),
            sum(obs$dvid == "cp"), sum(obs$dvid == "pca")))

doses <- vapply(subjects, function(s) {
  rec <- warfarin[warfarin$id == s & warfarin$evid == 1, , drop = FALSE][1, ]
  rec$amt
}, 0.0)

# Build the dMod long-format data frame with two observation names:
#   "y_cp"  - log-transformed cp
#   "y_pca" - linear PCA
# sigma = NA so both residual SDs are estimated jointly with the structural
# parameters and the omega Cholesky parameters.
data_df <- rbind(
  data.frame(name      = "y_cp",
             time      = obs$time[obs$dvid == "cp"],
             value     = log(obs$dv[obs$dvid == "cp"]),
             sigma     = NA_real_,
             condition = obs$id[obs$dvid == "cp"],
             stringsAsFactors = FALSE),
  data.frame(name      = "y_pca",
             time      = obs$time[obs$dvid == "pca"],
             value     = obs$dv[obs$dvid == "pca"],
             sigma     = NA_real_,
             condition = obs$id[obs$dvid == "pca"],
             stringsAsFactors = FALSE))
dlist <- as.datalist(data_df)

## ----------------------------------------------------------------------------
## 2. Coupled PK + PD ODE model
## ----------------------------------------------------------------------------
reactions <- eqnlist()
reactions <- addReaction(reactions, "Gut",   "",      "Ka * Gut",
                         "absorption")
reactions <- addReaction(reactions, "",     "Center", "Ka * Gut",
                         "appearance")
reactions <- addReaction(reactions, "Center","",     "(Cl/V) * Center",
                         "elimination")
reactions <- addReaction(reactions, "",     "PCA",
                         "kin * IC50 / (IC50 + Center/V)",
                         "PCA_synthesis")
reactions <- addReaction(reactions, "PCA",  "",      "kout * PCA",
                         "PCA_degradation")
m <- odemodel(reactions, modelname = "warf_ode", compile = F,
              solver = "CppODE", deriv2 = T)
x <- Xs(m, optionsSens = list(atol = 1e-4, rtol = 1e-4))

# Two observables in one Y() call. log(Cc + 1e-9) keeps the spurious t == 0
# prediction row CppODE emits (Center(0) = 0 -> cp(0) = 0) finite, which
# dMod's prediction sanity check requires. The offset is far below the
# residual SD so it does not affect the fit.
g <- Y(c(y_cp  = "log(Center/V + 1e-3)",
         y_pca = "PCA"),
       x, modelname = "warf_obs", compile = F, deriv2 = T)

# Two independent additive residual SDs, one per observable.
e <- Y(eqnvec(y_cp  = "sigma_cp",
              y_pca = "sigma_pca"),
       g, attach.input = FALSE, deriv2 = T,
       compile = F, modelname = "warf_err")

## ----------------------------------------------------------------------------
## 3. Per-subject parameter trafos
##
## kin is derived as pca0 * kout in the trafo so the random-effect-bearing
## parameters are the interpretable ones (pca0, kout). Initial conditions:
## Gut(0) = dose_i, Center(0) = 0, PCA(0) = pca0_i (start at PD baseline).
## ----------------------------------------------------------------------------
trafo <- eqnvec(
  Ka        = "exp(tka)",
  V         = "exp(tv  + eta_V)",
  Cl        = "exp(tcl + eta_Cl)",
  kout      = "exp(tkout + eta_kout)",
  IC50      = "exp(tic50)",
  kin       = "exp(tpca0 + eta_pca0) * exp(tkout + eta_kout)",
  Gut       = "Gut_init",
  Center    = "0",
  PCA       = "exp(tpca0 + eta_pca0)",
  sigma_cp  = "exp(log_sigma_cp)",
  sigma_pca = "exp(log_sigma_pca)")
subj_table <- data.frame(
  eta_V    = paste0("eta_V_",    subjects),
  eta_Cl   = paste0("eta_Cl_",   subjects),
  eta_pca0 = paste0("eta_pca0_", subjects),
  eta_kout = paste0("eta_kout_", subjects),
  Gut_init = doses,
  row.names = subjects, stringsAsFactors = F)
trafos <- branch(trafo, table = subj_table, apply = "insert")
p <- P(trafos, method = "explicit", compile = F, modelname = "warf_p",
       deriv2 = T)
prd <- g * x * p

compile(prd, e, cores = 8)

om <- omega(eta = c("eta_V", "eta_Cl", "eta_pca0", "eta_kout"),
            subjects = subjects)
obj <- normL2(dlist, prd, errmodel = e, use.bessel = FALSE) +
         constraintL2(mu = 0, Omega = om)

## Starting values (literature defaults rounded to a convenient log-grid):
##   Ka_pop   = exp(0.0)   = 1.0    /h
##   V_pop    = exp(2.0)   = 7.4    L
##   Cl_pop   = exp(-2.0)  = 0.14   L/h
##   pca0_pop = exp(4.6)   = 100
##   kout_pop = exp(-3.0)  = 0.05   /h
##   IC50_pop = exp(0.7)   = 2.0    mg/L
##   sigma_cp  = exp(-1.6) = 0.20    (log-scale, ~ 20% CV)
##   sigma_pca = exp(2.0)  = 7.4     (linear PCA units)
## eta-SDs at 0.3 -> log(0.3) per cholPars (diagonal omega).
init <- c(tka            = 0.0,
          tv             = 2.0,
          tcl            = -2.0,
          tpca0          = 4.6,
          tkout          = -3.0,
          tic50          = 0.7,
          log_sigma_cp   = log(0.2),
          log_sigma_pca  = log(7.4))
init[om$cholPars] <- rep(log(0.3), length(om$cholPars))

## ----------------------------------------------------------------------------
## 4. Single FOCEI fit (reference run)
## ----------------------------------------------------------------------------
cat("\n== nlmeFit(method = 'focei') ==\n")
fit_focei <- nlmeFit(obj, om, init,
                     prdfn   = prd,
                     data    = dlist,
                     errfn   = e,
                     method  = "focei",
                     control = list(focei = list(
                       innerControl = list(rtol = 1e-7, maxit = 50),
                       trustControl = list(iterlim = 200))))
print(fit_focei)

## ----------------------------------------------------------------------------
## 5. Multistart FOCEI + waterfall plot
##
## msnlmeFit() forks workers via parallel::mclapply on Unix (PSOCK + foreach
## on Windows) and returns a parlist of fits. The classical dMod multistart
## workflow converts the parlist to a parframe and feeds it to plotValues()
## for the waterfall (OFV vs sorted-fit index). Plateaus / steps in the
## waterfall identify distinct local minima; a single flat plateau at the
## lowest OFV is the signature of a globally identified fit.
## ----------------------------------------------------------------------------
cat("\n== msnlmeFit (20 starts, perturbed around init) ==\n")
ms_fits <- msnlmeFit(obj, om, init,
                     prdfn  = prd, data = dlist, errfn = e,
                     method = "focei",
                     control = list(focei = list(
                       innerControl = list(rtol = 1e-7, maxit = 50),
                       trustControl = list(iterlim = 200))),
                     fits  = 50L,
                     cores = 10L,
                     sd    = 2)


summary(ms_fits)

# parlist -> parframe -> waterfall via plotValues().
ms_pframe   <- as.parframe(ms_fits)
p_waterfall <- plotValues(ms_pframe, tol = 1)

## ----------------------------------------------------------------------------
## 6. Population summary (linear-scale pharma units + %CV + OFV)
## ----------------------------------------------------------------------------
cv_pct <- function(omega2) 100 * sqrt(exp(omega2) - 1)
Omega <- fit_focei$Omega
arg   <- fit_focei$argument
pop <- c(
  Ka_pop     = exp(as.numeric(arg["tka"])),
  V_pop      = exp(as.numeric(arg["tv"])),
  Cl_pop     = exp(as.numeric(arg["tcl"])),
  pca0_pop   = exp(as.numeric(arg["tpca0"])),
  kout_pop   = exp(as.numeric(arg["tkout"])),
  IC50_pop   = exp(as.numeric(arg["tic50"])),
  sigma_cp   = exp(as.numeric(arg["log_sigma_cp"])),
  sigma_pca  = exp(as.numeric(arg["log_sigma_pca"])),
  `CV%_V`    = cv_pct(Omega["eta_V",    "eta_V"   ]),
  `CV%_Cl`   = cv_pct(Omega["eta_Cl",   "eta_Cl"  ]),
  `CV%_pca0` = cv_pct(Omega["eta_pca0", "eta_pca0"]),
  `CV%_kout` = cv_pct(Omega["eta_kout", "eta_kout"]),
  OFV        = fit_focei$value)
cat("\n== Population summary ==\n")
print(round(pop, 3))

## ----------------------------------------------------------------------------
## 7. Diagnostic plots (returned as ggplot objects)
##
## Two observables triggers the facet_grid(name ~ condition) layout in
## plotIndivs and the name-as-row layout in plot.nlmeFit / plotResiduals.
## Print any of these to view interactively, or wrap with ggplot2::ggsave()
## if you want to write them out. cowplot is optional; without it,
## plotHistIndivs() returns list(hist, qq) instead of one stacked panel.
## ----------------------------------------------------------------------------
p_indivs   <- plotIndivs   (fit_focei)
p_obs_pred <- plot         (fit_focei)
p_resid    <- plotResiduals(fit_focei)
p_etas     <- plotHistIndivs(fit_focei)

print(p_waterfall)
print(p_indivs)
print(p_obs_pred)
print(p_resid)
print(p_etas)
