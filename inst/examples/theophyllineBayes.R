## ============================================================================
## Theophylline NLME, Bayesian variant: MAP-Laplace + fixed-metric SMC
## ============================================================================
##
## Companion to inst/examples/theophylline.R. Same 1-compartment oral PK
## chemistry, same lognormal residual, same per-subject trafo. The frequentist
## FOCEI fit (nlmeFit(method = "focei")) is the warm-start for the Bayesian
## run; on top of it we draw a posterior over
##
##     (tka, tv, tcl, log_sigma_add, omega_chol_*)
##
## using sequential SMC. Per-subject etas are integrated out via FOCEI's
## Laplace approximation, so the SMC samples directly from the marginal
## posterior p(theta_struct, omega_chol | y).
##
## Tuning recipe (the one setting that actually gives clean posteriors):
##
##   1. FOCEI MAP first  -> point estimate + Hessian (= -2 log L curvature).
##   2. priorSample      -> Laplace cloud N(MAP, H^{-1}) so particles start
##                          where the likelihood is already concentrated.
##                          One pass of adaptive-ESS bisect then reaches
##                          beta = 1 in ~5-8 levels instead of ~20-30.
##   3. metric = "fixed" -> langevinControl(GFixed = mapfit$hessian) freezes
##                          the preconditioner at MAP curvature. For
##                          well-conditioned posteriors (which Theoph-PK
##                          near the MAP is) this beats both "euclidean"
##                          (ignores correlation) and "riemannFisher"
##                          (recomputes Fisher per step, fragile in flat
##                          regions). dMod stores Hessians in -2 log L
##                          units; the kernel divides by 2 internally.
##   4. essThreshold = 0.7 + malaSteps = 5 -> tight ESS trigger means small
##                          beta steps, more SMC levels, less particle
##                          depletion; 5 inner moves per level lets dual
##                          averaging settle on a good step size.
##
## Wallclock guidance (n_cores = detectCores() - 1):
##   - SCALE_UP = FALSE  ~3-6 min (4 subjects, 200 particles)
##   - SCALE_UP = TRUE   ~15-30 min (12 subjects, 500 particles)
##
## Returns chain, fit_focei, and the ggplot objects p_marginals / p_trace /
## p_pairs in the global env.
## ============================================================================

rm(list = ls(all.names = TRUE))
setwd(tempdir())
set.seed(3333)

library(dMod)

## ----------------------------------------------------------------------------
## 0. Configuration
## ----------------------------------------------------------------------------
SCALE_UP <- FALSE
n_subjects_use <- if (SCALE_UP) 12L else 4L
n_particles    <- if (SCALE_UP) 500L else 200L
mala_steps     <- 5L
ess_threshold  <- 0.7
n_cores        <- max(1L, min(n_particles,
                              parallel::detectCores() - 1L))

## ----------------------------------------------------------------------------
## 1. Data (same as theophylline.R, restricted to the first n_subjects_use)
## ----------------------------------------------------------------------------
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
Theoph_pos <- Theoph[Theoph$Time > 0, , drop = FALSE]
all_subjects <- sort(unique(Theoph_pos$Subject))
subjects <- all_subjects[seq_len(n_subjects_use)]
Theoph_pos <- Theoph_pos[Theoph_pos$Subject %in% subjects, , drop = FALSE]
N <- length(subjects)
cat(sprintf("Theoph (Bayes): N=%d subjects, %d observations\n",
            N, nrow(Theoph_pos)))

doses <- vapply(subjects, function(s) {
  rec <- Theoph[Theoph$Subject == s, ][1, ]
  rec$Dose * rec$Wt
}, 0.0)

data_df <- data.frame(
  name      = "y",         time      = Theoph_pos$Time,
  value     = log(Theoph_pos$conc), sigma     = NA_real_,
  condition = Theoph_pos$Subject,
  stringsAsFactors = FALSE)
dlist <- as.datalist(data_df)

plot(dlist)

## ----------------------------------------------------------------------------
## 2. Prediction chain: 1-cmt oral PK + log observable + additive errfn
## ----------------------------------------------------------------------------
reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")
m <- odemodel(reactions, modelname = "theoph_bayes_ode", compile = F,
              solver = "CppODE", deriv2 = T)
x <- Xs(m)
g <- Y(c(y = "log(Cc + 1e-9)"), x, modelname = "theoph_bayes_obs",
       compile = F, deriv2 = T)
e <- Y(eqnvec(y = "sigma_add"), g, attach.input = FALSE, deriv2 = T,
       compile = F, modelname = "theoph_bayes_err")

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
p <- P(trafos, method = "explicit", compile = F, modelname = "theoph_bayes_p",
       deriv2 = T)
prd <- g * x * p

compile(prd, e, cores = n_cores)

om  <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
obj <- normL2(dlist, prd, errmodel = e, use.bessel = FALSE) +
         constraintL2(mu = 0, Omega = om)

init_outer <- c(tka           = 0.0,
                tv            = 3.0,
                tcl           = 0.7,
                log_sigma_add = log(0.2))
init_outer[om$cholPars] <- rep(log(0.3), length(om$cholPars))

## ----------------------------------------------------------------------------
## 3. FOCEI MAP -- warm-start for the Bayesian run
## ----------------------------------------------------------------------------
cat("\n== nlmeFit(method = 'focei')  (MAP warm-start) ==\n")
fit_focei <- nlmeFit(obj, om, init_outer,
                     prdfn   = prd, data = dlist, errfn = e,
                     method  = "focei",
                     control = list(focei = list(
                       innerControl = list(rtol = 1e-7, maxit = 30),
                       trustControl = list(iterlim = 100))))
print(fit_focei)

## ----------------------------------------------------------------------------
## 4. Laplace cloud + priors centered at the MAP
##
## priorSample draws from N(MAP, H^{-1}) where H is the FOCEI Hessian on the
## -2 log L scale; the covariance is therefore 2 * H^{-1} in log-posterior
## units. The half-Cauchy / LKJ priorOmega still applies independently as the
## hyperprior on omega_chol; we only use the Laplace cloud to *seed* the
## particles, not as a hard prior.
## ----------------------------------------------------------------------------
mapArg    <- fit_focei$argument
mapCov    <- 2 * solve(fit_focei$hessian)         # -2 log L -> log p
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

tgt <- bayesNLMEMarginal(obj           = obj,
                         omegaSpec     = om,
                         prdfn         = prd, data = dlist, errfn = e,
                         init          = init_outer,
                         priorTheta    = priorTheta,
                         priorOmegaObj = priorOmegaObj,
                         priorSample   = priorSample)

## ----------------------------------------------------------------------------
## 5. SMC with fixed (MAP-curvature) Langevin preconditioner
## ----------------------------------------------------------------------------
cat(sprintf("\n== mcmc(SMC, langevin + fixed-MAP-metric) -- %d particles, %d mala/level ==\n",
            n_particles, mala_steps))
t_smc <- system.time(
  chain <- mcmc(target           = tgt,
                 sequenceType     = "sequential",
                 sequenceSchedule = "adaptiveEss",
                 moveType         = "langevin",
                 metric           = "fixed",
                 populationSize   = n_particles,
                 moveControl      = langevinControl(GFixed = fit_focei$hessian),
                 sequenceControl  = smcControl(malaSteps    = mala_steps,
                                                essThreshold = ess_threshold,
                                                verbose      = TRUE),
                 cores            = n_cores)
)
cat(sprintf("[SMC] elapsed = %.1fs, levels = %d, logZ = %.3f\n",
            t_smc[["elapsed"]], chain$nLevels, chain$logEvidence))

## ----------------------------------------------------------------------------
## 6. MAP / posterior comparison
## ----------------------------------------------------------------------------
shared    <- intersect(names(mapArg), colnames(chain$samples))
postMean  <- colMeans(chain$samples)
postSD    <- apply(chain$samples, 2L, sd)

summary_tbl <- data.frame(
  parameter     = shared,
  MAP           = round(mapArg[shared], 4),
  posteriorMean = round(postMean[shared], 4),
  posteriorSD   = round(postSD[shared], 4),
  diffOverSD    = round((postMean[shared] - mapArg[shared]) /
                        pmax(postSD[shared], 1e-8), 2),
  row.names     = NULL)
cat("\n== MAP vs posterior summary ==\n")
print(summary_tbl)
cat(sprintf("\n  logEvidence = %+8.3f\n", chain$logEvidence))
cat(sprintf("  max |Bayes mean - MAP| / SD = %.2f\n",
            max(abs(summary_tbl$diffOverSD))))

## ----------------------------------------------------------------------------
## 7. Diagnostic plots (ggplot objects)
## ----------------------------------------------------------------------------
p_marginals <- plot(chain)
p_trace     <- plotTrace(chain)
p_pairs     <- plotPairs(chain)

print(p_marginals)
print(p_trace)
print(p_pairs)
