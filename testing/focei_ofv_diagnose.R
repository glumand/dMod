## OFV breakdown diagnostic to localise FOCEI/nlmixr2 discrepancy.
## Setup mirrors inst/examples/theophylline.R.

rm(list = ls(all.names = TRUE))
devtools::load_all(quiet = TRUE)

oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
       force = TRUE)

set.seed(1)
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
subjects <- sort(unique(Theoph$Subject))
N <- length(subjects)
doses <- vapply(subjects, function(s) {
  rec <- Theoph[Theoph$Subject == s, ][1, ]; rec$Dose * rec$Wt }, 0.0)

data_df <- data.frame(
  name = "y", time = Theoph$Time, value = Theoph$conc,
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

trafo <- eqnvec(Ka = "Ka_pop * exp(eta_Ka)",
                V  = "V_pop  * exp(eta_V)",
                Cl = "Cl_pop * exp(eta_Cl)",
                Ag = "Ag_init", Cc = "0",
                sigma_add = "sigma_add")
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
joint <- normL2(dlist, model, errmodel = err) +
         constraintL2(mu = 0, Omega = om)

init <- c(Ka_pop = 1.5, V_pop = 32, Cl_pop = 3.0, sigma_add = 0.5)
init[om$cholPars] <- c(log(0.5), log(0.1), log(0.3))

em_lap <- emObjfn(joint, om, model = model, data = dlist, errmodel = err,
                  method = "laplace",
                  control = list(innerControl = list(rtol = 1e-7, maxit = 30)))
fit_focei <- focei(em_lap, init, trustControl = list(iterlim = 100))
cat("\n== dMod FOCEI converged ==\n")
cat(sprintf("  OFV (reported by ecmFit)         : %.6f\n", fit_focei$value))

## At convergence: recompute the OFV breakdown manually.
psi <- fit_focei$argument
eta_modes <- fit_focei$etaModes   # N x K, rows = subjects, cols = base eta
all_eta_names <- as.vector(om$subjectEtas)
eta_full <- as.vector(eta_modes); names(eta_full) <- all_eta_names
full_pars <- c(psi, eta_full)

## Evaluate components separately by reaching into the joint structure.
norm_part <- normL2(dlist, model, errmodel = err)
prior_part <- constraintL2(mu = 0, Omega = om)

normL2_at_modes  <- norm_part(pars = full_pars, deriv = FALSE)
priorL2_at_modes <- prior_part(pars = full_pars, deriv = FALSE)

cat(sprintf("  normL2 value at eta*               : %.6f\n", normL2_at_modes$value))
cat(sprintf("  constraintL2_mvn value at eta*     : %.6f\n", priorL2_at_modes$value))
cat(sprintf("  -> joint value sum                 : %.6f\n",
            normL2_at_modes$value + priorL2_at_modes$value))

## Now compute log|H_i| per subject from the converged inner fit.
final_diag <- fit_focei$emDiag
cat(sprintf("  sum_i log|H_i| (eigen-floored)     : %.6f\n",
            sum(final_diag$logdet)))

## Quick reality check on the inner mode-finding solver convergence
cat(sprintf("  inner solver converged for all i?  : %s\n",
            all(final_diag$converged)))
cat(sprintf("  inner solver max iter              : %d\n",
            max(final_diag$iter, na.rm = TRUE)))

cat(sprintf("\n  TOTAL OFV (joint + sum_logH)       : %.6f\n",
            normL2_at_modes$value + priorL2_at_modes$value +
              sum(final_diag$logdet)))

## Now: count of data points, K, N
T_total <- nrow(Theoph)
K <- om$K
cat(sprintf("\n  N=%d subjects, K=%d random effects, T_total=%d obs\n",
            N, K, T_total))

## What does -2 log p(eta_i | Omega) include?
##   K log(2pi) + log|Omega| + (eta-mu)^T Omega^-1 (eta-mu)
## constraintL2_mvn returns: sum_i quad_i + N log|Omega|.
## So PER SUBJECT it is MISSING K * log(2 pi).
missing_per_subj_prior <- K * log(2 * pi)
cat(sprintf("\n  Per-subject log(2pi)^(K/2) missing from constraintL2_mvn : %.4f\n",
            missing_per_subj_prior))
cat(sprintf("  Total over %d subjects                                    : %.4f\n",
            N, N * missing_per_subj_prior))

## And the Laplace -2 log L formula in our convention:
##   joint(eta*)$value is -2 log p(y, eta*) MISSING N*K*log(2pi).
##   H_i (dMod) = 2 * Hess(-log p)|_eta*, so |H_i|^(1/2) =
##                2^(K/2) * |Hess(-log p)|^(1/2). The standard Laplace formula
##   -2 log L ~ -2 log p(y, eta*) - K log(2pi) + log|Hess(-log p)|
##            = -2 log p(y, eta*) - K log(2pi) + log|H_dMod/2|
##            = -2 log p(y, eta*) - K log(2pi) + log|H_dMod| - K log(2)
##            = -2 log p(y, eta*) + log|H_dMod| - K log(4 pi).
## With dMod's joint missing N*K*log(2pi) for the prior, the correct OFV is:
##   OFV_correct = joint(eta*)$value + N*K*log(2pi)
##                + sum_i log|H_dMod,i| - N*K*log(4pi)
##              = joint(eta*)$value + sum_i log|H_dMod,i| - N*K*log(2)
cat(sprintf("\n  Predicted constant gap N*K*log(2) = %.4f\n", N*K*log(2)))
cat(sprintf("  -> Corrected OFV should be         : %.4f\n",
            fit_focei$value - N*K*log(2)))

cat("\n== Gradient & Hessian at dMod FOCEI optimum ==\n")
out_at_opt <- em_lap(fit_focei$argument, deriv = TRUE)
cat("Gradient (should be near zero at optimum):\n")
print(out_at_opt$gradient)
cat(sprintf("\nMax |gradient|: %.4e\n", max(abs(out_at_opt$gradient))))
cat("\nHessian eigenvalues (positive => local min):\n")
ev_outer <- eigen(out_at_opt$hessian, symmetric = TRUE, only.values = TRUE)$values
print(ev_outer)
cat(sprintf("\nCondition number of Hessian: %.2e\n",
            max(abs(ev_outer)) / min(abs(ev_outer))))

## Try forcing further trust steps from dMod optimum
cat("\n== Try restarting trust from dMod optimum with stricter tols ==\n")
fit2 <- focei(em_lap, fit_focei$argument,
              trustControl = list(iterlim = 200, fterm = 1e-10, mterm = 1e-10))
cat(sprintf("  Re-converged OFV : %.6f  (was %.6f)\n", fit2$value, fit_focei$value))
print(fit2$argument)

## And from nlmixr2's optimum parameter values (back-transformed):
cat("\n== Evaluate dMod em_lap at nlmixr2's optimum ==\n")
psi_nl <- c(Ka_pop = 1.6129351, V_pop = 31.5113788, Cl_pop = 2.7845595,
            sigma_add = 0.7022597)
# nlmixr2 reports BSV(CV%): convert CV% -> omega^2 -> chol-diag
# BSV(CV%) = 100 * sqrt(exp(omega^2) - 1)  =>  omega^2 = log(1 + (CV/100)^2)
om2_nl <- c(log(1 + (68.48155/100)^2),
            log(1 + (10.46610/100)^2),
            log(1 + (28.52731/100)^2))
chol_diag <- 0.5 * log(om2_nl)
psi_nl[om$cholPars] <- chol_diag
out_at_nl <- em_lap(psi_nl, deriv = FALSE)
cat(sprintf("  OFV(dMod_em_lap at nlmixr2 pars): %.6f\n", out_at_nl$value))
cat(sprintf("  Difference vs dMod-optimum     : %.6f\n",
            out_at_nl$value - fit_focei$value))
cat(sprintf("  (if positive: dMod has the better-likelihood optimum,\n",
            "   if negative: dMod missed the global optimum)\n"))

## Now nlmixr2
if (requireNamespace("nlmixr2", quietly = TRUE) &&
    requireNamespace("rxode2",  quietly = TRUE) &&
    requireNamespace("nlmixr2data", quietly = TRUE)) {
  one.cmt <- function() {
    ini({
      tka <- log(1.5); tcl <- log(3.0); tv <- log(32)
      eta.ka ~ 0.25; eta.cl ~ 0.09; eta.v ~ 0.01
      add.sd <- 0.5 })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl/v * center
      cp <- center / v
      cp ~ add(add.sd) }) }
  data(theo_sd, package = "nlmixr2data")
  fit_fn <- get("nlmixr", asNamespace("nlmixr2est"))
  fit_nl <- suppressMessages(fit_fn(one.cmt, theo_sd, est = "focei"))

  cat("\n== nlmixr2 FOCEI ==\n")
  cat(sprintf("  -2 * logLik(fit_nlmixr2)           : %.6f\n",
              -2 * as.numeric(stats::logLik(fit_nl))))
  cat(sprintf("  AIC                                : %.6f\n", AIC(fit_nl)))
  cat(sprintf("  BIC                                : %.6f\n", BIC(fit_nl)))
  if (!is.null(fit_nl$objf))
    cat(sprintf("  NONMEM-style objf                  : %.6f\n", fit_nl$objf))
  if (!is.null(fit_nl$ofv))
    cat(sprintf("  ofv                                : %.6f\n", fit_nl$ofv))
  cat("\n  Parameter estimates (nlmixr2):\n")
  fe <- if (!is.null(fit_nl$parFixedDf)) fit_nl$parFixedDf else fit_nl$parFixed
  print(fe)
  cat("\n  dMod parameter estimates:\n")
  print(fit_focei$argument)

  cat(sprintf("\n  OFV diff (dMod - nlmixr2)          : %.6f\n",
              fit_focei$value - (-2 * as.numeric(stats::logLik(fit_nl)))))
  cat(sprintf("  After -N*K*log(2) correction       : %.6f\n",
              (fit_focei$value - N*K*log(2)) - (-2 * as.numeric(stats::logLik(fit_nl)))))
}
