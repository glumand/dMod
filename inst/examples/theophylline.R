## ============================================================================
## Theophylline NLME example: FOCEI in dMod vs nlmixr2
## ============================================================================
##
## End-to-end NLME walkthrough on the classic Theophylline dataset
## (12 subjects, 1-compartment oral PK). Demonstrates:
##
##   1. Building the joint objfn (normL2 + constraintL2 MVN prior).
##   2. Running FOCEI via the unified API: nlmeFit(method = "focei").
##   3. Diagnostic plots from R/plots.R (camelCase S3 generics):
##      - plotIndivs(fit)      per-subject IPRED/PRED curves + observed dots
##      - plot(fit)            DV vs IPRED + DV vs PRED scatter (plot.nlmeFit)
##      - plotResiduals(fit)   IWRES vs IPRED + IWRES vs TIME
##      - plotHistIndivs(fit)  eta histogram + QQ vs N(0, Omega_kk)
##   4. Cross-tool comparison against nlmixr2 (FOCEI in R), if installed:
##      reported on scale-invariant interpretable metrics, population pars
##      in pharma units (Ka, V, Cl), %CV per random effect, OFV, and the
##      per-subject Cl_i table.
##
## Parameterisation matches nlmixr2 exactly (and the NONMEM/Monolix
## convention): structural parameters live on the natural-log scale,
##
##     Ka_i = exp(tka + eta_Ka_i),   V_i = exp(tv  + eta_V_i),
##     Cl_i = exp(tcl + eta_Cl_i),   sigma_add = exp(log_sigma_add).
##
## The log trafo enforces positivity of Ka_pop, V_pop, Cl_pop, and the
## additive residual SD without any optimiser-side bounds, and keeps omega
## directly comparable to nlmixr2 (same natural-log scale).
##
## Intended use: source() interactively. Generates compiled C/C++ artefacts
## in tempdir(); outputs all plots to <tempdir>/theoph_plots/.
##
## Reference values (Pinheiro & Bates 2000):
##   Ka_pop  ~ 1.50 (1/h)   sd(eta_Ka) ~ 0.5   (% CV ~ 53)
##   V_pop   ~ 32   (L)     sd(eta_V)  ~ 0.1   (% CV ~ 10)
##   Cl_pop  ~ 3.0  (L/h)   sd(eta_Cl) ~ 0.3   (% CV ~ 31)
## ============================================================================

rm(list = ls(all.names = TRUE))

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
} else {
  library(dMod)
}

oldwd <- setwd(tempdir())
on.exit(setwd(oldwd), add = TRUE)
plot_dir <- file.path(tempdir(), "theoph_plots")
dir.create(plot_dir, showWarnings = FALSE)
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
       force = TRUE)
cat("Plots will be saved to: ", plot_dir, "\n")

set.seed(1)

## ----------------------------------------------------------------------------
## 1. Data
## ----------------------------------------------------------------------------
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
subjects <- sort(unique(Theoph$Subject))
N <- length(subjects)
cat(sprintf("Theoph: N=%d subjects, %d observations\n", N, nrow(Theoph)))

doses <- vapply(subjects, function(s) {
  rec <- Theoph[Theoph$Subject == s, ][1, ]
  rec$Dose * rec$Wt
}, 0.0)

# NB: sigma = NA so the errmodel below estimates the additive residual SD
# (matches the nlmixr2 / NONMEM convention). Hard-coding sigma in the data
# table would freeze it and make the cross-tool OFV comparison apples-to-
# oranges.
data_df <- data.frame(
  name      = "y",        time      = Theoph$Time,
  value     = Theoph$conc, sigma     = NA_real_,
  condition = Theoph$Subject,
  stringsAsFactors = FALSE)
dlist <- as.datalist(data_df)

## ----------------------------------------------------------------------------
## 2. ODE model: 1-compartment oral PK
##      dAg/dt = -Ka * Ag                  (gut depot)
##      dCc/dt =  Ka * Ag / V - Cl/V * Cc  (central concentration)
## ----------------------------------------------------------------------------
reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")
m <- odemodel(reactions, modelname = "theoph_ode", compile = TRUE,
              solver = "CppODE", deriv2 = TRUE)
x <- Xs(m)
g <- Y(c(y = "Cc"), x, modelname = "theoph_obs", compile = TRUE, deriv2 = TRUE)

# Additive residual error model: sigma_add is the SD of an additive Gaussian
# residual on the observed concentration. Estimated alongside the structural
# parameters and omegas, matching the nlmixr2 `add.sd` convention.
err <- Y(eqnvec(y = "sigma_add"), g, attach.input = FALSE,
         compile = TRUE, modelname = "theoph_err")

## ----------------------------------------------------------------------------
## 3. Per-subject parameter trafos (log-scale, positivity guaranteed)
##      theta_i = exp(t_theta + eta_i)
##      sigma_add = exp(log_sigma_add)
##    Matches nlmixr2 / NONMEM / Monolix exactly so omega is directly
##    comparable on natural-log scale.
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
  Ag_init = doses, row.names = subjects, stringsAsFactors = FALSE)
trafos <- branch(trafo, table = subj_table, apply = "insert")
p <- P(trafos, method = "explicit", compile = TRUE, modelname = "theoph_p",
       deriv2 = TRUE)
model <- g * x * p

om <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
# use.bessel = FALSE: Bessel correction is a finite-sample n/(n-p) inflation
# meant for unbiased sigma estimation in *fixed-effects* regression. In NLME
# the marginal likelihood is the right objective (no Bessel), and the effec-
# tive degrees of freedom shrink with the random-effects prior. NONMEM /
# Monolix / nlmixr2 all run plain ML; keeping Bessel on here would shift the
# OFV systematically and break the cross-tool comparison.
joint <- normL2(dlist, model, errmodel = err, use.bessel = FALSE) +
         constraintL2(mu = 0, Omega = om)

## Starting values close to (but not at) the Pinheiro-Bates reference:
##   Ka_pop  = exp(0.0)  = 1.00  (ref 1.5)
##   V_pop   = exp(3.0)  = 20.1  (ref 32)
##   Cl_pop  = exp(0.7)  = 2.01  (ref 3.0)
##   sigma_add = exp(-0.7) = 0.50  (ref ~0.7)
## eta-SDs at 0.3 -> cholPars (diag of L, exp-parametrised) = log(0.3).
init <- c(tka           = 0.0,
          tv            = 3.0,
          tcl           = 0.7,
          log_sigma_add = log(0.5))
init[om$cholPars] <- rep(log(0.3), length(om$cholPars))

## ----------------------------------------------------------------------------
## 4. FOCEI fit
## ----------------------------------------------------------------------------
cat("\n== nlmeFit(method = 'focei') ==\n")
t_focei <- system.time(
  fit_focei <- nlmeFit(joint, om, init,
                       model    = model,
                       data     = dlist,
                       errmodel = err,
                       method   = "focei",
                       control  = list(focei = list(
                         innerControl = list(rtol = 1e-7, maxit = 30),
                         trustControl = list(iterlim = 100))))
)
print(t_focei); print(fit_focei)

## ----------------------------------------------------------------------------
## 5. Cross-tool comparison: nlmixr2 FOCEI run
##
## We prefer the canonical dataset from `nlmixr2data::theo_sd` (pre-formatted
## with EVID/CMT/AMT/DV columns) over building events manually. The model
## syntax follows the canonical nlmixr2 vignette: log-transformed fixed
## effects, eta variances in `~`, additive residual. Starting values match
## the dMod fit above for an apples-to-apples comparison.
##
## If nlmixr2 isn't installed (e.g. symengine/MPFR system deps missing), the
## section is skipped and the final comparison table shows only the dMod
## column.
## ----------------------------------------------------------------------------
fit_nlmixr2 <- NULL
have_nlmixr2 <- requireNamespace("nlmixr2",     quietly = TRUE) &&
                requireNamespace("rxode2",      quietly = TRUE) &&
                requireNamespace("nlmixr2data", quietly = TRUE)
if (have_nlmixr2) {
  cat("\n== nlmixr2 FOCEI ==\n")
  one.cmt <- function() {
    ini({
      # Same starting values as the dMod fit above (log-scale fixed effects,
      # eta variance = 0.3^2 = 0.09, residual SD = 0.5). Cross-tool
      # comparison is apples-to-apples.
      tka <- 0.0
      tv  <- 3.0
      tcl <- 0.7
      eta.ka ~ 0.09
      eta.v  ~ 0.09
      eta.cl ~ 0.09
      add.sd <- 0.5
    })
    model({
      ka <- exp(tka + eta.ka)
      v  <- exp(tv  + eta.v )
      cl <- exp(tcl + eta.cl)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl/v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  data(theo_sd, package = "nlmixr2data")
  # theo_sd already has columns: ID, TIME, DV, AMT, EVID, CMT, WT (canonical).
  # The fit function exported by package `nlmixr2` is named `nlmixr` (legacy
  # , the package keeps the original nlmixr function name). Fall back to
  # `nlmixr2est::nlmixr2` for older / forked layouts.
  fit_fn <- if (exists("nlmixr", asNamespace("nlmixr2"), inherits = FALSE))
              get("nlmixr", asNamespace("nlmixr2"))
            else if (requireNamespace("nlmixr2est", quietly = TRUE) &&
                     exists("nlmixr2", asNamespace("nlmixr2est"), inherits = FALSE))
              get("nlmixr2", asNamespace("nlmixr2est"))
            else stop("Could not find a nlmixr/nlmixr2 fitting function.")
  t_nlmixr2 <- system.time(
    fit_nlmixr2 <- try(
      suppressMessages(fit_fn(one.cmt, theo_sd, est = "focei")),
      silent = TRUE)
  )
  if (inherits(fit_nlmixr2, "try-error")) {
    cat("nlmixr2 fit failed:\n"); print(attr(fit_nlmixr2, "condition")$message)
    fit_nlmixr2 <- NULL
  } else {
    cat("nlmixr2 elapsed (s):\n"); print(t_nlmixr2)
    cat("nlmixr2 fit OK. parFixed:\n"); print(fit_nlmixr2$parFixed)
  }
} else {
  cat("\n(nlmixr2 stack not installed; cross-tool comparison column skipped.)\n")
}

## ----------------------------------------------------------------------------
## 6. Per-subject Cl extraction (interpretable, scale-invariant)
## ----------------------------------------------------------------------------
# dMod fits carry etaModes (N x K) on natural-log scale. Backtransform via
#   Cl_i = exp(tcl + eta_Cl_i).
dmod_ind_Cl <- function(fit) {
  tcl    <- as.numeric(fit$argument["tcl"])
  eta_Cl <- fit$etaModes[, "eta_Cl"]
  setNames(exp(tcl + eta_Cl), rownames(fit$etaModes))
}
Cl_ind <- data.frame(
  Subject  = rownames(fit_focei$etaModes),
  Cl_FOCEI = round(dmod_ind_Cl(fit_focei), 3)
)
if (!is.null(fit_nlmixr2)) {
  # nlmixr2 stores etas in fit_nlmixr2$eta (per-ID, columns eta.ka/eta.v/eta.cl).
  # tcl is the fixed-effect on log scale; Cl_pop = exp(tcl). Column name for
  # the point estimate is "Est" in current nlmixr2; "Estimate" in some older
  # vintages, fall back gracefully.
  # Prefer parFixedDf (numeric) over parFixed (formatted character cols).
  fe_num <- if (!is.null(fit_nlmixr2$parFixedDf)) fit_nlmixr2$parFixedDf
            else fit_nlmixr2$parFixed
  est_candidates <- c("Est.", "Est", "Estimate")
  est_col <- est_candidates[est_candidates %in% colnames(fe_num)][1]
  if (is.na(est_col))
    stop("Unknown parFixed column layout: ",
         paste(colnames(fe_num), collapse = ", "))
  tcl_nl <- as.numeric(fe_num["tcl", est_col])
  eta_df <- fit_nlmixr2$eta
  eta_df <- eta_df[order(as.character(eta_df$ID)), , drop = FALSE]
  Cl_ind$Cl_nlmixr2 <- round(exp(tcl_nl + eta_df$eta.cl), 3)
}
cat("\n== Per-subject Cl (L/h) on linear scale ==\n")
print(Cl_ind, row.names = FALSE)

## ----------------------------------------------------------------------------
## 7. Side-by-side population summary (linear pop pars + %CV + OFV)
## ----------------------------------------------------------------------------
# %CV for lognormal eta with variance omega2 (natural-log):
#   CV = sqrt(exp(omega2) - 1) * 100%.
# Scale-invariant: any tool reports CV the same way regardless of internal
# log10 vs natural-log parameterisation.
cv_pct <- function(omega2) 100 * sqrt(exp(omega2) - 1)

dmod_summary <- function(fit) {
  Omega <- fit$omega
  omegas2 <- diag(Omega)
  c(Ka_pop    = exp(as.numeric(fit$argument["tka"])),
    V_pop     = exp(as.numeric(fit$argument["tv"])),
    Cl_pop    = exp(as.numeric(fit$argument["tcl"])),
    sigma_add = exp(as.numeric(fit$argument["log_sigma_add"])),
    `CV%_Ka`  = cv_pct(omegas2[1]),
    `CV%_V`   = cv_pct(omegas2[2]),
    `CV%_Cl`  = cv_pct(omegas2[3]),
    OFV       = fit$value)
}
tab <- cbind(FOCEI = dmod_summary(fit_focei))
if (!is.null(fit_nlmixr2)) {
  fe_num <- if (!is.null(fit_nlmixr2$parFixedDf)) fit_nlmixr2$parFixedDf
            else fit_nlmixr2$parFixed
  est_candidates <- c("Est.", "Est", "Estimate")
  est_col <- est_candidates[est_candidates %in% colnames(fe_num)][1]
  if (is.na(est_col))
    stop("Unknown parFixed column layout: ",
         paste(colnames(fe_num), collapse = ", "))
  om_nl <- fit_nlmixr2$omega
  nl_row <- c(Ka_pop    = exp(as.numeric(fe_num["tka", est_col])),
              V_pop     = exp(as.numeric(fe_num["tv",  est_col])),
              Cl_pop    = exp(as.numeric(fe_num["tcl", est_col])),
              sigma_add = as.numeric(fe_num["add.sd", est_col]),
              `CV%_Ka`  = cv_pct(om_nl["eta.ka","eta.ka"]),
              `CV%_V`   = cv_pct(om_nl["eta.v", "eta.v" ]),
              `CV%_Cl`  = cv_pct(om_nl["eta.cl","eta.cl"]),
              # Use -2 * logLik for the maximum-likelihood value. This INCLUDES
              # the n*log(2*pi) Gaussian normalisation that nlmixr2's OBJF
              # (NONMEM convention) drops. Now directly comparable to dMod's
              # OFV without any additive-constant correction.
              OFV       = as.numeric(-2 * stats::logLik(fit_nlmixr2)))
  tab <- cbind(tab, nlmixr2 = nl_row)
}
cat("\n== Population summary (linear scale, %CV, OFV) ==\n")
print(round(tab, 3))
cat("\nNote on the OFV column:\n",
    "  Both values are the full -2 log L of the marginal likelihood,\n",
    "  including the n*log(2*pi*sigma^2) Gaussian normalisation. For nlmixr2\n",
    "  we use -2 * logLik(fit) (NOT the NONMEM-style objf, which drops the\n",
    "  n*log(2*pi) constant and would be ~242 OFV units lower).\n",
    sep = "")

## ----------------------------------------------------------------------------
## 8. Diagnostic plots
## ----------------------------------------------------------------------------
cat("\n== Diagnostic plots ==\n")
save_plot <- function(p, filename, w = 8, h = 5) {
  if (inherits(p, "list") && !inherits(p, "ggplot")) {
    fp1 <- file.path(plot_dir, sub("\\.png$", "_hist.png", filename))
    fp2 <- file.path(plot_dir, sub("\\.png$", "_qq.png",   filename))
    ggplot2::ggsave(fp1, p$hist, width = w, height = h, dpi = 110)
    ggplot2::ggsave(fp2, p$qq,   width = w, height = h, dpi = 110)
    cat("  wrote", fp1, "and", fp2, "\n")
  } else {
    fp <- file.path(plot_dir, filename)
    ggplot2::ggsave(fp, p, width = w, height = h, dpi = 110)
    cat("  wrote", fp, "\n")
  }
}

save_plot(plotIndivs    (fit_focei), "focei_individuals.png", w = 10, h = 7)
save_plot(plot          (fit_focei), "focei_obs_vs_pred.png", w = 8,  h = 5)
save_plot(plotResiduals (fit_focei), "focei_residuals.png",   w = 9,  h = 5)
save_plot(plotHistIndivs(fit_focei), "focei_etas.png",        w = 9,  h = 7)

cat(sprintf("\nDone. dMod FOCEI elapsed = %.1fs.\n", t_focei["elapsed"]))
if (!is.null(fit_nlmixr2))
  cat(sprintf("       nlmixr2 FOCEI elapsed = %.1fs   (ratio dMod/nlmixr2 = %.2f).\n",
              t_nlmixr2["elapsed"], t_focei["elapsed"] / t_nlmixr2["elapsed"]))
cat("Open plots from: ", plot_dir, "\n")
