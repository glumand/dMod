## Re-fit Theophylline FOCEI with correction = "eager" (the proper FOCEI
## that includes the d log|H|/d theta term in the gradient).

rm(list = ls(all.names = TRUE))
devtools::load_all(quiet = TRUE)
oldwd <- setwd(tempdir()); on.exit(setwd(oldwd), add = TRUE)
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE), force = TRUE)

set.seed(1)
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
subjects <- sort(unique(Theoph$Subject)); N <- length(subjects)
doses <- vapply(subjects, function(s){ r <- Theoph[Theoph$Subject == s,][1,]; r$Dose*r$Wt }, 0.0)
dlist <- as.datalist(data.frame(
  name="y", time=Theoph$Time, value=Theoph$conc, sigma=NA_real_,
  condition=Theoph$Subject, stringsAsFactors=FALSE))

reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")
m <- odemodel(reactions, modelname = "th_ode", compile = TRUE,
              solver = "CppODE", deriv2 = TRUE)
x <- Xs(m)
g <- Y(c(y = "Cc"), x, modelname = "th_obs", compile = TRUE, deriv2 = TRUE)
err <- Y(eqnvec(y = "sigma_add"), g, attach.input = FALSE,
         compile = TRUE, modelname = "th_err")
trafo <- eqnvec(Ka = "Ka_pop * exp(eta_Ka)",
                V  = "V_pop  * exp(eta_V)",
                Cl = "Cl_pop * exp(eta_Cl)",
                Ag = "Ag_init", Cc = "0", sigma_add = "sigma_add")
subj_table <- data.frame(
  eta_Ka  = paste0("eta_Ka_", subjects),
  eta_V   = paste0("eta_V_",  subjects),
  eta_Cl  = paste0("eta_Cl_", subjects),
  Ag_init = doses, row.names = subjects, stringsAsFactors = FALSE)
trafos <- branch(trafo, table = subj_table, apply = "insert")
p <- P(trafos, method = "explicit", compile = TRUE, modelname = "th_p", deriv2 = TRUE)
model <- g * x * p
om <- omega(eta = c("eta_Ka","eta_V","eta_Cl"), subjects = subjects)
joint <- normL2(dlist, model, errmodel = err) + constraintL2(mu = 0, Omega = om)
init <- c(Ka_pop = 1.5, V_pop = 32, Cl_pop = 3.0, sigma_add = 0.5)
init[om$cholPars] <- c(log(0.5), log(0.1), log(0.3))

cat("\n== Try FOCEI with correction = 'eager' (full d log|H|/d theta) ==\n")
em_eager <- emObjfn(joint, om, model=model, data=dlist, errmodel=err,
                    method="laplace",
                    control = list(innerControl = list(rtol = 1e-7, maxit = 30),
                                   correction = "eager"))
t1 <- system.time(fit_eager <- focei(em_eager, init,
                                     trustControl = list(iterlim = 500,
                                                         fterm=1e-10, mterm=1e-10)))
print(t1)
cat(sprintf("  OFV : %.6f\n", fit_eager$value))
cat(sprintf("  converged : %s\n", fit_eager$converged))
cat(sprintf("  iterations : %d\n", fit_eager$iterations))
cat("  pars:\n"); print(fit_eager$argument)
out <- em_eager(fit_eager$argument, deriv = TRUE)
cat(sprintf("  max|grad| : %.4e\n", max(abs(out$gradient))))
cat("  gradient:\n"); print(out$gradient)
cat("  Hessian eigenvalues:\n")
print(eigen(out$hessian, symmetric=TRUE, only.values=TRUE)$values)

cat("\n  Constant correction (-N*K*log(2)):  ", N*om$K*log(2), "\n")
cat(sprintf("  -> OFV - N*K*log(2) = %.6f  (cf nlmixr2's -2logLik=360.35)\n",
            fit_eager$value - N*om$K*log(2)))
