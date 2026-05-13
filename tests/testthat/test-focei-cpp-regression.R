context("C++ FOCEI kernel regression vs pre-rewrite baseline")


# Anchor against the Phase 0 baseline recorded before the C++ kernel landed.
# The fixture stores ($value, $argument) at convergence with eager Stage-2
# correction; the consolidated kernel must reproduce them within Schur/eigen
# tolerance.
fixture_path <- "fixtures/focei_theoph_reference.rds"


test_that("nlmeFit(method='focei') matches the pre-rewrite Theoph baseline", {
  skip_on_cran()
  if (!requireNamespace("CppODE", quietly = TRUE))
    skip("CppODE not installed")
  if (!file.exists(fixture_path))
    skip(paste0("regression fixture not found at ", fixture_path))

  ref <- readRDS(fixture_path)

  set.seed(1)
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
         force = TRUE)

  data(Theoph, package = "datasets")
  Theoph$Subject <- as.character(Theoph$Subject)
  subjects <- sort(unique(Theoph$Subject))

  doses <- vapply(subjects, function(s) {
    rec <- Theoph[Theoph$Subject == s, ][1, ]
    rec$Dose * rec$Wt
  }, 0.0)
  dlist <- as.datalist(data.frame(
    name = "y", time = Theoph$Time, value = Theoph$conc,
    sigma = NA_real_, condition = Theoph$Subject,
    stringsAsFactors = FALSE))

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
  reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
  reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")
  m <- odemodel(reactions, modelname = "theoph_cppreg", compile = TRUE,
                solver = "CppODE", deriv2 = TRUE)
  x <- Xs(m)
  g <- Y(c(y = "Cc"), x, modelname = "theoph_cppreg_obs",
         compile = TRUE, deriv2 = TRUE)
  err <- Y(eqnvec(y = "sigma_add"), g, attach.input = FALSE,
           compile = TRUE, modelname = "theoph_cppreg_err")

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
  p <- P(trafos, method = "explicit", compile = TRUE,
         modelname = "theoph_cppreg_p", deriv2 = TRUE)
  prdfn <- g * x * p

  om <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
  obj <- normL2(dlist, prdfn, errmodel = err, use.bessel = FALSE) +
           constraintL2(mu = 0, Omega = om)

  fit <- nlmeFit(obj, om, ref$init,
                 prdfn    = prdfn,
                 data     = dlist,
                 errfn = err,
                 method   = "focei",
                 control  = list(focei = list(
                   innerControl = list(rtol = 1e-7, maxit = 30),
                   trustControl = list(iterlim = 100))),
                 verbose = FALSE)

  expect_true(fit$converged)
  expect_lt(abs(fit$value - ref$value), 1e-2)
  expect_lt(max(abs(fit$argument[names(ref$argument)] - ref$argument)), 5e-3)
})
