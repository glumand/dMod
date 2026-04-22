## End-to-end smoke test for the new derivMode arg in Pexpl/Y.
## Builds the README enzyme-kinetics example twice (default + swapped derivModes)
## and compares the objective gradient to confirm both paths agree.

setwd(tempdir())  # so generated C++ files don't pollute the repo
library(dMod)

run_once <- function(pexpl_mode, y_mode, label) {
  cat("\n==== ", label, "  (Pexpl=", pexpl_mode, ", Y=", y_mode, ") ====\n", sep = "")

  f <- NULL
  f <- dMod::addReaction(f, from = "Enz + Sub", to = "Compl",
                         rate = "k1*Enz*Sub", description = "production of complex")
  f <- dMod::addReaction(f, from = "Compl", to = "Enz + Sub",
                         rate = "k2*Compl", description = "decay of complex")
  f <- dMod::addReaction(f, from = "Compl", to = "Enz + Prod",
                         rate = "k3*Compl", description = "production of product")
  f <- dMod::addReaction(f, from = "Enz", to = "",
                         rate = "k4*Enz", description = "enzyme degradation")

  model <- dMod::odemodel(f, modelname = paste0("ek_", label))
  x <- dMod::Xs(model)

  observables <- dMod::eqnvec(
    product   = "Prod",
    substrate = "(Sub + Compl)",
    enzyme    = "(Enz + Compl)"
  )
  g <- dMod::Y(observables, x, compile = TRUE,
               modelname = paste0("ek_obs_", label),
               attach.input = FALSE,
               derivMode = y_mode)

  innerpars <- dMod::getParameters(g * x)
  trafo <- dMod::repar("x~x", x = innerpars)
  trafo <- dMod::repar("x~0", x = c("Compl", "Prod"), trafo)
  trafo <- dMod::repar("x~exp(x)", x = innerpars, trafo)
  trafo1 <- trafo2 <- trafo
  trafo1["k4"] <- "0"

  p <- NULL
  p <- p + dMod::Pexpl(trafo1, condition = "noDegradation",
                   compile = TRUE, modelname = paste0("ek_p1_", label), derivMode = pexpl_mode)
  p <- p + dMod::Pexpl(trafo2, condition = "withDegradation",
                   compile = TRUE, modelname = paste0("ek_p2_", label), derivMode = pexpl_mode)

  set.seed(1)
  outerpars <- dMod::getParameters(p)
  pouter <- structure(rnorm(length(outerpars), -2, .5), names = outerpars)

  data <- dMod::datalist(
    noDegradation = data.frame(
      name  = c("product","product","product","substrate","substrate","substrate"),
      time  = c(0,25,100,0,25,100),
      value = c(0.0025,0.2012,0.3080,0.3372,0.1662,0.0166),
      sigma = 0.02),
    withDegradation = data.frame(
      name  = c("product","product","product","substrate","substrate","substrate","enzyme","enzyme","enzyme"),
      time  = c(0,25,100,0,25,100,0,25,100),
      value = c(-0.0301,0.1512,0.2403,0.3013,0.1635,0.0411,0.4701,0.2001,0.0383),
      sigma = 0.02)
  )

  prior <- structure(rep(0, length(pouter)), names = names(pouter))
  obj <- dMod::normL2(data, g * x * p) + dMod::constraintL2(mu = prior, sigma = 10)

  runtime <- system.time({o <- obj(pouter)})
  cat("value:    ", o$value, "\n")
  cat("run time: ", runtime["elapsed"], "\n")
  cat("gradient: ", paste(format(o$gradient, digits = 8), collapse = " "), "\n")
  invisible(o)
}

a <- run_once("symbolic", "ad",       "default")
b <- run_once("symbolic", "symbolic", "Y_symb")
c <- run_once("ad",       "ad",       "Pexpl_ad")
d <- run_once("ad",       "symbolic", "all_swapped")

cat("\n==== diffs vs default ====\n")
cat("|val(b)-val(a)|: ", abs(b$value - a$value), "\n")
cat("|val(c)-val(a)|: ", abs(c$value - a$value), "\n")
cat("|val(d)-val(a)|: ", abs(d$value - a$value), "\n")
cat("max|grad(b)-grad(a)|: ", max(abs(b$gradient - a$gradient)), "\n")
cat("max|grad(c)-grad(a)|: ", max(abs(c$gradient - a$gradient)), "\n")
cat("max|grad(d)-grad(a)|: ", max(abs(d$gradient - a$gradient)), "\n")

