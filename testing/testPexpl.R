rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:1, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE), force = TRUE)

set.seed(5555)

library(dMod)
library(ggplot2)
library(dplyr)


# Set up reactions
r <- eqnvec() %>%
  addReaction("", "Stim", "0") %>%
  addReaction("", "R", "k_act_R_bas + k_act_R * Stim", "stim act R") %>%
  addReaction("R", "", "k_deact_R * R", "deact R") %>% 
  addReaction("A", "pA", "k1 * A * R / (Km + A)",  "phos A") %>%
  addReaction("pA", "A", "k2 * pA / (Km2 + pA)",   "dephos pA")

mysteadies <- steadyStates(r, forcings = "Stim")

trafo <- eqnvec() %>% 
  define("x~x", x = getParameters(r)) %>% 
  define("Stim~1") %>% 
  insert("x~y", x = names(mysteadies), y = mysteadies) %>% 
  insert("A~y", y = "totA/(1+k1*k_act_R_bas/(k2*k_deact_R))") %>% 
  insert("x~10^x", x = .currentSymbols[!grepl("Stim", .currentSymbols)]) %>%
  {.}


p <- P(trafo, compile = T, condition = "cond1")

outerpars <- getParameters(p)
pars <- structure(rep(-1, length(outerpars)), names = outerpars)

# debugonce(p)
pout <- p(pars)
p(pars) %>% getDerivs()


x <- odemodel(r, modelname = "test") %>% Xs()

times <- seq(0,10,len = 100)

prd <- x*p
debugonce(x)
out <- prd(times, pars)
getDerivs(out) %>% plot()
