rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:1, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)
set.seed(5555)

library(dMod)
library(ggplot2)
library(dplyr)


# Set up reactions
r <- eqnlist() %>%
  addReaction("", "Stim", "0") %>%
  addReaction("", "R", "k_act_R_bas + k_act_R * Stim", "stim act R") %>%
  addReaction("R", "", "k_deact_R * R", "deact R") %>% 
  addReaction("A", "pA", "k1 * A * R / (Km + A)",  "phos A") %>%
  addReaction("pA", "A", "k2 * pA / (Km2 + pA)",   "dephos pA")

mysteadies <- steadyStates(r, forcings = "Stim")

p.expl <- eqnvec() %>% 
  define("x~x", x = getParameters(r)) %>% 
  insert("x~y", x = names(mysteadies), y = mysteadies) %>% 
  insert("x~1", x = "Stim") %>% 
  P(compile = T, modelname = "parfn", condition = "C1")


outerpars <- getParameters(p.expl)
pouter <- structure(rep(0.1, length(outerpars)), names = outerpars)

p.expl(pouter)
p.expl(pouter) %>% getDerivs()
getEquations(p.expl)


r.ss <- eqnlist() %>%
  addReaction("", "R", "k_act_R_bas", "stim act R") %>%
  addReaction("R", "", "k_deact_R * R", "deact R") %>% 
  addReaction("A", "pA", "k1 * A * R / (Km + A)",  "phos A") %>%
  addReaction("pA", "A", "k2 * pA / (Km2 + pA)",   "dephos pA")

conservedQuantities(r.ss$smatrix)

f.ss <- as.eqnvec(r.ss)

replacement <- c(pA = "A + pA - totA")
f.ss[names(replacement)] <- replacement

p.impl <- P(f.ss, parameters = "totA", method = "implicit", compile = T, condition = "C1")
outerpars <- getParameters(p.impl)
pouter <- structure(rep(0.1, length(outerpars)), names = outerpars)
# pouter["Stim"] = 0
p.impl(pouter)
p.impl(pouter[-1], fixed = c(k_act_R_bas = 0.1)) %>% getDerivs()


p.log <- eqnvec() %>% 
  define("x~x", x = outerpars) %>% 
  insert("x~exp10(x)", x = .currentSymbols) %>% 
  P(compile = T, modelname = "parfnlog", condition = "C1")


p <- p.impl*p.log
outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)

p(pouter[-1], fixed = c(k_act_R_bas = 0.1)) %>% getDerivs()
