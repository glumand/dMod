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
r <- eqnvec() %>%
  addReaction("", "Stim", "0") %>%
  addReaction("", "R", "k_act_R_bas + k_act_R * Stim", "stim act R") %>%
  addReaction("R", "", "k_deact_R * R", "deact R") %>% 
  addReaction("A", "pA", "k1 * A * R / (Km + A)",  "phos A") %>%
  addReaction("pA", "A", "k2 * pA / (Km2 + pA)",   "dephos pA")

mysteadies <- steadyStates(r, forcings = "Stim")

p.expl <- eqnvec() %>% 
  define("x~x", x = getParameters(r)) %>% 
  insert("x~y", x = names(mysteadies), y = mysteadies) %>% 
  insert("A~y", y = "totA/(1+k1*k_act_R_bas/(k2*k_deact_R))") %>% 
  insert("x~1", x = "Stim") %>% 
  P(compile = T, modelname = "parfn", condition = "C1")

# p.impl <- P(r, compile = T)

outerpars <- getParameters(p.expl)
pouter <- structure(rep(0.1, length(outerpars)), names = outerpars)
pouter["Stim"] = 0

p.expl(pouter)
p.expl(pouter) %>% getDerivs()


getEquations(p.expl)
