rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:1, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)
set.seed(5555)

## Load Libraries ---
library(dMod)
library(dplyr)
library(ggplot2)


reactions <- eqnlist() %>% 
  addReaction("Prt", "pPrt", "vm1 * Prt/(Km1 + Prt)") %>% 
  addReaction("pPrt", "Prt", "vm2 * pPrt/(Km2 + pPrt)")


mysteadies <- steadyStates(reactions)

totExpress <- c(Prt = "alpha * totPrt", pPrt = "(1 - alpha)*totPrt")

trafo <- eqnvec() %>% 
  define("x~x", x = getParameters(reactions)) %>% 
  define("x~y", x = names(mysteadies), y = mysteadies) %>% 
  insert("x~y", x = names(totExpress), y = totExpress) %>% 
  insert("alpha~1/2 + 1/pi * arctan(z)")


trafo_eqns <- unclass(trafo)
derivs <- CppODE::derivSymb(trafo_eqns, real = T)
