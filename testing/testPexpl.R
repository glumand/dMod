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

mysteadies <- steadyStates(r, forcings = "Stim", resolve = T)

trafo <- eqnvec() %>% 
  define("x~x", x = getParameters(r)) %>% 
  define("Stim~1") %>% 
  insert("x~y", x = names(mysteadies), y = mysteadies) %>% 
  insert("x~10^x", x = .currentSymbols[!grepl("Stim", .currentSymbols)]) %>%
  {.}


p <- P(trafo, compile = T, condition = "cond1", attach.input = T)

outerpars <- setdiff(getParameters(p), "k2")
pars <- structure(rep(-1, length(outerpars)), names = outerpars)
fixed = c(k2 = -1)

pout <- p(pars, fixed = fixed)
pout
p(pars, fixed = fixed) %>% getDerivs()

x.boost <- odemodel(r, modelname = "test_boost", solver = "boost", compile = F) %>% Xs()
x.dS <- odemodel(r, modelname = "test_dS", compile = F) %>% Xs()
compile(x.boost, x.dS, cores = 4)
times <- seq(0,100,len = 300)

prd.bs <- x.boost*p
prd.dS <- x.dS*p
# debugonce(x.dS)
out.bs <- prd.bs(times, pars, fixed = fixed)
out.dS <- prd.dS(times, pars, fixed = fixed)

out.bs %>% plot()
out.dS %>% plot()
getDerivs(out.bs) %>% plot()
getDerivs(out.dS) %>% plot()
