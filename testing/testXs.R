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
eqns <- eqnvec(x = "-k*x")

events <- eventlist() %>% 
  addEvent(var = "x", time = NA, value = "v", root = "x-xcrit", method = "add")

x <- odemodel(eqns, events = events, modelname = "test_Xs", compile = F, solver = "boost", useDenseOutput = F) %>% Xs()

p <- eqnvec() %>% 
  define("x~x", x = getParameters(x)) %>% 
  insert("x~1") %>% P(compile = F)

getEquations(p)

compile(x,p, cores = 10)
# loadDLL(x,p)

getParameters(p)

pars <- c(k=1, xcrit = 0.25, v=1)
# debugonce(p)
p(pars)
p(pars) %>% getDerivs()


times <- seq(0,10,len = 300)
prd <- x*p
# 
out <- prd(times,pars)
out %>% plot()
# debugonce(getDerivs)
out %>% getDerivs() %>% plot()
