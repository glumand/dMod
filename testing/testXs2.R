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
  addEvent(var = "x", time = "te", value = 1, root = NA, method = "add")

x.ds <- odemodel(eqns, events = events, modelname = "test_Xs_ds", compile = F, solver = "deSolve") %>% Xs()
x.bs <- odemodel(eqns, events = events, modelname = "test_Xs_bs", compile = F, solver = "boost") %>% Xs()

p <- eqnvec() %>% 
  define("x~x", x = getParameters(x.ds)) %>% 
  insert("x~1") %>% P(compile = F)

getEquations(p)

compile(x.ds, x.bs ,p, verbose = F)
# loadDLL(x,p)

getParameters(p)

pars <- c(k=1, te = 1)
# debugonce(p)
p(pars)
p(pars) %>% getDerivs()


times <- seq(0,10,len = 300)
prd.ds <- x.ds*p
prd.bs <- x.bs*p
out.ds <- prd.ds(times,pars)
out.bs <- prd.bs(times,pars)
out.ds %>% plot()
out.bs %>% plot()
# debugonce(getDerivs)

derivs <- (getDerivs(out.bs))[[1]] 

out.ds %>% getDerivs() %>% plot()
out.bs %>% getDerivs() %>% plot()
