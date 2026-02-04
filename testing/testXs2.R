rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:1, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)
# unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE), force = TRUE)

set.seed(5555)

library(dMod)
library(ggplot2)
library(dplyr)


# Set up reactions
eqns <- eqnvec(x = "-k*x")

events <- eventlist() %>% 
  addEvent(var = "x", time = "te", value = "v", root = NA, method = "add")

x.ds <- odemodel(eqns, events = events, modelname = "test_Xs_ds", compile = F, solver = "deSolve") %>% Xs()
x.bs <- odemodel(eqns, events = events, modelname = "test_Xs_bs", compile = F, solver = "boost") %>% Xs()

observables <- eqnvec(obs = "offset-x")

g <- Y(observables, f=x.bs, modelname = "obsfun_test", compile = F)

p <- eqnvec() %>% 
  define("x~x", x = getParameters(x.ds,g)) %>% 
  insert("x~1") %>% P(compile = F)

getEquations(p)
getEquations(x.bs)
getEquations(g)


compile(g, x.ds, x.bs ,p, cores = 10)
# loadDLL(g, x.ds, x.bs ,p)

getParameters(p)

pars <- c(k=1, te = 1, v=0, offset = 1)
# debugonce(p)
p(pars)
p(pars) %>% getDerivs()


times <- seq(0,10,len = 300)
prd.ds <- x.ds*p
prd.bs <- x.bs*p
# debugonce(x.bs)
out.ds <- prd.ds(times,pars)
out.bs <- prd.bs(times,pars)
out.ds %>% plot()
out.bs %>% plot()
# debugonce(getDerivs)

out.ds %>% getDerivs() %>% plot()
out.bs %>% getDerivs() %>% plot()


prd.ds <- g*x.ds*p
# debugonce(g)
out.ds <- prd.ds(times,pars)

out.ds %>% plot()
out.ds %>% getDerivs() %>% plot()

