rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:1, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)
set.seed(5555)

library(dMod)
library(dplyr)
library(microbenchmark)

eqns <- c(A = "-k*A^2")
# events <- eventlist() %>% 
#   addEvent(var = "A", time = NA, value = "v", root = "A-Acrit", method = "add")

events <- NULL

x_ds <- odemodel(eqns, events = events, solver = "deSolve", 
                 modelname = "m_deSolve", compile = F) %>% Xs(condition = "lsodes", optionsSens = list(method = "lsodes"))
x_tsit <- odemodel(eqns, events = events, solver = "CppODE", 
                   modelname = "m_tsit5", method = "tsit5", compile = F) %>% Xs(condition = "tsit5") 
x_rb4 <- odemodel(eqns, events = events, solver = "CppODE", 
                  modelname = "m_rn4", method = "rb4", compile = F) %>% Xs(condition = "rb4") 
x_ndf <- odemodel(eqns, events = events, solver = "CppODE", 
                  modelname = "m_ndf", method = "bdf", compile = F) %>% Xs(condition = "ndf") 
x_cvode <- odemodel(eqns, events = events, solver = "Sundials", 
                    modelname = "m_covde", method = "bdf", compile = F) %>% Xs(condition = "cvode")


compile(x_ds, x_tsit, x_rb4, x_ndf, x_cvode, cores = 10, verbose = T)

pars = c(A = 1, k = 0.5)
times = seq(0,50, len = 3000)

x_ds(times, pars) %>% plot()
out_ds <- x_ds(times, pars) %>% getDerivs()
out_tsit <- x_tsit(times, pars) %>% getDerivs() 
out_rb4 <- x_rb4(times, pars) %>% getDerivs() 
out_ndf <- x_ndf(times, pars) %>% getDerivs() 
out_cvode <- x_cvode(times, pars) %>% getDerivs() 

c(out_ds, out_tsit, out_rb4, out_ndf, out_cvode) %>% plot()

microbenchmark(
  lsodes = x_ds(times, pars),
  Tsit5 = x_tsit(times, pars),
  Rb4 = x_rb4(times, pars),
  NDF = x_ndf(times, pars),
  CVODE = x_ndf(times, pars),
  times = 100
)

profvis::profvis({
  replicate(1000, x_ds(times, pars))
})
profvis::profvis({
  replicate(1000, x_rb4(times, pars))
})
