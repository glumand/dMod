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


eqns <- c(phi = "v", v = "(-cN * phi - d*v - cR * swtich * (phi - phiR)) / tau", switch = "0")


events <- eventlist() %>% 
  addEvent(var = "phi", time = "s0", value = "F0/tau", mehtod = "add") %>%
  addEvent(var = "switch", time = NA, value = "1", root = "phi - phiR", method = "replace")


mymodel <- odemodel(eqns, events = events, solver = "boost", modelname = "msd", compile = T)

x <- Xs(mymodel)

times <- seq(0, 10, len = 300)

innerpars <- getParameters(x)

steadies <- c(v = 0, phi = 0)

trafo <- eqnvec() %>% 
  define("x~x", x = innerpars) %>% 
  define("switch~0") %>% # initial of switch
  define("x~y", x = names(steadies), y = steadies) %>% 
  insert("x~10^x", x = setdiff(.currentSymbols, "tau")) %>% 
  branch(table = NULL, conditions = c("groc", "health")) %>% 
  insert("tau~1", conditionMatch = "groc") %>% 
  insert("tau~2", conditionMatch = "health")

p <- P(trafo, modelname = "msd", compile = T)

outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)

pouter["s0"] <- log10(2)

p(pouter)


(x*p)(times, pouter) %>% plot()

timesD <- c(0.1, 0.5, 1, 2, 3, 5, 7, 9)

data <- (x*p)(timesD, pouter) %>% as.data.frame() %>% 
  filter(name == "phi") %>% 
  mutate(sigma = 0.1) %>% 
  as.datalist()


obj <- normL2(data, x*p) + constraintL2(pouter, sigma = 2)


obj(pouter)


myfit <- trust(obj, pouter, rinit = 0.1, rmax = 10)


bestfit <- myfit$argument

plot((x*p)(times, bestfit), data)


msout <- mstrust(obj, pouter, sd = 4, fits = 100, cores = 10, studyname = "msd_fitting")

outframe <- as.parframe(msout)
bestfit <- as.parvec(outframe)


plotValues(outframe, tol = 0.1, value < 1e3)


profiles <- profile(obj, bestfit, names(bestfit), method = "optimize", cores = 10)

plotProfile(profiles, mode %in% c("data", "prior"))
