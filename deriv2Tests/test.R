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

## Load data ---
data(badata)
data <- badata %>% subset(condition == "closed") %>% as.datalist()
plot(data)

## Set up the ODE model ---

# Define via reactions
reactions <- eqnlist() %>%
  addReaction("TCA_buffer", "TCA_cell",  rate = "k_import*TCA_buffer", description = "Uptake") %>%
  addReaction("TCA_cell", "TCA_buffer",  rate = "k_export_sinus*TCA_cell", description = "Sinusoidal export") %>%
  addReaction("TCA_cell", "TCA_cana",    rate = "k_export_cana*TCA_cell", description = "Canalicular export") %>%
  addReaction("TCA_cana", "TCA_buffer",  rate = "k_reflux*TCA_cana", description = "Reflux into the buffer")

# Define via ODE system (equivalent)

# f <- eqnvec(TCA_buffer = "-k_import * TCA_buffer + k_export_sinus * TCA_cell + k_reflux * TCA_cana",
#             TCA_cana = "k_export_cana * TCA_cell - k_reflux * TCA_cana",
#             TCA_cell = "k_import * TCA_buffer - k_export_sinus * TCA_cell - k_export_cana * TCA_cell")

# Translate reactions into ODE model object
mymodel <- odemodel(reactions, modelname = "bamodel", compile = F, solver = "boost")
# Generate trajectories for the default condition
x <- Xs(mymodel)

# Define observables buffer and cellular
observables <- eqnvec(buffer = "s*TCA_buffer", cellular = "s*(TCA_cana + TCA_cell)")
g <- Y(observables, f = x, condition = NULL, compile = F, modelname = "obsfn_bamodel", attach.input = T)

# Define parameter transformations using define(), insert() and branch(). Old function repar also avaiable!
innerpars <- getParameters(x,g)
trafo <- NULL %>%
  define("x~x", x = innerpars) %>% # identity
  define("TCA_buffer~0") %>%
  insert("x~10^y", x = .currentSymbols, y = toupper(.currentSymbols))


# # # Explicit trafo (this is equivalent to the lines above)
# trafo <- eqnvec(TCA_buffer = "0",
#                 TCA_cell = "exp(log(10)*TCA_cell)",
#                 TCA_cana = "exp(log(10)*TCA_cana)",
#                 k_import = "exp(log(10)*k_import)",
#                 k_export_sinus = "exp(log(10)*k_export_sinus)",
#                 k_export_cana = "exp(log(10)*k_export_cana)",
#                 k_reflux = "exp(log(10)*k_reflux)",
#                 s = "exp(log(10)*s)")

p <- P(trafo, condition = "closed", compile = F)


# Compile the objects
compile(g, x, p, output = "bamodelSO", cores = 8) # Compile C/C++ output of odemodel in parallel

## Use simulate data to calibrate outer model parameters ---
outerpars <- getParameters(p)
pouter <- structure(rep(-1,length(outerpars)), names = outerpars)

prd <- g*x*p
# debugonce(x)
times <- seq(0, 45, len = 300)
# debugonce(g)
out <- prd(times, pouter)
out %>% plot()
out %>% getDerivs() %>% plot()

myderivs <- attr(out$closed, "deriv")
# Define objective function
obj <- normL2(data, g * x * p)

# Test objective function with and without explicit calculation of second derivatives
obj(pouter)

# Fit on time (starting from pouter)
myfit <- trust(obj, pouter, rinit = 0.1, rmax = 10, iterlim = 500, printIter = T)
times <- seq(0, 45, len = 300)
mypred <- (g * x * p)(times, myfit$argument)
plot(mypred, data)

obj(myfit$argument)

system.time({obj(myfit$argument)})
## Handling different experimental conditions

# Parameter Trafo, usage of "+" operator for trafo functions (output of P())
data <- as.datalist(badata) # full data
plot(data)

trafo <- getEquations(p, conditions = "closed") %>% 
  insert("K_REFLUX~K_REFLUX_OPEN")
p <- p + P(trafo, condition = "open", compile = T)

outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
p(pouter)
(g*x*p)(times, pouter) %>% plot(data)

# Objective function
obj <- normL2(data, g * x * p) + constraintL2(pouter, sigma = 4)

# Evaluation of obj at pouter
system.time({obj(pouter)})

myfit <- trust(obj, pouter, rinit = 0.1, rmax = 5, iterlim = 500, printIter = T)

# Fit 50 times, sample with sd=4 around pouter
out_frame <- mstrust(obj, pouter, sd = 4, studyname = "bamodel", cores=detectFreeCores(), fits=100, iterlim = 1e3)

outknecht <- runbg({
  mstrust(obj, pouter, sd = 4, studyname = "bamodel", cores=10, fits=100, iterlim = 1e3)
}, machine = "knecht3", filename = "testJoschi")


out_frame <- as.parframe(out_frame)
plotValues(out_frame) # Show "Waterfall" plot
plotPars(out_frame) # Show parameter plot
bestfit <- as.parvec(out_frame)


# Plot predictions along data
plot((g*x*p)(times, bestfit), data)
# 
# # Plot sensis
plot(getDerivs((g*x*p)(times, bestfit)))
# 
# Calculate Parameter Profiles and plot different contributions (for identifiablility only "data" is of interest)
profiles_integrate <- profile(obj, bestfit, whichPar = names(bestfit), method = "integrate", cores = 10, limits = c(lower = -5, upper = 5), 
                              stepControl = list(stop = "data"))

profiles_optimize <- profile(obj, bestfit, whichPar = names(bestfit), method = "optimize", cores = 10, limits = c(lower = -5, upper = 5), 
                             stepControl = list(stepsize = 1e-4, min = 1e-4, max = Inf, atol = 1e-2, rtol = 1e-2, limit = 200, stop = "data"))

proflist <- list(integrate = profiles_integrate, optimize = profiles_optimize) # The best tactic is to use method = "integrate" with reoptimize = TRUE in algoControl
plotProfile(proflist, mode %in% c("data", "prior"))

# plotProfile(profiles_integrate, mode %in% c("data", "prior"))
# plotProfile(profiles_optimize, mode %in% c("data", "prior"))
# plotPaths(profiles, whichPar = "TCA_CANA")
# plotPaths(profiles, whichPar = "K_EXPORT_CANA")
# plotPaths(profiles, whichPar = "K_REFLUX_OPEN")
# plotPaths(profiles, whichPar = "S") # The triple S, TCE_CELL, TCA_CANA compensate by structure of the model


# Tighten model assumptions with steady state constraint
mysteadies <- steadyStates(reactions, forcings = "k_import")

trafo <- eqnvec() %>%
  define("x~x", x = innerpars) %>% # identity
  define("TCA_buffer~0") %>%
  define("x~y", x = names(mysteadies), y = mysteadies) %>% 
  insert("x~10^y", x = .currentSymbols, y = toupper(.currentSymbols)) %>% 
  branch(conditions = c("closed", "open")) %>% 
  define("k_reflux~10^K_REFLUX_OPEN", conditionMatch = "open") %>% 
  insert("S~0") # fixed structural non identifiablility
  
# debugonce(P)
p <- P(trafo, modelname = "bamodel_SS", compile = T)

outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
p(pouter)
plot((g*x*p)(times, pouter),data)

# Objective function
obj <- normL2(data, g * x * p, attr.name = "data") + constraintL2(pouter, sigma = 20, attr.name = "prior")

# Multistart fit
outms <- mstrust(obj, pouter, sd = 4, iterlim = 1e3, studyname = "bamodel_ss", cores=20, fits=100)
outframe <- as.parframe(outms)
plotValues(outframe) # Show "Waterfall" plot
plotPars(outframe) # Show parameter plot
bestfit <- as.parvec(outframe)
pred <- (g * x * p)(times, bestfit)
plot(pred, data)


# Calculate Parameter Profiles
profiles <- profile(obj, bestfit, whichPar = names(bestfit), method = "integrate", cores = 10, limits = c(lower = -5, upper = 5), 
                    stepControl = list(stop = "data"),
                    algoControl = list(gamma = 1, reoptimize = T))
plotProfile(profiles,mode %in% c("data", "prior"))
plotPaths(profiles, whichPar = "K_REFLUX_OPEN")

# Note: The parameter "reflux_open" is still practical non-identifiable
# The profile is open to the left -> A possible reduction of the model would be the assumption of an immediate reflux, i.e. the limit case reflux_open -> infinity
# An pragmatic but ugly way to circumvent the reformulation of the ODE system is to fix the reflux_open parameter to a high value. E.g. 1e3

# One could also check the models ability to produce reliable predictions, by the calculation of prediction uncertainty with profile likelihood
# The calculation of prediction confidence intervals is done at next

## Prediction uncertainty taken from validation profile --------------------------------------------------------------------------

# choose sigma below 1 percent of the prediction in order to pull the prediction strongly towards d1
obj.validation <- normL2(data, g * x * p, times = c(20), attr.name = "data") +
  datapointL2(name = "TCA_cell", time = 20, value = "v", sigma = 1, attr.name = "validation", condition = "closed")

# If sigma is not known, and you therefore decide to calculate prediction confidence intervals, just choose a very small sigma, in order to "pull strongly" on the trajectory
obj.validation(c(v = 180, bestfit))

# refit
myfit <- trust(obj.validation, parinit = c(v = 190, bestfit), rinit = 1, rmax = 10, iterlim = 1000)

# Calculate profile
validation_profile <- profile(obj.validation, myfit$argument, "v", cores = 4, method = "integrate",
                              stepControl = list(stop = "data"),
                              algoControl = list(gamma = 1, reoptimize = T),
                              optControl = list(rinit = .1, rmax = 10, iterlim = 100, fterm = 1e-5, mterm = 1e-5))


# plotProfile(validation_profile) # This also plos the prediction colums, which is a bug in the code.
plotProfile(validation_profile, mode %in% c("validation", "data")) # Plots only the two contributions validation and data, along with the sum (total)
# Is the contribution of validation small?? # If yes: The total aligns with a prediction profile

# Confidence Interval of the prediction
confint(validation_profile, val.column = "value")



## Prediction band (prediction uncertainty for several time points) --------------------------------------------------------------
# Here we calculate a prediction CI for different timepoints. In the end we interpolate to a "prediction band"
library(parallel)
predprofs <- list()
prediction_band <- do.call(rbind, mclapply(seq(10, 50, 10), function(t) {

  cat("Computing prediction profile for t =", t, "\n")

  obj.validation <- normL2(data, g * x * p, times = c(t), attr.name = "data") +
    datapointL2(name = "TCA_cell", time = t, value = "v", sigma = 1, attr.name = "validation", condition = "closed")

  refit <- trust(obj.validation, parinit = c(v = 190, bestfit), rinit = 1, rmax = 10, iterlim = 1000)

  profile_prediction <- profile(obj.validation, refit$argument, "v", cores = 1, method = "integrate",
                                stepControl = list(stop = "data"),
                                algoControl = list(gamma = 1, reoptimize = T),
                                optControl = list(rinit = .1, rmax = 10, iterlim = 100, fterm = 1e-5, mterm = 1e-5))
  
  proflist <- c(predprofs, list(profile = profile_prediction, time = t))

  d1 <- confint(profile_prediction, val.column = "value")

  # Output
  data.frame(time = t, condition = "closed", name = "TCA_cell",  d1[-1])

}, mc.cores = 10))

times <- seq(0,51,len=300)
prediction <- (g * x * p)(times, bestfit) %>%
  as.data.frame()


prediction_band_spline <- data.frame(
  time = prediction$time[prediction$time>=10],
  value = prediction$value[prediction$time>=10],
  condition = "closed",
  name = "TCA_cell",
  lower = spline(prediction_band$time, prediction_band$lower, xout = prediction$time[prediction$time>=10])$y,
  upper = spline(prediction_band$time, prediction_band$upper, xout = prediction$time[prediction$time>=10])$y
)

# Create the ggplot
ggplot(prediction, aes(x = time, y = value, color = condition)) +
  geom_line() +  # Line connecting the points for each condition
  geom_ribbon(data = prediction_band_spline, aes(x = time, ymin = lower, ymax = upper, fill = condition),
              lty = 0, alpha = .3, show.legend = F) +  # Show ribbon in the legend
  geom_point(data = prediction_band, aes(x = time, y = lower, color = condition), shape = 4, show.legend = F) +
  geom_point(data = prediction_band, aes(x = time, y = upper, color = condition), shape = 4, show.legend = F) +
  facet_wrap(~ name, scales = "free_y") +  # Facet by 'name' column
  labs(
    x = "Time",
    y = "Value",
    color = "Condition"
  ) +
  dMod::theme_dMod() +  # Apply dMod theme
  dMod::scale_color_dMod() +  # Apply dMod color scale to lines
  dMod::scale_fill_dMod()   # Apply the same color scale to the fill


## Alternative implementation of the Steady state via implicit parameter transformation of the steady state  -----------------------------------------------

# Redefine reactions in order to control the standard and open condition by events
reactions <- eqnlist() %>%
  addReaction("TCA_buffer", "TCA_cell", rate = "import*TCA_buffer", description = "Uptake")%>%
  addReaction("TCA_cell", "TCA_buffer", rate = "export_sinus*TCA_cell", description = "Sinusoidal export")%>%
  addReaction("TCA_cell", "TCA_cana", rate = "export_cana*TCA_cell", description = "Canalicular export")%>%
  addReaction("TCA_cana", "TCA_buffer", rate = "(reflux*(1-switch) + reflux_open*switch)*TCA_cana", description = "Reflux into the buffer")%>%
  addReaction("0", "switch", rate = "0", description = "Create a switch")

events <- NULL
events <- addEvent(events, var = "TCA_buffer", time = 0, value = 0)
events <- addEvent(events, var = "switch" , time = 0, value = "OnOff")
mymodel <- odemodel(reactions, modelname = "bamodel2", events = events)
x <- Xs(mymodel)



# Replace one reaction with a analytical expression for the conserved quantity: TCR_tot
f <- as.eqnvec(reactions)[c("TCA_buffer", "TCA_cana", "TCA_cell")]
f["TCA_cell"] <- "TCA_buffer + TCA_cana + TCA_cell - TCA_tot"
pSS <- P(f, method = "implicit", compile = TRUE, modelname = "pfn")

observables <- eqnvec(buffer = "s*TCA_buffer", cellular = "s*(TCA_cana + TCA_cell)")

innerpars <- unique(c(getParameters(mymodel), getSymbols(observables), getSymbols(f)))
trafo <- repar("x~x" , x = innerpars)
trafo <- repar("x~0" , x = reactions$states, trafo)

trafo <- repar("x~exp(log(10)*x)", x = setdiff(innerpars, "OnOff"), trafo)
p <- P(repar("OnOff~0", trafo), condition = "closed") + P(repar("OnOff~1", trafo), condition = "open")
g <- Y(observables, f = x, compile = TRUE, modelname = "obsfn2")


