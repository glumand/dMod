rm(list = ls(all.names = TRUE))
setwd(tempdir())
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
mymodel <- odemodel(reactions, modelname = "bamodel", compile = F)
x <- Xs(mymodel)

# Define observables buffer and cellular
observables <- eqnvec(buffer = "s*TCA_buffer", cellular = "s*(TCA_cana + TCA_cell)")
g <- Y(observables, f = x, condition = NULL, compile = F, modelname = "obsfn_bamodel", attach.input = T)

# Define parameter transformations using define(), insert() and branch(). Old function repar also avaiable!
innerpars <- getParameters(x,g)
trafo <- NULL %>%
  define("x~x", x = innerpars) %>% # identity
  define("TCA_buffer~0") %>%
  insert("x~exp10(y)", x = .currentSymbols, y = toupper(.currentSymbols))


# # Explicit trafo (this is equivalent to the lines above)
# trafo <- eqnvec(TCA_buffer = "0",
#                 TCA_cell = "10^TCA_cell",
#                 TCA_cana = "10^TCA_cana",
#                 k_import = "10^k_import",
#                 k_export_sinus = "10^k_export_sinus",
#                 k_export_cana = "10^k_export_cana",
#                 k_reflux = "10^k_reflux",
#                 s = "10^s")

p <- P(trafo, condition = "closed", compile = F)


# Compile the objects
compile(g, x, p, output = "bamodel", cores = 4) # Compile C/C++ output of odemodel in parallel

## Use simulate data to calibrate outer model parameters ---
outerpars <- getParameters(p)
pouter <- structure(runif(length(outerpars), min = -1, max = 0), names = outerpars)

prd <- g*x*p
# debugonce(x)
times <- seq(0, 45, len = 300)
# debugonce(g)
out <- prd(times, pouter)
plot(out, data)
plot(getDerivs(out))

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

# # Fit 50 times, sample with sd=4 around pouter
outms <- mstrust(obj, pouter, sd = 4, studyname = "bamodel", cores=detectFreeCores(), fits=50, iterlim = 1e3)

# ## Later: Fitting on Knecht machines
# outknecht <- runbg({
#   mstrust(obj, pouter, sd = 4, studyname = "bamodelms", cores=detectFreeCores(), fits=100, iterlim = 1e3)
# }, machine = "knecht1", filename = "bamodelms", link = T)
# outknecht$check()
# outknecht$get()
# 
# outms <- .runbgOutput$knecht1

out_frame <- as.parframe(outms)
plotValues(out_frame) # Show "Waterfall" plot
plotPars(out_frame) # Show parameter plot
bestfit <- as.parvec(out_frame)


# Plot predictions along data
plot((g*x*p)(times, bestfit), data)
# 
# Calculate Parameter Profiles and plot different contributions (for identifiablility only "data" is of interest)
profiles_integrate <- profile(obj, bestfit, whichPar = names(bestfit), method = "integrate", cores = detectFreeCores(), 
                              limits = c(lower = -5, upper = 5), stepControl = list(stop = "data"))

profiles_optimize <- profile(obj, bestfit, whichPar = names(bestfit), method = "optimize", cores = detectFreeCores(), 
                             limits = c(lower = -5, upper = 5), 
                             stepControl = list(stepsize = 1e-4, min = 1e-4, max = Inf, 
                                                atol = 1e-2, rtol = 1e-2, limit = 200, stop = "data"))

proflist <- list(integrate = profiles_integrate, optimize = profiles_optimize) 
# Integration based profiles fast but not exakt
# Best practice: use method = "integrate" with reoptimize = TRUE in algoControl, then the integrate step is already close to the new optimum

plotProfile(proflist, mode %in% c("data", "prior"))

plotProfile(profiles_integrate, mode %in% c("data", "prior"))
plotProfile(profiles_optimize, mode %in% c("data", "prior"))

# The triple S, TCE_CELL, TCA_CANA compensate by structure of the model
plotPaths(profiles_optimize, whichPar = "TCA_CANA")
plotPaths(profiles_optimize, whichPar = "K_EXPORT_CANA")
plotPaths(profiles_optimize, whichPar = "K_REFLUX_OPEN")
plotPaths(profiles_optimize, whichPar = "S") 

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
p <- P(trafo, modelname = "bamodel_SS", compile = TRUE)

outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
p(pouter)
plot((g*x*p)(times, pouter),data)

# Objective function
obj <- normL2(data, g * x * p, attr.name = "data") + constraintL2(pouter, sigma = 20, attr.name = "prior")

# Multistart fit
outms <- mstrust(obj, pouter, sd = 4, iterlim = 1e3, studyname = "bamodel_ss", cores = detectFreeCores(), fits = 50)
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
validation_profile <- profile(obj.validation, myfit$argument, "v", cores = 4, method = "optimize",
                              stepControl = list(stepsize = 1e-4, min = 1e-5, max = 1e3, atol = 1e-2, rtol = 1e-2, limit = 500, stop = "data"),
                              algoControl = list(reoptimize = T),
                              optControl = list(rinit = .1, rmax = 5, iterlim = 200, fterm = 1e-5, mterm = 1e-5),
                              cautiousMode = TRUE)

# plotProfile(validation_profile) # This also plos the prediction colums, which is a bug in the code.
plotProfile(validation_profile, mode %in% c("validation", "data")) # Plots only the two contributions validation and data, along with the sum (total)
# Is the contribution of validation small?? # If yes: The total aligns with a prediction profile

# Confidence Interval of the prediction
confint(validation_profile, val.column = "value")



## Prediction band (prediction uncertainty for several time points) --------------------------------------------------------------
# Here we calculate a prediction CI for different timepoints. In the end we interpolate to a "prediction band"
library(parallel)
predprofs <- list()
prediction_band <- do.call(rbind, mclapply(c(0,1,2,3,4,seq(5, 50, 5)), function(t) {
  
  cat("Computing prediction profile for t =", t, "\n")
  
  obj.validation <- normL2(data, g * x * p, times = c(t), attr.name = "data") +
    datapointL2(name = "TCA_cell", time = t, value = "v", sigma = 0.001, attr.name = "validation", condition = "closed")
  
  refit <- trust(obj.validation, parinit = c(v = 190, bestfit), rinit = 1, rmax = 10, iterlim = 1000)
  
  profile_prediction <- profile(obj.validation, refit$argument, "v", cores = 1, method = "integrate",
                                stepControl = list(stop = "data"),
                                algoControl = list(gamma = 1, reoptimize = T),
                                optControl = list(rinit = .1, rmax = 10, iterlim = 100, fterm = 1e-5, mterm = 1e-5))
  
  proflist <- c(predprofs, list(profile = profile_prediction, time = t))
  
  d1 <- confint(profile_prediction, val.column = "data")
  
  # Output
  data.frame(time = t, condition = "closed", name = "TCA_cell",  d1[-1])
  
}, mc.cores = detectFreeCores()-1))

times <- seq(0,50,len=300)
prediction <- (g * x * p)(times, bestfit) %>%
  as.data.frame()


prediction_band_spline <- data.frame(
  time = prediction$time,
  value = prediction$value,
  condition = "closed",
  name = "TCA_cell",
  lower = spline(prediction_band$time, prediction_band$lower, xout = prediction$time)$y,
  upper = spline(prediction_band$time, prediction_band$upper, xout = prediction$time)$y
)

# Create the ggplot using dMod theme and color
ggplot(prediction, aes(x = time, y = value, color = condition)) +
  geom_line() +  # Line connecting the points for each condition
  geom_point(data = badata, aes(x = time, y = value, color = condition)) + 
  geom_errorbar(data = badata, aes(x = time, y = value,
                                   ymin = value - sigma, ymax = value + sigma,
                                   colour = condition)) +
  geom_ribbon(data = prediction_band_spline, aes(x = time, ymin = lower, ymax = upper, fill = condition), lty = 0, alpha = .3, show.legend = F) +  # Show ribbon in the legend
  geom_point(data = prediction_band, aes(x = time, y = lower, color = condition), shape = 4, show.legend = F) +
  geom_point(data = prediction_band, aes(x = time, y = upper, color = condition), shape = 4, show.legend = F) +
  facet_wrap(~ name, scales = "free_y") +  # Facet by 'name' column
  labs(
    x = "Time",
    y = "Value",
    color = "Condition"
  ) +
  dMod::theme_dMod() +        # Apply dMod theme
  dMod::scale_color_dMod() +  # Apply dMod color scale to lines
  dMod::scale_fill_dMod()     # Apply the same color scale to the fill



