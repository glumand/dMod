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
f <- eqnvec() %>%
  addReaction("A", "B", "k1*A", "Production of B") %>%
  addReaction("B", "C", "k2*B", "Production of C")

# Define observables and error model
observables <- eqnvec(B_obs = "B + off_B")
errors <- eqnvec(B_obs = "sigma_rel*B_obs")

# Generate dMod objects
model <- odemodel(f, modelname = "errtest", compile = F, solver = "boost")
x     <- Xs(model)
g     <- Y(observables, x, 
           compile = F, modelname = "obsfn")
e     <- Y(errors, g, attach.input = FALSE,
           compile = F, modelname = "errfn")

# Generate parameter transformation
innerpars <- getParameters(model, g, e)
covariates <- data.frame(Aini = c("C1", "C2"), row.names = c("C.1", "C.2"))

trafo <- eqnvec() %>% 
  define("x~x", x = innerpars) %>%
  define("x~0", x = c("B", "C")) %>%
  branch(table = covariates) %>%
  insert("A~Aini", Aini = Aini) %>%
  insert("x~exp(x)", x = .currentSymbols) %>%
  {.}

p <- P(trafo, modelname = "parfn", compile = F)

compile(g, x, e, p, output = "errtest_total", cores = 8)
#compile(g, x, e, p, cores = 4)


## Simulate data
ptrue <- c(C1 = 1, C2 = 2, k1 = -2, k2 = -3, off_B = -3, sigma_rel = log(.1))
times <- seq(0, 50, 1)
prediction <- (g*x*p)(times, ptrue, deriv = TRUE)
datasheet <- subset(as.data.frame(prediction, errfn = e), name == "B_obs")
datasheet$value <- datasheet$value + rnorm(length(datasheet$value), sd = datasheet$sigma)
data <- as.datalist(datasheet)

plotData(data)

# Remove sigma
datasheet$sigma <- NA
data <- as.datalist(datasheet)
plotData(data)

## Fit data with error model
obj <- normL2(data, g*x*p, e)
myfit <- trust(obj, ptrue, rinit = 1, rmax = 10, printIter = TRUE)
fits <- mstrust(obj, center = ptrue, sd = 3, fits = 100, cores = 20, studyname = "msrun1")

outframe <- fits %>% as.parframe()
plotValues(outframe)
bestfit <- as.parvec(outframe)

obj(bestfit)

bestprediction <- (g*x*p)(times, bestfit, deriv = F)
pred <- subset(as.data.frame(bestprediction, errfn = e), name == "B_obs")

ggplot(data = datasheet, aes(time, value, color = condition)) + 
  geom_point() + 
  geom_line(data = pred, aes(time, value, color = condition)) + 
  geom_ribbon(data = pred, aes(ymin = value - sigma, ymax = value + sigma, fill = condition), alpha = 0.2, linewidth = 0) + 
  dMod::theme_dMod() + 
  dMod::scale_color_dMod() + 
  dMod::scale_fill_dMod()

profiles <- profile(obj, 
                    bestfit, names(bestfit), 
                    limits = c(-5, 5), 
                    cores = length(bestfit),
                    method = "optimize",
                    stepControl = list(stop = "neg2ll"))

plotProfile(profiles, mode %in% c("data", "neg2ll"))
