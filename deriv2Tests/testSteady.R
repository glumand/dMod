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
r <- eqnlist() %>%
  # addReaction("P0", "P1", "2*ka*P0") %>% 
  # addReaction("P1", "P0", "ki*P1") %>% 
  # addReaction("P1", "P2", "ka*P1") %>% 
  # addReaction("P2", "P1", "2*ki*P2") %>% 
  
  addReaction("TCR", "y10", "kon*P1*TCR") %>% 
  addReaction("y10", "TCR", "koff*y10") %>% 
  addReaction("y10", "TCR", "ki*y10") %>% 
  
  addReaction("TCR", "y11", "2*kon*P2*TCR") %>% 
  addReaction("y11", "TCR", "koff*y11") %>% 
  addReaction("y11", "TCR", "ki*y11") %>% 
  
  addReaction("y10", "y11", "ka*y10") %>% 
  addReaction("y11", "y10", "ki*y11") %>%
  
  addReaction("y11 + TCR", "y22", "qon*y11*TCR") %>% 
  addReaction("y22", "y11 + TCR", "2*qoff*y22") %>% 
  addReaction("y22", "y10 + TCR", "2*ki*y22")

r$rates <- replaceSymbols("ka","4*ki", r$rates)
r$rates <- replaceSymbols("ki","ci*Intsty", r$rates)
# r$rates <- replaceSymbols("qoff","koff", r$rates)

r$rates <- replaceSymbols("P1","2*K_L/(1+K_L)^2 * Ptot", r$rates)
r$rates <- replaceSymbols("P2","(K_L/(1+K_L))^2 * Ptot", r$rates)

mysteadies <- steadyStates(r, neglect = "Intsty")


innerpars <- getParameters(r)
trafo <- eqnvec() %>% 
  define("x~x", x=innerpars) %>% 
  insert("x~y", x=names(mysteadies$equations), y=mysteadies$equations) %>% 
  # insert("x~y", x=c("TCR","K_L", "Intsty"), y=c(1,4,1))
  {.}
  
p <- P(trafo, compile = T)


x <- odemodel(r, modelname = "testAlyssa_KPR", compile = T, solver = "boost") %>% Xs()

# compile(x,p, output = "AlyssaTest", cores = 4)


times <- seq(0,400,len = 300)
outerpars <- getParameters(p)
pouter <- structure(runif(length(outerpars)), names = outerpars)
(x*p)(times, pouter) %>% plot()

p(pouter) %>% getDerivs()
