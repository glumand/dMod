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
  addReaction("", "Stim", "0") %>%
  addReaction("", "R", "k_R_act*Stim+k_R_act_basal") %>%
  addReaction("R", "", "k_R_deact * R") %>%
  addReaction("K1","pK1","Vmax_K1_act*R*K1/(Km_K1_act+K1)") %>%
  addReaction("pK1","K1","Vmax_K1_deact*pK1 * Phos/(Km_Phos+pK1)") %>%
  addReaction("K2","pK2","Vmax_K2_act*pK1*K2/(Km_K2_act+K2)") %>%
  addReaction("pK2","K2","Vmax_K2_deact*pK2/(Km_Phos+pK2)") %>%
  addReaction("K3","pK3","Vmax_K3_act*pK2*K3/(Km_K3_act+K3)") %>%
  addReaction("pK3","K3","Vmax_K3_deact*pK3/(Km_Phos+pK3)") %>%
  addReaction("TF","pTF","Vmax_TF_act*pK3*TF/(Km_TF_act+TF)") %>%
  addReaction("pTF","TF","Vmax_TF_deact*pTF/(Km_TF_deact+pTF)") %>%
  addReaction("","mRNA_Phos","Vmax_tx*pTF^n_hill/(K_hill^n_hill+pTF^n_hill)") %>%
  addReaction("mRNA_Phos","","k_mRNA_deg*mRNA_Phos") %>%
  addReaction("","Phos","k_translation*mRNA_Phos") %>%
  addReaction("Phos","","k_phos_deg*Phos")


mysteadies <- steadyStates(r, forcings = "Stim")


p.eql <- P(r, forcings = "Stim", compile = F, keep.root = F)

outerpars <- getParameters(p.eql)
pouter <- structure(runif(length(outerpars)), names = outerpars)
pouter
p.eql(pouter)

cond.grid <- data.frame(Stim = c(1, 0))
rownames(cond.grid) <- c("Stim", "CTRL")
innpers_equil <- getParameters(p.eql)
p.log <- eqnvec() %>% 
  define("x~x", x = innpers_equil) %>%
  branch(cond.grid, apply = "insert") %>%
  insert("x~exp10(x)", x = .currentSymbols[!grepl("n_hill", .currentSymbols)]) %>% 
  P(compile = F)


x <- odemodel(r, compile = F) %>% Xs()

# debugonce(compile)
compile(p.eql, p.log , x, cores = 8)

times <- seq(0,100,len = 300)
p <- p.eql*p.log
outerpars <- getParameters(p)
set.seed(5555)
pouter <- structure(runif(length(outerpars)), names = outerpars)
pouter[!grepl("n_hill", names(pouter))] <- -1*pouter[!grepl("n_hill", names(pouter))]

system.time({p(pouter)})
system.time({p(pouter)})

set.seed(5555)
pouter <- structure(runif(length(outerpars)), names = outerpars)
pouter[!grepl("n_hill", names(pouter))] <- -1*pouter[!grepl("n_hill", names(pouter))]
p(pouter)
(x*p)(times, pouter) %>% plot()
