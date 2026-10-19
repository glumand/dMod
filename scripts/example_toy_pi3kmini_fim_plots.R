## ---------------------------------------------------------------------------
## Plots for example_toy_pi3kmini_fim_fit.R. Reads
## scripts/plots/pi3kmini_fim_fit/pi3kmini_results.rds (parframes, best fits,
## and precomputed fine-grid prediction curves) and does NOT touch dMod's
## compiled model at all -- run the fitting script first, then this one, as
## often as you like to restyle plots without re-fitting anything.
##
## Produces figures in scripts/plots/pi3kmini_fim_fit/:
##   pi3kmini_waterfalls.png            -- multistart waterfall, all 3 methods
##   pi3kmini_group_scatter.png         -- raw group-parameter (U_cond1)
##                                          scatter across near-optimal fits
##                                          (boundedness)
##   pi3kmini_fit_vs_data.png           -- fitted curves (each method) over
##                                          the simulated data, faceted by
##                                          observable x condition
##   pi3kmini_parscatter_<method>.png   -- dMod's own plotPars() boxplot-per-
##                                          parameter scatter (same tool as
##                                          plotValues(), just for parameters
##                                          instead of objective values), one
##                                          per method, run directly on the
##                                          saved parframe
## ---------------------------------------------------------------------------

KNOWN_REPO_ROOT <- "/home/mio/Phd/dMod"
repo_root <- local({
  full <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", full, value = TRUE)
  ofile <- if (!is.null(sys.frames()) && length(sys.frames()) >= 1L) sys.frame(1)$ofile else NULL
  candidates <- c(if (length(hit)) dirname(sub("^--file=", "", hit[1])),
                   if (!is.null(ofile)) dirname(ofile), getwd(), KNOWN_REPO_ROOT)
  for (cand in candidates) {
    if (is.null(cand) || is.na(cand)) next
    here <- normalizePath(cand, mustWork = FALSE, winslash = "/")
    for (i in seq_len(8)) {
      if (file.exists(file.path(here, "DESCRIPTION"))) return(here)
      parent <- dirname(here); if (parent == here) break; here <- parent
    }
  }
  KNOWN_REPO_ROOT
})
repo_lib <- file.path(repo_root, ".Rlib")

## plotPars()/theme_dMod() need the dMod namespace (S3 dispatch on the saved
## parframe objects) -- this is package-level R code only, no compiled model
## is touched (no odemodel()/compile()/mstrust() call here).
library(dMod, lib.loc = repo_lib)
suppressMessages(library(ggplot2))


sigma <- 0.05
out_dir <- file.path(repo_root, "scripts", "plots", "pi3kmini_fim_fit")
res <- readRDS(file.path(out_dir, paste0("pi3kmini_results_sigma", sigma, ".rds")))

method_levels <- c("no regularization", "FIM null-space projector", "L2 (uniform)")
method_colors <- setNames(c("#E15759", "#59A14F", "#4E79A7"), method_levels)

## --- 1. Waterfall: penalized-objective values, ranked, per method ----------
mk_wf <- function(pf, label) {
  d <- pf[order(pf$value), ]
  data.frame(index = seq_len(nrow(d)), value = d$value - min(d$value), method = label)
}
wf <- rbind(mk_wf(res$pf_unreg, "no regularization"),
            mk_wf(res$pf_l2, "L2 (uniform)"),
            mk_wf(res$pf_fim, "FIM null-space projector"))
wf$method <- factor(wf$method, levels = method_levels)

## y-axis in units of 1e-3: the shared exponent goes once into the axis
## title (expression(...~(10^-3))) and the tick labels stay plain numbers,
## rather than repeating "1e-03" etc. on every tick.
# wf_scale <- 1e-1
p_wf <- ggplot(wf, aes(index, value / wf_scale, color = method)) +
 # scale_y_continuous(limits = c(0, 5), breaks = seq(0, 1, 10)) +
  ylim(0, 0.01) +
  geom_point(size = 1) +
  scale_color_manual(values = method_colors) +
  labs(x = "start, ranked", 
       # y = expression(Delta*"objective value"~(10^-3)),
       y = expression(Delta*"objective value"),
       title = "PI3K-mini: multistart waterfalls (own objective)") +
  theme_bw() + 
  theme(legend.position = "top",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
     #   panel.border = element_blank(),
        panel.background = element_blank()) 
ggsave(file.path(out_dir, paste0("pi3kmini_waterfalls_sigma", sigma, ".pdf")), p_wf, width = 6, height = 4)

## --- 2. Group parameter (U_cond1) scatter across near-optimal fits --------
top <- function(pf) pf[pf$value < min(pf$value) + 0.5, ]
scat <- rbind(
  data.frame(U_cond1 = top(res$pf_unreg)$U_cond1, method = "no regularization"),
  data.frame(U_cond1 = top(res$pf_l2)$U_cond1, method = "L2 (uniform)"),
  data.frame(U_cond1 = top(res$pf_fim)$U_cond1, method = "FIM null-space projector")
)
scat$method <- factor(scat$method, levels = method_levels)
truth_line <- log(res$true_pars_inner["U_cond1"])
center_line <- res$mu_center["U_cond1"]

p_scat <- ggplot(scat, aes(method, U_cond1, color = method)) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  geom_hline(aes(yintercept = truth_line, linetype = "true value"), color = "black") +
  geom_hline(aes(yintercept = center_line, linetype = "regularization center"), color = "grey40") +
  scale_linetype_manual(name = NULL, values = c("true value" = "dashed", "regularization center" = "dotted")) +
  scale_color_manual(values = method_colors, guide = "none") +
  labs(y = expression(log(U[cond1])~"(raw group parameter)"), x = NULL,
       title = "Individual group-parameter estimates across near-optimal fits",
       subtitle = "L2 and FIM now share the same regularization center (mu = -1)") +
  theme_bw() + theme(legend.position = "top")
ggsave(file.path(out_dir, "pi3kmini_group_scatter.png"), p_scat, width = 6, height = 4.5, dpi = 130)

## --- 3. Fit shown on the data -----------------------------------------------
obs_names <- c("obsA", "obsB", "obsMpp")
mk_pred <- function(df, label) {
  d <- df[df$name %in% obs_names, c("time", "name", "condition", "value")]
  d$method <- label
  d
}
pred_all <- rbind(mk_pred(res$pred_unreg, "no regularization"),
                   mk_pred(res$pred_l2, "L2 (uniform)"),
                   mk_pred(res$pred_fim, "FIM null-space projector"))
pred_all$method <- factor(pred_all$method, levels = method_levels)

data_obs <- res$sim[res$sim$name %in% obs_names, c("time", "name", "condition", "value", "sigma")]

p_fit <- ggplot() +
  geom_errorbar(data = data_obs, aes(time, ymin = value - sigma, ymax = value + sigma),
                color = "grey30", width = 0, alpha = 0.6) +
  geom_point(data = data_obs, aes(time, value), color = "grey30", size = 0.7, alpha = 0.6) +
  geom_line(data = pred_all, aes(time, value, color = method), linewidth = 0.6) +
  scale_color_manual(values = method_colors) +
  facet_grid(name ~ condition, scales = "free_y") +
  labs(x = "time", y = "observable value",
       title = "PI3K-mini: fitted curves vs. simulated data",
       subtitle = "grey points +- sigma = simulated data; all three methods reproduce the data essentially identically") +
  theme_bw() + 
  theme(legend.position = "top",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
       # panel.border = element_blank(),
        panel.background = element_blank()) 
ggsave(file.path(out_dir, paste0("pi3kmini_fit_vs_data_sigma", sigma, ".pdf")), p_fit, width = 9, height = 6, dpi = 130)

## --- 4. Standard dMod parameter scatter (plotPars), one per method ---------
## Same idiom as fim-portugal/ab_fit.R's `parscatter <- plotPars(pf)`: a
## boxplot of every outer parameter's value, grouped/colored by which
## objective-value "step" (stepDetect(), same grouping plotValues() uses for
## the waterfall) each fit belongs to. Run directly on the saved parframes.
method_pf <- list(
  "no regularization"        = res$pf_unreg,
  "L2 (uniform)"              = res$pf_l2,
  "FIM null-space projector"  = res$pf_fim
)
method_file <- c("no regularization" = "unreg", "L2 (uniform)" = "l2",
                  "FIM null-space projector" = "fim")

for (m in names(method_pf)) {
  p_pars <- plotPars(method_pf[[m]]) +
    labs(title = paste0("PI3K-mini: parameter scatter (plotPars) -- ", m))
  fname <- file.path(out_dir, paste0("pi3kmini_parscatter_", method_file[m], "_sigma", sigma, ".pdf"))
  ggsave(fname, p_pars, width = 8, height = 4.5, dpi = 130)
  cat("Saved", fname, "\n")
}

## --- 5. All-parameter scatter clouds, all three methods overlaid ----------
## Same idea as the U_cond1-only scatter above, extended to every outer
## parameter and faceted, with all three methods' point clouds overlaid in
## each panel (rather than dMod's own plotPars(), which only ever shows one
## parframe/method at a time as boxplots).
all_scatter <- do.call(rbind, lapply(names(method_pf), function(m) {
  d <- top(method_pf[[m]])
  do.call(rbind, lapply(res$outerpars, function(par) {
    data.frame(name = par, value = d[[par]], method = m)
  }))
}))
all_scatter$method <- factor(all_scatter$method, levels = method_levels)
all_scatter$name <- factor(all_scatter$name, levels = res$outerpars)

truth_df  <- data.frame(name = factor(names(res$pouter_true), levels = res$outerpars),
                         truth = as.numeric(res$pouter_true))
center_df <- data.frame(name = factor(names(res$mu_center), levels = res$outerpars),
                         center = as.numeric(res$mu_center))

p_clouds <- ggplot(all_scatter, aes(method, value, color = method)) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  geom_hline(data = truth_df, aes(yintercept = truth), linetype = "dashed", color = "black") +
  geom_hline(data = center_df, aes(yintercept = center), linetype = "dotted", color = "grey40") +
  facet_wrap(~ name, scales = "free_y", ncol = 5) +
  scale_color_manual(values = method_colors) +
  labs(x = NULL, y = "log-parameter value",
       title = "PI3K-mini: all-parameter scatter clouds, all 3 methods overlaid",
       subtitle = "dashed = true value, dotted = regularization center; near-optimal fits only") +
  theme_bw() +
  theme(legend.position = "top",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
ggsave(file.path(out_dir, paste0("pi3kmini_parscatter_clouds_sigma", sigma, ".pdf")), p_clouds, width = 11, height = 6, dpi = 130)
cat("Saved", file.path(out_dir, paste0("pi3kmini_parscatter_clouds_sigma", sigma, ".pdf")), "\n")

cat("\nPlots saved to", out_dir, "\n")
cat(" -", file.path(out_dir, "pi3kmini_waterfalls.png"), "\n")
cat(" -", file.path(out_dir, "pi3kmini_group_scatter.png"), "\n")
cat(" -", file.path(out_dir, "pi3kmini_fit_vs_data.png"), "\n")
cat(" -", file.path(out_dir, "pi3kmini_parscatter_{unreg,l2,fim}.png"), "\n")
