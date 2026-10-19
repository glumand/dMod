## ---------------------------------------------------------------------------
## Reads every saved pi3kmini_results*.rds (one per sigma/data setting run
## via example_toy_pi3kmini_fim_fit.R) and reports the Gauss-Newton/FIM
## eigenvalue spectrum at the INITIAL reference point (mu_center) that was
## already computed and saved there (`res$eig`) -- no model rebuild, no
## dMod loading, nothing recomputed. Intended to let you track how the
## spectral gap structure (the true exact-symmetry gap vs. any secondary/
## emergent gaps) depends on sigma/data density.
##
## This script PRINTS a CSV to stdout -- it does not write the CSV to a file
## itself. Redirect it yourself, e.g.:
##
##   Rscript scripts/example_toy_pi3kmini_fim_eigenvalue_report.R \
##     > scripts/plots/pi3kmini_fim_fit/pi3kmini_eigenvalue_report.csv
##
## It ALSO saves a plot of log10(eigenvalue) vs. sigma, one line per rank --
## the gap between adjacent rank-lines at a given sigma is exactly the
## spectral gap; watching how those gaps open/close across sigma is the
## point -- to scripts/plots/pi3kmini_fim_fit/pi3kmini_eigenvalues_vs_sigma.pdf
## (this file write happens regardless of stdout redirection).
## ---------------------------------------------------------------------------

out_dir <- "/home/mio/Phd/dMod/scripts/plots/pi3kmini_fim_fit"
rds_files <- list.files(out_dir, pattern = "^pi3kmini_results.*\\.rds$", full.names = TRUE)

cat("file,sigma,n_timepoints,rank,eigenvalue_initial,log10_eigenvalue_initial\n", sep = "")

report <- data.frame(file = character(0), sigma = numeric(0), n_timepoints = integer(0),
                      rank = integer(0), eigenvalue_initial = numeric(0),
                      log10_eigenvalue_initial = numeric(0))

for (rds in rds_files) {
  res <- readRDS(rds)
  eig_initial <- res$eig$values

  for (rank in seq_along(eig_initial)) {
    log10_val <- log10(abs(eig_initial[rank]) + 1e-300)
    cat(basename(rds), ",",
        res$sigma, ",",
        length(res$times), ",",
        rank, ",",
        format(eig_initial[rank], scientific = TRUE, digits = 6), ",",
        format(log10_val, digits = 6),
        "\n", sep = "")
    report <- rbind(report, data.frame(file = basename(rds), sigma = res$sigma,
                                        n_timepoints = length(res$times), rank = rank,
                                        eigenvalue_initial = eig_initial[rank],
                                        log10_eigenvalue_initial = log10_val))
  }
}

suppressMessages(library(ggplot2))
report$rank <- factor(report$rank, levels = sort(unique(report$rank)))

p_eig <- ggplot(report, aes(sigma, log10_eigenvalue_initial, color = rank, group = rank)) +
  geom_line() +
  geom_point(size = 1.3) +
  labs(x = expression(sigma), y = expression(log[10]*"(eigenvalue)"),
       color = "rank",
       title = "PI3K-mini: initial FIM eigenvalue spectrum vs. sigma",
       subtitle = "one line per eigenvalue rank; vertical distance between lines = spectral gap") +
  theme_bw() + theme(legend.position = "right")
ggsave(file.path(out_dir, "pi3kmini_eigenvalues_vs_sigma.pdf"), p_eig, width = 7, height = 5)

