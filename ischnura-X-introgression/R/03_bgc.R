# Bayesian genomic cline (BGC) analysis: X-linked versus autosomal SNPs
#
# Input: bgc output for the diagnostic SNP set (females only; 5 chains of 50,000
# steps, burn-in 25,000, thinning 20), combined with ClineHelpR.
#   alpha: direction of introgression (> 0: alleles move from I. graellsii into I. elegans)
#   beta:  steepness of the cline (higher = stronger barrier to gene flow)
#
# Tests whether X-linked SNPs have higher beta and alpha than autosomal SNPs with
# a permutation test on the medians (paper: figure 2, table S3).
#
# Output: results/R/bgc/

library(ClineHelpR)   # devtools::install_github("btmartin721/ClineHelpR")
library(ggplot2)
source("R/functions.R")

# --- Settings -----------------------------------------------------------------
bgc_dir   <- "data/bgc/long_runs_females"          # bgc output files with prefix "eatt"
popmap    <- "data/bgc/eatt.bgc.popmap_final.txt"
loci_file <- file.path(bgc_dir, "eatt_bgc_loci.txt")
out_dir   <- "results/R/bgc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Combine chains, check convergence and get per-SNP alpha and beta --------
bgc <- combine_bgc_output(results.dir = bgc_dir, prefix = "eatt")
plot_traces(df.list = bgc, prefix = "eatt", plotDIR = out_dir)

outliers <- get_bgc_outliers(df.list = bgc, admix.pop = "EATT", popmap = popmap,
                             loci.file = loci_file, qn = 0.975)
write.table(outliers[[1]], file.path(out_dir, "bgc_snp_parameters.txt"), quote = FALSE)
write.table(outliers[[3]], file.path(out_dir, "bgc_hybrid_index.txt"),    quote = FALSE)

# --- 2. X-linked versus autosomal SNPs -------------------------------------------
# Per-SNP alpha and beta, split into X-linked and autosomal SNPs
ab_X   <- read.table("data/bgc/alpha_beta_females_longruns_X.txt",   header = TRUE)
ab_aut <- read.table("data/bgc/alpha_beta_females_longruns_aut.txt", header = TRUE)

set.seed(1)
for (param in c("beta", "alpha")) {
  res <- perm_median_test(ab_X[[param]], ab_aut[[param]], alternative = "greater")
  cat(sprintf("%-5s  median X = %.3f  median autosomes = %.3f  permutation p = %.4f\n",
              param, res$median_X, res$median_aut, res$p))
  cat(sprintf("       positive values: X %.0f%%, autosomes %.0f%%\n",
              100 * mean(ab_X[[param]] > 0), 100 * mean(ab_aut[[param]] > 0)))
}
# Paper: beta  median X = 0.166, autosomes = 0.022, p < 0.001
#        alpha median X = 0.015, autosomes = -0.008, p = 0.034

# --- 3. Figure 2: distributions of beta and alpha --------------------------------
ab <- rbind(data.frame(ab_X[, c("alpha", "beta")],   Xaut = "X-linked"),
            data.frame(ab_aut[, c("alpha", "beta")], Xaut = "autosomal"))

for (param in c("beta", "alpha")) {
  p <- ggplot(ab, aes(x = .data[[param]], fill = Xaut)) +
    geom_histogram(colour = "#e9ecef", alpha = 0.7, position = "identity", binwidth = 0.03) +
    scale_fill_manual(values = c("autosomal" = "#6E825A", "X-linked" = "#1F3500")) +
    labs(x = param, y = "count", fill = "") +
    theme_bw()
  ggsave(file.path(out_dir, paste0("figure2_", param, ".pdf")), p, width = 6, height = 4)
}

# --- 4. Same test including males and females (table S3) -------------------------
ab_all <- read.table("data/bgc/beta_alpha_long_run.txt", header = TRUE)   # column Xaut: "X" / "aut"
for (param in c("beta", "alpha")) {
  res <- perm_median_test(ab_all[ab_all$Xaut == "X",   param],
                          ab_all[ab_all$Xaut == "aut", param], alternative = "greater")
  cat(sprintf("males + females, %-5s: median X = %.3f, autosomes = %.3f, p = %.4f\n",
              param, res$median_X, res$median_aut, res$p))
}
