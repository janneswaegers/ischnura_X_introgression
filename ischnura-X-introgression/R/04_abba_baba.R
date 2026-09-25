# ABBA-BABA analysis: introgression at X-linked versus autosomal windows
#
# Input: Dsuite Dinvestigate output (scripts/08_dsuite.sh), windows of 50 SNPs,
# step 25, with an added column Xaut ("X" / "aut").
#
# Six analyses (females only): sympatric I. elegans and sympatric I. graellsii,
# each with three autosomal admixture (Q) cut-offs used to decide which sympatric
# individuals count as genomically I. elegans or I. graellsii.
#
# For each analysis, a permutation test asks whether X-linked windows show lower
# f_dM (less introgression) than autosomal windows (paper: figure 3), with the same
# test for D and f_d (table S5). A Wilcoxon test compares overall introgression into
# sympatric I. elegans and into sympatric I. graellsii.
#
# Output: results/R/abba_baba/

library(ggplot2)
source("R/functions.R")

# --- Settings -----------------------------------------------------------------
dsuite_dir <- "data/dsuite"
sample_set <- "females"        # "malesfemales" repeats the analysis with both sexes (table S4)
out_dir    <- "results/R/abba_baba"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analyses <- data.frame(
  species  = rep(c("I. elegans", "I. graellsii"), each = 3),
  Q_cutoff = c("Q = 0", "Q < 0.1", "Q < 0.25", "Q = 1", "Q > 0.9", "Q > 0.75"),
  folder   = paste0("ALL_filtered_", sample_set,
                    rep(c("", "_graellsii"), each = 3), c(1, 10, 25))
)

read_dsuite <- function(folder) {
  file <- list.files(file.path(dsuite_dir, folder), pattern = "localFstats.*Xaut\\.txt$",
                     full.names = TRUE)[1]
  read.table(file, header = TRUE)
}

# --- 1. Permutation tests: X-linked versus autosomal windows --------------------
set.seed(1)
results <- list()
windows <- list()

for (i in seq_len(nrow(analyses))) {
  d <- read_dsuite(analyses$folder[i])
  d$analysis <- paste(analyses$species[i], analyses$Q_cutoff[i], sep = ", ")
  windows[[i]] <- d

  for (stat in c("f_dM", "D", "f_d")) {
    res <- perm_median_test(d[d$Xaut == "X", stat], d[d$Xaut == "aut", stat],
                            alternative = "less")
    results[[length(results) + 1]] <- data.frame(
      analyses[i, c("species", "Q_cutoff")], statistic = stat,
      median_X = res$median_X, median_aut = res$median_aut, p_permutation = res$p)
  }
}

results <- do.call(rbind, results)
windows <- do.call(rbind, windows)
print(results[results$statistic == "f_dM", ])
write.table(results, file.path(out_dir, "permutation_tests.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

# --- 2. Figure 3: f_dM distributions --------------------------------------------
windows$analysis <- factor(windows$analysis, levels = unique(windows$analysis))
p <- ggplot(windows, aes(x = f_dM, fill = Xaut)) +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_histogram(colour = "#e9ecef", alpha = 0.7, position = "identity", binwidth = 0.008) +
  scale_fill_manual(values = c(aut = "#6E825A", X = "#1F3500"),
                    labels = c(aut = "autosomal", X = "X-linked"), name = "") +
  facet_wrap(~ analysis, ncol = 2, dir = "v") +
  xlim(-0.3, 0.3) +
  theme_bw()
ggsave(file.path(out_dir, "figure3_fdM.pdf"), p, width = 8, height = 9)

# --- 3. Introgression into sympatric I. elegans versus sympatric I. graellsii ----
for (q in 1:3) {
  el <- windows$f_dM[windows$analysis == levels(windows$analysis)[q]]
  gr <- windows$f_dM[windows$analysis == levels(windows$analysis)[q + 3]]
  w  <- wilcox.test(el, gr, alternative = "two.sided")
  cat(sprintf("%s vs %s: Wilcoxon p = %.2g\n",
              levels(windows$analysis)[q], levels(windows$analysis)[q + 3], w$p.value))
}

# --- 4. Candidate 'reproductive barrier' windows on the X --------------------------
# X-linked windows in the lowest 2.5% of the f_dM distribution, per analysis;
# their scaffolds were compared with the scaffolds of the top BGC beta outliers.
candidates <- do.call(rbind, lapply(split(windows[windows$Xaut == "X", ], ~ analysis), function(d) {
  d[d$f_dM < quantile(d$f_dM, 0.025, na.rm = TRUE), c("analysis", "chr", "windowStart", "windowEnd", "f_dM")]
}))
write.table(candidates, file.path(out_dir, "X_windows_lowest_fdM.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")
