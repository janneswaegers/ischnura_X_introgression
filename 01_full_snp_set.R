# Genotypic sexing and filtering of the full SNP set
#
# 1. Calculates per-individual homozygosity at X-linked SNPs. Males are hemizygous
#    for the X (X0), so their homozygosity is ~1; samples with homozygosity > 0.96
#    were classified as males, samples below as females.
# 2. Removes X-linked SNPs where any male is heterozygous.
# 3. Keeps one random SNP per RAD-tag and writes X-linked / autosomal SNP lists.
#
# Input:  Stacks VCF of the full SNP set (scripts/04_populations.slurm),
#         lists of X-linked and autosomal scaffolds, sample information.
# Output: results/R/full/

library(adegenet)
library(vcfR)
library(dplyr)
library(ggplot2)
source("R/functions.R")

# --- Settings -----------------------------------------------------------------
vcf_file         <- "data/populations_full/populations.snps.vcf"
sample_info_file <- "data/sample_info.txt"    # columns: ID, PopID, Sex_reclassified
x_scaffolds      <- read.table("data/X_chrom_ids")$V1
aut_scaffolds    <- read.table("data/Autosome_ids")$V1
out_dir          <- "results/R/full"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gl <- read_genlight(vcf_file)

# --- 1. Genotypic sexing ------------------------------------------------------
chrom <- as.character(gl@chromosome)
sexing <- data.frame(
  ID      = indNames(gl),
  hom_X   = prop_hom_individual(gl[, chrom %in% x_scaffolds]),
  hom_aut = prop_hom_individual(gl[, chrom %in% aut_scaffolds])
)
sexing$sex_genotypic <- ifelse(sexing$hom_X > 0.96, "male", "female")
write.table(sexing, file.path(out_dir, "homozygosity_per_individual.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

# Homozygosity at X-linked SNPs is bimodal; the cut-off sits in the valley
# between the female and male peaks (electronic supplementary material, figure S1)
p <- ggplot(sexing, aes(x = hom_X)) +
  geom_histogram(binwidth = 0.006) +
  geom_vline(xintercept = 0.96, colour = "red") +
  labs(x = "Proportion of homozygous X-linked SNPs", y = "Number of individuals") +
  theme_bw()
ggsave(file.path(out_dir, "homozygosity_X_histogram.pdf"), p, width = 6, height = 4)

# The final sex assignment used below (column Sex_reclassified) combines this
# genotypic classification with the phenotypic sexing (see paper, Methods 2b).

# --- 2. Remove excluded samples and filter X-linked SNPs ------------------------
excluded <- c("M5", "M6", "I8", "K7", "No4", "N7", "K10", "K5", "S7", "N6")
writeLines(excluded, file.path(out_dir, "excluded_samples.txt"))
gl <- drop_samples(gl, excluded)

sample_info <- read.table(sample_info_file, header = TRUE)
sample_info <- sample_info[match(indNames(gl), sample_info$ID), ]   # same order as the VCF
is_male <- tolower(sample_info$Sex_reclassified) == "male"
is_male[is.na(is_male)] <- FALSE

gl_filtered <- filter_X_hom_males(gl, is_male, x_scaffolds)
# Paper: 50,733 SNPs, of which 2,469 X-linked
writeLines(locNames(gl_filtered), file.path(out_dir, "loc_filtered.txt"))

# --- 3. One random SNP per RAD-tag (unlinked SNP set) ---------------------------
# Note: the published locus lists were drawn without a fixed seed;
# set.seed makes future runs reproducible but gives a different random subset.
set.seed(1)
gl_unlinked <- one_snp_per_radtag(gl_filtered)
write_X_aut_loci(gl_unlinked, x_scaffolds, aut_scaffolds,
                 file.path(out_dir, "full_1SNPperRADtag"))
