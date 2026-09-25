# Filtering of the diagnostic SNP set
#
# Same steps as R/01_full_snp_set.R, applied to the diagnostic SNPs (loci with fixed
# differences between allopatric I. elegans and I. graellsii, scripts/05_diagnostic_snps.slurm):
# remove X-linked SNPs that are heterozygous in any male, keep one random SNP per
# RAD-tag and write X-linked / autosomal SNP lists.
# Also writes the sex file needed for ADMIXTURE on the X (scripts/07_pca_admixture.sh).
#
# Output: results/R/diagnostic/

library(adegenet)
library(vcfR)
library(dplyr)
source("R/functions.R")

# --- Settings -----------------------------------------------------------------
vcf_file         <- "data/populations_diagnostic/populations.snps.vcf"
sample_info_file <- "data/sample_info_diagnostic.txt"   # columns: ID, PopID, Sex
x_scaffolds      <- read.table("data/X_chrom_ids")$V1
aut_scaffolds    <- read.table("data/Autosome_ids")$V1
out_dir          <- "results/R/diagnostic"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gl <- read_genlight(vcf_file)

# --- Remove excluded samples ----------------------------------------------------
excluded <- c("No6", "S8", "BEL-2", "VMM.12", "Pop1_1", "Pop1_4", "Pop1_9", "LSA.04",
              "RGC.07", "RGC.08", "Pop1_6", "LSA.10", "M5", "M6", "I8", "K7", "No4",
              "K10", "K5")
writeLines(excluded, file.path(out_dir, "excluded_samples.txt"))
gl <- drop_samples(gl, excluded)

sample_info <- read.table(sample_info_file, header = TRUE)
sample_info <- sample_info[match(indNames(gl), sample_info$ID), ]
is_male <- tolower(sample_info$Sex) == "male"
is_male[is.na(is_male)] <- FALSE

# PLINK sex file (FID, IID, sex: 1 = male, 2 = female) for haploid males in ADMIXTURE
write.table(data.frame(indNames(gl), indNames(gl), ifelse(is_male, 1, 2)),
            file.path(out_dir, "plink_sex.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# --- Filter X-linked SNPs and keep one SNP per RAD-tag --------------------------
gl_filtered <- filter_X_hom_males(gl, is_male, x_scaffolds)
writeLines(locNames(gl_filtered), file.path(out_dir, "loc_filtered.txt"))

set.seed(1)   # see note in R/01_full_snp_set.R
gl_unlinked <- one_snp_per_radtag(gl_filtered)
# Paper: 1,931 diagnostic SNPs, of which 111 X-linked
write_X_aut_loci(gl_unlinked, x_scaffolds, aut_scaffolds,
                 file.path(out_dir, "diagnostic_1SNPperRADtag"))
