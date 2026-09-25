# Helper functions shared by the R scripts.
#
# Genotypes in a genlight object are coded as the number of alternative alleles:
# 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternative, NA = missing.

# Read a Stacks VCF into a genlight object --------------------------------------
read_genlight <- function(vcf_file) {
  vcfR::vcfR2genlight(vcfR::read.vcfR(vcf_file, verbose = FALSE))
}

# Remove individuals by ID ------------------------------------------------------
drop_samples <- function(gl, ids) {
  gl[!(adegenet::indNames(gl) %in% ids)]
}

# Proportion of homozygous genotypes -------------------------------------------
# As in the original analysis, "homozygous" is counted as genotype 0, i.e. the
# expression `geno == 0 & 2` used there, which R evaluates as `geno == 0`.
prop_hom_individual <- function(gl) {          # per individual (rows)
  g <- as.matrix(gl)
  hom <- rowSums(g == 0, na.rm = TRUE)
  het <- rowSums(g == 1, na.rm = TRUE)
  hom / (hom + het)
}

prop_hom_snp <- function(gl) {                 # per SNP (columns)
  g <- as.matrix(gl)
  hom <- colSums(g == 0, na.rm = TRUE)
  het <- colSums(g == 1, na.rm = TRUE)
  hom / (hom + het)
}

# Remove X-linked SNPs where any male is heterozygous ---------------------------
# Males are hemizygous for the X (X0), so a heterozygous call in a male at an
# X-linked SNP indicates a genotyping error or a mis-assigned scaffold.
filter_X_hom_males <- function(gl, is_male, x_scaffolds) {
  is_X <- as.character(gl@chromosome) %in% x_scaffolds
  all_males_hom <- prop_hom_snp(gl[is_male, ]) == 1
  keep <- ifelse(!all_males_hom & is_X, FALSE, TRUE)
  message(sprintf("X-linked SNPs: %d; removed (heterozygous in >= 1 male): %d; SNPs kept: %d",
                  sum(is_X), sum(!keep, na.rm = TRUE), sum(keep, na.rm = TRUE)))
  gl[, keep]
}

# Keep one random SNP per RAD-tag -----------------------------------------------
# Stacks SNP IDs look like "<RAD-tag ID>:<position>:<strand>".
one_snp_per_radtag <- function(gl) {
  loci <- data.frame(locus = adegenet::locNames(gl))
  loci$radtag <- sapply(strsplit(loci$locus, ":"), `[`, 1)
  chosen <- loci |>
    dplyr::group_by(radtag) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::pull(locus)
  gl[, adegenet::locNames(gl) %in% chosen]
}

# Write X-linked and autosomal SNP ID lists (input for vcftools --snps) ---------
write_X_aut_loci <- function(gl, x_scaffolds, aut_scaffolds, out_prefix) {
  chrom <- as.character(gl@chromosome)
  loci  <- adegenet::locNames(gl)
  writeLines(loci[chrom %in% x_scaffolds],   paste0(out_prefix, "_X_loci.txt"))
  writeLines(loci[chrom %in% aut_scaffolds], paste0(out_prefix, "_aut_loci.txt"))
  message(sprintf("Written: %d X-linked and %d autosomal SNPs",
                  sum(chrom %in% x_scaffolds), sum(chrom %in% aut_scaffolds)))
}

# Permutation test: X-linked versus autosomal median ----------------------------
# Draws n_perm random sets of autosomal values (without replacement), each the same
# size as the X-linked set, and compares their medians with the observed X median.
#   alternative = "greater": p = proportion of permuted medians > X median
#                            (X higher, e.g. BGC beta: stronger barrier to gene flow)
#   alternative = "less":    p = proportion of permuted medians < X median
#                            (X lower, e.g. f_dM: less introgression)
perm_median_test <- function(x_values, aut_values,
                             alternative = c("greater", "less"), n_perm = 10000) {
  alternative <- match.arg(alternative)
  x_values   <- x_values[!is.na(x_values)]
  aut_values <- aut_values[!is.na(aut_values)]
  obs  <- median(x_values)
  perm <- replicate(n_perm, median(sample(aut_values, length(x_values), replace = FALSE)))
  p <- if (alternative == "greater") mean(perm > obs) else mean(perm < obs)
  list(median_X = obs, median_aut = median(aut_values), p = p, perm = perm)
}
