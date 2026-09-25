#!/bin/bash
# Population structure:
#   - PCA on autosomal SNPs of the full SNP set (PLINK 1.9)
#   - supervised ADMIXTURE (K = 2) on autosomal and X-linked diagnostic SNPs,
#     with allopatric I. elegans and I. graellsii as reference samples

set -euo pipefail
source config.sh

module load bioinfo-tools vcftools plink/1.90b4.9 ADMIXTURE/1.3.0

cd "$ANALYSIS_DIR"

# --- PCA -------------------------------------------------------------------
vcftools --vcf full_aut.recode.vcf --plink --out full_aut
plink --file full_aut --make-bed --allow-extra-chr --out full_aut
plink --bfile full_aut --pca --allow-extra-chr --out full_aut_pca

# --- Supervised ADMIXTURE ----------------------------------------------------
# Needs a <prefix>.pop file next to each .bed: one line per individual (same order as the .fam),
# "elegans" / "graellsii" for allopatric reference individuals and "-" for sympatric ones.
for part in aut X; do
    vcftools --vcf diagnostic_${part}.recode.vcf --plink --out diagnostic_${part}
    plink --file diagnostic_${part} --make-bed --allow-extra-chr --out diagnostic_${part}
    # ADMIXTURE does not accept scaffold names as chromosome codes: set them to 0
    awk '{$1=0; print $0}' diagnostic_${part}.bim > tmp.bim && mv tmp.bim diagnostic_${part}.bim
done

admixture -C 0.000001 --supervised diagnostic_aut.bed 2 | tee admixture_aut.log

# X chromosome: add genotypic sex (written by R/02_diagnostic_snp_set.R) and
# treat all males as haploid, since they carry a single X (X0 system)
plink --bfile diagnostic_X --update-sex "$R_RESULTS/diagnostic/plink_sex.txt" \
    --make-bed --allow-extra-chr --out diagnostic_X_sex
cp diagnostic_X.pop diagnostic_X_sex.pop
admixture -C 0.000001 --supervised --haploid="male:*" diagnostic_X_sex.bed 2 | tee admixture_X.log
