#!/bin/bash
# Make separate VCFs for X-linked and autosomal SNPs.
# Run after R/01_full_snp_set.R and R/02_diagnostic_snp_set.R, which write the locus lists
# (one random SNP per RAD-tag, X-linked SNPs heterozygous in any male removed).

set -euo pipefail
source config.sh

module load bioinfo-tools vcftools

mkdir -p "$ANALYSIS_DIR"

split_vcf () {
    local vcf=$1 set=$2
    for part in X aut; do
        vcftools --vcf "$vcf" \
            --snps "$R_RESULTS/${set}/${set}_1SNPperRADtag_${part}_loci.txt" \
            --remove "$R_RESULTS/${set}/excluded_samples.txt" \
            --recode --recode-INFO-all \
            --out "$ANALYSIS_DIR/${set}_${part}"
    done
}

split_vcf "$POP_FULL/populations.snps.vcf"       full
split_vcf "$POP_DIAGNOSTIC/populations.snps.vcf" diagnostic
