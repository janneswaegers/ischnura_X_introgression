#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 22:00:00
#SBATCH -J stacks_populations

# Filter SNPs with Stacks populations to create the full and outgroup SNP sets.
# Filters: locus present in 80% of individuals of a population (-r 0.80) and in
# 20 populations (-p 20), minor allele frequency > 0.05, observed heterozygosity < 0.7.
# Samples with mean depth < 20x were removed from the popmaps beforehand.
# --fstats gives the diversity and differentiation statistics (Ho, He, pi, FST, dxy).

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh

module load bioinfo-tools Stacks

FILTERS="-r 0.80 -p 20 --min_maf 0.05 --max_obs_het 0.70"

# Full SNP set: I. elegans and I. graellsii (25 populations)
populations -P "$STACKS_DIR" -M "$POPMAP" -O "$POP_FULL" -t "$THREADS" \
    $FILTERS --vcf --plink --fstats

# Outgroup SNP set: including I. genei, I. fountaineae and I. saharensis (28 populations)
populations -P "$STACKS_DIR_OUTGROUP" -M "$POPMAP_OUTGROUP" -O "$POP_OUTGROUP" -t "$THREADS" \
    $FILTERS --vcf
