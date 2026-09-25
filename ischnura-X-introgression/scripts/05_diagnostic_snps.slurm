#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 22:00:00
#SBATCH -J diagnostic_snps

# Diagnostic SNP set: loci that are fixed for different alleles in allopatric
# I. elegans and allopatric I. graellsii, genotyped in all samples.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh

module load bioinfo-tools Stacks

mkdir -p "$POP_FIXED" "$POP_DIAGNOSTIC"

# 1. Find fixed differences between the two allopatric groups.
#    The popmap contains only allopatric samples, grouped by species.
populations -P "$STACKS_DIR" -M "$POPMAP_ALLOPATRIC" -O "$POP_FIXED" -t "$THREADS" --phylip

# 2. Turn the list of fixed loci into a whitelist (locus ID + SNP column), dropping the header line.
cut -f 2,3 "$POP_FIXED/populations.fixed.phylip.log" | sed '1d' > "$POP_FIXED/whitelist_fixed.txt"

# 3. Genotype only the whitelisted SNPs in all samples.
populations -P "$STACKS_DIR" -M "$POPMAP" -O "$POP_DIAGNOSTIC" -t "$THREADS" \
    -W "$POP_FIXED/whitelist_fixed.txt" \
    -r 0.80 --min_maf 0.05 --max_obs_het 0.70 --vcf
