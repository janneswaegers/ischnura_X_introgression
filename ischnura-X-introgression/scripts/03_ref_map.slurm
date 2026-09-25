#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -J ref_map

# Build loci and call SNPs from the aligned reads with Stacks ref_map (default parameters).
# Run twice: once for I. elegans + I. graellsii only, and once including the outgroup species.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh

module load bioinfo-tools Stacks

mkdir -p "$STACKS_DIR" "$STACKS_DIR_OUTGROUP"

ref_map.pl -T "$THREADS" --samples "$BAM_DIR" --popmap "$POPMAP"          -o "$STACKS_DIR"
ref_map.pl -T "$THREADS" --samples "$BAM_DIR" --popmap "$POPMAP_OUTGROUP" -o "$STACKS_DIR_OUTGROUP"
