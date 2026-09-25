#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 20:00:00
#SBATCH -J clone_filter

# Remove PCR duplicates from the paired-end RAD reads with Stacks clone_filter.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh

module load bioinfo-tools Stacks   # module names as on UPPMAX; adjust for other clusters

mkdir -p "$DECLONED_DIR"

while read -r sample; do
    clone_filter \
        -1 "$CLEAN_READS_DIR/${sample}.1.fq.gz" \
        -2 "$CLEAN_READS_DIR/${sample}.2.fq.gz" \
        -i gzfastq \
        -o "$DECLONED_DIR"

    # clone_filter names its output <sample>.1.1.fq.gz / <sample>.2.2.fq.gz;
    # rename so the next steps find <sample>.1.fq.gz / <sample>.2.fq.gz
    mv "$DECLONED_DIR/${sample}.1.1.fq.gz" "$DECLONED_DIR/${sample}.1.fq.gz"
    mv "$DECLONED_DIR/${sample}.2.2.fq.gz" "$DECLONED_DIR/${sample}.2.fq.gz"

    echo "${sample} processed"
done < "$SAMPLES"
