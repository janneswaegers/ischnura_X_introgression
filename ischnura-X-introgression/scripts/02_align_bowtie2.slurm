#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -J bowtie2_align

# Align decloned reads to the I. elegans draft genome with Bowtie2 and sort the alignments.
# Settings: 1 mismatch per seed (-N 1), mismatch penalty max 6 / min 2 (--mp 6,2),
# fragment length 100-1000 bp (-I 100 -X 1000).

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source config.sh

module load bioinfo-tools bowtie2 samtools

mkdir -p "$BAM_DIR/logs"

while read -r sample; do
    bowtie2 -N 1 --mp 6,2 -I 100 -X 1000 -p "$THREADS" \
        -x "$REF_INDEX" \
        -1 "$DECLONED_DIR/${sample}.1.fq.gz" \
        -2 "$DECLONED_DIR/${sample}.2.fq.gz" \
        2> "$BAM_DIR/logs/${sample}.bowtie2.log" \
    | samtools sort -@ 2 -o "$BAM_DIR/${sample}.bam" -

    echo "${sample} aligned"
done < "$SAMPLES"
