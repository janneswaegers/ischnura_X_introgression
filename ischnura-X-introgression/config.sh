#!/bin/bash
# Project settings shared by all scripts in scripts/.
# Edit these paths before running. All scripts are submitted from the repository root,
# e.g.:  sbatch -A <your_account> scripts/01_clone_filter.slurm

PROJECT_DIR="/path/to/project"

# Input
CLEAN_READS_DIR="$PROJECT_DIR/data/clean_reads"        # demultiplexed reads (process_radtags): <sample>.1.fq.gz / <sample>.2.fq.gz
REF_INDEX="$PROJECT_DIR/reference/ischnura_index"      # Bowtie2 index of the I. elegans genome assembly (Chauhan et al. 2021)
X_SCAFFOLDS="$PROJECT_DIR/reference/X_chrom_ids"       # one X-linked scaffold name per line
SAMPLES="$PROJECT_DIR/info/samples.txt"                # one sample ID per line
POPMAP="$PROJECT_DIR/info/popmap_all.txt"              # I. elegans + I. graellsii samples (25 populations)
POPMAP_OUTGROUP="$PROJECT_DIR/info/popmap_outgroup.txt"  # as above + I. genei, I. fountaineae, I. saharensis (28 populations)
POPMAP_ALLOPATRIC="$PROJECT_DIR/info/popmap_allopatric_species.txt"  # allopatric samples only, grouped by species

# Output
DECLONED_DIR="$PROJECT_DIR/results/decloned_reads"
BAM_DIR="$PROJECT_DIR/results/sorted_bam"
STACKS_DIR="$PROJECT_DIR/results/stacks_ref"
STACKS_DIR_OUTGROUP="$PROJECT_DIR/results/stacks_ref_outgroup"
POP_FULL="$PROJECT_DIR/results/populations_full"
POP_OUTGROUP="$PROJECT_DIR/results/populations_outgroup"
POP_FIXED="$PROJECT_DIR/results/populations_fixed"
POP_DIAGNOSTIC="$PROJECT_DIR/results/populations_diagnostic"
R_RESULTS="$PROJECT_DIR/results/R"                     # locus lists written by the R filtering scripts
ANALYSIS_DIR="$PROJECT_DIR/results/analyses"
DSUITE_DIR="$PROJECT_DIR/results/dsuite"               # one sub-folder per ABBA-BABA analysis

THREADS=16
