#!/bin/bash
# ABBA-BABA statistics (D, f_d, f_dM) in sliding windows with Dsuite Dinvestigate,
# using the outgroup SNP set and females only.
#
# Each sub-folder of $DSUITE_DIR is one analysis (sympatric I. elegans or I. graellsii,
# at three admixture (Q) cut-offs) and contains:
#   SETS.txt        sample <tab> group (P1/P2/P3 groups, "Outgroup", or "xxx" to ignore a sample, e.g. males)
#   test_trios.txt  P1 <tab> P2 <tab> P3, e.g. elegans_allo  elegans_sym  graellsii_allo
#
# Windows of 50 informative SNPs, moving in steps of 25 SNPs.

set -euo pipefail
source config.sh

WINDOW="50,25"
VCF="$POP_OUTGROUP/populations.snps.vcf"

for dir in "$DSUITE_DIR"/*/; do
    echo "Running Dinvestigate in $dir"
    (cd "$dir" && Dsuite Dinvestigate -w "$WINDOW" "$VCF" SETS.txt test_trios.txt)

    # Label each window as X-linked or autosomal (column "Xaut") for the R analysis
    for f in "$dir"/*localFstats_*_50_25.txt; do
        awk 'NR == FNR { x[$1]; next }
             FNR == 1  { print $0 "\tXaut"; next }
                       { print $0 "\t" (($1 in x) ? "X" : "aut") }' \
            "$X_SCAFFOLDS" "$f" > "${f%.txt}_Xaut.txt"
    done
done
