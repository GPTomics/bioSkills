#!/usr/bin/env bash
# Reference: GISTIC 2.0.23 | Verify API if version differs
#
# Cohort-level recurrent CNV detection with GISTIC2.
# GISTIC2 is a MATLAB-compiled binary (needs the MATLAB Compiler Runtime); there is no
# conda/pip package. The segment file MUST be diploid-centered before pooling -- a
# mis-centered profile inverts amplification/deletion calls before any statistics run.

set -euo pipefail

SEG=cohort.seg              # 6 cols: sample, chrom, start, end, num_markers, seg.mean (log2)
REFGENE=hg38.refgene.mat    # reference .mat -- MUST match the genome build of the seg file
OUTDIR=gistic_output
mkdir -p "$OUTDIR"

# -conf 0.99   wider, conservative peaks (default 0.75) -> higher confidence the driver is inside
# -brlen 0.7   events longer than 70% of a chromosome arm are treated as broad
# -armpeel 1   peel arm-level events before focal testing (Ziggurat deconstruction)
# -genegistic 1 run the gene-level recurrence test
# -rx 0        keep sex chromosomes
gistic2 \
    -b "$OUTDIR" \
    -seg "$SEG" \
    -refgene "$REFGENE" \
    -genegistic 1 \
    -broad 1 \
    -brlen 0.7 \
    -conf 0.99 \
    -armpeel 1 \
    -savegene 1 \
    -gcm extreme \
    -rx 0

echo "Done. Key outputs in $OUTDIR:"
echo "  amp_genes.txt / del_genes.txt  - peak regions, q-values, genes"
echo "  all_lesions.txt                - per-sample focal/broad calls per peak"
echo "  broad_values_by_arm.txt        - arm-level event matrix"
echo
echo "Interpretation reminders:"
echo "  - q-values are cohort-size dependent; compare recurrence FREQUENCY across cohorts."
echo "  - A wide peak localizes a region, not a gene; cross-check known drivers."
echo "  - Spurious narrow peaks usually mean the input segmentation was over-fragmented."
