#!/usr/bin/env bash
# Reference: AmpliconSuite-pipeline 1.3+, AmpliconArchitect 1.3+, CNVkit 0.9.10+ | Verify API if version differs
#
# Reconstruct focal-amplification architecture (ecDNA / BFB / HSR / linear) from WGS.
# AmpliconArchitect does not call amplifications from scratch -- it reconstructs the
# breakpoint graph of the high-copy focal SEEDS it is given. Garbage seeds produce
# garbage amplicons, so seeds must be vetted focal high-CN regions.
#
# Requires $AA_DATA_REPO (reference data) and a free academic Mosek license configured.

set -euo pipefail

SAMPLE=sample_id
TUMOR_BAM=tumor.bam        # whole-genome sequencing BAM (AA is not for panels/WES)
REF=GRCh38                 # must match the BAM build and the installed $AA_DATA_REPO
THREADS=8
SEED_CN=4.5                # minimum copy number for a focal amplicon seed
OUTDIR=ampliconsuite_out
mkdir -p "$OUTDIR"

# Step 1: call copy number and derive focal high-CN seeds. A pooled panel of normals is
# strongly preferred -- a flat-reference callset produces false high-CN seeds.
# cnvkit.py batch names its output from the BAM basename, then `call` adds the integer
# `cn` column.
BN=$(basename "$TUMOR_BAM" .bam)
cnvkit.py batch "$TUMOR_BAM" --method wgs --reference pooled_reference.cnn \
    --output-dir "$OUTDIR/cnvkit" --drop-low-coverage -p "$THREADS"
cnvkit.py call "$OUTDIR/cnvkit/${BN}.cns" -o "$OUTDIR/cnvkit/${BN}.call.cns"

# Keep focal segments above the seed copy-number threshold. The < 10 Mb size cutoff
# keeps seeds focal (sub-chromosome-arm); arm-level seeds dilute the focal amplicon
# signal and produce spurious amplicons. .call.cns columns: chrom/start/end/gene/log2/cn/...
# AmpliconSuite --cnv_bed expects 4 columns: chrom, start, end, copy_number.
FOCAL_MAX_BP=10000000   # ~sub-arm; segments larger than this are treated as arm-level
awk -v cn="$SEED_CN" -v maxbp="$FOCAL_MAX_BP" \
    'NR>1 && $6 >= cn && ($3-$2) < maxbp {print $1"\t"$2"\t"$3"\t"$6}' \
    "$OUTDIR/cnvkit/${BN}.call.cns" > "$OUTDIR/focal_seeds.bed"

# Step 2: run AmpliconArchitect + AmpliconClassifier on the vetted seeds.
AmpliconSuite-pipeline.py \
    -s "$SAMPLE" \
    -t "$THREADS" \
    --bam "$TUMOR_BAM" \
    --ref "$REF" \
    --cnv_bed "$OUTDIR/focal_seeds.bed" \
    --run_AA --run_AC \
    -o "$OUTDIR"

echo "Done. Key outputs in $OUTDIR:"
echo "  *_amplicon*_graph.txt   - breakpoint graphs (closed cycle = circular evidence)"
echo "  *_amplicon*_cycles.txt  - reconstructed amplicon cycles"
echo "  *_classification*.tsv   - AmpliconClassifier ecDNA / BFB / HSR / linear calls"
echo
echo "An ecDNA call needs a closed-cycle graph AND a classifier ecDNA label. Confirm"
echo "high-stakes calls with FISH, single-cell copy number, or optical mapping."
