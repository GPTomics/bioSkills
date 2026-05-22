#!/usr/bin/env bash
# Reference: CNVkit 0.9.10+, samtools 1.19+ | Verify API if version differs
#
# CNVkit tumor-normal CNV calling with a pooled panel of normals.
# Demonstrates the decision-grade defaults: pooled reference over flat reference,
# accessible-genome restriction, low-coverage dropout handling, and purity-aware
# integer calling. Edit the variables below for the local dataset.

set -euo pipefail

PANEL_BED=panel.bed                 # capture / target regions
FASTA=reference.fa                  # same reference used for alignment
REFFLAT=refFlat.txt                 # gene annotation (optional but recommended)
NORMAL_GLOB="normals/*.bam"         # process-matched normals for the PoN (>= 5)
TUMOR_BAM=tumor.bam
TUMOR_PURITY=0.65                   # from pathology or an allele-specific caller
OUTDIR=cnvkit_results
THREADS=8

mkdir -p "$OUTDIR"

# Step 1: define the accessible genome once (mappable, non-gap regions).
# Antitarget bins are drawn only from accessible regions; skipping this inflates noise.
cnvkit.py access "$FASTA" -o "$OUTDIR/access.bed"

# Step 2: build a pooled reference from process-matched normals.
# A pooled PoN averages out per-normal capture noise. A flat reference (no --normal)
# corrects only GC content and leaves capture bias as recurrent false focal calls.
cnvkit.py batch --normal $NORMAL_GLOB \
    --targets "$PANEL_BED" \
    --annotate "$REFFLAT" \
    --fasta "$FASTA" \
    --access "$OUTDIR/access.bed" \
    --output-reference "$OUTDIR/pooled_reference.cnn" \
    -p "$THREADS"

# Step 3: call CNVs on the tumor against the pooled reference.
# --drop-low-coverage prevents FFPE/low-input zero-coverage bins from being read as
# homozygous deletions.
cnvkit.py batch "$TUMOR_BAM" \
    --reference "$OUTDIR/pooled_reference.cnn" \
    --output-dir "$OUTDIR" \
    --drop-low-coverage \
    --scatter --diagram \
    -p "$THREADS"

SAMPLE=$(basename "$TUMOR_BAM" .bam)

# Step 4: re-segment with an explicit method. `batch` already produced a default-CBS
# .cns; this step overrides it for explicit method control. CBS is precise on focal
# events; switch to hmm-tumor for impure or heterogeneous tumors where CBS over-fragments.
cnvkit.py segment "$OUTDIR/${SAMPLE}.cnr" -m cbs --drop-low-coverage \
    -o "$OUTDIR/${SAMPLE}.cns"

# Step 5: purity-aware integer calling. The clonal method rescales by purity before
# rounding; without it an impure tumor's true CN=4 rounds down to CN=2 or 3.
cnvkit.py call "$OUTDIR/${SAMPLE}.cns" \
    -m clonal --purity "$TUMOR_PURITY" --ploidy 2 \
    -o "$OUTDIR/${SAMPLE}.call.cns"

# Step 6: QC. MAD < 0.5 is the minimum bar for confident focal calls; > 0.5 is unreliable.
cnvkit.py metrics "$OUTDIR/${SAMPLE}.cnr" -s "$OUTDIR/${SAMPLE}.cns"
cnvkit.py sex "$OUTDIR/${SAMPLE}.cnr"

# Step 7: gene-level report with bootstrap confidence intervals.
cnvkit.py genemetrics "$OUTDIR/${SAMPLE}.cnr" -s "$OUTDIR/${SAMPLE}.cns" \
    -t 0.2 --ci --bootstrap 100 -o "$OUTDIR/${SAMPLE}.genemetrics.tsv"

echo "Done. Purity-corrected calls: $OUTDIR/${SAMPLE}.call.cns"
echo "If absolute CN, LOH, or whole-genome-doubling status is needed, escalate to an"
echo "allele-specific caller (see copy-number/allele-specific-copy-number)."
