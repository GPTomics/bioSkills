#!/bin/bash
# Reference: NanoPlot 1.42+, cramino 0.14+, chopper 0.7+, seqkit 2.5+ | Verify API if version differs
# Long-read QC: FASTQ overview (length/posterior-Q), real identity from a BAM, and an
# intent-conditioned filter. Read-only Qscore overstates accuracy - identity needs alignment.
set -euo pipefail

READS=${1:?Usage: $0 <reads.fq.gz> [aligned.bam] [output_dir] [goal]}
BAM=${2:-}                # optional: aligned BAM for REAL percent identity
OUTPUT_DIR=${3:-qc_output}
GOAL=${4:-variant}        # variant | assembly (sets the filter strategy)
mkdir -p "$OUTPUT_DIR"

echo 'Length / yield overview (FASTQ Qscore is a posterior, not real accuracy):'
seqkit stats -a "$READS"
NanoPlot --fastq "$READS" -o "${OUTPUT_DIR}/nanoplot" --N50 --plots dot

# Real accuracy needs alignment - gap-compressed identity from the BAM.
if [ -n "$BAM" ]; then
    echo 'Real gap-compressed identity from the BAM:'
    cramino "$BAM"
fi

# Filter conditioned on intent. Assembly preserves long reads + small replicons (subsample
# by QUALITY); variant calling filters almost nothing.
if [ "$GOAL" = assembly ]; then
    echo 'Assembly: light pass, then quality-subsample (no hard length floor that erases plasmids)...'
    chopper -q 10 -l 1000 -i "$READS" | gzip > "${OUTPUT_DIR}/clean.fq.gz"
    # Set --target_bases to ~100x your genome size. Filtlong scores length AND quality equally
    # by default; --mean_q_weight 10 tilts the subsampling toward quality (the assembly goal).
    filtlong --target_bases 500000000 --mean_q_weight 10 "${OUTPUT_DIR}/clean.fq.gz" \
        | gzip > "${OUTPUT_DIR}/filtered.fq.gz"
else
    echo 'Variant calling: light quality only, keep depth...'
    chopper -q 10 -i "$READS" | gzip > "${OUTPUT_DIR}/filtered.fq.gz"
fi

echo 'Filtered stats:'
seqkit stats -a "${OUTPUT_DIR}/filtered.fq.gz"
echo 'For run health (pore death, speed drift), run: pycoQC -f sequencing_summary.txt -o run.html'
