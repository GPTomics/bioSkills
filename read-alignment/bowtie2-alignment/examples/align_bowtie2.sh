#!/usr/bin/env bash
# Reference: bowtie2 2.5+, samtools 1.19+ | Verify API if version differs
# Bowtie2 alignment for ChIP-seq / ATAC-seq: mode + fragment-geometry flags drive the peak coordinates.

set -euo pipefail

R1=${1:-reads_1.fastq.gz}
R2=${2:-reads_2.fastq.gz}
INDEX=${3:-bt2_index}            # index BASENAME, not a .bt2 file
OUTPUT=${4:-aligned.bam}
ASSAY=${5:-atac}                 # atac | chip
THREADS=${6:-8}

echo "=== Bowtie2 alignment ($ASSAY) ==="

# ATAC: local mode soft-clips adapter read-through; --dovetail + wide -X admit short nucleosome fragments.
# ChIP (trimmed): end-to-end is fine. Both drop singletons/discordants for clean fragment-level signal.
if [[ "$ASSAY" == "atac" ]]; then
    MODE=(--local --dovetail -X 2000)
else
    MODE=(--end-to-end)
fi

# MAPQ 30 drops multimappers on the Bowtie2 scale (max 42 e2e / 44 local -- NOT a BWA-style 60).
bowtie2 -p "$THREADS" --very-sensitive --no-mixed --no-discordant "${MODE[@]}" \
    -x "$INDEX" -1 "$R1" -2 "$R2" \
    2> "${OUTPUT%.bam}_stats.txt" | \
    samtools view -@ "$THREADS" -bS -q 30 -F 1804 - | \
    samtools sort -@ "$THREADS" -o "$OUTPUT" -

samtools index "$OUTPUT"

echo "Output: $OUTPUT"
cat "${OUTPUT%.bam}_stats.txt"
# The Tn5 +4/-5 cut-site shift for ATAC is a downstream signal-track transform -> atac-seq.
