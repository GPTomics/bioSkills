#!/usr/bin/env bash
# Reference: bwa-mem2 2.2.1+, samtools 1.19+ | Verify API if version differs
# DNA alignment with bwa-mem2: read groups + streaming dedup (collate -> fixmate -m -> sort -> markdup).

set -euo pipefail

R1=${1:-reads_1.fastq.gz}
R2=${2:-reads_2.fastq.gz}
REFERENCE=${3:-reference.fa}
OUTPUT=${4:-aligned.markdup.bam}
SAMPLE=${5:-sample}
THREADS=${6:-8}

echo "=== bwa-mem2 alignment + duplicate marking ==="

# Index reference if needed (bwa-mem2 emits .bwt.2bit.64, NOT interchangeable with `bwa index`).
if [[ ! -f "${REFERENCE}.bwt.2bit.64" ]]; then
    echo "Indexing reference..."
    bwa-mem2 index "$REFERENCE"
fi

# Align with read groups (SM/ID/PL/LB are a hard GATK requirement), then the strict dedup order.
# -K 100000000 pins per-batch insert-size estimation -> thread-count-invariant output.
# Do NOT run markdup on amplicon/PCR data (identical ends are by design); use UMIs there instead.
bwa-mem2 mem -t "$THREADS" -K 100000000 \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1" \
    "$REFERENCE" "$R1" "$R2" | \
    samtools collate -@ "$THREADS" -O -u - | \
    samtools fixmate -m -@ "$THREADS" -u - - | \
    samtools sort -@ "$THREADS" -u - | \
    samtools markdup -@ "$THREADS" - "$OUTPUT"

samtools index "$OUTPUT"

# QC gate: confirm reads aligned and are informative before any variant calling.
echo "=== flagstat (alignment + properly-paired + duplicate rate) ==="
samtools flagstat "$OUTPUT"
echo "Output: $OUTPUT"
