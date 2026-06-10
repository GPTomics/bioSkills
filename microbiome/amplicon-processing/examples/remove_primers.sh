#!/usr/bin/env bash
# Reference: cutadapt 4.6+ | Verify API if version differs
# Primers OFF FIRST: this runs BEFORE DADA2 (see dada2_workflow.R). Leftover primer sequence is synthetic,
# often degenerate, and not template - left on, it corrupts the per-run error model, shifts the truncLen
# frame, and masquerades as chimeras. The order primers -> filter -> learnErrors is non-negotiable.
set -euo pipefail

# 16S V4 EMP primers (degenerate IUPAC bases are expected); replace with the primers used for the library.
FWD='GTGYCAGCMGCCGCGGTAA'    # 515F
REV='GGACTACNVGGGTWTCTAAT'   # 806R

raw_dir='raw_reads'
out_dir='trimmed'
mkdir -p "$out_dir"

for r1 in "$raw_dir"/*_R1_001.fastq.gz; do
    r2="${r1/_R1_/_R2_}"
    base=$(basename "$r1" _R1_001.fastq.gz)
    # -g/-G match the forward/reverse primer as 5' adapters on R1/R2;
    # --discard-untrimmed drops pairs lacking the primer (a primerless read is suspect).
    cutadapt \
        -g "$FWD" \
        -G "$REV" \
        --discard-untrimmed \
        -o "$out_dir/${base}_R1_001.fastq.gz" \
        -p "$out_dir/${base}_R2_001.fastq.gz" \
        "$r1" "$r2"
done
