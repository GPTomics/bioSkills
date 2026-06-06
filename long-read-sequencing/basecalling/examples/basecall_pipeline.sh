#!/bin/bash
# Reference: Dorado 1.0+, pod5 0.3+, samtools 1.19+, chopper 0.7+ | Verify API if version differs
# Simplex basecalling pipeline: POD5 -> reads -> filtered FASTQ -> QC.
# The MODEL string is propagated to downstream tools; pin it for reproducibility.
set -euo pipefail

INPUT=${1:?Usage: $0 <input_dir> [output_dir] [tier] [min_qual] [min_len]}
OUTPUT=${2:-basecall_output}
TIER=${3:-sup}          # sup for analysis; hac for routine; fast for previews only
MIN_QUAL=${4:-10}       # Q10 ~ 90% nominal; a permissive QC floor, not a hard rule
MIN_LEN=${5:-500}       # drop sub-500 bp reads; below typical long-read utility

mkdir -p "$OUTPUT"

# FAST5 is legacy and slow to basecall directly; convert to POD5 first.
if ls "$INPUT"/*.fast5 >/dev/null 2>&1; then
    echo 'Converting FAST5 to POD5...'
    mkdir -p "$OUTPUT/pod5"
    pod5 convert fast5 "$INPUT"/*.fast5 --output "$OUTPUT/pod5/"
    POD5_DIR="$OUTPUT/pod5"
else
    POD5_DIR="$INPUT"
fi

# Bare tier auto-detects chemistry from POD5 metadata. BAM is the default output.
echo "Basecalling at $TIER..."
dorado basecaller "$TIER" "$POD5_DIR" > "$OUTPUT/calls.bam"

# Record the model used so downstream medaka/Clair3 models can match it. The model
# string lives in the BAM @RG header (DS:basecall_model=...), not in `dorado summary`.
samtools view -H "$OUTPUT/calls.bam" | grep -oE 'basecall_model=[^[:space:]]+' \
    > "$OUTPUT/basecaller_model.txt" || true

echo 'Converting to FASTQ and filtering...'
samtools fastq "$OUTPUT/calls.bam" | gzip > "$OUTPUT/calls.fastq.gz"
gunzip -c "$OUTPUT/calls.fastq.gz" \
    | chopper -q "$MIN_QUAL" -l "$MIN_LEN" \
    | gzip > "$OUTPUT/filtered.fastq.gz"

echo 'QC report...'
NanoPlot --fastq "$OUTPUT/filtered.fastq.gz" -o "$OUTPUT/qc/" --plots dot

echo "Done. Reads: $(samtools view -c "$OUTPUT/calls.bam"); QC: $OUTPUT/qc/NanoPlot-report.html"
