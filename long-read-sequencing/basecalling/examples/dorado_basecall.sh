#!/bin/bash
# Reference: Dorado 1.0+, samtools 1.19+ | Verify API if version differs
# Methylation-aware basecalling and barcode demultiplexing - the two operations
# whose order/flags most often go wrong.
set -euo pipefail

POD5_DIR=${1:?Usage: $0 <pod5_dir> [output_dir] [kit_name]}
OUTPUT=${2:-dorado_output}
KIT=${3:-}              # e.g. SQK-NBD114-24; leave empty for a non-barcoded run

mkdir -p "$OUTPUT"

# Methylation MUST be requested at basecall time - it cannot be recovered from a plain
# BAM. 5mCG_5hmCG calls CpG 5mC and 5hmC; the result carries MM/ML tags. KEEP the POD5.
echo 'Basecalling with CpG methylation...'
if [ -n "$KIT" ]; then
    # Barcoded: basecall WITHOUT trimming, then demux (demux trims barcodes itself).
    # Trimming during basecalling would strip barcodes before demux can read them.
    dorado basecaller sup,5mCG_5hmCG "$POD5_DIR" --no-trim > "$OUTPUT/calls.bam"
    dorado demux --kit-name "$KIT" --output-dir "$OUTPUT/demux/" "$OUTPUT/calls.bam"
    echo "Per-barcode BAMs (with MM/ML tags preserved): $OUTPUT/demux/"
else
    dorado basecaller sup,5mCG_5hmCG "$POD5_DIR" > "$OUTPUT/calls.bam"
fi

# Confirm the modification tags survived - their absence is the silent failure mode.
echo 'Checking for MM/ML methylation tags...'
# Capture a count (grep -cm1); `head | grep -q` would return 141 under `set -o pipefail`
# (samtools killed by SIGPIPE) and invert the check.
if [ "$(samtools view "$OUTPUT/calls.bam" | grep -cm1 'MM:Z' || true)" -ge 1 ]; then
    echo 'MM/ML tags present. Carry them through alignment with minimap2 -y -Y (see nanopore-methylation).'
else
    echo 'WARNING: no MM/ML tags found - was a modification model requested?'
fi

echo "Done: $OUTPUT/calls.bam"
