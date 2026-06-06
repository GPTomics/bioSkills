#!/bin/bash
# Reference: minimap2 2.28+, samtools 1.19+ | Verify API if version differs
# Align a Dorado modified-basecalling BAM WITHOUT losing the MM/ML methylation tags.
# This is the single most common silent failure in long-read methylation pipelines:
# samtools fastq strips tags unless -T, and minimap2 ignores them unless -y.
set -euo pipefail

MOD_BAM=${1:?Usage: $0 <dorado_mod.bam> <ref.fa> [output.bam] [threads]}
REFERENCE=${2:?reference required}
OUTPUT=${3:-meth.aligned.bam}
THREADS=${4:-16}

# Confirm the input actually carries modification tags before spending compute. Capture a
# count (grep -cm1) rather than piping through `head | grep -q`, which would return 141 under
# `set -o pipefail` (samtools killed by SIGPIPE) and invert the check.
if [ "$(samtools view "$MOD_BAM" | grep -cm1 'MM:Z' || true)" -eq 0 ]; then
    echo 'ERROR: input BAM has no MM/ML tags - it was not basecalled with a mods model.'
    exit 1
fi

# -T MM,ML keeps the tags through FASTQ conversion; -y copies them into the alignment; -Y
# soft-clips supplementary records so hard-clipping does not break the MM per-base skip counting.
samtools fastq -T MM,ML "$MOD_BAM" \
    | minimap2 -ax lr:hq -y -Y --MD -t "$THREADS" "$REFERENCE" - \
    | samtools sort -@4 -o "$OUTPUT"
samtools index "$OUTPUT"

# Verify the tags survived alignment - their absence here is the failure to catch.
if [ "$(samtools view "$OUTPUT" | grep -cm1 'MM:Z' || true)" -ge 1 ]; then
    echo "MM/ML tags preserved. Ready for: modkit pileup $OUTPUT ..."
else
    echo 'WARNING: MM/ML tags lost during alignment - check the -T/-y flags.'
fi
