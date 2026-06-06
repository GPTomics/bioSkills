#!/bin/bash
# Reference: minimap2 2.28+, samtools 1.19+ | Verify API if version differs
# Map long reads with the ERROR-RATE-matched preset, keeping the tags downstream
# variant/SV callers need (--MD, -Y soft-clipped supplementaries).
set -euo pipefail

REFERENCE=${1:?Usage: $0 <ref.fa> <reads.fq.gz> <platform> [output_prefix] [threads]}
READS=${2:?reads required}
PLATFORM=${3:?platform required: ont-r9 | ont-r10 | hifi | clr}
OUTPUT=${4:-aligned}
THREADS=${5:-16}
SAMPLE=$(basename "$OUTPUT")

# Preset is a statement about read error rate, not just platform.
case "$PLATFORM" in
    ont-r9)  PRESET=map-ont ;;   # noisy R9 / fast / hac
    ont-r10) PRESET=lr:hq ;;     # accurate R10 sup / Q20+ / duplex (2.27+)
    hifi)    PRESET=map-hifi ;;  # PacBio HiFi/CCS
    clr)     PRESET=map-pb ;;    # PacBio CLR (legacy); never use for HiFi
    *) echo "Unknown platform: $PLATFORM"; exit 1 ;;
esac

# Bake the SAME preset's k/w into the index, or the preset's k/w is silently ignored.
minimap2 -x "$PRESET" -d "${REFERENCE}.${PRESET//:/_}.mmi" "$REFERENCE"

# -Y keeps full SEQ on supplementary (split) reads so SV callers recover breakpoints.
# --MD gives mismatch positions for small-variant callers (needs minimap2 >= 2.28).
minimap2 -ax "$PRESET" -t "$THREADS" --MD -Y \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
    "${REFERENCE}.${PRESET//:/_}.mmi" "$READS" \
    | samtools sort -@4 -o "${OUTPUT}.bam"
samtools index "${OUTPUT}.bam"

echo 'Alignment statistics:'
samtools flagstat "${OUTPUT}.bam"
