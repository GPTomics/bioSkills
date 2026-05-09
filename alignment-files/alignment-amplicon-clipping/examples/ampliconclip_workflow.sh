#!/bin/bash
# Reference: samtools 1.19+ | Verify API if version differs
# Standard amplicon primer-clipping workflow:
#   ampliconclip -> sort by name -> fixmate -> sort by coord -> calmd -> index
# Inputs: input BAM (coord-sorted, indexed), primer BED (with strand in col 6),
#         reference FASTA (indexed with samtools faidx)

set -euo pipefail

BAM=${1:?usage: $0 input.bam primers.bed reference.fa output.bam}
PRIMERS=${2:?primer BED required}
REF=${3:?reference FASTA required}
OUT=${4:?output BAM required}
THREADS=${THREADS:-4}

WORK=$(mktemp -d)
trap "rm -rf $WORK" EXIT

# 1. Soft-clip primers from BED.
#    --both-ends:   primers may appear at 5' or 3' end of read
#    --strand:      respects column-6 strand in BED (only clips matching strand)
#    --soft-clip:   reversible (CIGAR S, bases retained); preferred over --hard-clip
samtools ampliconclip \
    --both-ends \
    --strand \
    --soft-clip \
    -b "$PRIMERS" \
    "$BAM" \
    -o "$WORK/clipped.bam"

# 2. Repair mate info (CIGARs changed by clipping invalidate MC, MS, TLEN).
#    collate (faster) -> fixmate -m -> coord sort.
samtools collate -O -u "$WORK/clipped.bam" "$WORK/collate" | \
    samtools fixmate -m -u - - | \
    samtools sort -@ "$THREADS" -o "$WORK/sorted.bam" -

# 3. Repair MD/NM tags (clipping invalidates them; bcftools mpileup BAQ
#    and IGV mismatch coloring depend on these).
samtools calmd -@ "$THREADS" -b "$WORK/sorted.bam" "$REF" > "$OUT" 2>/dev/null

# 4. Index for region access.
samtools index -@ "$THREADS" "$OUT"

echo "Clipped BAM: $OUT"
echo "Diagnostic: regions covered by primers should now show CIGAR S ops at read 5' ends:"
echo "  samtools view $OUT | awk '\$6 ~ /^[0-9]+S/' | head"
