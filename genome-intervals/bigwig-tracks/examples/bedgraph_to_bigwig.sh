#!/bin/bash
# Reference: ucsc-bedgraphtobigwig 469+, ucsc-bigwiginfo 469+ | Verify API if version differs
# Build a valid bigWig from a coverage bedGraph, then sanity-check it.
# Usage: ./bedgraph_to_bigwig.sh input.bedGraph chrom.sizes output.bw

BEDGRAPH="${1:-coverage.bedGraph}"
CHROM_SIZES="${2:-chrom.sizes}"
OUTPUT="${3:-output.bw}"

if [[ ! -f "$BEDGRAPH" ]]; then
    echo "bedGraph not found: $BEDGRAPH"
    echo "Usage: $0 <bedGraph> <chrom.sizes> <output.bw>"
    exit 1
fi

if [[ ! -f "$CHROM_SIZES" ]]; then
    echo "chrom.sizes not found: $CHROM_SIZES"
    echo "Create from FASTA index: cut -f1,2 reference.fa.fai > chrom.sizes"
    echo "Or fetch from UCSC:       fetchChromSizes hg38 > chrom.sizes"
    exit 1
fi

# bedGraphToBigWig REQUIRES coordinate-sorted, non-overlapping input (signal is one value per base).
# -c verifies sort order without rewriting; only sort if needed.
echo "=== Checking sort order ==="
if sort -k1,1 -k2,2n -c "$BEDGRAPH" 2>/dev/null; then
    echo "already sorted"
    SORTED="$BEDGRAPH"
else
    SORTED="${BEDGRAPH%.bedGraph}.sorted.bedGraph"
    echo "sorting -> $SORTED"
    sort -k1,1 -k2,2n "$BEDGRAPH" > "$SORTED"
fi

echo "=== Converting to bigWig ==="
# chrom.sizes names/lengths must match the bedGraph exactly (chr1 vs 1 mismatch -> dropped intervals or error).
bedGraphToBigWig "$SORTED" "$CHROM_SIZES" "$OUTPUT"

echo "=== bigWigInfo sanity check ==="
# Confirm zoom levels exist (browser zoom-out reads them) and basesCovered vs genome size
# (low coverage means most positions are NaN -> the mean-vs-mean0 choice downstream matters).
bigWigInfo "$OUTPUT"

[[ "$SORTED" != "$BEDGRAPH" ]] && rm -f "$SORTED"
echo "=== Done: $OUTPUT ==="
