#!/bin/bash
# Convert bedGraph to bigWig

# Usage: ./bedgraph_to_bigwig.sh input.bedGraph chrom.sizes output.bw

BEDGRAPH="${1:-coverage.bedGraph}"
CHROM_SIZES="${2:-hg38.chrom.sizes}"
OUTPUT="${3:-output.bw}"

# Check inputs
if [[ ! -f "$BEDGRAPH" ]]; then
    echo "bedGraph file not found: $BEDGRAPH"
    echo "Usage: $0 <bedGraph> <chrom.sizes> <output.bw>"
    exit 1
fi

if [[ ! -f "$CHROM_SIZES" ]]; then
    echo "Chromosome sizes file not found: $CHROM_SIZES"
    echo ""
    echo "Create from FASTA: cut -f1,2 reference.fa.fai > chrom.sizes"
    echo "Download hg38:     wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
    exit 1
fi

# Check for bedGraphToBigWig
if ! command -v bedGraphToBigWig &> /dev/null; then
    echo "bedGraphToBigWig not found"
    echo "Install with: conda install -c bioconda ucsc-bedgraphtobigwig"
    exit 1
fi

echo "Input:        $BEDGRAPH"
echo "Chrom sizes:  $CHROM_SIZES"
echo "Output:       $OUTPUT"

# Check if already sorted
echo -e "\n=== Checking if bedGraph is sorted ==="
SORTED_CHECK=$(sort -k1,1 -k2,2n -c "$BEDGRAPH" 2>&1)
if [[ -z "$SORTED_CHECK" ]]; then
    echo "bedGraph is sorted"
    SORTED_FILE="$BEDGRAPH"
else
    echo "bedGraph needs sorting"
    SORTED_FILE="${BEDGRAPH%.bedGraph}.sorted.bedGraph"
    echo "Sorting to: $SORTED_FILE"
    sort -k1,1 -k2,2n "$BEDGRAPH" > "$SORTED_FILE"
fi

# Convert to bigWig
echo -e "\n=== Converting to bigWig ==="
bedGraphToBigWig "$SORTED_FILE" "$CHROM_SIZES" "$OUTPUT"

if [[ $? -eq 0 ]]; then
    echo "Success!"
    echo ""
    ls -lh "$OUTPUT"

    # Show summary
    echo -e "\n=== BigWig Summary ==="
    if command -v bigWigInfo &> /dev/null; then
        bigWigInfo "$OUTPUT" 2>/dev/null | head -10
    fi
else
    echo "Conversion failed"
    exit 1
fi

# Cleanup sorted file if we created it
if [[ "$SORTED_FILE" != "$BEDGRAPH" ]]; then
    echo -e "\nCleaning up sorted intermediate file"
    rm -f "$SORTED_FILE"
fi

echo -e "\n=== Done ==="
