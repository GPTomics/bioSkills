#!/bin/bash
# miRDeep2 workflow for novel miRNA discovery and quantification

COLLAPSED_READS=$1
GENOME_FA=$2
SPECIES=${3:-"Human"}

# miRBase files (download if not present)
MIRBASE_DIR="mirbase"
mkdir -p $MIRBASE_DIR

if [ ! -f "$MIRBASE_DIR/mature.fa" ]; then
    echo "Downloading miRBase..."
    wget -q -O $MIRBASE_DIR/mature.fa https://www.mirbase.org/download/mature.fa
    wget -q -O $MIRBASE_DIR/hairpin.fa https://www.mirbase.org/download/hairpin.fa
fi

# Extract species-specific sequences
# hsa = human, mmu = mouse, rno = rat
case $SPECIES in
    Human) PREFIX="hsa" ;;
    Mouse) PREFIX="mmu" ;;
    Rat)   PREFIX="rno" ;;
    *)     PREFIX="hsa" ;;
esac

grep -A1 "^>$PREFIX-" $MIRBASE_DIR/mature.fa | grep -v "^--$" > ${PREFIX}_mature.fa
grep -A1 "^>$PREFIX-" $MIRBASE_DIR/hairpin.fa | grep -v "^--$" > ${PREFIX}_hairpin.fa

# Get mature miRNAs from other species for conservation scoring
grep -v "^>$PREFIX-" $MIRBASE_DIR/mature.fa | grep "^>" -A1 | grep -v "^--$" > other_mature.fa

echo "Using $PREFIX miRBase references"

# Build bowtie index if needed
GENOME_INDEX="${GENOME_FA%.fa}_index"
if [ ! -f "${GENOME_INDEX}.1.ebwt" ]; then
    echo "Building bowtie index..."
    bowtie-build $GENOME_FA $GENOME_INDEX
fi

# Step 1: Map reads (if not already mapped)
ARF_FILE="${COLLAPSED_READS%.fa}_vs_genome.arf"
if [ ! -f "$ARF_FILE" ]; then
    echo "Mapping reads to genome..."
    mapper.pl $COLLAPSED_READS \
        -c \
        -p $GENOME_INDEX \
        -t $ARF_FILE \
        -o 4 \
        2> mapper.log
fi

# Step 2: Run miRDeep2 discovery
# This is the main analysis step - discovers novel miRNAs and quantifies known ones
echo "Running miRDeep2..."
miRDeep2.pl \
    $COLLAPSED_READS \
    $GENOME_FA \
    $ARF_FILE \
    ${PREFIX}_mature.fa \
    other_mature.fa \
    ${PREFIX}_hairpin.fa \
    -t $SPECIES \
    2> mirdeep2_report.log

# Results are in result_*.html and result_*.csv
echo "miRDeep2 complete!"

# Parse high-confidence novel miRNAs
# Score thresholds:
# >10: High confidence - likely real novel miRNA
# 5-10: Medium confidence - may need validation
# <5: Low confidence - likely false positive
echo ""
echo "High-confidence novel miRNAs (score > 10):"
RESULT_CSV=$(ls -t result_*.csv 2>/dev/null | head -1)
if [ -f "$RESULT_CSV" ]; then
    awk -F'\t' 'NR>1 && $2>10 {print $1, "score="$2, "reads="$5}' "$RESULT_CSV" | head -20
fi

echo ""
echo "Output files:"
ls -lh result_* miRNAs_expressed* 2>/dev/null
