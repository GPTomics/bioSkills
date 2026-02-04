#!/bin/bash
# Differential splicing analysis with rMATS-turbo
# Compares splicing between two conditions

set -e

# Configuration
GTF="annotation.gtf"
READ_LENGTH=150
THREADS=8
OUTPUT_DIR="rmats_output"
TMP_DIR="rmats_tmp"

# Create sample list files
# Each file contains comma-separated BAM paths for one condition
# Example: /path/to/sample1.bam,/path/to/sample2.bam,/path/to/sample3.bam
CONDITION1_BAMS="condition1_bams.txt"
CONDITION2_BAMS="condition2_bams.txt"

# Run rMATS-turbo
# Performs both quantification and differential testing
rmats.py \
    --b1 "$CONDITION1_BAMS" \
    --b2 "$CONDITION2_BAMS" \
    --gtf "$GTF" \
    -t paired \
    --readLength "$READ_LENGTH" \
    --nthread "$THREADS" \
    --od "$OUTPUT_DIR" \
    --tmp "$TMP_DIR" \
    --cstat 0.01  # Chi-square cutoff for significance

# Output files:
# SE.MATS.JC.txt - Skipped exons (junction counts only)
# SE.MATS.JCEC.txt - Skipped exons (junction + exon body counts)
# A5SS, A3SS, MXE, RI - Alternative splice sites, mutually exclusive, retained introns

echo "rMATS analysis complete. Results in $OUTPUT_DIR"
echo ""
echo "Key output columns:"
echo "  - IncLevelDifference: deltaPSI (positive = higher in condition1)"
echo "  - FDR: Benjamini-Hochberg corrected p-value"
echo "  - IJC/SJC: Inclusion/Skipping junction counts"

# Filter significant events (example)
echo ""
echo "Filtering significant SE events (|deltaPSI| > 0.1, FDR < 0.05):"
awk -F'\t' 'NR==1 || ($20 < 0.05 && ($23 > 0.1 || $23 < -0.1))' \
    "$OUTPUT_DIR/SE.MATS.JC.txt" | head -20
