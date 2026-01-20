#!/bin/bash
# CRISPResso2 editing analysis example

set -e

# Configuration
AMPLICON="AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT"
GUIDE="CTGCCCACAGGTGAGGAGGT"
OUTPUT_DIR="crispresso_results"

mkdir -p "$OUTPUT_DIR"

# Analyze control sample
echo "Analyzing control sample..."
CRISPResso \
    --fastq_r1 control_R1.fastq.gz \
    --fastq_r2 control_R2.fastq.gz \
    --amplicon_seq "$AMPLICON" \
    --guide_seq "$GUIDE" \
    --output_folder "$OUTPUT_DIR" \
    --name control \
    --quantification_window_size 10

# Analyze edited sample
echo "Analyzing edited sample..."
CRISPResso \
    --fastq_r1 edited_R1.fastq.gz \
    --fastq_r2 edited_R2.fastq.gz \
    --amplicon_seq "$AMPLICON" \
    --guide_seq "$GUIDE" \
    --output_folder "$OUTPUT_DIR" \
    --name edited \
    --quantification_window_size 10

# Compare samples
echo "Comparing samples..."
CRISPRessoCompare \
    --crispresso_output_folder_1 "$OUTPUT_DIR/CRISPResso_on_control" \
    --crispresso_output_folder_2 "$OUTPUT_DIR/CRISPResso_on_edited" \
    --output_folder "$OUTPUT_DIR/comparison"

echo "Analysis complete!"
echo "Results in: $OUTPUT_DIR"
