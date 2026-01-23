#!/bin/bash
# Complete TCR repertoire analysis pipeline

R1=$1
R2=$2
SAMPLE_NAME=$3
OUTPUT_DIR=${4:-"tcr_results"}
SPECIES=${5:-"hsa"}  # hsa for human, mmu for mouse

mkdir -p ${OUTPUT_DIR}/{mixcr,vdjtools,plots}

echo "=== Step 1: MiXCR Alignment ==="
# Align reads to V(D)J reference
# -s: species (hsa=human, mmu=mouse)
# -p: protocol preset (rna-seq for bulk RNA-seq derived)
mixcr align \
    -s $SPECIES \
    -p rna-seq \
    --report ${OUTPUT_DIR}/mixcr/align_report.txt \
    $R1 $R2 \
    ${OUTPUT_DIR}/mixcr/${SAMPLE_NAME}.vdjca

# QC check: alignment rate
grep "Successfully aligned" ${OUTPUT_DIR}/mixcr/align_report.txt

echo "=== Step 2: MiXCR Assembly ==="
# Assemble clonotypes from aligned reads
mixcr assemble \
    --report ${OUTPUT_DIR}/mixcr/assemble_report.txt \
    ${OUTPUT_DIR}/mixcr/${SAMPLE_NAME}.vdjca \
    ${OUTPUT_DIR}/mixcr/${SAMPLE_NAME}.clns

# QC check: clonotype count
grep "Clonotypes" ${OUTPUT_DIR}/mixcr/assemble_report.txt

echo "=== Step 3: Export Clonotypes ==="
mixcr exportClones \
    ${OUTPUT_DIR}/mixcr/${SAMPLE_NAME}.clns \
    ${OUTPUT_DIR}/mixcr/${SAMPLE_NAME}_clones.txt

echo "=== Step 4: VDJtools Analysis ==="
# Convert MiXCR output to VDJtools format
vdjtools Convert \
    -S mixcr \
    ${OUTPUT_DIR}/mixcr/${SAMPLE_NAME}_clones.txt \
    ${OUTPUT_DIR}/vdjtools/

# Calculate diversity metrics
# Outputs: Shannon entropy, Simpson diversity, Chao1, etc.
vdjtools CalcDiversityStats \
    ${OUTPUT_DIR}/vdjtools/${SAMPLE_NAME}_clones.txt \
    ${OUTPUT_DIR}/vdjtools/diversity

echo "=== Step 5: Visualization ==="
# V-J usage heatmap
vdjtools PlotFancyVJUsage \
    ${OUTPUT_DIR}/vdjtools/${SAMPLE_NAME}_clones.txt \
    ${OUTPUT_DIR}/plots/vj_usage

# CDR3 spectratype
vdjtools PlotFancySpectratype \
    ${OUTPUT_DIR}/vdjtools/${SAMPLE_NAME}_clones.txt \
    ${OUTPUT_DIR}/plots/spectratype

echo "=== Pipeline Complete ==="
echo "Results in: $OUTPUT_DIR"
echo "Key outputs:"
echo "  - Clonotypes: ${OUTPUT_DIR}/mixcr/${SAMPLE_NAME}_clones.txt"
echo "  - Diversity: ${OUTPUT_DIR}/vdjtools/diversity.txt"
echo "  - Plots: ${OUTPUT_DIR}/plots/"
