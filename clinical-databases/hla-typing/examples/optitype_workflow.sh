#!/bin/bash
# OptiType HLA Class I typing from WES/WGS
#
# Prerequisites:
#   conda install -c bioconda optitype razers3 samtools
#
# OptiType types HLA-A, HLA-B, HLA-C at 4-field resolution
# For Class II (HLA-DR, DQ, DP), use HLA-HD or arcasHLA

set -e

# Input BAM (WES or WGS)
INPUT_BAM="${1:-input.bam}"
SAMPLE_NAME="${2:-sample}"
OUTPUT_DIR="${3:-optitype_output}"
THREADS=4

# HLA region on chromosome 6
HLA_REGION="chr6:28000000-34000000"

echo "=== OptiType HLA Typing ==="
echo "Input: $INPUT_BAM"
echo "Sample: $SAMPLE_NAME"
echo "Output: $OUTPUT_DIR"
echo ""

# Check input exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM not found: $INPUT_BAM"
    echo ""
    echo "Usage: $0 input.bam sample_name output_dir"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Step 1: Extract HLA region reads
# This reduces input size significantly and speeds up typing
echo "Extracting HLA region reads..."
samtools view -h "$INPUT_BAM" "$HLA_REGION" | \
    samtools fastq -@ "$THREADS" \
    -1 "${OUTPUT_DIR}/${SAMPLE_NAME}_hla_R1.fq.gz" \
    -2 "${OUTPUT_DIR}/${SAMPLE_NAME}_hla_R2.fq.gz" \
    -

echo "Extracted reads: ${OUTPUT_DIR}/${SAMPLE_NAME}_hla_R*.fq.gz"

# Step 2: Create OptiType config if not exists
CONFIG_FILE="${OUTPUT_DIR}/config.ini"
if [ ! -f "$CONFIG_FILE" ]; then
    cat > "$CONFIG_FILE" << 'EOF'
[mapping]
razers3=razers3
threads=4

[ilp]
solver=glpk
threads=4

[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
EOF
fi

# Step 3: Run OptiType
# -d for DNA mode (use -r for RNA-seq)
echo ""
echo "Running OptiType (DNA mode)..."
OptiTypePipeline.py \
    -i "${OUTPUT_DIR}/${SAMPLE_NAME}_hla_R1.fq.gz" \
       "${OUTPUT_DIR}/${SAMPLE_NAME}_hla_R2.fq.gz" \
    -d \
    -o "${OUTPUT_DIR}/${SAMPLE_NAME}" \
    -c "$CONFIG_FILE"

# Step 4: Parse and display results
RESULT_FILE=$(find "${OUTPUT_DIR}/${SAMPLE_NAME}" -name "*_result.tsv" | head -1)

if [ -f "$RESULT_FILE" ]; then
    echo ""
    echo "=== HLA Typing Results ==="
    cat "$RESULT_FILE"
    echo ""
    echo "Results saved to: $RESULT_FILE"

    # Parse for clinical report format
    echo ""
    echo "=== Clinical Report Format ==="
    awk -F'\t' 'NR==2 {
        printf "HLA-A: %s / %s\n", $2, $3
        printf "HLA-B: %s / %s\n", $4, $5
        printf "HLA-C: %s / %s\n", $6, $7
    }' "$RESULT_FILE"
else
    echo "Error: Results file not found"
    exit 1
fi

# Step 5: Check pharmacogenomic alleles
echo ""
echo "=== Pharmacogenomic Screening ==="

# Extract alleles for checking
ALLELES=$(awk -F'\t' 'NR==2 {print $2, $3, $4, $5, $6, $7}' "$RESULT_FILE")

# Check for high-risk alleles
# B*57:01: Abacavir hypersensitivity
if echo "$ALLELES" | grep -q "B\*57:01"; then
    echo "WARNING: HLA-B*57:01 detected - Abacavir contraindicated"
fi

# B*15:02: Carbamazepine SJS/TEN (high risk in Asian populations)
if echo "$ALLELES" | grep -q "B\*15:02"; then
    echo "WARNING: HLA-B*15:02 detected - Carbamazepine SJS/TEN risk"
fi

# B*58:01: Allopurinol SJS/TEN
if echo "$ALLELES" | grep -q "B\*58:01"; then
    echo "WARNING: HLA-B*58:01 detected - Allopurinol SJS/TEN risk"
fi

# A*31:01: Carbamazepine DRESS
if echo "$ALLELES" | grep -q "A\*31:01"; then
    echo "WARNING: HLA-A*31:01 detected - Carbamazepine DRESS risk"
fi

echo ""
echo "HLA typing complete."
