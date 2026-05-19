#!/bin/bash
# Reference: OptiType 1.3.5+, samtools 1.19+, bwa-mem 0.7.17+ | Verify CLI flags if version differs
#
# OptiType class I HLA typing from WES/WGS (Szolek 2014 Bioinformatics).
# OptiType is the TCGA convention and the class-I 4-digit benchmark anchor.
# For class II + KIR use T1K (t1k_workflow.sh) or HLA-LA.

set -euo pipefail

INPUT_BAM="${1:-input.bam}"
SAMPLE_NAME="${2:-sample}"
OUTPUT_DIR="${3:-optitype_output}"
THREADS="${THREADS:-4}"

HLA_REGION="chr6:28000000-34000000"
HLA_ALT_CONTIGS=(
    chr6_GL000250v2_alt chr6_GL000251v2_alt chr6_GL000252v2_alt chr6_GL000253v2_alt
    chr6_GL000254v2_alt chr6_GL000255v2_alt chr6_GL000256v2_alt
)

echo "=== OptiType class I HLA typing ==="
mkdir -p "$OUTPUT_DIR"

# Alt-aware extraction (critical for GRCh38 accuracy)
echo "Extracting MHC reads (region + alt contigs)..."
samtools view -bh -@ "$THREADS" "$INPUT_BAM" "$HLA_REGION" "${HLA_ALT_CONTIGS[@]}" \
    > "$OUTPUT_DIR/hla_region.bam"
samtools sort -n -@ "$THREADS" "$OUTPUT_DIR/hla_region.bam" -o "$OUTPUT_DIR/hla_sorted.bam"
samtools fastq -@ "$THREADS" \
    -1 "$OUTPUT_DIR/${SAMPLE_NAME}_hla_R1.fq.gz" \
    -2 "$OUTPUT_DIR/${SAMPLE_NAME}_hla_R2.fq.gz" \
    -0 /dev/null -s "$OUTPUT_DIR/${SAMPLE_NAME}_singletons.fq.gz" \
    "$OUTPUT_DIR/hla_sorted.bam"

# Config -- ILP solver = GLPK (open-source); use CPLEX for production speedup
CONFIG="$OUTPUT_DIR/config.ini"
cat > "$CONFIG" << EOF
[mapping]
razers3=razers3
threads=$THREADS

[ilp]
solver=glpk
threads=$THREADS

[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
EOF

# DNA mode (-d). Use -r for RNA-seq.
echo "Running OptiType..."
OptiTypePipeline.py \
    -i "$OUTPUT_DIR/${SAMPLE_NAME}_hla_R1.fq.gz" \
       "$OUTPUT_DIR/${SAMPLE_NAME}_hla_R2.fq.gz" \
    -d -o "$OUTPUT_DIR/${SAMPLE_NAME}" -c "$CONFIG"

RESULT_FILE=$(find "$OUTPUT_DIR/${SAMPLE_NAME}" -name "*_result.tsv" | head -1)
if [ ! -f "$RESULT_FILE" ]; then
    echo "Error: result.tsv not found"
    exit 1
fi

echo
echo "=== OptiType Class I Result ==="
cat "$RESULT_FILE"

# Format for clinical report
echo
echo "=== Clinical Report Format ==="
awk -F'\t' 'NR==2 {
    printf "HLA-A: %s / %s\n", $2, $3
    printf "HLA-B: %s / %s\n", $4, $5
    printf "HLA-C: %s / %s\n", $6, $7
}' "$RESULT_FILE"

# PGx screen -- 4-field specificity required
echo
echo "=== Pharmacogenomic Screen (4-field) ==="
ALLELES=$(awk -F'\t' 'NR==2 {print $2, $3, $4, $5, $6, $7}' "$RESULT_FILE")

echo "$ALLELES" | grep -q "B\*57:01" && \
    echo "WARNING: HLA-B*57:01 detected -- Abacavir contraindicated (Mallal 2008 NEJM; OR ~100)"
echo "$ALLELES" | grep -q "B\*15:02" && \
    echo "WARNING: HLA-B*15:02 detected -- Carbamazepine SJS/TEN risk (Chung 2004 Nature; OR ~2500)"
echo "$ALLELES" | grep -q "B\*58:01" && \
    echo "WARNING: HLA-B*58:01 detected -- Allopurinol SJS/TEN risk (Hung 2005 PNAS; OR ~580)"
echo "$ALLELES" | grep -q "A\*31:01" && \
    echo "WARNING: HLA-A*31:01 detected -- Carbamazepine DRESS risk (McCormack 2011 NEJM)"
echo "$ALLELES" | grep -q "B\*13:01" && \
    echo "WARNING: HLA-B*13:01 detected -- Dapsone hypersensitivity (Zhang 2013 NEJM)"

echo
echo "Note: OptiType covers class I only. For class II (B*15:02 + DPB1/DRB1 transplant matching),"
echo "      run t1k_workflow.sh or HLA-LA."
