#!/bin/bash
# Complete Mutect2 somatic variant calling pipeline

set -euo pipefail

TUMOR_BAM=${1:-tumor.bam}
NORMAL_BAM=${2:-normal.bam}
REFERENCE=${3:-reference.fa}
OUTPUT_PREFIX=${4:-somatic}
GNOMAD=${5:-gnomad.vcf.gz}
PON=${6:-pon.vcf.gz}

echo "=== Somatic Variant Calling Pipeline ==="
echo "Tumor: $TUMOR_BAM"
echo "Normal: $NORMAL_BAM"

# Get normal sample name
NORMAL_NAME=$(samtools view -H $NORMAL_BAM | grep '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/' | head -1)

echo "Normal sample name: $NORMAL_NAME"
echo ""

# Step 1: Call variants
echo "=== Step 1: Calling variants with Mutect2 ==="
gatk Mutect2 \
    -R $REFERENCE \
    -I $TUMOR_BAM \
    -I $NORMAL_BAM \
    -normal $NORMAL_NAME \
    --germline-resource $GNOMAD \
    --panel-of-normals $PON \
    --f1r2-tar-gz ${OUTPUT_PREFIX}_f1r2.tar.gz \
    -O ${OUTPUT_PREFIX}_raw.vcf.gz

# Step 2: Learn read orientation model
echo "=== Step 2: Learning read orientation model ==="
gatk LearnReadOrientationModel \
    -I ${OUTPUT_PREFIX}_f1r2.tar.gz \
    -O ${OUTPUT_PREFIX}_read_orientation.tar.gz

# Step 3: Calculate contamination
echo "=== Step 3: Calculating contamination ==="
gatk GetPileupSummaries \
    -I $TUMOR_BAM \
    -V $GNOMAD \
    -L $GNOMAD \
    -O ${OUTPUT_PREFIX}_tumor_pileups.table

gatk GetPileupSummaries \
    -I $NORMAL_BAM \
    -V $GNOMAD \
    -L $GNOMAD \
    -O ${OUTPUT_PREFIX}_normal_pileups.table

gatk CalculateContamination \
    -I ${OUTPUT_PREFIX}_tumor_pileups.table \
    -matched ${OUTPUT_PREFIX}_normal_pileups.table \
    -O ${OUTPUT_PREFIX}_contamination.table \
    --tumor-segmentation ${OUTPUT_PREFIX}_segments.table

# Step 4: Filter variants
echo "=== Step 4: Filtering variants ==="
gatk FilterMutectCalls \
    -R $REFERENCE \
    -V ${OUTPUT_PREFIX}_raw.vcf.gz \
    --contamination-table ${OUTPUT_PREFIX}_contamination.table \
    --tumor-segmentation ${OUTPUT_PREFIX}_segments.table \
    --ob-priors ${OUTPUT_PREFIX}_read_orientation.tar.gz \
    -O ${OUTPUT_PREFIX}_filtered.vcf.gz

# Extract PASS variants
bcftools view -f PASS ${OUTPUT_PREFIX}_filtered.vcf.gz -Oz -o ${OUTPUT_PREFIX}_pass.vcf.gz
bcftools index -t ${OUTPUT_PREFIX}_pass.vcf.gz

# Statistics
echo ""
echo "=== Results Summary ==="
bcftools stats ${OUTPUT_PREFIX}_pass.vcf.gz | grep -E "^SN"

echo ""
echo "Output files:"
echo "  Raw calls: ${OUTPUT_PREFIX}_raw.vcf.gz"
echo "  Filtered: ${OUTPUT_PREFIX}_filtered.vcf.gz"
echo "  PASS only: ${OUTPUT_PREFIX}_pass.vcf.gz"
