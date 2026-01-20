#!/bin/bash
# GATK joint calling workflow for multiple samples

set -euo pipefail

REFERENCE=${1:-reference.fa}
OUTPUT_DIR=${2:-joint_calling_output}
INTERVALS=${3:-intervals.bed}

mkdir -p ${OUTPUT_DIR}

echo "=== Step 1: Generate per-sample GVCFs ==="
for bam in *.bam; do
    sample=$(basename $bam .bam)
    echo "Processing $sample..."

    gatk HaplotypeCaller \
        -R ${REFERENCE} \
        -I ${bam} \
        -O ${OUTPUT_DIR}/${sample}.g.vcf.gz \
        -ERC GVCF \
        -L ${INTERVALS}
done

echo "=== Step 2: Combine GVCFs ==="
GVCF_ARGS=""
for gvcf in ${OUTPUT_DIR}/*.g.vcf.gz; do
    GVCF_ARGS="${GVCF_ARGS} -V ${gvcf}"
done

gatk CombineGVCFs \
    -R ${REFERENCE} \
    ${GVCF_ARGS} \
    -O ${OUTPUT_DIR}/combined.g.vcf.gz

echo "=== Step 3: Joint Genotyping ==="
gatk GenotypeGVCFs \
    -R ${REFERENCE} \
    -V ${OUTPUT_DIR}/combined.g.vcf.gz \
    -O ${OUTPUT_DIR}/cohort.vcf.gz

echo "=== Step 4: Basic Filtering ==="
gatk VariantFiltration \
    -R ${REFERENCE} \
    -V ${OUTPUT_DIR}/cohort.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "HighFS" \
    --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
    -O ${OUTPUT_DIR}/cohort_filtered.vcf.gz

# Index and stats
bcftools index -t ${OUTPUT_DIR}/cohort_filtered.vcf.gz
bcftools stats ${OUTPUT_DIR}/cohort_filtered.vcf.gz > ${OUTPUT_DIR}/cohort_stats.txt

echo "=== Complete ==="
echo "Output: ${OUTPUT_DIR}/cohort_filtered.vcf.gz"
