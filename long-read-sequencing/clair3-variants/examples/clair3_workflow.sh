#!/bin/bash
# Clair3 variant calling workflow for long reads

set -euo pipefail

SAMPLE="sample"
BAM="sample.bam"
REFERENCE="reference.fasta"
PLATFORM="ont"  # or "hifi"
THREADS=32

# Model path (adjust for your installation)
MODEL_PATH="${CONDA_PREFIX}/bin/models/${PLATFORM}"

echo "=== Step 1: Input validation ==="
samtools quickcheck ${BAM} && echo "BAM OK" || { echo "BAM invalid"; exit 1; }
samtools view -H ${BAM} | head -5

echo "=== Step 2: Run Clair3 ==="
run_clair3.sh \
    --bam_fn=${BAM} \
    --ref_fn=${REFERENCE} \
    --threads=${THREADS} \
    --platform=${PLATFORM} \
    --model_path=${MODEL_PATH} \
    --output=${SAMPLE}_clair3

VCF="${SAMPLE}_clair3/merge_output.vcf.gz"

echo "=== Step 3: Index VCF ==="
bcftools index -t ${VCF}

echo "=== Step 4: Variant statistics ==="
echo "Total variants:"
bcftools stats ${VCF} | grep "^SN"

echo "=== Step 5: Filter high-quality variants ==="
bcftools view -i 'QUAL>20 && GQ>30' ${VCF} -Oz -o ${SAMPLE}_hq.vcf.gz
bcftools index -t ${SAMPLE}_hq.vcf.gz

echo "High-quality variants:"
bcftools stats ${SAMPLE}_hq.vcf.gz | grep "^SN"

echo "=== Step 6: Separate SNPs and indels ==="
bcftools view -v snps ${SAMPLE}_hq.vcf.gz -Oz -o ${SAMPLE}_snps.vcf.gz
bcftools view -v indels ${SAMPLE}_hq.vcf.gz -Oz -o ${SAMPLE}_indels.vcf.gz

bcftools index -t ${SAMPLE}_snps.vcf.gz
bcftools index -t ${SAMPLE}_indels.vcf.gz

echo "SNPs: $(bcftools view -H ${SAMPLE}_snps.vcf.gz | wc -l)"
echo "Indels: $(bcftools view -H ${SAMPLE}_indels.vcf.gz | wc -l)"

echo "=== Workflow complete ==="
echo "Raw VCF: ${VCF}"
echo "High-quality VCF: ${SAMPLE}_hq.vcf.gz"
echo "SNPs: ${SAMPLE}_snps.vcf.gz"
echo "Indels: ${SAMPLE}_indels.vcf.gz"
