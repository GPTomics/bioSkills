#!/bin/bash

BAM="aligned.bam"
REFERENCE="reference.fa"
OUTPUT_DIR="sv_calls"
THREADS=8

mkdir -p $OUTPUT_DIR

echo "Calling SVs with Sniffles2..."
sniffles --input $BAM \
    --vcf ${OUTPUT_DIR}/sniffles.vcf \
    --reference $REFERENCE \
    --threads $THREADS \
    --minsupport 3 \
    --minsvlen 50 \
    --output-rnames

echo ""
echo "Calling SVs with cuteSV..."
mkdir -p ${OUTPUT_DIR}/cutesv_work
cuteSV $BAM $REFERENCE ${OUTPUT_DIR}/cutesv.vcf ${OUTPUT_DIR}/cutesv_work/ \
    --threads $THREADS \
    --min_support 3 \
    --min_size 50 \
    --genotype \
    --max_cluster_bias_INS 100 \
    --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 100 \
    --diff_ratio_merging_DEL 0.3

echo ""
echo "Sniffles2 results:"
bcftools stats ${OUTPUT_DIR}/sniffles.vcf | grep "^SN"
echo ""
for svtype in DEL INS INV DUP BND; do
    count=$(grep -c "SVTYPE=${svtype}" ${OUTPUT_DIR}/sniffles.vcf || echo "0")
    echo "${svtype}: ${count}"
done

echo ""
echo "cuteSV results:"
bcftools stats ${OUTPUT_DIR}/cutesv.vcf | grep "^SN"
echo ""
for svtype in DEL INS INV DUP BND; do
    count=$(grep -c "SVTYPE=${svtype}" ${OUTPUT_DIR}/cutesv.vcf || echo "0")
    echo "${svtype}: ${count}"
done

bcftools filter -i 'QUAL>=20' ${OUTPUT_DIR}/sniffles.vcf > ${OUTPUT_DIR}/sniffles.filtered.vcf
bcftools filter -i 'QUAL>=20' ${OUTPUT_DIR}/cutesv.vcf > ${OUTPUT_DIR}/cutesv.filtered.vcf

echo ""
echo "Results saved to ${OUTPUT_DIR}/"
