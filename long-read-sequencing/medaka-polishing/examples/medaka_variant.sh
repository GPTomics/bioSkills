#!/bin/bash

READS="reads.fastq.gz"
REFERENCE="reference.fa"
OUTPUT_DIR="medaka_variants"
THREADS=4
MODEL="r1041_e82_400bps_sup_v5.0.0"

# medaka v2.0+ uses medaka_variant for haploid samples
medaka_variant \
    -i $READS \
    -r $REFERENCE \
    -o $OUTPUT_DIR \
    -m $MODEL \
    -t $THREADS

if [ -f "${OUTPUT_DIR}/medaka.annotated.vcf" ]; then
    echo "Variant calling complete!"

    echo ""
    echo "Variant summary:"
    bcftools stats ${OUTPUT_DIR}/medaka.annotated.vcf | grep "^SN"

    bcftools filter -i 'QUAL>20' ${OUTPUT_DIR}/medaka.annotated.vcf \
        > ${OUTPUT_DIR}/medaka.filtered.vcf

    echo ""
    echo "Filtered variants (QUAL>20):"
    bcftools stats ${OUTPUT_DIR}/medaka.filtered.vcf | grep "^SN"
else
    echo "Error: Variant calling failed"
    exit 1
fi
