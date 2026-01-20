#!/bin/bash
# Single-sample GATK variant calling pipeline

SAMPLE=$1
REF=$2

if [ -z "$SAMPLE" ] || [ -z "$REF" ]; then
    echo "Usage: $0 <sample_bam> <reference_fasta>"
    exit 1
fi

NAME=$(basename $SAMPLE .bam)

echo "Calling variants for $NAME..."

gatk HaplotypeCaller \
    -R $REF \
    -I $SAMPLE \
    -O ${NAME}.g.vcf.gz \
    -ERC GVCF

gatk GenotypeGVCFs \
    -R $REF \
    -V ${NAME}.g.vcf.gz \
    -O ${NAME}.vcf.gz

# GATK hard filter thresholds per GATK Best Practices documentation
gatk VariantFiltration \
    -R $REF \
    -V ${NAME}.vcf.gz \
    -O ${NAME}.filtered.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "HighFS" \
    --filter-expression "MQ < 40.0" --filter-name "LowMQ"
# QD<2.0: Quality by depth; low values indicate false positives
# FS>60.0: Fisher strand bias; high values suggest strand artifacts
# MQ<40.0: Mapping quality; low values indicate ambiguous alignments

echo "Done: ${NAME}.filtered.vcf.gz"
