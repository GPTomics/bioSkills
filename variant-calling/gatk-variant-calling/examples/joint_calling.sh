#!/bin/bash
# Reference: GATK 4.5+, bcftools 1.19+ | Verify API if version differs
# Joint genotyping from GVCFs

REF=$1
INTERVAL_LIST=$2
OUTPUT=$3
shift 3
GVCFS=("$@")   # array, not a string: a path with a space must stay one argument

if [ -z "$REF" ] || [ -z "$INTERVAL_LIST" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: $0 <reference> <intervals> <output_prefix> <gvcf1> [gvcf2] ..."
    exit 1
fi

# GenomicsDBImport --sample-name-map expects tab-separated sampleName<TAB>path lines, NO header
: > sample_map.txt
for gvcf in "${GVCFS[@]}"; do
    name=$(basename "$gvcf" .g.vcf.gz)
    printf '%s\t%s\n' "$name" "$gvcf" >> sample_map.txt
done

echo "Importing GVCFs..."
gatk GenomicsDBImport \
    --genomicsdb-workspace-path genomicsdb_${OUTPUT} \
    --sample-name-map sample_map.txt \
    -L $INTERVAL_LIST

echo "Joint genotyping..."
gatk GenotypeGVCFs \
    -R $REF \
    -V gendb://genomicsdb_${OUTPUT} \
    -O ${OUTPUT}.vcf.gz

echo "Done: ${OUTPUT}.vcf.gz"
