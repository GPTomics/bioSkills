#!/bin/bash
# Joint genotyping from GVCFs

REF=$1
INTERVAL_LIST=$2
OUTPUT=$3
shift 3
GVCFS="$@"

if [ -z "$REF" ] || [ -z "$INTERVAL_LIST" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: $0 <reference> <intervals> <output_prefix> <gvcf1> [gvcf2] ..."
    exit 1
fi

echo "sample\tgvcf" > sample_map.txt
i=1
for gvcf in $GVCFS; do
    name=$(basename $gvcf .g.vcf.gz)
    echo "${name}\t${gvcf}" >> sample_map.txt
    i=$((i+1))
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
