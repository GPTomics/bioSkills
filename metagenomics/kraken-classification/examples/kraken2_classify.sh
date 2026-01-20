#!/bin/bash

KRAKEN_DB="/path/to/kraken2_db"
READS_R1="sample_R1.fastq.gz"
READS_R2="sample_R2.fastq.gz"
OUTPUT_DIR="kraken_output"
SAMPLE="sample"

mkdir -p $OUTPUT_DIR

kraken2 --db $KRAKEN_DB \
    --threads 8 \
    --paired \
    --gzip-compressed \
    --output ${OUTPUT_DIR}/${SAMPLE}.kraken \
    --report ${OUTPUT_DIR}/${SAMPLE}_report.txt \
    --use-names \  # outputs human-readable taxon names instead of numeric NCBI IDs
    $READS_R1 $READS_R2

echo "Classification summary:"
head -20 ${OUTPUT_DIR}/${SAMPLE}_report.txt

echo ""
echo "Top 10 species:"
awk '$4 == "S"' ${OUTPUT_DIR}/${SAMPLE}_report.txt | sort -k1 -nr | head -10

classified=$(grep -c "^C" ${OUTPUT_DIR}/${SAMPLE}.kraken || echo "0")
unclassified=$(grep -c "^U" ${OUTPUT_DIR}/${SAMPLE}.kraken || echo "0")
total=$((classified + unclassified))
pct=$(echo "scale=2; $classified * 100 / $total" | bc)
echo ""
# Classification rate >70% typical for well-characterized samples; lower rates suggest novel taxa or contamination
echo "Classification rate: ${pct}% (${classified}/${total})"
