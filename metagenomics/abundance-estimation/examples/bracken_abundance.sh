#!/bin/bash

KRAKEN_DB="/path/to/kraken2_db"
REPORTS_DIR="kraken_reports"
OUTPUT_DIR="bracken_output"
READ_LENGTH=150

mkdir -p $OUTPUT_DIR

for report in ${REPORTS_DIR}/*_report.txt; do
    sample=$(basename $report _report.txt)
    echo "Processing ${sample}..."

    bracken -d $KRAKEN_DB \
        -i $report \
        -o ${OUTPUT_DIR}/${sample}_species.txt \
        -w ${OUTPUT_DIR}/${sample}_bracken_report.txt \
        -r $READ_LENGTH \
        -l S \
        -t 10
done

echo "Combining samples..."
combine_bracken_outputs.py \
    --files ${OUTPUT_DIR}/*_species.txt \
    -o ${OUTPUT_DIR}/combined_species_abundance.txt

echo ""
echo "Summary of top species:"
head -20 ${OUTPUT_DIR}/combined_species_abundance.txt

echo ""
echo "Results saved to ${OUTPUT_DIR}/"
