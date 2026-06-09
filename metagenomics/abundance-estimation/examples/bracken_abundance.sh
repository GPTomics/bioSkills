#!/bin/bash
# Reference: Bracken 2.9+, Kraken2 2.1.3+ | Verify API if version differs
# Re-estimate species abundance from Kraken2 reports and combine into a matrix.
# Bracken is step one: the output is a read-fraction composition, not a final answer.
set -euo pipefail

KRAKEN_DB="/path/to/kraken2_db"
REPORTS_DIR="kraken_reports"
OUTPUT_DIR="bracken_output"
READ_LENGTH=150          # MUST equal the bracken-build -l AND the actual post-trim read length (not auto-detected)
THRESHOLD=10             # redistribution floor: taxa with fewer than 10 clade-level reads (strict <) are dropped before redistribution

mkdir -p "$OUTPUT_DIR"

for report in "${REPORTS_DIR}"/*_report.txt; do
    sample=$(basename "$report" _report.txt)
    echo "Processing ${sample}..."
    bracken -d "$KRAKEN_DB" \
        -i "$report" \
        -o "${OUTPUT_DIR}/${sample}_species.txt" \
        -w "${OUTPUT_DIR}/${sample}_bracken_report.txt" \
        -r "$READ_LENGTH" \
        -l S \
        -t "$THRESHOLD"
done

# Taxa-by-sample matrix (read fractions). Do NOT merge with a MetaPhlAn table - different estimand.
echo "Combining samples..."
combine_bracken_outputs.py \
    --files "${OUTPUT_DIR}"/*_species.txt \
    -o "${OUTPUT_DIR}/combined_species_abundance.txt"

# Species dominated by added_reads (vs kraken_assigned_reads) may be redistribution artifacts -
# fabricated relatives of a database-absent taxon. Inspect before trusting the tail.
echo "Top species in the combined matrix:"
head -20 "${OUTPUT_DIR}/combined_species_abundance.txt"
echo "Results saved to ${OUTPUT_DIR}/"
