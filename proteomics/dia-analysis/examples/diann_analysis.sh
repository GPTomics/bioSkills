#!/bin/bash
# DIA-NN library-free proteomics analysis

set -e

# Configuration
FASTA="uniprot_human_reviewed.fasta"
OUTPUT_DIR="diann_results"
THREADS=8

mkdir -p "$OUTPUT_DIR"

# Find all mzML files
MZML_FILES=$(ls *.mzML | tr '\n' ' ')

echo "Running DIA-NN library-free analysis..."
diann \
    --f $MZML_FILES \
    --lib "" \
    --threads $THREADS \
    --verbose 1 \
    --out "$OUTPUT_DIR/report.tsv" \
    --qvalue 0.01 \
    --matrices \
    --out-lib "$OUTPUT_DIR/generated_lib.tsv" \
    --gen-spec-lib \
    --predictor \
    --fasta "$FASTA" \
    --fasta-search \
    --min-fr-mz 200 \
    --max-fr-mz 1800 \
    --met-excision \
    --cut K*,R* \
    --missed-cleavages 1 \
    --min-pep-len 7 \
    --max-pep-len 30 \
    --min-pr-mz 300 \
    --max-pr-mz 1800 \
    --min-pr-charge 1 \
    --max-pr-charge 4 \
    --unimod4 \
    --var-mods 1 \
    --var-mod UniMod:35,15.994915,M \
    --reanalyse \
    --smart-profiling

echo "Analysis complete!"
echo "Protein matrix: $OUTPUT_DIR/report.pg_matrix.tsv"
echo "Statistics: $OUTPUT_DIR/report.stats.tsv"
