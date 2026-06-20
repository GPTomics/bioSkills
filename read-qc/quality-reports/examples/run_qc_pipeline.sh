#!/bin/bash
# Reference: FastQC 0.12+, MultiQC 1.21+, falco 1.2+ | Verify API if version differs
# Per-file QC then cross-sample aggregation. FastQC fails are calibrated to WGS DNA;
# read the resulting MultiQC plots against the assay (see SKILL.md), do not trust the
# traffic light alone. For long reads use NanoPlot instead of FastQC.
# Usage: ./run_qc_pipeline.sh <input_dir> <output_dir> <threads>

set -euo pipefail

INPUT_DIR="${1:-.}"
OUTPUT_DIR="${2:-qc_results}"
THREADS="${3:-4}"

# falco is a drop-in FastQC replacement (~3x faster); prefer it on large cohorts
QC_TOOL=fastqc
command -v falco >/dev/null 2>&1 && QC_TOOL=falco

mkdir -p "$OUTPUT_DIR/perfile" "$OUTPUT_DIR/multiqc"

echo "=== Per-file QC ($QC_TOOL) ==="
"$QC_TOOL" -t "$THREADS" -o "$OUTPUT_DIR/perfile" "$INPUT_DIR"/*.fastq.gz

echo "=== Aggregating with MultiQC ==="
multiqc "$OUTPUT_DIR/perfile" -o "$OUTPUT_DIR/multiqc" -f

# Count input files (falco and fastqc name their HTML differently, so count the inputs)
N_SAMPLES=$(ls "$INPUT_DIR"/*.fastq.gz 2>/dev/null | wc -l)
echo "Report: $OUTPUT_DIR/multiqc/multiqc_report.html  (inputs: $N_SAMPLES)"
echo "Read the General Statistics table first; investigate outliers, not absolute pass/fail."
