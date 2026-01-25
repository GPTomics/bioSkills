#!/bin/bash
# MultiQC Report Generation Examples
# Demonstrates various QC report configurations

# =============================================================================
# Basic Usage
# =============================================================================

# Generate report from all QC outputs in directory
multiqc results/ -o qc_report/

# Specify custom report name
multiqc results/ -n "project_xyz_qc" -o qc_report/

# =============================================================================
# Selective Module Inclusion
# =============================================================================

# RNA-seq specific modules
multiqc results/ \
    --module fastqc \
    --module star \
    --module featurecounts \
    --module salmon \
    -o qc_report/

# Variant calling modules
multiqc results/ \
    --module fastqc \
    --module bwa \
    --module samtools \
    --module bcftools \
    -o qc_report/

# =============================================================================
# Configuration File
# =============================================================================

# Create configuration file for reproducible reports
cat > multiqc_config.yaml << 'EOF'
title: "RNA-seq QC Report"
subtitle: "Project XYZ - Batch 1"
intro_text: "Quality control metrics for all samples in batch 1."

# Report options
report_comment: "Generated automatically by pipeline"
show_analysis_paths: False
show_analysis_time: True

# Sample name cleaning
extra_fn_clean_exts:
  - '.sorted'
  - '.dedup'
  - '.trimmed'
  - '_R1'
  - '_R2'

# Module order (top to bottom in report)
module_order:
  - fastqc:
      name: "Read Quality (FastQC)"
  - star:
      name: "Alignment (STAR)"
  - featurecounts:
      name: "Quantification"

# Highlight problematic samples
table_cond_formatting_rules:
  percent_gc:
    warn: [{lt: 35}, {gt: 65}]
  percent_duplicates:
    warn: [{gt: 50}]
  percent_mapped:
    fail: [{lt: 50}]
    warn: [{lt: 70}]

# Custom color scheme
custom_plot_config:
  percent_mapped:
    min: 0
    max: 100
EOF

# Run with config
multiqc results/ -c multiqc_config.yaml -o qc_report/

# =============================================================================
# Comparing Multiple Runs
# =============================================================================

# Compare pre and post-trimming FastQC
multiqc \
    raw_fastqc/ \
    trimmed_fastqc/ \
    -o comparison_report/ \
    -n "trim_comparison"

# =============================================================================
# Export Data
# =============================================================================

# Export parsed data as TSV
multiqc results/ -o qc_report/ --export

# Flat data format (easier to parse)
multiqc results/ -o qc_report/ --data-format tsv --flat

# =============================================================================
# Integration with Workflows
# =============================================================================

# Function for Snakemake/Nextflow integration
generate_qc_report() {
    local results_dir=$1
    local output_dir=$2
    local project_name=$3

    multiqc "$results_dir" \
        -o "$output_dir" \
        -n "$project_name" \
        -c multiqc_config.yaml \
        --force \
        --quiet

    echo "QC report generated: $output_dir/$project_name.html"
}

# Example call
# generate_qc_report "results/" "qc_reports/" "rnaseq_batch1"

echo "MultiQC examples complete"
