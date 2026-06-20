#!/bin/bash
# Reference: MultiQC 1.21+ | Verify API if version differs
# MultiQC report generation, scoping, and turning the report into a QC gate.
# MultiQC aggregates metrics other tools wrote; it does not gate. Gating is the
# separate parse-and-exit step at the bottom.

set -euo pipefail

# =============================================================================
# Basic report
# =============================================================================

multiqc results/ -o qc_report/                 # writes qc_report/multiqc_report.html + multiqc_data/
multiqc results/ -n project_xyz_qc -o qc_report/

# =============================================================================
# Scope detection (avoid phantom samples / false module matches)
# =============================================================================

# Run ONLY the modules for an RNA-seq run, and ignore scratch/work dirs.
multiqc results/ \
    --ignore "work/" --ignore "*_tmp/" \
    -m fastqc -m star -m featurecounts -m salmon \
    -o qc_report/

# =============================================================================
# Reproducible config (pin title, order, sample-name cleaning, thresholds)
# =============================================================================

cat > multiqc_config.yaml << 'EOF'
title: "RNA-seq QC Report"
subtitle: "Project XYZ - Batch 1"
intro_text: "Quality control metrics for all samples in batch 1."

show_analysis_paths: False

# Append suffixes to strip from sample names. extra_fn_clean_exts APPENDS to the
# ~100 defaults; overriding fn_clean_exts would REPLACE them and break cleaning.
extra_fn_clean_exts:
  - '.sorted'
  - '.dedup'
  - '.trimmed'

module_order:
  - fastqc:
      name: "Read Quality (FastQC)"
  - star:
      name: "Alignment (STAR)"
  - featurecounts:
      name: "Quantification"

# Conditional formatting is a CONFIGURED convenience, not a biological verdict.
# These thresholds color cells; they do not stop the pipeline.
table_cond_formatting_rules:
  percent_mapped:
    fail: [{lt: 50}]
    warn: [{lt: 70}]
  percent_duplicates:
    warn: [{gt: 50}]

# Disable AI summaries so no metrics or sample names leave the network.
ai_summary: False
EOF

multiqc results/ -c multiqc_config.yaml -o qc_report/ --force --quiet

# =============================================================================
# Pre vs post-trimming comparison (one report, two FastQC passes)
# =============================================================================

multiqc raw_fastqc/ trimmed_fastqc/ -o comparison_report/ -n trim_comparison

# =============================================================================
# QC GATE - the decision MultiQC does not make
# =============================================================================
# MultiQC always exits 0. To actually fail on bad QC, emit the machine-readable
# multiqc_data.json and apply thresholds in a separate step.

multiqc results/ -o qc_report/ --data-format json --force --quiet

gate_on_mapping_rate() {
    local data_json=$1
    local min_mapped=$2
    python3 - "$data_json" "$min_mapped" << 'PY'
import json, sys
data = json.load(open(sys.argv[1]))
min_mapped = float(sys.argv[2])
stats = data.get('report_general_stats_data', [])
failed = []
for block in stats:
    for sample, metrics in block.items():
        # general-stats keys are module-specific: STAR=uniquely_mapped_percent, samtools=reads_mapped_percent
        mapped = (metrics.get('uniquely_mapped_percent')
                  or metrics.get('mapped_percent')
                  or metrics.get('reads_mapped_percent'))
        if mapped is not None and mapped < min_mapped:
            failed.append((sample, mapped))
if failed:
    for sample, mapped in failed:
        print(f'FAIL {sample}: {mapped:.1f}% mapped (< {min_mapped}%)', file=sys.stderr)
    sys.exit(1)
print('QC gate passed')
PY
}

# gate_on_mapping_rate qc_report/multiqc_data/multiqc_data.json 70

echo "MultiQC examples complete"
