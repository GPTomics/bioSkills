# Quality Reports - Usage Guide

## Overview

Quality reports are the first step in any NGS analysis. FastQC generates per-sample reports showing quality scores, adapter content, GC bias, and duplication levels. MultiQC aggregates multiple FastQC reports into a single interactive summary.

## When to Use This Skill

- Initial QC check of raw FASTQ files
- Comparing quality before/after trimming
- Identifying adapter contamination
- Detecting sample issues (low quality, contamination)
- Generating project-wide QC summaries

## Installation

```bash
# Conda (recommended)
conda install -c bioconda fastqc multiqc

# pip (MultiQC only)
pip install multiqc
```

## Basic Workflow

### 1. Generate FastQC Reports

```bash
# Create output directory
mkdir -p qc_reports

# Run FastQC on all samples
fastqc -t 4 -o qc_reports/ *.fastq.gz
```

### 2. Aggregate with MultiQC

```bash
# Generate summary report
multiqc qc_reports/ -o multiqc_output/
```

### 3. Review Reports

Open `multiqc_output/multiqc_report.html` in a browser.

## Key Quality Metrics

| Metric | Good | Warning | Action |
|--------|------|---------|--------|
| Per base quality | >Q30 throughout | <Q20 at ends | Trim |
| Adapter content | <1% | >5% | Trim adapters |
| Duplication | <20% | >50% | Dedup |
| GC content | Normal curve | Secondary peak | Investigate |

## Common Issues

### Low Quality at 3' End

Normal for Illumina sequencing. Use quality trimming:
```bash
fastp -i sample.fastq.gz -o trimmed.fastq.gz -q 20
```

### Adapter Contamination

Adapters visible in "Overrepresented sequences" or "Adapter Content":
```bash
cutadapt -a AGATCGGAAGAGC -o trimmed.fastq.gz sample.fastq.gz
```

### High Duplication

May indicate low library complexity or excessive PCR:
- Mark duplicates after alignment
- Consider deeper sequencing

### Unusual GC Distribution

Secondary peaks suggest contamination:
- Run FastQ Screen to identify source
- Check for adapter dimers

## MultiQC Configuration

Create `multiqc_config.yaml` for custom settings:

```yaml
title: "My Project QC Report"
report_header_info:
  - Project: "RNA-seq Analysis"
  - Sequencing: "NovaSeq 6000"

# Custom module order
module_order:
  - fastqc
  - fastp
  - star
  - picard

# Sample name cleaning
fn_clean_exts:
  - ".fastq.gz"
  - "_R1"
  - "_R2"
```

## Resources

- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC Documentation](https://multiqc.info/)
- [FastQC Example Reports](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
