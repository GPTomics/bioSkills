# fastp Workflow - Usage Guide

## Overview

fastp is a modern, all-in-one FASTQ preprocessing tool that replaces the need for separate adapter trimming (Cutadapt), quality filtering (Trimmomatic), and QC reporting (FastQC) steps. It's fast, memory-efficient, and produces comprehensive HTML reports.

## When to Use This Skill

- Starting a new NGS project (recommended default)
- Need fast, single-step preprocessing
- Processing NovaSeq/NextSeq data (poly-G handling)
- Want integrated QC reports
- Deduplication without alignment

## Installation

```bash
conda install -c bioconda fastp
```

## Basic Workflow

### Minimal Command

```bash
fastp -i R1.fq.gz -I R2.fq.gz -o R1_clean.fq.gz -O R2_clean.fq.gz
```

This automatically:
- Detects and removes adapters
- Filters low-quality reads (Q15)
- Removes poly-G tails (if detected)
- Generates HTML and JSON reports

### Recommended Settings

```bash
fastp \
    -i R1.fq.gz -I R2.fq.gz \
    -o R1_clean.fq.gz -O R2_clean.fq.gz \
    --cut_right -q 20 -l 36 \
    --thread 8 \
    -h report.html -j report.json
```

## Comparison with Traditional Pipeline

| Task | Traditional | fastp |
|------|-------------|-------|
| QC report | FastQC | Built-in |
| Adapter trim | Cutadapt | Built-in |
| Quality trim | Trimmomatic | Built-in |
| Poly-G | Manual | Auto |
| Dedup | After align | Optional |
| Time | 3 steps | 1 step |

## Key Features

### Auto Adapter Detection

fastp identifies Illumina adapters automatically:
```bash
# Report shows detected adapters
fastp -i in.fq -o out.fq
# Check HTML report for "Adapter" section
```

### Poly-G Trimming

Essential for NovaSeq/NextSeq (two-color chemistry):
```bash
fastp -i in.fq -o out.fq --trim_poly_g
```

### Sliding Window Quality Trim

`--cut_right` is similar to Trimmomatic's SLIDINGWINDOW:
```bash
fastp -i in.fq -o out.fq \
      --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20
```

### Read Merging

For paired reads with overlap (e.g., amplicons):
```bash
fastp -i R1.fq -I R2.fq \
      --merge --merged_out merged.fq
```

## MultiQC Integration

fastp JSON reports are compatible with MultiQC:
```bash
# Process multiple samples
for sample in sample1 sample2 sample3; do
    fastp -i ${sample}_R1.fq.gz -I ${sample}_R2.fq.gz \
          -o ${sample}_R1_clean.fq.gz -O ${sample}_R2_clean.fq.gz \
          -j ${sample}_fastp.json -h ${sample}_fastp.html
done

# Aggregate reports
multiqc .
```

## Resources

- [fastp GitHub](https://github.com/OpenGene/fastp)
- [fastp Publication](https://doi.org/10.1093/bioinformatics/bty560)
