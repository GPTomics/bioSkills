# fastp Workflow - Usage Guide

## Overview
fastp is a modern, all-in-one FASTQ preprocessing tool that replaces the need for separate adapter trimming (Cutadapt), quality filtering (Trimmomatic), and QC reporting (FastQC) steps. It is fast, memory-efficient, and produces comprehensive HTML reports.

## Prerequisites
```bash
conda install -c bioconda fastp
```

## Quick Start
Tell your AI agent what you want to do:
- "Process my FASTQ files with fastp"
- "Run all-in-one preprocessing on my paired-end reads"
- "Clean up my sequencing data with automatic adapter detection"

## Example Prompts

### Basic Processing
> "Run fastp on my paired-end FASTQ files with default settings"

> "Process all samples in my data directory with fastp and generate reports"

### Custom Parameters
> "Use fastp with Q20 quality threshold and minimum length of 36bp"

> "Run fastp with sliding window trimming from the 3' end"

### Special Features
> "Merge overlapping paired-end reads using fastp"

> "Deduplicate my reads at the FASTQ level with fastp"

> "Process my NovaSeq data with poly-G trimming enabled"

### Batch Processing
> "Run fastp on all my samples and aggregate reports with MultiQC"

## What the Agent Will Do
1. Run fastp with automatic adapter detection
2. Apply quality filtering and trimming
3. Handle platform-specific issues (poly-G for NovaSeq)
4. Generate HTML and JSON reports
5. Optionally merge overlapping reads or deduplicate

## Comparison with Traditional Pipeline

| Task | Traditional | fastp |
|------|-------------|-------|
| QC report | FastQC | Built-in |
| Adapter trim | Cutadapt | Built-in |
| Quality trim | Trimmomatic | Built-in |
| Poly-G | Manual | Auto |
| Dedup | After align | Optional |
| Time | 3 steps | 1 step |

## Tips
- fastp auto-detects Illumina adapters, no need to specify them manually
- Use `--cut_right` for sliding window trimming similar to Trimmomatic
- JSON reports integrate seamlessly with MultiQC
- Read merging is useful for amplicon sequencing with overlapping pairs
- FASTQ-level deduplication can reduce file sizes before alignment
- Always review the HTML report for QC metrics and filtering statistics

## Resources
- [fastp GitHub](https://github.com/OpenGene/fastp)
- [fastp Publication](https://doi.org/10.1093/bioinformatics/bty560)
