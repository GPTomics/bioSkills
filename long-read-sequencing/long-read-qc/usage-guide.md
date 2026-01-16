# Long-Read Quality Control - Usage Guide

## Overview

Quality control for long-read data involves checking read length distribution, quality scores, and throughput. Tools like NanoPlot visualize these metrics, while chopper/NanoFilt filter reads.

## When to Use This Skill

- You received raw long-read data and want to assess quality
- You need to filter reads by length or quality before analysis
- You want to generate QC reports for publication

## Installation

```bash
conda install -c bioconda nanoplot nanostat chopper seqkit
```

## Quick QC Workflow

```bash
# 1. Get statistics
NanoStat --fastq reads.fastq.gz

# 2. Generate plots
NanoPlot --fastq reads.fastq.gz -o qc_output

# 3. Filter if needed
gunzip -c reads.fastq.gz | chopper -q 10 -l 1000 | gzip > filtered.fastq.gz
```

## Interpreting Results

### Good Quality Indicators

- Mean quality Q15+ (ONT SUP) or Q30+ (HiFi)
- N50 > 10kb for most applications
- High percentage of reads passing filters

### Warning Signs

- Bimodal quality distribution
- Very short N50 compared to expected
- Large fraction of low-quality reads

## Filtering Guidelines

| Application | Min Quality | Min Length |
|-------------|-------------|------------|
| Assembly | Q10 | 1000bp |
| Variant calling | Q15 | 500bp |
| Structural variants | Q10 | 1000bp |
| Transcriptome | Q10 | 300bp |

## Resources

- [NanoPlot GitHub](https://github.com/wdecoster/NanoPlot)
- [chopper GitHub](https://github.com/wdecoster/chopper)
