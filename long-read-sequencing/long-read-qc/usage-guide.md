# Long-Read Quality Control - Usage Guide

## Overview
Quality control for long-read data involves checking read length distribution, quality scores, and throughput. NanoPlot visualizes these metrics, while chopper filters reads by quality and length.

## Prerequisites
```bash
conda install -c bioconda nanoplot nanostat chopper seqkit
```

## Quick Start
Tell your AI agent what you want to do:
- "Generate QC report for my Nanopore reads"
- "Filter reads below Q10 and shorter than 1kb"

## Example Prompts

### Generate QC Report
> "Run NanoPlot on reads.fastq.gz and generate quality plots"

> "Get read statistics for my ONT run using NanoStat"

### Filtering
> "Filter my reads to keep only those with Q10+ quality and 1000bp+ length"

> "Remove low-quality reads below Q15 from my HiFi dataset"

### Combined Workflow
> "Run QC on my Nanopore data, show me the statistics, then filter for assembly-quality reads"

## What the Agent Will Do
1. Run NanoStat to get summary statistics (mean quality, N50, total bases)
2. Generate NanoPlot visualizations (length histogram, quality distribution)
3. Apply filters based on your application requirements
4. Report before/after statistics

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

## Tips
- Run QC before any downstream analysis to catch issues early
- N50 is more informative than mean length for long reads
- Bimodal quality distribution may indicate mixed flowcell performance
- Use chopper for fast filtering; it's the successor to NanoFilt
- Keep original reads - filtering is lossy

## Resources
- [NanoPlot GitHub](https://github.com/wdecoster/NanoPlot)
- [chopper GitHub](https://github.com/wdecoster/chopper)
