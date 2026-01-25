# Methylation Pipeline - Usage Guide

## Overview

This workflow processes bisulfite sequencing data from FASTQ to differentially methylated regions (DMRs), covering alignment, methylation calling, and statistical analysis.

## Prerequisites

```bash
# CLI tools
conda install -c bioconda bismark bowtie2 trim-galore samtools

# R packages
BiocManager::install(c('methylKit', 'genomation'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the methylation pipeline on my bisulfite-seq data"
- "Find differentially methylated regions between treatment and control"
- "Align my WGBS data with Bismark and call methylation"

## Example Prompts

### Starting from FASTQ
> "Process my RRBS data through methylation calling"

> "Align bisulfite sequencing reads and extract methylation"

> "Run Bismark on my whole-genome bisulfite data"

### Analysis
> "Find DMRs with methylKit"

> "Compare methylation between tumor and normal"

> "Annotate my DMRs with gene features"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| FASTQ files | .fastq.gz | Paired-end bisulfite-treated reads |
| Reference | FASTA | Genome (Bismark will prepare) |

## What the Workflow Does

1. **Quality Control** - Trim adapters and low-quality bases
2. **Alignment** - Map bisulfite-converted reads with Bismark
3. **Deduplication** - Remove PCR duplicates
4. **Methylation Calling** - Extract methylation status per CpG
5. **Analysis** - Statistical analysis with methylKit
6. **DMR Detection** - Find differentially methylated regions

## Tips

- **Coverage**: WGBS needs 10-30x coverage; RRBS can work with less
- **Conversion rate**: Should be >99%; check with spike-in controls
- **M-bias**: Check for position bias and trim if needed
- **Replicates**: Minimum 2-3 per condition for reliable DMR calling
