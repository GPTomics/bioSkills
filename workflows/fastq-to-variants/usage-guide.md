# FASTQ to Variants - Usage Guide

## Overview

This workflow takes you from raw DNA sequencing FASTQ files to a filtered set of variant calls (SNPs and indels). It covers the entire process from quality control through alignment and variant calling.

## Prerequisites

```bash
# CLI tools
conda install -c bioconda fastp bwa-mem2 samtools bcftools

# For GATK path
conda install -c bioconda gatk4
```

## Quick Start

Tell your AI agent what you want to do:
- "Call variants from my whole genome sequencing FASTQ files"
- "Run the FASTQ to variants pipeline on my exome data"
- "I have paired-end DNA reads, help me find variants"

## Example Prompts

### Starting from FASTQ
> "I have FASTQ files for 5 samples, call variants jointly"

> "Process my WGS data from raw reads to VCF"

> "Use GATK HaplotypeCaller instead of bcftools"

### Customizing the workflow
> "Add BQSR to my variant calling pipeline"

> "Call variants only in the exome target regions"

> "Run the pipeline on a specific chromosome"

### From alignment to variants
> "I already have BAM files, just call variants"

> "My BAMs are not duplicate-marked, process and call variants"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| FASTQ files | .fastq.gz | Paired-end reads (R1 and R2 per sample) |
| Reference | FASTA | Reference genome (indexed for bwa-mem2) |
| Targets (optional) | BED | For exome/targeted sequencing |
| Known sites (GATK) | VCF | dbSNP for BQSR |

## What the Workflow Does

1. **Quality Control** - Trim adapters and low-quality bases
2. **Alignment** - Map reads to reference genome
3. **BAM Processing** - Sort, mark duplicates, index
4. **Variant Calling** - Identify SNPs and indels
5. **Filtering** - Remove low-quality calls

## Choosing Between bcftools and GATK

| Use bcftools when | Use GATK when |
|-------------------|---------------|
| Speed is important | Quality is paramount |
| Germline variants | Somatic variants |
| Small cohort | Large cohort with VQSR |
| Resource-limited | Resources available |

## Tips

- **Read groups**: Always add read group information during alignment
- **Duplicates**: Mark duplicates for PCR-based libraries
- **Depth**: Check coverage before variant calling (aim for 30x WGS, 100x exome)
- **Joint calling**: Improves sensitivity, especially for rare variants
- **Filtering**: Start with default filters, adjust based on Ti/Tv ratio and known site overlap
