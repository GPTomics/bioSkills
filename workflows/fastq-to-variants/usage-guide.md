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

## What the Agent Will Do

1. Trim adapters and low-quality bases with fastp
2. Align reads to reference genome with bwa-mem2
3. Sort, mark duplicates, and index BAM files
4. Call SNPs and indels with bcftools or GATK HaplotypeCaller
5. Filter variants by quality metrics and validate with Ti/Tv ratio

## Choosing Between bcftools and GATK

| Use bcftools when | Use GATK when |
|-------------------|---------------|
| Speed is important | Quality is paramount |
| Germline variants | Somatic variants |
| Small cohort | Large cohort with VQSR |
| Resource-limited | Resources available |

## Tips

- Always add read group information during alignment -- required for GATK
- Mark duplicates for PCR-based libraries (PCR-free libraries still benefit)
- Check coverage before variant calling (aim for 30x WGS, 100x exome)
- Joint calling improves sensitivity, especially for rare variants
- Start with default filters, adjust based on Ti/Tv ratio and known site overlap

## Related Skills

- variant-calling/variant-calling - bcftools calling details
- variant-calling/gatk-variant-calling - GATK HaplotypeCaller and DRAGEN mode
- variant-calling/filtering-best-practices - VQSR and hard filtering
- read-qc/fastp-workflow - Read QC details
- read-alignment/bwa-alignment - Alignment details
