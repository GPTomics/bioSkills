# GATK Variant Calling - Usage Guide

## Overview

Germline variant calling using GATK HaplotypeCaller, the gold standard for production-quality SNP and indel detection following GATK Best Practices.

## Prerequisites

```bash
# GATK 4.x
conda install -c bioconda gatk4

# Or download from Broad
# https://github.com/broadinstitute/gatk/releases
```

## Quick Start

Tell your AI agent what you want to do:
- "Call germline variants from my aligned BAM file"
- "Run HaplotypeCaller in GVCF mode for joint calling"
- "Apply BQSR before variant calling"
- "Set up a variant calling pipeline for whole genome data"

## Example Prompts

### Single Sample Calling
> "Call variants from sample.bam using GATK HaplotypeCaller"

> "Run HaplotypeCaller on my exome data with the capture regions BED file"

### GVCF for Joint Calling
> "Generate GVCFs from all BAM files in my cohort for joint genotyping"

> "Run HaplotypeCaller in GVCF mode so I can combine samples later"

### With BQSR
> "Apply base quality score recalibration before variant calling"

> "Run the full GATK best practices pipeline with BQSR and variant calling"

### Filtering
> "Apply VQSR to my whole genome variant calls"

> "Use hard filtering for my small exome cohort"

## What the Agent Will Do

1. Verify input BAM is sorted, indexed, and has read groups
2. Check for required reference files (FASTA, dbSNP, known sites)
3. Run BQSR if requested (BaseRecalibrator + ApplyBQSR)
4. Execute HaplotypeCaller with appropriate settings (VCF or GVCF mode)
5. Apply filtering (VQSR for large WGS cohorts, hard filters otherwise)
6. Generate variant statistics for quality assessment

## Tips

- Always use `-ERC GVCF` when planning to do joint genotyping later
- VQSR requires many variants to train the model; use hard filters for exomes or small cohorts
- Mark duplicates before variant calling (Picard MarkDuplicates or samtools markdup)
- BQSR improves accuracy but requires known sites files (dbSNP, Mills indels)
- For WGS, expect Ti/Tv ratio around 2.0-2.1 for high-quality SNP calls

## When to Use GATK

- Production-quality variant calls
- Cohort analysis with joint genotyping
- When VQSR is possible (many variants)
- Following GATK Best Practices

## Workflow Options

### Single Sample
1. Mark duplicates (Picard/samtools)
2. BQSR (optional but recommended)
3. HaplotypeCaller
4. Hard filtering

### Cohort (Recommended)
1. Mark duplicates per sample
2. BQSR per sample
3. HaplotypeCaller with -ERC GVCF per sample
4. GenomicsDBImport or CombineGVCFs
5. GenotypeGVCFs
6. VQSR or hard filtering

## GVCF vs VCF

- **VCF**: Final variant calls, reference sites not included
- **GVCF**: Includes reference blocks for joint genotyping
- Always use GVCF for cohorts to capture reference confidence

## VQSR vs Hard Filtering

### Use VQSR when:
- Whole genome sequencing
- Large cohort (many variants)
- Resource files available

### Use Hard Filtering when:
- Exome/panel sequencing
- Single sample or small cohort
- Non-model organism
- VQSR fails (not enough variants)

## Resource Files

Download from GATK Resource Bundle:
- gs://genomics-public-data/resources/broad/hg38/v0/

Key files:
- Homo_sapiens_assembly38.fasta
- Homo_sapiens_assembly38.dbsnp138.vcf
- hapmap_3.3.hg38.vcf.gz
- 1000G_omni2.5.hg38.vcf.gz
- Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
