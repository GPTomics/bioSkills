# FASTQ to Variants - Usage Guide

## Overview

This workflow takes raw DNA sequencing FASTQ files to a filtered, normalized, and benchmarked set of germline variant calls (SNPs and indels). It is a chaining skill: it owns the pipeline-level decisions that only arise between steps -- the pipeline-wide reference-genome commitment, the defensible step order (normalize before annotate; filter site-level then genotype-level), engine and single-sample-vs-cohort selection, filtering strategy by cohort size, and stratified benchmarking -- and hands the mechanism of each step off to the variant-calling and read-alignment component skills.

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
- "Which reference genome should I commit to before I start calling?"
- "Benchmark my callset against GIAB truth"

## Example Prompts

### Starting from FASTQ
> "I have FASTQ files for 5 samples, call variants jointly"

> "Process my WGS data from raw reads to a normalized, filtered VCF"

> "Use GATK HaplotypeCaller instead of bcftools"

### Pipeline-level decisions
> "Should I use the GRCh38 analysis set or T2T-CHM13 for this cohort?"

> "Do I still need BQSR on NovaSeq data?"

> "My cohort is 4 exomes -- should I use VQSR or hard filters?"

> "In what order do I normalize, filter, and annotate?"

### Benchmarking and QC
> "Run hap.py against the HG002 truth set stratified by region"

> "My Ti/Tv is 1.8 on an exome -- what is wrong?"

### From alignment to variants
> "I already have BAM files, just call and normalize variants"

> "My BAMs are not duplicate-marked, process and call variants"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| FASTQ files | .fastq.gz | Paired-end reads (R1 and R2 per sample) |
| Reference | FASTA | Reference genome (indexed for bwa-mem2) |
| Targets (optional) | BED | For exome/targeted sequencing |
| Known sites (GATK) | VCF | dbSNP for BQSR |
| Truth set (benchmarking) | VCF + BED | GIAB truth VCF and confident-region BED, only when the sample is a GIAB genome |

## What the Agent Will Do

1. Confirm the pipeline-wide reference-genome commitment (GRCh38 analysis set vs T2T, ALT/decoy handling) before aligning anything
2. Trim adapters and low-quality bases with fastp
3. Align reads to the committed reference with bwa-mem2, adding read groups
4. Sort, mark duplicates, and index BAM files (skipping BQSR where the instrument makes it optional)
5. Call SNPs and indels per-sample, or per-sample gVCFs joint-genotyped for a cohort, with the engine chosen for the data
6. Normalize the VCF (left-align, split multiallelics) against the same reference, before any annotation or comparison
7. Filter site-level then genotype-level, with a strategy chosen by cohort size
8. Annotate on the normalized representation, then benchmark stratified with hap.py/vcfeval and sanity-check Ti/Tv

## Choosing an Engine and Cohort Mode

Engine selection is a pipeline-level decision handed to variant-calling/variant-calling: bcftools for exploratory/non-model/resource-limited work, GATK HaplotypeCaller for auditable open-source cohorts, DeepVariant for best indel/difficult-region accuracy, DRAGEN for maximum throughput. Single-sample vs cohort is a separate chaining decision: for a cohort, emit per-sample gVCFs and joint-genotype them (variant-calling/joint-calling) rather than merging single-sample VCFs.

## Tips

- Commit the reference genome before aligning a single read; every downstream coordinate (annotation DB, benchmark truth, cohort) inherits it
- Always add read group information during alignment -- required for GATK
- BQSR is honestly optional on modern binned-quality instruments (NovaSeq); DeepVariant recommends skipping it
- Normalize BEFORE annotating -- an un-left-aligned indel silently misses its database record
- Never build a cohort by bcftools merge of single-sample VCFs -- absence becomes a fabricated hom-ref genotype
- Filter site-level first, then genotype-level, then recompute cohort QC on the filtered matrix
- Do not trust a single genome-wide F1 -- benchmark stratified within the GIAB confident region, SNPs and indels separately

## Related Skills

- database-access/sra-data - Pull public WGS FASTQ for reanalysis
- database-access/ncbi-datasets-cli - Pull reference genome assembly via Datasets v2 CLI
- variant-calling/variant-calling - Engine selection and bcftools calling details
- variant-calling/gatk-variant-calling - GATK HaplotypeCaller and DRAGEN mode
- variant-calling/deepvariant - Deep-learning calling and GLnexus cohorts
- variant-calling/joint-calling - Cohort joint genotyping and scaling
- variant-calling/variant-normalization - Normalize before annotate/compare
- variant-calling/filtering-best-practices - VQSR and hard filtering
- variant-calling/variant-annotation - Annotate on the normalized representation
- variant-calling/vcf-statistics - Ti/Tv, het/hom, and identity QC
- read-qc/fastp-workflow - Read QC details
- sequence-io/fastq-quality - Confirm the FASTQ quality encoding before trimming public or legacy data
- sequence-io/paired-end-fastq - Keep R1/R2 mates synchronized so filtering does not corrupt mapping
- read-alignment/bwa-alignment - Alignment details, read groups, ALT/decoy analysis set, the dedup ordering
