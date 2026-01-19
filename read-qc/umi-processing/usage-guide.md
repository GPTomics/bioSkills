# UMI Processing - Usage Guide

## Overview

Unique Molecular Identifiers (UMIs) are random sequences added to molecules before PCR amplification. They enable distinguishing PCR duplicates from biological duplicates, crucial for accurate quantification in RNA-seq, targeted sequencing, and single-cell applications.

## When to Use This Skill

- Processing libraries prepared with UMI-containing adapters
- Single-cell RNA-seq (10X Genomics, Drop-seq, etc.)
- Low-input RNA-seq requiring accurate molecule counting
- Targeted sequencing panels
- Any workflow where PCR duplicate rates are high

## Installation

```bash
conda install -c bioconda umi_tools
```

## Basic Workflow

### 1. Extract UMIs

Move UMIs from read sequence to read header:

```bash
umi_tools extract \
    --stdin=R1.fastq.gz \
    --read2-in=R2.fastq.gz \
    --stdout=R1_umi.fastq.gz \
    --read2-out=R2_umi.fastq.gz \
    --bc-pattern=NNNNNNNN
```

### 2. Align Reads

Use your preferred aligner (STAR, BWA, etc.):

```bash
STAR --genomeDir ref --readFilesIn R1_umi.fq.gz R2_umi.fq.gz ...
samtools sort -o aligned_sorted.bam aligned.bam
samtools index aligned_sorted.bam
```

### 3. Deduplicate

Remove PCR duplicates based on UMI + alignment position:

```bash
umi_tools dedup -I aligned_sorted.bam -S deduplicated.bam
```

## UMI Pattern Syntax

| Pattern | Description |
|---------|-------------|
| `N` | UMI base |
| `C` | Cell barcode base |
| `X` | Base to discard |
| `NNNNNNNN` | 8bp UMI |
| `CCCCCCCCNNNNNNNN` | Cell barcode + UMI |

## Common Library Types

| Library | Pattern |
|---------|---------|
| NEBNext | `NNNNNNNN` (8bp in R1) |
| 10X 3' v3 | `CCCCCCCCCCCCCCCCNNNNNNNNNNNN` (16bp CB + 12bp UMI) |
| TruSeq UMI | `NNNNNNNNN` (9bp in index) |

## Deduplication Methods

| Method | Best For |
|--------|----------|
| `directional` | Most RNA-seq (default) |
| `unique` | High diversity libraries |
| `cluster` | High error rates |

## QC Metrics

Good libraries typically show:
- 20-60% deduplication rate (RNA-seq)
- >70% may indicate over-amplification
- <10% may indicate under-sequencing

## Resources

- [umi_tools documentation](https://umi-tools.readthedocs.io/)
- [umi_tools publication](https://doi.org/10.1101/gr.209601.116)
