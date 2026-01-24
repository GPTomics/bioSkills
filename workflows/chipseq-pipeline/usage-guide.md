# ChIP-seq Pipeline - Usage Guide

## Overview

This workflow processes ChIP-seq data from raw FASTQ files to annotated peaks. It handles both narrow peaks (transcription factors) and broad peaks (histone modifications).

## Prerequisites

```bash
# CLI tools
conda install -c bioconda fastp bowtie2 samtools macs3 deeptools bedtools

# R packages
BiocManager::install(c('ChIPseeker', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the ChIP-seq pipeline on my IP and Input samples"
- "Call peaks from my H3K4me3 ChIP-seq data"
- "Process my transcription factor ChIP-seq experiment"

## Example Prompts

### Starting from FASTQ
> "I have ChIP-seq FASTQ files for 2 IP replicates and 2 input controls"

> "Call broad peaks for my H3K27me3 ChIP-seq"

> "Run ChIP-seq analysis with mouse genome"

### Annotation and downstream
> "Annotate my peaks with nearby genes"

> "Find peaks in promoter regions"

> "Compare peaks between treatment and control"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| IP FASTQ | .fastq.gz | ChIP samples (immunoprecipitated) |
| Input FASTQ | .fastq.gz | Input controls (critical for peak calling) |
| Reference | FASTA | Reference genome + Bowtie2 index |

## What the Workflow Does

1. **Quality Control** - Trim adapters and low-quality bases
2. **Alignment** - Map reads to reference with Bowtie2
3. **BAM Processing** - Remove duplicates and low-quality alignments
4. **Peak Calling** - Identify enriched regions with MACS3
5. **QC** - Calculate FRiP and enrichment metrics
6. **Annotation** - Link peaks to genes with ChIPseeker

## Narrow vs Broad Peaks

| Type | Examples | MACS3 Flag |
|------|----------|------------|
| Narrow | TFs, H3K4me3, H3K27ac | Default |
| Broad | H3K36me3, H3K27me3, H3K9me3 | --broad |

## Tips

- **Input controls**: Always include matched input controls
- **Replicates**: Use biological replicates for reproducibility (check with IDR)
- **Sequencing depth**: 20-50M reads for TFs, 50-100M for histone marks
- **Blacklist**: Filter peaks against ENCODE blacklist regions
- **FRiP**: Low FRiP (<1%) indicates poor enrichment
