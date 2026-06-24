# RNA-seq to Differential Expression - Usage Guide

## Overview

This workflow takes you from raw RNA-seq FASTQ files to a list of differentially expressed genes. It combines quality control, quantification, and statistical analysis into a complete pipeline.

## Prerequisites

```bash
# CLI tools
conda install -c bioconda fastp salmon star subread

# R packages
install.packages('BiocManager')
BiocManager::install(c('DESeq2', 'tximport', 'apeglm'))
install.packages(c('ggplot2', 'pheatmap', 'ggrepel'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the RNA-seq to DE workflow on my FASTQ files"
- "I have paired-end RNA-seq data, help me find differentially expressed genes"
- "Process my RNA-seq from raw reads through DESeq2"

## Example Prompts

### Starting from FASTQ
> "I have FASTQ files for 3 control and 3 treated samples, run the full RNA-seq pipeline"

> "Quantify my RNA-seq with Salmon and run DESeq2"

> "Use STAR alignment instead of Salmon for my RNA-seq"

### Customizing the workflow
> "Run the pipeline but use a custom GTF file"

> "Add batch correction to my RNA-seq analysis"

> "Skip the QC step, my reads are already trimmed"

### From quantification to DE
> "I already have Salmon quant files, just run DESeq2"

> "Import my featureCounts output into DESeq2"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| FASTQ files | .fastq.gz | Paired-end reads (R1 and R2 per sample) |
| Sample metadata | CSV/TSV | Sample names and experimental conditions |
| Reference | FASTA + GTF | Transcriptome for Salmon, genome + GTF for STAR |
| tx2gene | CSV | Transcript-to-gene mapping (for Salmon) |

## What the Workflow Does

1. **Quality Control** - Trim adapters and low-quality bases with fastp
2. **Quantification** - Count reads per gene/transcript (Salmon or STAR+featureCounts)
3. **Import** - Load counts into R with proper normalization offsets
4. **Differential Expression** - Run DESeq2 statistical analysis
5. **Results** - Export significant genes and visualizations

## Choosing Between Paths

| Criterion | Salmon Path | STAR Path |
|-----------|-------------|-----------|
| Speed | Faster (no alignment) | Slower |
| Storage | Lower (no BAM files) | Higher (BAM files) |
| Use case | Standard DE analysis | Need BAMs for other analyses |
| Accuracy | Excellent for DE | Excellent for DE |
| Novel junctions | No | Yes (with 2-pass) |

## Tips

- **Replicates**: Minimum 3 biological replicates per condition (more is better)
- **Sequencing depth**: 20-30M reads per sample is typical for DE analysis
- **Decoy index**: Build the Salmon index with genome decoys; a transcriptome-only index misassigns intron and pseudogene reads
- **Library type**: Salmon auto-detects with `-l A`; verify the call by reading `lib_format_counts.json`, and set featureCounts `-s` to match (dUTP/TruSeq is reverse, `-s 2`)
- **Batch effects**: If samples were processed in batches, include batch in the design formula, unless batch is confounded with condition (then the design is unfixable)
- **Outliers**: Check PCA; let DESeq2 handle single-gene outliers via Cook's distance, and remove a whole-sample outlier only with a documented technical cause
- **Transcript-level**: For differential transcript expression or usage, generate Salmon Gibbs samples and use the transcript-aware tools in alternative-splicing/isoform-switching
