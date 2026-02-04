# cfDNA Preprocessing - Usage Guide

## Overview
Preprocess cell-free DNA sequencing data with UMI-aware duplicate removal. Essential for accurate mutation detection at low variant allele fractions in liquid biopsy samples.

## Prerequisites
```bash
# fgbio (recommended for UMI handling)
conda install -c bioconda fgbio

# Alignment
conda install -c bioconda bwa samtools

# Python dependencies
pip install pysam pandas numpy matplotlib
```

## Quick Start
Tell your AI agent what you want to do:
- "Preprocess my cfDNA FASTQ files with UMI handling"
- "Create UMI consensus reads from my plasma DNA sequencing"
- "Analyze fragment size distribution from my cfDNA BAM"
- "Set up a cfDNA preprocessing pipeline with fgbio"

## Example Prompts

### UMI Processing
> "Extract UMIs from my reads and create consensus sequences using fgbio."

> "Process my cfDNA BAM with UMI-aware duplicate removal requiring at least 2 reads per family."

### Fragment Analysis
> "Analyze the fragment size distribution to verify this is good quality cfDNA."

> "Check if my sample has the expected ~167bp mononucleosomal peak."

### Quality Control
> "Generate cfDNA quality metrics for my preprocessed samples."

## What the Agent Will Do
1. Extract UMIs from read sequences
2. Align reads with soft-clipping preserved
3. Group reads by UMI
4. Generate consensus sequences
5. Filter low-quality consensus reads
6. Analyze fragment size distribution for QC

## Tips
- Use fgbio 3.0+ for best UMI handling
- Modal cfDNA fragment size should be ~167 bp (mononucleosome)
- Pre-analytical factors matter: use Streck tubes for delayed processing
- Minimum 2 reads per UMI family for consensus
- Avoid cfDNApipe (uncertain maintenance) - use fgbio + individual tools

## Related Skills
- fragment-analysis - Analyze fragmentomics after preprocessing
- tumor-fraction-estimation - Estimate ctDNA from sWGS
- ctdna-mutation-detection - Detect mutations from panel data
