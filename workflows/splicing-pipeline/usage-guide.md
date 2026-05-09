# Splicing Pipeline - Usage Guide

## Overview
Complete bulk short-read alternative splicing analysis workflow from raw RNA-seq FASTQ files to differential splicing results and visualizations. Includes QC checkpoints and best practices. For variant-driven splice prediction use splice-variant-prediction; for rare-disease single-patient outlier detection use outlier-splicing-detection; for PacBio/ONT full-isoform analysis use long-read-splicing.

## Prerequisites
```bash
# Core tools
conda install -c bioconda star rmats-turbo ggsashimi fastp rseqc

# Python dependencies
pip install pandas numpy matplotlib

# R packages (optional, for isoform switching)
BiocManager::install('IsoformSwitchAnalyzeR')
```

## Quick Start
Tell your AI agent what you want to do:
- "Run a complete splicing analysis from my FASTQ files"
- "Set up a splicing pipeline for my RNA-seq experiment"
- "Analyze differential splicing between my conditions"
- "Generate sashimi plots for my significant splicing events"

## Example Prompts

### Full Pipeline
> "I have RNA-seq FASTQs from 3 control and 3 treatment samples. Run a complete splicing analysis."

> "Set up STAR 2-pass alignment and rMATS differential splicing for my samples."

### Individual Steps
> "Align my samples with STAR 2-pass mode for optimal junction detection."

> "Run junction saturation QC to check if I have enough sequencing depth."

> "Visualize the top 20 differential splicing events with sashimi plots."

## What the Agent Will Do
1. Perform read QC and adapter trimming (fastp)
2. Align with STAR 2-pass mode for junction detection
3. Run junction saturation QC checkpoint
4. Perform differential splicing analysis (rMATS-turbo)
5. Filter significant events (FDR < 0.05, |deltaPSI| > 0.1)
6. Optionally analyze isoform switches
7. Generate sashimi plots for visualization

## Tips
- STAR 2-pass mode is essential for novel junction discovery
- Check junction saturation curves before trusting results
- Use `--outSJfilterOverhangMin 8 8 8 8` for stringent junction filtering
- rMATS-turbo combines quantification and differential testing
- Standard thresholds: |deltaPSI| > 0.1, FDR < 0.05
- Always require minimum junction reads (>= 10) for reliability

## Related Skills

- alternative-splicing/splicing-quantification - PSI computation, event taxonomy
- alternative-splicing/differential-splicing - Tool selection, MAJIQ V3, Shiba, leafcutter
- alternative-splicing/isoform-switching - DTU and NMD/domain consequences
- alternative-splicing/sashimi-plots - Visualization tools
- alternative-splicing/splicing-qc - QC prerequisites
- alternative-splicing/single-cell-splicing - Chemistry-first single-cell decision
- alternative-splicing/splice-variant-prediction - Variant impact (SpliceAI/Pangolin)
- alternative-splicing/outlier-splicing-detection - Rare-disease single-patient (FRASER 2.0)
- alternative-splicing/long-read-splicing - Full-isoform PacBio/ONT analysis
- read-alignment/star-alignment - STAR 2-pass cohort-style configuration
