# Methylation-Based Detection - Usage Guide

## Overview
Analyze cfDNA methylation patterns for cancer detection and tissue-of-origin analysis. Uses bisulfite sequencing or cfMeDIP-seq for methylation profiling.

## Prerequisites
```bash
# MethylDackel
conda install -c bioconda methyldackel

# Dependencies
pip install pandas numpy scipy matplotlib
```

## Quick Start
Tell your AI agent what you want to do:
- "Extract methylation from my bisulfite-seq cfDNA BAM"
- "Find differentially methylated regions between cancer and normal"
- "Perform tissue-of-origin deconvolution from cfDNA methylation"
- "Analyze MCED-style regions for cancer detection"

## Example Prompts

### Methylation Extraction
> "Extract CpG methylation levels from my bisulfite BAM using MethylDackel."

> "Calculate beta values from methylation calls."

### DMR Analysis
> "Find differentially methylated regions between my cancer and normal samples."

> "Identify hypermethylated regions in cancer cfDNA."

### Tissue Deconvolution
> "Deconvolve tissue composition from cfDNA methylation using reference atlas."

## What the Agent Will Do
1. Extract methylation from bisulfite-seq BAM
2. Calculate beta values for CpG sites
3. Identify differentially methylated regions
4. Perform tissue deconvolution if reference available
5. Score cancer-specific methylation signatures

## Tips
- MethylDackel is actively maintained and nf-core integrated
- cfMeDIP-seq is good for low input (>= 5ng)
- Bisulfite-seq gives single-base resolution but needs >= 10ng
- Tissue deconvolution requires reference methylomes
- Limitation: Often requires >= 10ng cfDNA input

## Related Skills
- cfdna-preprocessing - Preprocess before methylation analysis
- fragment-analysis - Complement with fragmentomics
- methylation-analysis/bismark-alignment - General methylation processing
