# Tumor Fraction Estimation - Usage Guide

## Overview
Estimate circulating tumor DNA fraction from shallow whole-genome sequencing using ichorCNA. Detects copy number alterations and calculates ctDNA percentage for treatment monitoring.

## Prerequisites
```r
# R packages
install.packages('devtools')
devtools::install_github('GavinHaLab/ichorCNA')

# HMMcopy for read counting
# BiocManager::install('HMMcopy')
```

```bash
# For read counting
conda install -c bioconda hmmcopy
```

## Quick Start
Tell your AI agent what you want to do:
- "Estimate tumor fraction from my shallow WGS data"
- "Run ichorCNA on my sWGS BAM files"
- "Calculate ctDNA percentage for treatment monitoring"
- "Set up a tumor fraction estimation pipeline"

## Example Prompts

### Single Sample
> "Run ichorCNA to estimate tumor fraction from my 0.5x sWGS BAM."

> "Generate read count WIG files and run ichorCNA analysis."

### Batch Processing
> "Process all my sWGS samples through ichorCNA in parallel."

> "Parse ichorCNA results from all samples into a summary table."

### Interpretation
> "What tumor fraction does my sample have and is it reliable?"

## What the Agent Will Do
1. Generate read counts in genomic bins (readCounter)
2. Apply GC and mappability correction
3. Run HMM segmentation for CNAs
4. Estimate tumor fraction and ploidy
5. Output segments and summary statistics

## Tips
- Use GavinHaLab/ichorCNA fork (v0.5.1+, actively maintained)
- Requires sWGS data (0.1-1x), NOT targeted panels
- Sensitivity: 97-100% at >= 3% tumor fraction
- 0.5x coverage is optimal for most applications
- Tumor fractions < 3% are near detection limit

## Related Skills
- cfdna-preprocessing - Preprocess BAMs before ichorCNA
- fragment-analysis - Complementary fragmentomics analysis
- ctdna-mutation-detection - Mutation detection from panel data
- copy-number/cnvkit-analysis - CNV concepts
