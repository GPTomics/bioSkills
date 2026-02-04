# Differential Splicing - Usage Guide

## Overview
Detect differential alternative splicing between experimental conditions. Identifies splicing events that change significantly between groups, reporting both statistical significance and effect size (delta PSI).

## Prerequisites
```bash
# rMATS-turbo
pip install rmats-turbo

# SUPPA2
pip install suppa

# leafcutter (R)
# install.packages('devtools')
# devtools::install_github('davidaknowles/leafcutter/leafcutter')

# Python dependencies
pip install pandas numpy scipy
```

## Quick Start
Tell your AI agent what you want to do:
- "Find differential splicing between tumor and normal samples"
- "Compare splicing patterns between treatment and control groups"
- "Identify genes with significant splicing changes"
- "Run differential splicing analysis on my RNA-seq data"

## Example Prompts

### rMATS Analysis
> "I have BAM files from control and treatment conditions. Run rMATS-turbo to find differential splicing events."

> "Compare exon skipping between my two sample groups using rMATS."

### SUPPA2 Analysis
> "Use SUPPA2 diffSplice to compare splicing between conditions from my TPM data."

> "Find differential splicing events with delta PSI greater than 0.2."

### leafcutter Analysis
> "Run leafcutter differential splicing analysis to find novel differential intron usage."

> "Identify differential splicing at intron clusters between my conditions."

## What the Agent Will Do
1. Organize samples by condition/group
2. Run differential splicing analysis with selected tool
3. Calculate FDR-corrected p-values and delta PSI
4. Filter significant events (default: FDR < 0.05, |deltaPSI| > 0.1)
5. Generate ranked list of differential events by significance and effect size

## Tips
- Use |deltaPSI| > 0.1 for discovery, > 0.2 for high-confidence sets
- rMATS-turbo combines quantification and differential testing in one run
- SUPPA2 tends to be more stringent, returning fewer significant events
- leafcutter is good for novel junction discovery without annotation bias
- Always check junction read support for significant events

## Related Skills
- splicing-quantification - Calculate PSI values first
- isoform-switching - Analyze functional consequences
- sashimi-plots - Visualize significant events
- read-alignment/star-alignment - STAR 2-pass alignment required
