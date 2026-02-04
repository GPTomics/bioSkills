# Isoform Switching - Usage Guide

## Overview
Analyze isoform switching events and predict their functional consequences. Identifies cases where the dominant isoform changes between conditions and predicts impacts on protein domains, NMD sensitivity, and coding potential.

## Prerequisites
```r
# R packages
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('IsoformSwitchAnalyzeR')

# External tools for functional annotation:
# - CPC2: http://cpc2.gao-lab.org/
# - Pfam: https://www.ebi.ac.uk/Tools/hmmer/
# - SignalP: https://services.healthtech.dtu.dk/service.php?SignalP
# - IUPred2A: https://iupred2a.elte.hu/
```

## Quick Start
Tell your AI agent what you want to do:
- "Analyze isoform switches between my conditions"
- "Find genes where the dominant isoform changes"
- "Predict functional consequences of splicing changes"
- "Identify isoform switches that affect protein domains"

## Example Prompts

### Basic Analysis
> "Import my Salmon quantification and identify isoform switches between control and treatment."

> "Find significant isoform switching events with dIF greater than 0.1."

### Functional Annotation
> "Analyze the functional consequences of isoform switches including domain changes and NMD sensitivity."

> "Identify switches that result in loss of protein domains."

### Visualization
> "Generate switch plots for my top differentially used isoforms."

> "Summarize consequence types across all switching genes."

## What the Agent Will Do
1. Import transcript quantification from Salmon/kallisto
2. Filter lowly expressed isoforms
3. Test for isoform switches using DEXSeq
4. Run external annotation tools (CPC2, Pfam, etc.)
5. Predict functional consequences
6. Generate visualizations and summary statistics

## Tips
- IsoformSwitchAnalyzeR version 2.x supports long-read and single-cell data
- Switch q-value < 0.05 and dIF > 0.1 are typical thresholds
- External tools (CPC2, Pfam) require running separately and importing results
- Focus on switches with functional consequences for biological interpretation
- NMD-sensitive switches may reduce protein levels rather than change function

## Related Skills
- differential-splicing - Identify differential events first
- splicing-quantification - PSI-level analysis
- pathway-analysis/go-enrichment - Pathway enrichment of switching genes
