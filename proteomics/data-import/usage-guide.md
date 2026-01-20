# Data Import - Usage Guide

## Overview
Load and parse mass spectrometry data from various formats (mzML, MaxQuant, DIA-NN) for proteomics analysis.

## Prerequisites
```bash
pip install pyopenms pandas numpy
# R alternative: BiocManager::install("MSnbase")
```

## Quick Start
Tell your AI agent what you want to do:
- "Load my MaxQuant proteinGroups.txt file and filter contaminants"
- "Read mzML files from my experiment folder"
- "Import DIA-NN output and prepare for statistical analysis"

## Example Prompts

### Loading Search Engine Output
> "Load MaxQuant proteinGroups.txt, remove contaminants and reverse hits, and log2-transform the LFQ intensities"

> "Import the DIA-NN report.tsv and create a protein abundance matrix"

### Loading Raw MS Data
> "Read all mzML files in my data folder using pyOpenMS"

> "Parse the MS1 spectra from my mzML file and extract precursor information"

### Data Preparation
> "Check the missing value pattern in my MaxQuant output"

> "Filter to proteins with at least 2 unique peptides and valid values in 70% of samples"

## What the Agent Will Do
1. Load data from specified format (mzML, MaxQuant, DIA-NN)
2. Apply standard filtering (contaminants, decoys, only-by-site)
3. Log2-transform intensities for downstream analysis
4. Report missing value statistics

## Supported Formats

| Format | Description | Tool |
|--------|-------------|------|
| mzML | Open standard for MS data | pyOpenMS, MSnbase |
| mzXML | Legacy open format | pyOpenMS |
| proteinGroups.txt | MaxQuant protein output | pandas |
| evidence.txt | MaxQuant peptide output | pandas |
| report.tsv | DIA-NN output | pandas |

## Missing Value Patterns
- **MCAR**: Missing completely at random (rare in proteomics)
- **MAR**: Missing at random (can impute)
- **MNAR**: Missing not at random (low abundance) - use left-censored imputation

## Tips
- Use `low_memory=False` when loading large MaxQuant files
- Check for batch effects in missing value patterns
- Log2-transform intensities before statistical analysis
- Filter contaminants, reverse, and only-by-site columns before analysis
