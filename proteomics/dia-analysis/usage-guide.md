# DIA Analysis - Usage Guide

## Overview
Process data-independent acquisition (DIA) mass spectrometry data for comprehensive proteome quantification with fewer missing values than DDA.

## Prerequisites
```bash
pip install pandas numpy
# CLI: DIA-NN (recommended), MSFragger-DIA, OpenSWATH
# Commercial: Spectronaut
```

## Quick Start
Tell your AI agent what you want to do:
- "Run DIA-NN on my mzML files in library-free mode"
- "Process DIA data using a spectral library I built"
- "Load DIA-NN results and prepare for statistical analysis"

## Example Prompts

### Library-Free Analysis
> "Run DIA-NN in library-free mode against the UniProt human FASTA with 1% FDR"

> "Set up library-free DIA analysis with trypsin digestion and standard modifications"

> "Process my DIA mzML files without a spectral library using deep learning prediction"

### Library-Based Analysis
> "Search my DIA data against the spectral library from my previous DDA experiments"

> "Run DIA-NN with my Prosit-predicted library for targeted analysis"

> "Use match-between-runs with my spectral library for improved coverage"

### Parameter Configuration
> "Configure DIA-NN with 10 ppm mass accuracy and 1 missed cleavage"

> "Set up DIA analysis with phosphorylation as a variable modification"

> "Enable two-pass analysis (reanalyse) for improved quantification"

### Results Processing
> "Load the DIA-NN report.tsv and create a protein abundance matrix"

> "Filter DIA results to precursor q-value < 0.01 and protein q-value < 0.01"

> "Compare protein identifications between library-free and library-based modes"

## What the Agent Will Do
1. Configure DIA-NN parameters (FASTA, library, tolerances)
2. Run search in library-free or library-based mode
3. Apply FDR filtering at precursor and protein level
4. Export protein/precursor matrices
5. Load results for downstream analysis

## Library-Free vs Library-Based

| Mode | Description | Use When |
|------|-------------|----------|
| Library-free | Deep learning predicts spectra | Quick analysis, no prior data |
| Library-based | Match to experimental spectra | Higher sensitivity for known targets |
| Hybrid | Predicted + empirical library | Best coverage for large studies |

## Key Parameters
- `--qvalue 0.01` - 1% FDR at precursor and protein level
- `--reanalyse` - Two-pass analysis for match-between-runs
- `--smart-profiling` - Improved quantification accuracy
- `--min-fr-mz 200 --max-fr-mz 1800` - Fragment ion range

## Tips
- Library-free mode is often sufficient for discovery proteomics
- Use --reanalyse for better quantification across many samples
- DIA typically has fewer missing values than DDA
- Filter to 1% FDR at both precursor and protein levels
- For large cohorts, consider building a project-specific library
