# Peptide Identification - Usage Guide

## Overview
Match MS/MS spectra to peptide sequences using database search or spectral library matching.

## Prerequisites
```bash
pip install pyopenms
# CLI tools: MSFragger, Comet, or X!Tandem
# R alternative: BiocManager::install(c("mzID", "MSnbase"))
```

## Quick Start
Tell your AI agent what you want to do:
- "Run a database search on my mzML files against UniProt human"
- "Set up peptide identification with trypsin digestion and 10 ppm tolerance"
- "Parse search results and filter to 1% FDR"

## Example Prompts

### Database Search Setup
> "Configure a database search with trypsin, 2 missed cleavages, carbamidomethyl C as fixed mod, and oxidation M as variable"

> "Set up MSFragger search with 10 ppm precursor tolerance and 0.02 Da fragment tolerance"

### Running Searches
> "Run a peptide search against the reviewed human proteome from UniProt"

> "Perform an open modification search to identify unknown PTMs"

### Results Processing
> "Parse the mzIdentML output and filter to 1% peptide FDR"

> "Convert pepXML results to a pandas DataFrame with PSM scores"

### Spectral Library Search
> "Search my DIA data against a spectral library for targeted quantification"

## What the Agent Will Do
1. Configure search parameters (enzyme, tolerances, modifications)
2. Set up target-decoy strategy for FDR control
3. Run database or spectral library search
4. Filter results to specified FDR threshold
5. Report identification statistics

## Key Parameters

| Parameter | Typical Value | Effect |
|-----------|---------------|--------|
| Precursor tolerance | 10-20 ppm | Match window for precursor m/z |
| Fragment tolerance | 0.02-0.05 Da | Match window for fragment ions |
| Missed cleavages | 2 | Allow incomplete digestion |
| Fixed modifications | Carbamidomethyl (C) | Always present |
| Variable modifications | Oxidation (M), Phospho (STY) | May be present |

## Common Enzymes

| Enzyme | Cleavage Rule |
|--------|---------------|
| Trypsin | After K, R (not before P) |
| Lys-C | After K |
| Chymotrypsin | After F, W, Y, L |
| Asp-N | Before D |

## Tips
- Always use target-decoy FDR for quality control
- 1% FDR at peptide level is standard; protein FDR is separate
- Use >= 2 unique peptides per protein for confident identification
- Variable modifications increase search space exponentially
