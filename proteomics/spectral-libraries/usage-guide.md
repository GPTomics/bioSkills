# Spectral Libraries - Usage Guide

## Overview
Build and use spectral libraries containing reference MS2 spectra for faster and more sensitive peptide identification in DIA and targeted proteomics.

## Prerequisites
```bash
pip install pandas numpy
# CLI: EasyPQP, SpectraST, DIA-NN
# Deep learning: Prosit (web), MS2PIP, DeepLC
```

## Quick Start
Tell your AI agent what you want to do:
- "Build a spectral library from my DDA search results"
- "Generate a predicted library using Prosit for my protein list"
- "Convert my library to DIA-NN format"

## Example Prompts

### Building Empirical Libraries
> "Create a spectral library from my MaxQuant msms.txt using EasyPQP"

> "Build a SpectraST library from my pepXML search results"

> "Combine libraries from multiple DDA experiments into a consensus library"

### Generating Predicted Libraries
> "Use Prosit to predict spectra for the human proteome"

> "Generate a predicted library with MS2PIP for my target protein list"

> "Create a DeepLC retention time library for my FASTA sequences"

### Hybrid Libraries
> "Merge my empirical library with Prosit predictions for missing peptides"

> "Supplement my DDA library with predicted spectra for low-abundance proteins"

### Format Conversion
> "Convert my SpectraST .splib to DIA-NN .tsv format"

> "Export my library in PQP format for OpenSWATH"

> "Transform the Spectronaut .kit library to a format compatible with DIA-NN"

### Quality Assessment
> "Check the number of precursors and transitions in my library"

> "Analyze the retention time coverage and charge state distribution"

> "Identify peptides missing from my library compared to the proteome"

## What the Agent Will Do
1. Load source data (search results, FASTA, or existing library)
2. Build or generate spectral library
3. Filter for quality (RT, intensity, charge)
4. Convert to target format
5. Generate library statistics

## Library Types

| Type | Source | Pros | Cons |
|------|--------|------|------|
| Empirical | DDA experiments | Highest quality spectra | Limited coverage |
| Predicted | Deep learning (Prosit) | Complete proteome | Slightly lower accuracy |
| Hybrid | Both combined | Best coverage + quality | More complex to build |

## Key Formats

| Format | Tool | Extension |
|--------|------|-----------|
| SpectraST | TPP | .splib |
| PQP | OpenMS/Skyline | .pqp |
| DLIB | EncyclopeDIA | .dlib |
| TSV | DIA-NN/OpenSWATH | .tsv |
| SSL | Spectronaut | .kit |

## Tips
- Empirical libraries from same sample type give best results
- Predicted libraries enable analysis without prior DDA data
- Aim for 6-10 transitions per precursor
- Check retention time coverage spans your gradient
- Library-free DIA-NN is a good alternative when no library exists
