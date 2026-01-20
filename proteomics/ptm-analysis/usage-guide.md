# PTM Analysis - Usage Guide

## Overview
Identify, localize, and quantify post-translational modifications (phosphorylation, acetylation, ubiquitination, etc.) that regulate protein function.

## Prerequisites
```bash
pip install numpy pandas scipy
# R packages: BiocManager::install(c("MSstatsPTM", "PhosR"))
# CLI: MaxQuant (with PTM search), MSFragger
```

## Quick Start
Tell your AI agent what you want to do:
- "Analyze phosphorylation sites from my MaxQuant Phospho(STY)Sites.txt"
- "Find differentially regulated phosphosites between conditions"
- "Perform kinase enrichment analysis on my phosphoproteomics data"

## Example Prompts

### Site Identification
> "Load the MaxQuant Phospho(STY)Sites.txt and filter to class I sites (localization probability > 0.75)"

> "Extract modification sites with confident localization from my search results"

> "Summarize PTM sites per protein and identify multiply-modified proteins"

### Quantification
> "Normalize phosphosite intensities to total protein abundance"

> "Calculate site occupancy (modified / total) for each phosphorylation site"

> "Use MSstatsPTM to compare PTM levels adjusted for protein changes"

### Differential Analysis
> "Find phosphosites changing significantly after drug treatment"

> "Identify acetylation sites regulated by the histone deacetylase inhibitor"

> "Compare ubiquitination profiles between wild-type and mutant cells"

### Motif and Kinase Analysis
> "Run motif-x to find enriched sequence patterns around phosphosites"

> "Perform kinase enrichment analysis to identify active kinases"

> "Map my phosphosites to known kinase-substrate relationships"

## What the Agent Will Do
1. Load PTM site data from search engine output
2. Filter by localization probability (class I/II/III)
3. Normalize (optionally to protein level)
4. Perform differential analysis
5. Run motif/kinase enrichment
6. Generate site-level results and visualizations

## Common PTMs

| Modification | Sites | Mass Shift | Enrichment |
|--------------|-------|------------|------------|
| Phosphorylation | S, T, Y | +79.97 | TiO2, IMAC |
| Acetylation | K, N-term | +42.01 | Anti-acetyl antibody |
| Methylation | K, R | +14.02 | Anti-methyl antibody |
| Ubiquitination | K | +114.04 (GG) | Anti-K-GG antibody |
| Glycosylation | N, S, T | Variable | Lectin enrichment |

## Site Localization Classes
- **Class I (>0.75)**: Confident site assignment
- **Class II (0.50-0.75)**: Probable site
- **Class III (<0.50)**: Ambiguous site

## Tips
- Filter to class I sites for confident analysis
- Normalize to protein level to distinguish PTM changes from abundance changes
- Use MSstatsPTM for rigorous PTM vs protein statistical testing
- Check PhosphoSitePlus for known functions of regulated sites
