# Diversity Analysis - Usage Guide

## Overview

Diversity analysis characterizes microbial community structure through alpha (within-sample) and beta (between-sample) diversity metrics using phyloseq and vegan.

## Prerequisites

```bash
# R packages
BiocManager::install('phyloseq')
install.packages('vegan')
```

## Quick Start

Tell your AI agent what you want to do:
- "Calculate alpha and beta diversity for my microbiome data"
- "Run PERMANOVA to test for differences between treatment groups"

## Example Prompts

### Alpha Diversity
> "Calculate Shannon and observed richness for each sample"

> "Compare alpha diversity between control and treatment groups"

> "Generate rarefaction curves to assess sampling depth"

### Beta Diversity
> "Calculate Bray-Curtis distances and make a PCoA plot"

> "Run weighted UniFrac analysis on my phyloseq object"

> "Create an NMDS ordination colored by sample groups"

### Statistical Testing
> "Test if diversity differs between groups using Kruskal-Wallis"

> "Run PERMANOVA on my distance matrix with treatment as factor"

> "Check beta dispersion homogeneity before PERMANOVA"

## What the Agent Will Do

1. Optionally rarefy data to normalize sequencing depth
2. Calculate alpha diversity metrics per sample
3. Test alpha diversity differences between groups
4. Calculate beta diversity distance matrix
5. Generate ordination plots (PCoA, NMDS)
6. Run PERMANOVA for group differences
7. Test dispersion homogeneity

## Tips

- Rarefaction normalizes depth but discards data - consider alternatives
- UniFrac requires a phylogenetic tree
- NMDS stress <0.2 indicates good fit
- PERMANOVA R2 >0.2 indicates moderate effect
- Always check betadisper before interpreting PERMANOVA

## Alpha Diversity Metrics

| Metric | Measures | Interpretation |
|--------|----------|----------------|
| Observed | Richness | Number of taxa |
| Chao1 | Richness | Estimated total richness |
| Shannon | Richness + Evenness | Higher = more diverse |
| Simpson | Dominance | Higher = more even |
| Faith's PD | Phylogenetic diversity | Requires tree |

## Beta Diversity Metrics

| Metric | Type | When to Use |
|--------|------|-------------|
| Bray-Curtis | Abundance-based | General purpose |
| Jaccard | Presence/absence | Focus on composition |
| UniFrac (weighted) | Phylogenetic | Abundant taxa matter |
| UniFrac (unweighted) | Phylogenetic | Rare taxa matter |
| Aitchison | Compositional | Compositionally-aware |

## Statistical Tests

| Test | Use Case |
|------|----------|
| Kruskal-Wallis | Alpha diversity, >2 groups |
| Wilcoxon | Alpha diversity, 2 groups |
| PERMANOVA | Beta diversity group differences |
| betadisper | Test dispersion homogeneity |
