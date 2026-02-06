# Colocalization Analysis - Usage Guide

## Overview

Test whether two association signals at the same genomic locus are driven by the same causal variant using Bayesian colocalization. Distinguishes shared causality from coincidental LD overlap between GWAS and eQTL signals.

## Prerequisites

```r
install.packages('coloc')
install.packages('susieR')

# For multi-trait colocalization
install.packages('remotes')
remotes::install_github('jrs95/hyprcoloc')

# Visualization
install.packages(c('ggplot2', 'patchwork'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Test if my GWAS hit and the eQTL for gene X share a causal variant"
- "Run colocalization between GWAS and eQTL summary stats at this locus"
- "Check if three traits colocalize at the same locus using HyPrColoc"

## Example Prompts

### Standard Colocalization

> "I have GWAS and eQTL summary stats for chromosome 6 -- run coloc.abf and interpret the posterior probabilities"

> "Test colocalization between my disease GWAS and this gene's eQTL data"

### Multiple Causal Variants

> "There might be multiple signals at this locus -- run SuSiE-coloc instead of standard coloc"

> "Use coloc.susie with my LD matrix to test colocalization allowing multiple causal variants"

### Multi-Trait

> "Test whether a GWAS signal, eQTL, and sQTL all share the same causal variant at this locus"

### Visualization

> "Make a LocusCompare plot showing GWAS vs eQTL p-values at this region"

## What the Agent Will Do

1. Extract locus data (typically 1 Mb window around lead SNP)
2. Format summary statistics into coloc input lists
3. Run coloc.abf to compute posterior probabilities (H0-H4)
4. Interpret PP.H4 (shared variant) and PP.H3 (distinct variants)
5. Run sensitivity analysis across prior values
6. Optionally run SuSiE-coloc for multi-signal loci or HyPrColoc for multi-trait
7. Generate regional association and LocusCompare plots

## Tips

- **PP.H4 threshold** - PP.H4 > 0.8 is strong evidence for colocalization; 0.5-0.8 is suggestive
- **Sample size matters** - Small eQTL studies (N < 200) reduce power; larger is always better
- **Multiple signals** - If LD is complex, use SuSiE-coloc to avoid false H3 conclusions
- **Prior sensitivity** - Run sensitivity() to check if conclusions depend on prior choices
- **Both signals needed** - Both traits must have significant associations at the locus for coloc to be informative
- **LD matrix** - SuSiE-coloc requires LD from a matching ancestry reference panel

## Related Skills

- causal-genomics/mendelian-randomization - Causal effect estimation with genetic IVs
- causal-genomics/fine-mapping - Credible sets for prioritizing causal variants
- population-genetics/linkage-disequilibrium - LD reference panels
- differential-expression/deseq2-basics - eQTL data generation
