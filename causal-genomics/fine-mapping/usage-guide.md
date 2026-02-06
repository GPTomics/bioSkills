# Fine-Mapping - Usage Guide

## Overview

Identify likely causal variants within GWAS loci by computing posterior inclusion probabilities (PIPs) and constructing credible sets. SuSiE handles multiple causal variants and integrates with colocalization analysis.

## Prerequisites

```r
install.packages('susieR')
install.packages(c('ggplot2', 'patchwork'))

# FINEMAP (command-line tool)
# Download from http://www.christianbenner.com/
```

```bash
# For LD matrix generation
# plink v1.9+ required
conda install -c bioconda plink
```

## Quick Start

Tell your AI agent what you want to do:
- "Fine-map this GWAS locus to a 95% credible set"
- "Compute PIPs for all variants at my GWAS hit using SuSiE"
- "Run FINEMAP on this locus with the LD matrix from 1000 Genomes"

## Example Prompts

### SuSiE Fine-Mapping

> "I have GWAS summary stats and an LD matrix for a locus on chromosome 6 -- run SuSiE and report the credible sets"

> "Fine-map my GWAS signal using susie_rss with L = 10"

### FINEMAP

> "Set up FINEMAP input files for this locus and run the stochastic search"

> "Compare SuSiE and FINEMAP results at my top GWAS locus"

### Visualization

> "Make a PIP plot alongside the GWAS Manhattan for this locus"

> "Color the regional plot by credible set membership"

### LD Reference

> "Generate an LD matrix from 1000 Genomes Europeans for this 1 Mb region"

## What the Agent Will Do

1. Extract locus data (typically 1 Mb window around lead SNP)
2. Obtain or compute LD matrix from matched ancestry reference
3. Run SuSiE (susie_rss) with summary statistics and LD
4. Extract 95% credible sets and PIPs
5. Report credible set size, purity, and top PIP variants
6. Optionally run FINEMAP for comparison
7. Generate PIP and regional association plots

## Tips

- **L parameter** - Start with L = 10; SuSiE automatically prunes unused effects
- **LD ancestry match** - The LD reference must match the GWAS sample ancestry; mismatched LD causes false signals
- **Credible set purity** - Minimum absolute correlation > 0.5 indicates a well-resolved signal
- **PIP thresholds** - PIP > 0.95 is strong evidence; PIP > 0.5 is suggestive
- **FINEMAP comparison** - When SuSiE and FINEMAP agree, results are more reliable
- **Positive semi-definite** - Add a small ridge to the LD matrix diagonal if eigenvalues are negative

## Related Skills

- causal-genomics/colocalization-analysis - SuSiE-coloc for shared causal variants
- causal-genomics/mendelian-randomization - Fine-map instrument loci
- population-genetics/linkage-disequilibrium - LD reference panels
- variant-calling/variant-annotation - Annotate fine-mapped variants
