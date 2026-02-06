# Non-Coding RNA Annotation - Usage Guide

## Overview

Identify and annotate non-coding RNAs (tRNAs, rRNAs, snoRNAs, riboswitches, miRNAs, and other regulatory RNAs) in genome assemblies using Infernal covariance model searches against the Rfam database and tRNAscan-SE for specialized tRNA prediction.

## Prerequisites

```bash
# Infernal (covariance model search)
conda install -c bioconda infernal

# Download Rfam database
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz && cmpress Rfam.cm
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin

# tRNAscan-SE
conda install -c bioconda trnascan-se

# barrnap (rRNA, legacy)
conda install -c bioconda barrnap

# Python utilities
pip install pandas biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Find all non-coding RNAs in my genome assembly"
- "Annotate tRNAs in my bacterial genome"
- "Search my genome for riboswitches and snoRNAs"

## Example Prompts

### General ncRNA Annotation

> "Run Infernal cmscan against Rfam to find all ncRNAs in my assembly"

> "Annotate non-coding RNAs in my eukaryotic genome with Infernal and tRNAscan-SE"

### tRNA-Specific

> "Find all tRNA genes and their anticodons in my bacterial genome"

> "Run tRNAscan-SE in eukaryotic mode and report isotype counts"

### Combined Analysis

> "Combine Infernal and tRNAscan-SE results into a single GFF file"

> "How many tRNAs, rRNAs, and snoRNAs are in my genome?"

## What the Agent Will Do

1. Set up Rfam database if not present (cmpress)
2. Run Infernal cmscan with Rfam gathering thresholds
3. Run tRNAscan-SE with domain-appropriate mode
4. Combine results, preferring tRNAscan-SE for tRNA calls
5. Report ncRNA counts by category
6. Generate combined GFF3 output

## Tips

- **Use --cut_ga** - Rfam gathering thresholds are family-specific and well-calibrated; prefer over arbitrary E-value cutoffs
- **tRNAscan-SE for tRNAs** - Always use tRNAscan-SE over Infernal for tRNA detection; it is more sensitive and specific
- **Domain matters** - Use the correct tRNAscan-SE mode (-B, -A, -E) for accurate results
- **Clan resolution** - Supply `--clanin Rfam.clanin` to resolve overlapping hits from related Rfam families
- **Speed** - cmscan is slow on large genomes; use `--rfam` flag and consider parallelization
- **Pseudogenes** - Eukaryotic genomes have many tRNA pseudogenes; filter by score > 50 bits

## Related Skills

- genome-annotation/prokaryotic-annotation - Bakta includes ncRNA annotation
- genome-annotation/eukaryotic-gene-prediction - Gene prediction excludes ncRNAs
