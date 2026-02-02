# Outbreak Pipeline - Usage Guide

## Overview

This workflow investigates disease outbreaks using whole-genome sequencing data. It combines pathogen typing (MLST), antimicrobial resistance detection, phylodynamic analysis, and transmission inference to reconstruct outbreak dynamics and identify transmission chains.

## Prerequisites

```bash
conda install -c bioconda mlst abricate snippy iqtree fasttree gubbins

pip install treetime biopython pandas matplotlib seaborn

Rscript -e "install.packages('TransPhylo')"
```

**Required databases:**
- MLST schemes (auto-downloaded by mlst)
- AMR databases: NCBI AMRFinderPlus, CARD, ResFinder
- Reference genome for core alignment

## Quick Start

Tell your AI agent what you want to do:
- "Investigate this Salmonella outbreak with genomic data"
- "Build a transmission network from my isolate genomes"
- "Type my bacterial isolates and check for AMR genes"
- "Run phylodynamic analysis on my outbreak samples"

## Example Prompts

### Complete pipeline
> "Run outbreak investigation on my Klebsiella isolate genomes"

> "Build a transmission tree from my hospital outbreak samples"

### Typing and AMR
> "MLST type my E. coli isolates"

> "Screen my isolates for antibiotic resistance genes"

> "Compare AMR profiles across outbreak isolates"

### Phylodynamics
> "Build a time-scaled phylogeny from my outbreak with collection dates"

> "Estimate the clock rate for this outbreak"

### Transmission
> "Infer who infected whom in my outbreak"

> "Estimate R0 from my outbreak genomic data"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Isolate genomes | FASTA | Assembled genomes for each isolate |
| Collection dates | TSV | Sample names with dates (YYYY-MM-DD) |
| Reference genome | GenBank/FASTA | For core SNP alignment |
| Location data (optional) | TSV | Geographic coordinates for mapping |

## What the Agent Will Do

1. **MLST Typing** - Assigns sequence types to verify clonality
2. **AMR Detection** - Identifies resistance genes and mutations
3. **Core Alignment** - Extracts SNPs from conserved genome regions
4. **Phylogenetics** - Builds ML tree from core SNPs
5. **Phylodynamics** - Calibrates tree with collection dates
6. **Transmission Inference** - Reconstructs transmission network and R0

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Clock filter | 3 IQR | Remove temporal outliers |
| Generation time | 7 days | Mean serial interval (bacteria) |
| MCMC iterations | 10000 | TransPhylo sampling |
| Core coverage | 10x | Minimum for SNP calling |

## Pathogen-Specific Settings

| Pathogen | Generation Time | Clock Rate | Notes |
|----------|-----------------|------------|-------|
| Salmonella | 7-14 days | 1-3e-7 | Hospital outbreaks |
| E. coli | 7-10 days | 1-2e-7 | STEC outbreaks |
| Klebsiella | 7-14 days | 1-3e-7 | Hospital AMR |
| TB | 1-2 years | 0.5-1e-7 | Long generation time |
| SARS-CoV-2 | 4-5 days | 8e-4 | Rapid evolution |

## Tips

- **Verify clonality first**: MLST should show same or closely related STs
- **Remove recombination**: Use Gubbins for bacteria with high recombination
- **Check temporal signal**: Root-to-tip R2 > 0.5 indicates usable signal
- **Match generation time**: Use literature values for your pathogen
- **Integrate epi data**: Compare inferred transmission to known contacts
- **Multiple databases**: Run abricate with ncbi, card, and resfinder for complete AMR picture
