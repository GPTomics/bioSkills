# AMR Detection - Usage Guide

## Overview
Identify antimicrobial resistance genes in bacterial genomes and metagenomes using curated databases like NCBI AMRFinderPlus, ResFinder, and CARD.

## Prerequisites
```bash
# AMRFinderPlus (recommended)
conda install -c bioconda ncbi-amrfinderplus
amrfinder -u  # Update database

# ResFinder
pip install resfinder

# ABRicate (quick multi-database screening)
conda install -c bioconda abricate
abricate --setupdb
```

## Quick Start
Tell your AI agent what you want to do:
- "Screen my genome assembly for resistance genes"
- "Identify AMR genes in my metagenome contigs"
- "Compare resistance profiles across isolates"

## Example Prompts
### Basic Screening
> "Run AMRFinderPlus on contigs.fasta to find resistance genes"

> "Screen my assembly against the CARD database using ABRicate"

### Comprehensive Analysis
> "Run AMRFinderPlus with --plus to include virulence and stress genes"

> "Screen my genome against ResFinder, CARD, and NCBI databases and compare results"

### Metagenome Analysis
> "Identify AMR genes in my metagenomic assembly"

> "Find resistance genes in my assembled MAGs"

### Batch Processing
> "Run AMR detection on all fasta files in this directory and summarize results"

> "Create a presence/absence matrix of resistance genes across all isolates"

## What the Agent Will Do
1. Select appropriate tool and database for the analysis
2. Run resistance gene detection on assemblies or protein sequences
3. Parse output to identify gene classes and drug targets
4. Summarize findings across multiple samples if requested

## Tips
- AMRFinderPlus is most curated; use for clinical/publication work
- ABRicate is fast for screening multiple databases
- Always update databases before analysis (`amrfinder -u`, `abricate --setupdb`)
- Use `--plus` flag with AMRFinderPlus for virulence and stress genes
- Check % coverage and % identity to assess match quality

## Database Comparison

| Database | Focus | Best For |
|----------|-------|----------|
| NCBI (AMRFinderPlus) | Curated, comprehensive | General use |
| ResFinder | Acquired resistance | Clinical isolates |
| CARD | Mechanism-focused | Research |
| ARG-ANNOT | Acquired genes | Surveillance |

## Resistance Mechanisms

| Type | Example | Implication |
|------|---------|-------------|
| Acquired | blaCTX-M | Horizontal transfer |
| Mutational | gyrA | Chromosomal |
| Efflux | mexAB-oprM | Broad spectrum |

## Resources
- [AMRFinderPlus wiki](https://github.com/ncbi/amr/wiki)
- [CARD database](https://card.mcmaster.ca/)
- [ResFinder web server](https://cge.food.dtu.dk/services/ResFinder/)
