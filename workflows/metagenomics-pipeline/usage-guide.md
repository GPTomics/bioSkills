# Metagenomics Pipeline - Usage Guide

## Overview

This workflow processes metagenomic sequencing data to produce taxonomic and functional profiles. It supports both Kraken2/Bracken and MetaPhlAn approaches.

## Prerequisites

```bash
# CLI tools
conda install -c bioconda kraken2 bracken metaphlan humann bowtie2 fastp

# Download databases
kraken2-build --download-taxonomy --db kraken2_db
kraken2-build --download-library bacteria --db kraken2_db
kraken2-build --build --db kraken2_db
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the metagenomics pipeline on my shotgun metagenomic data"
- "Classify my metagenomic reads and get species abundances"
- "Profile pathways in my gut microbiome samples"

## Example Prompts

### Taxonomic profiling
> "Classify my metagenomic reads with Kraken2"

> "Use MetaPhlAn to profile my samples"

> "Get species-level abundances with Bracken"

### Functional profiling
> "Run HUMAnN on my metagenome"

> "Get pathway abundances for my samples"

### Visualization
> "Create a stacked bar plot of species abundances"

> "Make a heatmap of the top taxa"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| FASTQ files | .fastq.gz | Paired-end metagenomic reads |
| Kraken2 database | Directory | Pre-built or custom database |
| Host reference | FASTA | For host read removal |

## What the Workflow Does

1. **Quality Control** - Filter low-quality reads
2. **Host Removal** - Remove host contamination
3. **Classification** - Assign taxonomic labels
4. **Abundance** - Estimate species/genus abundances
5. **Functional** - Profile metabolic pathways

## Kraken2/Bracken vs MetaPhlAn

| Feature | Kraken2/Bracken | MetaPhlAn |
|---------|-----------------|-----------|
| Speed | Very fast | Moderate |
| Database | Large, customizable | Marker genes |
| Accuracy | Database-dependent | Standardized |
| Best for | Large-scale screening | Publication-ready profiles |

## Tips

- **Host removal**: Critical for human-associated samples
- **Database choice**: Standard for general use, custom for specific environments
- **Depth**: Aim for >5M reads per sample for reliable profiling
- **Replicates**: Important for differential abundance analysis
