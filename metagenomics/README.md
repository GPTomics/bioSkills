# metagenomics

## Overview

Taxonomic profiling of metagenomic data using Kraken2 for k-mer based classification and MetaPhlAn for marker gene profiling. Includes abundance estimation with Bracken and functional profiling with HUMAnN.

**Tool type:** cli | **Primary tools:** Kraken2, MetaPhlAn, Bracken, HUMAnN

## Skills

| Skill | Description |
|-------|-------------|
| kraken-classification | Taxonomic classification with Kraken2 |
| metaphlan-profiling | Marker gene profiling with MetaPhlAn |
| abundance-estimation | Species abundance with Bracken |
| metagenome-visualization | Visualize taxonomic profiles |
| functional-profiling | Pathway and gene family abundance with HUMAnN |
| amr-detection | Antimicrobial resistance with AMRFinderPlus, ResFinder |
| strain-tracking | Strain-level analysis with MASH, sourmash, inStrain |

## Example Prompts

- "Classify my metagenomic reads with Kraken2"
- "Run Kraken2 with the standard database"
- "What species are in my microbiome sample?"
- "Profile my metagenome with MetaPhlAn"
- "Get species-level abundances from my reads"
- "Run MetaPhlAn on paired-end reads"
- "Estimate species abundances with Bracken"
- "Convert Kraken2 output to abundance table"
- "Get genus-level abundance estimates"
- "Create a stacked bar chart of my samples"
- "Make a heatmap of species abundances"
- "Calculate alpha diversity metrics"
- "Profile functional potential with HUMAnN"
- "Get pathway abundances from my metagenome"

## Requirements

```bash
# Kraken2 and Bracken
conda install -c bioconda kraken2 bracken

# MetaPhlAn
conda install -c bioconda metaphlan

# HUMAnN
conda install -c bioconda humann
```

## Related Skills

- **read-qc** - Quality control before classification
- **sequence-io** - FASTQ handling
- **pathway-analysis** - Functional annotation
