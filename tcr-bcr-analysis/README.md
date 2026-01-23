# tcr-bcr-analysis

## Overview

Analyze T-cell receptor (TCR) and B-cell receptor (BCR) repertoires from bulk or single-cell sequencing data for immunology research, vaccine development, and cancer immunotherapy.

**Tool type:** mixed | **Primary tools:** MiXCR, VDJtools, Immcantation, scirpy

## Skills

| Skill | Description |
|-------|-------------|
| mixcr-analysis | V(D)J alignment and clonotype assembly from sequencing data |
| vdjtools-analysis | Repertoire diversity metrics and sample comparison |
| immcantation-analysis | BCR phylogenetics and somatic hypermutation analysis |
| scirpy-analysis | Single-cell TCR/BCR analysis with Scanpy integration |
| repertoire-visualization | Clone tracking, diversity plots, and network visualization |

## Example Prompts

- "Assemble TCR clonotypes from my FASTQ files"
- "Calculate Shannon diversity and clonality for my repertoire"
- "Find shared clonotypes between samples"
- "Analyze BCR somatic hypermutation patterns"
- "Integrate VDJ data with my scRNA-seq analysis"
- "Track clonal expansion between timepoints"
- "Create a circos plot of V-J gene usage"

## Requirements

```bash
# MiXCR
conda install -c bioconda mixcr

# VDJtools (requires Java)
wget https://github.com/mikessh/vdjtools/releases/download/1.2.1/vdjtools-1.2.1.zip
unzip vdjtools-1.2.1.zip

# Immcantation (R)
install.packages(c('alakazam', 'shazam', 'tigger', 'dowser'))

# scirpy (Python)
pip install scirpy mudata
```

## Related Skills

- **single-cell** - scRNA-seq analysis with VDJ enrichment
- **metagenomics** - Immune profiling in metagenomic samples
- **phylogenetics** - Tree reconstruction concepts
