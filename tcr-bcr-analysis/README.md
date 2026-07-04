# tcr-bcr-analysis

## Overview

Analyze T-cell receptor (TCR) and B-cell receptor (BCR) repertoires from bulk or single-cell sequencing data for immunology research, vaccine development, and cancer immunotherapy.

**Tool type:** mixed | **Primary tools:** MiXCR, VDJtools, immunarch, Immcantation, scirpy, tcrdist3

## Skills

| Skill | Description |
|-------|-------------|
| mixcr-analysis | Preset-driven V(D)J alignment and clonotype assembly, with library-chemistry matching and UMI/single-cell handling |
| vdjtools-analysis | Depth-normalized diversity (Hill profiles) and overlap metrics, with immunarch equivalents |
| immcantation-analysis | BCR clonal clustering at a data-derived threshold, somatic hypermutation, selection, and lineage trees |
| scirpy-analysis | Single-cell TCR/BCR clonotype definition, chain-pairing QC, and gene-expression integration |
| repertoire-visualization | Rarefaction, spectratype, clonal tracking, overlap and network figures with correct-comparison guidance |
| specificity-annotation | Antigen-specificity database annotation and sequence clustering with generation-probability nulls |

## Example Prompts

- "Assemble TCR clonotypes from my FASTQ files with the right preset for my kit"
- "Compare repertoire diversity between conditions after normalizing for depth"
- "Find shared clonotypes between samples and test whether sharing exceeds chance"
- "Cluster my BCR sequences into clonal lineages and quantify somatic hypermutation"
- "Integrate 10x VDJ data with my scRNA-seq analysis and map clonal expansion onto cell states"
- "Track clonal expansion between timepoints"
- "Annotate my TCRs against VDJdb and cluster them by likely antigen specificity"

## Requirements

```bash
# MiXCR (4.x requires an activated academic/commercial license)
conda install -c bioconda mixcr

# VDJtools (requires Java) and immunarch (R, actively maintained alternative)
wget https://github.com/mikessh/vdjtools/releases/download/1.2.1/vdjtools-1.2.1.zip
unzip vdjtools-1.2.1.zip

# Immcantation (R)
install.packages(c('alakazam', 'shazam', 'scoper', 'tigger', 'dowser', 'immunarch'))

# scirpy (Python)
pip install scirpy mudata

# specificity annotation (Python)
pip install tcrdist3 olga
```

## Related Skills

- **single-cell** - scRNA-seq analysis with VDJ enrichment
- **immunoinformatics** - MHC binding and epitope-directed specificity prediction
- **phylogenetics** - Tree reconstruction concepts for B-cell lineages
