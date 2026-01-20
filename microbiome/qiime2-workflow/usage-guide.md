# QIIME2 Workflow - Usage Guide

## Overview

QIIME2 is a comprehensive microbiome analysis platform with built-in provenance tracking, a plugin architecture, and visualization tools, providing an alternative to the DADA2/phyloseq R workflow.

## Prerequisites

```bash
# Install QIIME2 via conda
conda env create -n qiime2 --file https://data.qiime2.org/distro/core/qiime2-2024.5-py39-linux-conda.yml
conda activate qiime2
```

## Quick Start

Tell your AI agent what you want to do:
- "Process my 16S data through the full QIIME2 pipeline"
- "Run diversity analysis in QIIME2 with my metadata"

## Example Prompts

### Import and Denoising
> "Import my paired-end FASTQ files into QIIME2"

> "Run DADA2 denoising in QIIME2 with truncation at 240 and 160"

### Taxonomy
> "Classify features against the SILVA database in QIIME2"

> "Train a QIIME2 classifier for my V4 region primers"

### Diversity
> "Run core diversity analysis at sampling depth of 10000"

> "Generate alpha rarefaction curves to choose sampling depth"

### Visualization
> "Create an emperor PCoA plot colored by treatment group"

> "Generate a taxonomy barplot from my feature table"

## What the Agent Will Do

1. Import raw sequences as QIIME2 artifacts
2. Run quality filtering and denoising (DADA2/Deblur)
3. Assign taxonomy with trained classifier
4. Build phylogenetic tree
5. Calculate alpha and beta diversity metrics
6. Run statistical tests on diversity
7. Generate interactive visualizations

## Tips

- Artifacts (.qza) contain data + provenance tracking
- Visualizations (.qzv) viewable at view.qiime2.org
- Check rarefaction curves before choosing sampling depth
- Use ITSxpress plugin for fungal ITS data
- Export to phyloseq for custom R analysis

## When to Use QIIME2 vs R

| QIIME2 | DADA2/phyloseq |
|--------|----------------|
| Large studies (parallelization) | Custom R analyses |
| Reproducibility requirements | Bioconductor integration |
| Command-line preference | Interactive RStudio |
| Plugin ecosystem | Existing R workflows |

## Key Concepts

### Artifacts (.qza)
Binary files containing data + metadata + provenance. Can be visualized at view.qiime2.org

### Visualizations (.qzv)
Interactive HTML visualizations viewable in browser or at view.qiime2.org

### Plugins
Modular components: dada2, deblur, feature-classifier, diversity, etc.

## Common Workflows

### 16S V4 Region
```bash
qiime dada2 denoise-paired \
    --p-trunc-len-f 240 \
    --p-trunc-len-r 160 \
    ...
```

### ITS Fungal
```bash
# Extract ITS region first
qiime itsxpress trim-pair-output-unmerged \
    --p-region ITS2 --p-taxa F ...

# Then denoise without truncation
qiime dada2 denoise-paired \
    --p-trunc-len-f 0 --p-trunc-len-r 0 ...
```
