# Doublet Detection - Usage Guide

## Overview

Doublet detection identifies and removes droplets containing two or more cells, which can create artificial intermediate populations and lead to false biological conclusions.

## Prerequisites

```bash
# Python
pip install scrublet scanpy
```

```r
# R
install.packages('Seurat')
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
BiocManager::install('scDblFinder')
```

## Quick Start

Tell your AI agent what you want to do:
- "Detect doublets in my scRNA-seq data"
- "Remove doublets before clustering"
- "Run Scrublet on my AnnData object"

## Example Prompts

### Detection
> "Run Scrublet to identify doublets"
> "Use scDblFinder to detect doublets in my Seurat object"
> "Run DoubletFinder with optimized parameters"

### Filtering
> "Remove predicted doublets from my data"
> "Show me which cells are doublets on the UMAP"
> "What percentage of cells are doublets?"

### Troubleshooting
> "The doublet score distribution is not bimodal, what should I do?"
> "I'm getting too many doublets detected, how do I adjust?"
> "Run doublet detection on each sample separately"

## What the Agent Will Do

1. Simulate artificial doublets from the data
2. Train classifier to distinguish doublets from singlets
3. Score each cell for doublet probability
4. Identify threshold for calling doublets
5. Flag or remove predicted doublets
6. Visualize doublet scores on UMAP

## Method Selection

| Method | Strengths | Language |
|--------|-----------|----------|
| Scrublet | Fast, simple | Python |
| DoubletFinder | Most widely used | R |
| scDblFinder | Fastest, often most accurate | R |

## Expected Doublet Rates

Use the 10X formula: ~0.8% per 1,000 cells loaded

| Cells | Rate |
|-------|------|
| 5,000 | 4% |
| 10,000 | 8% |
| 15,000 | 12% |

## Tips

- **Run before normalization** - doublet detection works best on raw or minimally processed data
- **Run per sample** - if samples are pooled, detect doublets in each separately
- **Expect bimodal distribution** - doublet scores should show two peaks
- **High gene counts often indicate doublets** - filter these first if doublet detection fails
- **Check intermediate populations** - doublets often appear between cell types
- **Validate with markers** - doublets may express markers of multiple cell types
