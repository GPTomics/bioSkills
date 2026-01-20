# Automated Cell Annotation - Usage Guide

## Overview

Automated cell type annotation uses reference datasets or trained classifiers to consistently label cells, reducing manual effort and improving reproducibility compared to marker-based annotation.

## Prerequisites

```bash
# CellTypist (Python)
pip install celltypist
```

```r
# SingleR
BiocManager::install('SingleR')
BiocManager::install('celldex')

# Azimuth
remotes::install_github('satijalab/azimuth')

# scPred
devtools::install_github('powellgenomicslab/scPred')
```

## Quick Start

Tell your AI agent what you want to do:
- "Annotate my PBMC data with cell types"
- "Use a reference dataset to label my cells"
- "Train a classifier on my annotated data"

## Example Prompts

### Automated Annotation
> "Run CellTypist with the immune_all_low model"
> "Use SingleR with Human Primary Cell Atlas reference"
> "Annotate my lung data with Azimuth"

### Reference Selection
> "What reference datasets are available for my tissue?"
> "Which CellTypist model should I use for PBMCs?"
> "Download the celldex reference for mouse brain"

### Quality Assessment
> "Show annotation confidence scores"
> "Which cells have low confidence predictions?"
> "Compare automated labels to my manual annotations"

### Custom Training
> "Train a scPred model on my annotated reference"
> "Create a custom CellTypist model from my data"
> "Apply my trained classifier to new data"

### Refinement
> "Re-annotate cluster 5 with finer labels"
> "Merge similar cell type labels"
> "Transfer labels to a new dataset"

## What the Agent Will Do

1. Normalize query data appropriately
2. Select matching reference (tissue/species)
3. Run annotation algorithm
4. Assess prediction confidence
5. Filter low-quality predictions
6. Add labels to cell metadata
7. Validate against canonical markers

## Tool Selection

| Tool | Strengths | Language | Reference |
|------|-----------|----------|-----------|
| CellTypist | Fast, many immune models | Python | Pre-trained models |
| SingleR | Correlation-based, flexible | R | Any reference dataset |
| Azimuth | Seurat integration, mapping | R | Curated references |
| scPred | SVM classifier, trainable | R | Train your own |

## Tips

- **Reference quality matters** - annotations are only as good as the reference
- **Match tissue and species** - use appropriate reference for your data
- **Check confidence scores** - filter predictions below threshold (e.g., 0.5)
- **Validate with markers** - confirm that marker genes match predicted types
- **Fine vs coarse labels** - choose granularity appropriate for your question
- **Combine methods** - consensus across tools increases confidence
