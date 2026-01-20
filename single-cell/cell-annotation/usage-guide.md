# Automated Cell Annotation Usage Guide

## Overview

Automated cell type annotation uses reference datasets or trained classifiers to consistently label cells, reducing manual effort and improving reproducibility.

## Tool Selection

| Tool | Strengths | Language | Reference |
|------|-----------|----------|-----------|
| CellTypist | Fast, many immune models | Python | Pre-trained models |
| SingleR | Correlation-based, flexible | R | Any reference dataset |
| Azimuth | Seurat integration, mapping | R | Curated references |
| scPred | SVM classifier, trainable | R | Train your own |

## Quick Start Prompts

- "Annotate my PBMC data using CellTypist immune model"
- "Run SingleR with Human Primary Cell Atlas reference"
- "Use Azimuth to annotate my lung scRNA-seq data"
- "Train a custom annotation model on my reference dataset"

## Workflow

1. **Process query data** - Normalize, find variable genes
2. **Choose reference** - Match tissue/species
3. **Run annotation** - Apply chosen method
4. **Filter by confidence** - Remove low-quality predictions
5. **Validate** - Check canonical markers

## Requirements

```bash
# CellTypist
pip install celltypist

# SingleR
BiocManager::install('SingleR')
BiocManager::install('celldex')

# Azimuth
remotes::install_github('satijalab/azimuth')

# scPred
devtools::install_github('powellgenomicslab/scPred')
```

## Key Considerations

- **Reference quality** directly determines annotation quality
- **Species/tissue match** - Use appropriate reference
- **Granularity** - Fine vs coarse labels trade-off
- **Confidence thresholds** - Filter unreliable predictions

## Related Skills

- **single-cell/clustering-annotation** - Manual annotation
- **single-cell/cell-communication** - Downstream analysis
