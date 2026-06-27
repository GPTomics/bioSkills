# Doublet Detection - Usage Guide

## Overview

This skill identifies and removes droplets that captured two or more cells, which masquerade as artificial intermediate populations and corrupt clustering and trajectory inference. It covers scDblFinder (R, the current best-balance default), Scrublet (Python), and DoubletFinder (R, legacy Seurat), and the decisions around per-sample detection, expected-rate setting, and avoiding over-removal.

## Prerequisites

```bash
# Python
pip install scanpy        # provides sc.pp.scrublet (the maintained Scrublet path)
```

```r
# R
BiocManager::install('scDblFinder')
install.packages('Seurat')
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
```

## Quick Start

Tell the AI agent what is needed:
- "Run scDblFinder per sample and flag doublets before integration"
- "Detect doublets with Scrublet, setting the expected rate from my cell count"
- "My doublet score histogram is not bimodal - what should I do?"

## Example Prompts

### Detection
> "Run scDblFinder on my Seurat object, processing each sample separately"

> "Run Scrublet on this AnnData with the expected rate set from recovered cells"

> "Run DoubletFinder with a pK sweep and homotypic-adjusted nExp"

### Filtering and Interpretation
> "Flag doublets but keep the score so I can inspect suspect clusters"

> "Is this CD3+LYZ intermediate cluster real or a doublet artifact?"

> "Why might my trajectory have a bridge between two unrelated lineages?"

### Troubleshooting
> "I detected 0% doublets - is that plausible?"

> "Doublet removal deleted my megakaryocytes - how do I avoid that?"

## What the Agent Will Do

1. Confirm detection runs per sample on raw counts, before integration or clustering
2. Set the expected doublet rate from the recovered-cell count of each lane
3. Run the chosen method (scDblFinder by default; Scrublet for scanpy; DoubletFinder for legacy Seurat)
4. Inspect the score distribution and set a manual threshold if it is not bimodal
5. Flag rather than blindly delete; cross-check high-score cells against cell-cycle and activation signatures
6. Coordinate doublet removal with count-based QC to avoid double-penalizing high-RNA cells

## Method Selection

| Method | Strengths | Use when |
|--------|-----------|----------|
| scDblFinder | Fast, accurate, built-in per-sample handling | Default choice; R/Bioconductor pipelines |
| Scrublet | Simple, scanpy-native | Python workflows |
| DoubletFinder | Widely used historically | Legacy Seurat pipelines |

## Expected Doublet Rates

Set the rate from recovered cells, not a package default: rate ~= 0.008 x cells / 1000.

| Cells recovered | Rate |
|-----------------|------|
| 5,000 | ~3.9% |
| 10,000 | ~7.6-8% |
| 15,000 | ~12% |

## Tips

- **Detect per sample before integration** - cross-sample doublets are physically impossible and merging corrupts the scoring neighborhood.
- **Run on raw counts after basic QC** - detectors simulate doublets from counts.
- **Set the rate from recovered cells** - the flat 0.05 in Scrublet is a placeholder, not a recommendation.
- **scDblFinder is the current default** - the old "DoubletFinder is most accurate" ranking is superseded.
- **Homotypic doublets are invisible** - removal is never complete; do not claim a doublet-free dataset.
- **Flag, then inspect** - prefer keeping the score and examining high-score clusters over blind deletion.
- **Treat intermediate clusters as suspect** - any cluster co-expressing two lineage programs is doublet-suspect until ruled out.
- **Do not double-penalize** - doublet score correlates with total counts; coordinate with count-based QC to keep real high-RNA cells.
- **Use hashing as ground truth** - cell hashing and MULTI-seq call inter-sample doublets experimentally, regardless of expression similarity.

## Related Skills

- single-cell/preprocessing - QC and ambient-RNA handling before doublet detection
- single-cell/data-io - load raw per-sample matrices before processing
- single-cell/clustering - run clustering after doublet removal
- single-cell/batch-integration - integrate samples only after per-sample doublet calling
- single-cell/trajectory-inference - doublets create false bridges; remove them first
