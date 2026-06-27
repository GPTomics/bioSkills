# Automated Cell Annotation - Usage Guide

## Overview

Automated reference-based annotation transfers cell type labels from an annotated reference or pretrained model onto query cells. It is closed-world: every cell is forced toward the nearest reference label, so the output is a hypothesis to be triangulated with markers and triaged for artifacts, never a measurement. This skill covers CellTypist, SingleR, Azimuth, scANVI/scArches, and scmap; manual marker-based labeling lives in single-cell/markers-annotation.

## Prerequisites

```bash
# CellTypist and scANVI (Python)
pip install scanpy celltypist scvi-tools
```

```r
# SingleR, scmap (Bioconductor)
BiocManager::install(c('SingleR', 'celldex', 'scmap'))

# Azimuth
remotes::install_github('satijalab/azimuth')
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

### Atlas Mapping
> "Map my query onto a reference atlas with scANVI"
> "Create a custom CellTypist model from my annotated data"
> "Use scmap with an explicit unassigned category"

### Refinement
> "Re-annotate cluster 5 with finer labels"
> "Merge similar cell type labels"
> "Transfer labels to a new dataset"

## What the Agent Will Do

1. Normalize the query to each tool's required input (CellTypist=CP10K log1p, scANVI=raw counts)
2. Select a reference matching tissue, species, and platform
3. Run the annotation method and transfer labels
4. Assess per-cell confidence and apply a calibrated rejection threshold
5. Triage poorly-mapped clusters (doublet / low-quality / batch / ambient) before claiming novelty
6. Add labels to cell metadata and flag disagreements as ambiguous
7. Validate against canonical markers (triangulation, not proof)

## Tool Selection

| Method | Model | Language | Reference | World |
|--------|-------|----------|-----------|-------|
| CellTypist | Logistic regression (pretrained) | Python | Pretrained immune/cross-tissue models | Closed (+probability) |
| SingleR | Spearman correlation | R | celldex bulk or single-cell refs | Closed (+pruning) |
| Azimuth | Supervised PCA + anchors | R | Curated Seurat atlases | Closed (+mapping.score) |
| scANVI/scArches | Semi-supervised VAE | Python | Annotated atlas + raw counts | Closed (+latent uncertainty) |
| scmap | Nearest centroid/cell | R | Single-cell reference | Open (explicit unassigned) |
| LLM (GPTCelltype) | Prompted from markers | R/Python | None (uses marker list) | Open-ish |

## Tips

- **Annotation is a hypothesis** - reference-based methods are closed-world and force every cell to the nearest label, often confidently wrong for novel states.
- **Match the input normalization** - CellTypist needs CP10K log1p, scANVI needs raw counts; wrong input degrades labels silently with no error.
- **Triage before claiming novelty** - rule out doublets, low-quality, batch, and ambient RNA before naming an unexpected cluster a new cell type.
- **Calibrate rejection per dataset** - inspect score/delta distributions; a hard universal probability cutoff is not principled across models.
- **Marker confirmation is circular** - canonical markers are a prior, not independent proof; triangulate automated + markers + expert curation.
- **Annotate hierarchically** - report coarse compartments confidently, fine subtypes as hypotheses sensitive to reference and resolution.

## Related Skills

- markers-annotation - Manual marker discovery and hand-labeling
- clustering - Cluster cells before annotating
- batch-integration - Reference mapping vs de-novo integration
- differential-abundance - Test whether annotated proportions changed between conditions
- pathway-analysis/go-enrichment - Characterize a de-novo / novel population
