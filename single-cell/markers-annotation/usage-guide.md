# Marker Gene Detection and Manual Annotation - Usage Guide

## Overview

This skill covers discovering cluster marker genes and hand-labeling clusters from canonical markers in single-cell RNA-seq using Seurat (R) and Scanpy (Python). Marker detection is descriptive ranking, not inference: it labels clusters, it does not prove a cluster is a real cell type, and it is distinct from cross-condition DE (which needs pseudobulk + biological replicates). Automated reference-based label transfer lives in single-cell/cell-annotation.

## Prerequisites

**Python (Scanpy):**
```bash
pip install scanpy pandas matplotlib
```

**R (Seurat):**
```r
install.packages(c('Seurat', 'dplyr'))
# Optional for better DE:
BiocManager::install('MAST')
```

## Quick Start

Ask your AI agent:
- "Find marker genes for each cluster"
- "What cell types are in my data?"
- "Show a dot plot of marker expression"

## Example Prompts

### Finding Markers
> "Find differentially expressed genes for cluster 0"

> "What genes distinguish cluster 1 from cluster 2?"

> "Show top 10 markers for each cluster"

### Visualization
> "Create a dot plot of these marker genes"

> "Show a heatmap of the top markers"

> "Plot CD3D expression on UMAP"

### Annotation
> "Annotate clusters based on these markers"

> "Score cells for T cell signature genes"

> "Label the clusters with cell type names"

### Export
> "Export all markers to CSV"

> "Save the top 20 markers per cluster"

## What the Agent Will Do

1. Run a marker ranking test (Wilcoxon for scanpy/Seurat; pass it explicitly)
2. Filter markers by specificity (effect size + in/out fraction), not p-value alone
3. Visualize marker expression patterns (dot plot, heatmap, UMAP)
4. Map clusters to cell type labels from canonical markers (manual labeling)
5. Add labels to the object and flag unmapped clusters
6. For cross-condition DE, aggregate to pseudobulk and hand off to differential-expression/deseq2-basics

## Tips

- **Marker detection is not condition DE** - one-cluster-vs-rest ranking uses Wilcoxon; treatment-vs-control needs pseudobulk + DESeq2/edgeR across biological replicates.
- **Rank markers by specificity, not p-value** - require a large positive log fold change and a high in-group / low out-group fraction; tiny p-values just reflect large n.
- **Cluster-marker p-values are double-dipping** - they label clusters, they do not prove a cluster is a real cell type; significance-test a suspicious split (scSHC/ClusterDE) before naming it.
- **scanpy defaults to t-test** - pass `method='wilcoxon'` explicitly; Seurat v5 thresholds dropped to logfc 0.1 / min.pct 0.01, so refilter.
- **Install presto** - Seurat's Wilcoxon silently falls back to slow base-R without it.
- **Markers are context-dependent** - re-validate any panel ported across tissue or condition; aggregate raw counts (never normalized) for pseudobulk.

## Related Skills

- clustering - Cluster cells before finding markers
- cell-annotation - Automated reference-based label transfer
- differential-abundance - Test whether cell-type proportions changed between conditions
- differential-expression/deseq2-basics - Pseudobulk condition DE engine
- pathway-analysis/go-enrichment - Functional interpretation of marker / DE gene lists
