# scATAC-seq Analysis - Usage Guide

## Overview

Single-cell ATAC-seq measures chromatin accessibility at single-cell resolution, enabling identification of cell type-specific regulatory elements and transcription factor activity.

## Prerequisites

```r
# Signac
install.packages('Signac')
BiocManager::install(c('EnsDb.Hsapiens.v86', 'chromVAR', 'motifmatchr', 'JASPAR2020'))

# ArchR
devtools::install_github('GreenleafLab/ArchR')
```

```bash
# Python (SnapATAC2)
pip install snapatac2
```

## Quick Start

Tell your AI agent what you want to do:
- "Process my 10X scATAC-seq data"
- "Find cell type-specific peaks and motifs"
- "Integrate my scATAC with scRNA-seq"

## Example Prompts

### Data Loading
> "Load my 10X scATAC filtered_peak_bc_matrix.h5"
> "Create a ChromatinAssay with the fragments file"

### Quality Control
> "Calculate TSS enrichment and nucleosome signal"
> "Show QC violin plots for my scATAC data"
> "Filter cells with low TSS enrichment"

### Processing
> "Run TF-IDF normalization and LSI"
> "Cluster my scATAC data and generate UMAP"
> "Call peaks per cluster"

### Motif Analysis
> "Find enriched motifs in cluster-specific peaks"
> "Run chromVAR for transcription factor activity"
> "Which TFs are active in cluster 3?"

### Integration
> "Integrate scATAC with my scRNA-seq data"
> "Transfer cell type labels from RNA to ATAC"
> "Link peaks to genes"

## What the Agent Will Do

1. Load scATAC data with fragment files
2. Calculate QC metrics (TSS enrichment, nucleosome signal, fragments in peaks)
3. Filter low-quality cells
4. Normalize with TF-IDF
5. Run LSI for dimensionality reduction
6. Cluster and visualize with UMAP
7. Identify marker peaks and enriched motifs
8. Score transcription factor activity

## Tool Selection

| Tool | Best For |
|------|----------|
| Signac | Integration with Seurat scRNA-seq, familiar users |
| ArchR | Large datasets, memory efficiency |
| SnapATAC2 | Python users, very large datasets |

## Key Differences from scRNA-seq

| Aspect | scRNA-seq | scATAC-seq |
|--------|-----------|------------|
| Features | ~20,000 genes | ~100,000+ peaks |
| Sparsity | ~90% zeros | ~99% zeros |
| Normalization | Log normalize | TF-IDF |
| Dim reduction | PCA | LSI |
| Cell type markers | Gene expression | Accessibility, motifs |

## QC Metrics

| Metric | Good | Poor |
|--------|------|------|
| Unique fragments | > 3,000 | < 1,000 |
| TSS enrichment | > 4 | < 2 |
| Fraction in peaks | > 30% | < 15% |
| Nucleosome signal | < 2 | > 4 |

## Tips

- **Skip first LSI component** - it correlates with sequencing depth
- **Use gene activity scores** for initial cell type annotation
- **Fragment file required** - needed for coverage plots and peak calling
- **TF-IDF before LSI** - standard normalization for sparse chromatin data
- **Validate with RNA** - integrate with scRNA-seq to confirm cell types
- **Motifs need background** - use matched GC content for enrichment testing
