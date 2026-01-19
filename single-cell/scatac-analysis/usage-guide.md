# scATAC-seq Analysis - Usage Guide

## Overview

Single-cell ATAC-seq measures chromatin accessibility at single-cell resolution, enabling identification of cell type-specific regulatory elements and transcription factor activity.

## When to Use This Skill

- Processing 10X Genomics scATAC-seq data
- Identifying cell types by chromatin accessibility
- Finding cell type-specific regulatory elements
- Scoring transcription factor activity per cell
- Integrating scATAC with scRNA-seq data

## Installation

```r
# Signac
install.packages('Signac')
BiocManager::install(c('EnsDb.Hsapiens.v86', 'chromVAR', 'motifmatchr', 'JASPAR2020'))

# ArchR
devtools::install_github('GreenleafLab/ArchR')
```

## Tool Selection

| Tool | Best For |
|------|----------|
| Signac | Integration with Seurat scRNA-seq, familiar users |
| ArchR | Large datasets, memory efficiency |
| SnapATAC2 | Python users, very large datasets |

## Quick Start (Signac)

```r
library(Signac)
library(Seurat)

# Load data
counts <- Read10X_h5('filtered_peak_bc_matrix.h5')
obj <- CreateSeuratObject(counts = CreateChromatinAssay(counts = counts, fragments = 'fragments.tsv.gz'))

# QC
obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)
obj <- subset(obj, peak_region_fragments > 3000 & TSS.enrichment > 2)

# Process
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)
obj <- RunUMAP(obj, reduction = 'lsi', dims = 2:30)
obj <- FindNeighbors(obj, reduction = 'lsi', dims = 2:30)
obj <- FindClusters(obj, algorithm = 3)
```

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

## Resources

- [Signac Vignettes](https://satijalab.org/signac/)
- [ArchR Book](https://www.archrproject.com/bookdown/)
- [chromVAR Paper](https://doi.org/10.1038/nmeth.4401)
