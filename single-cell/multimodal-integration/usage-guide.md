# Multimodal Integration Usage Guide

Multi-modal single-cell technologies measure multiple biological layers per cell, enabling deeper insights into cell states.

## Technologies

| Technology | Measures | Cells/Run |
|------------|----------|-----------|
| CITE-seq | RNA + ~200 proteins | 10,000+ |
| 10X Multiome | RNA + chromatin accessibility | 10,000+ |
| SHARE-seq | RNA + ATAC | ~10,000 |
| Visium | RNA + spatial location | ~5,000 spots |

## Key Concepts

### Weighted Nearest Neighbors (WNN)
- Combines information from multiple modalities
- Weights each modality per cell based on informativeness
- Creates unified cell-cell graph for clustering

### CLR Normalization (for proteins)
- Centered Log-Ratio transformation
- Standard for ADT/protein data
- Handles compositional nature of antibody capture

## Quick Start (Seurat CITE-seq)

```r
library(Seurat)

# Create object with both modalities
obj <- CreateSeuratObject(counts = rna_counts)
obj[['ADT']] <- CreateAssayObject(counts = adt_counts)

# Normalize both
obj <- NormalizeData(obj, assay = 'RNA')
obj <- NormalizeData(obj, assay = 'ADT', normalization.method = 'CLR', margin = 2)

# Dimension reduction
obj <- ScaleData(obj) %>% RunPCA(reduction.name = 'pca')
DefaultAssay(obj) <- 'ADT'
obj <- ScaleData(obj) %>% RunPCA(reduction.name = 'apca')

# WNN clustering
obj <- FindMultiModalNeighbors(obj, reduction.list = list('pca', 'apca'),
                                dims.list = list(1:30, 1:18))
obj <- FindClusters(obj, graph.name = 'wsnn')
obj <- RunUMAP(obj, nn.name = 'weighted.nn')
```

## Package Requirements

```r
# R/Seurat
install.packages('Seurat')
BiocManager::install('Signac')  # For ATAC

# Python
# pip install muon scanpy anndata
```

## Best Practices

1. **QC each modality separately** before integration
2. **Normalize appropriately** - different methods for RNA vs protein
3. **Check modality weights** - ensure both contribute
4. **Validate with known biology** - protein markers should match cell types
