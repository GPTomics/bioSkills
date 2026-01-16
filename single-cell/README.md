# single-cell

## Overview

Single-cell RNA-seq analysis using Seurat (R) and Scanpy (Python). Covers the complete workflow from loading data through quality control, normalization, clustering, marker gene identification, and cell type annotation.

**Tool type:** mixed | **Primary tools:** Seurat (R), Scanpy/AnnData (Python)

## Skills

| Skill | Description |
|-------|-------------|
| data-io | Load 10X data, create Seurat/AnnData objects, read/write h5ad/RDS |
| preprocessing | QC metrics, filtering, normalization, HVGs, scaling |
| clustering | PCA, neighbors, Leiden/Louvain clustering, UMAP/tSNE |
| markers-annotation | Differential expression, marker genes, cell type annotation |
| multimodal-integration | CITE-seq, Multiome, WNN clustering for multi-modal data |

## Example Prompts

- "Load my 10X data into Scanpy"
- "Create a Seurat object from this count matrix"
- "Convert h5ad to Seurat format"
- "Run QC and filter cells with >20% mitochondrial"
- "Normalize and find highly variable genes"
- "Preprocess this data using SCTransform"
- "Run PCA and cluster at resolution 0.5"
- "Generate a UMAP colored by cluster"
- "Try different clustering resolutions"
- "Find marker genes for each cluster"
- "What cell types are these clusters?"
- "Show a dot plot of canonical PBMC markers"
- "Load my CITE-seq data with both RNA and ADT"
- "Run WNN clustering combining RNA and protein"
- "Analyze my 10X Multiome data"

## Requirements

```bash
# Python (Scanpy)
pip install scanpy anndata leidenalg matplotlib

# R (Seurat)
install.packages('Seurat')
```

## Related Skills

- **differential-expression** - Bulk RNA-seq DE analysis
- **pathway-analysis** - GO/KEGG enrichment of marker genes
- **rna-quantification** - Count matrix generation
