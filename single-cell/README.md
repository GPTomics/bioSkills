# Single-Cell RNA-seq

Single-cell RNA-seq analysis using Seurat (R) and Scanpy (Python) for preprocessing, clustering, and cell type annotation.

## Overview

This category covers the complete single-cell RNA-seq workflow: loading data, quality control, normalization, dimensionality reduction, clustering, finding marker genes, and annotating cell types. Both Seurat (R) and Scanpy (Python) are supported, allowing users to choose their preferred language.

**Tool type:** `mixed`
**Primary tools:** Seurat (R), Scanpy/AnnData (Python)

## Skills

| Skill | Description |
|-------|-------------|
| [data-io](data-io/) | Load 10X data, create Seurat/AnnData objects, read/write h5ad/RDS |
| [preprocessing](preprocessing/) | QC metrics, filtering, normalization, HVGs, scaling |
| [clustering](clustering/) | PCA, neighbors, Leiden/Louvain clustering, UMAP/tSNE |
| [markers-annotation](markers-annotation/) | Differential expression, marker genes, cell type annotation |

## Workflow

```
10X Genomics / Count Matrix
    |
    v
[data-io] ----------> Create Seurat or AnnData object
    |
    v
[preprocessing] ----> QC filtering, normalization, HVGs
    |
    v
[clustering] -------> PCA -> Neighbors -> Leiden -> UMAP
    |
    v
[markers-annotation] -> FindMarkers -> Annotate cell types
    |
    v
Annotated Dataset (h5ad / RDS)
```

## Tool Comparison

| Step | Scanpy (Python) | Seurat (R) |
|------|-----------------|------------|
| Data object | AnnData | Seurat object |
| Read 10X | `sc.read_10x_mtx()` | `Read10X()` |
| QC metrics | `calculate_qc_metrics()` | `PercentageFeatureSet()` |
| Normalize | `normalize_total()` + `log1p()` | `NormalizeData()` or `SCTransform()` |
| HVGs | `highly_variable_genes()` | `FindVariableFeatures()` |
| PCA | `sc.tl.pca()` | `RunPCA()` |
| Neighbors | `sc.pp.neighbors()` | `FindNeighbors()` |
| Cluster | `sc.tl.leiden()` | `FindClusters()` |
| UMAP | `sc.tl.umap()` | `RunUMAP()` |
| Markers | `rank_genes_groups()` | `FindAllMarkers()` |
| Native format | `.h5ad` | `.rds` |

## Example Prompts

### Data I/O
- "Load my 10X data into Scanpy"
- "Create a Seurat object from this count matrix"
- "Convert h5ad to Seurat format"

### Preprocessing
- "Run QC and filter cells with >20% mitochondrial"
- "Normalize and find highly variable genes"
- "Preprocess this data using SCTransform"

### Clustering
- "Run PCA and cluster at resolution 0.5"
- "Generate a UMAP colored by cluster"
- "Try different clustering resolutions"

### Markers & Annotation
- "Find marker genes for each cluster"
- "What cell types are these clusters?"
- "Show a dot plot of canonical PBMC markers"

## Requirements

**Python (Scanpy):**
```bash
pip install scanpy anndata leidenalg matplotlib
```

**R (Seurat):**
```r
install.packages('Seurat')
# For format conversion:
remotes::install_github('mojaveazure/seurat-disk')
```

## Notes

- **Seurat v5 uses layers** - Use `LayerData()` instead of `GetAssayData(slot=)`
- **SCTransform v2** is default in Seurat v5 - use `vst.flavor='v1'` for old behavior
- **Leiden preferred** over Louvain - better performance and consistency
- **Scanpy deprecations** - Use `highly_variable_genes()` not `filter_genes_dispersion()`
- **Store raw counts** before normalization for DE analysis

## Deprecations Avoided

### Seurat v5
- `slot` argument deprecated -> use `layer`
- `GetAssayData(slot='counts')` -> `LayerData(layer='counts')`
- SCTransform v1 requires explicit `vst.flavor='v1'`

### Scanpy
- `pp.filter_genes_dispersion` -> `pp.highly_variable_genes`
- `pp.normalize_per_cell` -> `pp.normalize_total`
- `pp.subsample` -> `pp.sample`
- `tl.louvain` -> `tl.leiden` (louvain still works but leiden preferred)

## Related Skills

- **sequence-io** - For working with FASTA/FASTQ before alignment
- **differential-expression** - Bulk RNA-seq DE analysis (DESeq2, edgeR)
- **pathway-analysis** - GO/KEGG enrichment of marker genes

## References

- [Scanpy documentation](https://scanpy.readthedocs.io/)
- [Seurat documentation](https://satijalab.org/seurat/)
- [Scanpy PBMC tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html)
- [Seurat PBMC tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
