---
name: bio-workflows-scrnaseq-pipeline
description: Orchestrates the end-to-end single-cell RNA-seq pipeline from 10x Cell Ranger output to annotated cell types, chaining ambient-RNA removal, doublet detection, MAD-adaptive QC, normalization, integration, clustering, marker annotation, and (separately) pseudobulk DE + differential abundance. Use when honoring the made-once counting commitments (reference build/Ensembl vintage, --include-introns, cell-calling, feature namespace), ordering correct-then-detect-then-normalize (ambient before doublet before normalize), running per-sample QC before merge, integrating-then-clustering (never testing on integrated values), aggregating to pseudobulk for condition DE instead of cells-as-replicates, or pairing DE with differential abundance. Hands mechanism to the single-cell component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: Seurat
goal_approach_exempt: true
workflow: true
depends_on:
  - single-cell/data-io
  - single-cell/preprocessing
  - single-cell/doublet-detection
  - single-cell/clustering
  - single-cell/markers-annotation
qc_checkpoints:
  - after_loading: "Expected cell count, reasonable UMI distribution"
  - after_qc: "Remove low-quality cells and doublets"
  - after_normalization: "No batch effects, HVGs look sensible"
  - after_clustering: "Clusters are biologically meaningful"
---

## Version Compatibility

Reference examples tested with: Cell Ranger 10.0+ (human/mouse refs on Ensembl v110), Seurat 5.1+, Scanpy 1.10+, scDblFinder 1.16+, SoupX 1.6+, CellBender 0.3+, harmony/scvi-tools current, ggplot2 3.5+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: the Cell Ranger `filtered_feature_bc_matrix` is cell-CALLING only, NOT ambient-corrected; recovering low-RNA cells or running SoupX/CellBender needs the RAW matrix. `--include-introns` is default TRUE since Cell Ranger 7 (essential for nuclei). `gex_only=True`/default silently drops Antibody Capture/CRISPR/HTO features. Confirm in-tool before quoting.

# Single-Cell RNA-seq Pipeline

**"Analyze my single-cell RNA-seq data from counts to cell types"** -> Orchestrate QC filtering, normalization (scanpy/Seurat), batch integration (scVI/Harmony), clustering, marker detection, cell type annotation, and trajectory inference.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## Made-once commitments (decided at `cellranger count`, inherited downstream)

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Reference build + Ensembl vintage (CR v10 = Ensembl v110) | Every gene ID/marker lookup/cross-dataset merge; join on `gene_ids` (stable), not symbols (lossy across releases) |
| `--include-introns` (default TRUE since CR 7) | Total UMI/cell and snRNA-vs-scRNA comparability; nuclei need introns; mixing intron-on/off samples is a hidden batch axis |
| Cell-vs-empty-droplet calling | The filtered matrix is cell-CALLING only, NOT ambient-corrected; SoupX/CellBender and low-RNA-cell rescue need the RAW matrix (filtered-only is irreversible loss) |
| Feature/barcode namespace (`gex_only`) | `gex_only=True` silently drops ADT/HTO/CRISPR features; set False and split by `feature_types` for CITE-seq/hashed pools |

## Pipeline orchestration: the ordering decisions that make or break the result

Stage order is not arbitrary; getting it wrong silently corrupts every downstream step. The cross-cutting decisions this pipeline must get right:
- Correct-then-detect-then-normalize, never reorder. Ambient-RNA removal (SoupX/CellBender) reads the RAW droplet matrix and MUST precede doublet detection and the final QC thresholds (ambient inflation distorts both doublet scores and per-cell QC); doublet detection precedes integration (doublets seed fake intermediate clusters). See single-cell/preprocessing.
- Multiplexed (hashed) pools are demultiplexed FIRST. When samples are pooled with cell hashing (HTO/MULTI-seq) or genotype multiplexing, assign each cell to its sample of origin and remove cross-sample doublets before any per-sample step; the per-sample QC and doublet operations below assume the pool is already split by sample. See single-cell/hashing-demultiplexing.
- Per-sample QC and doublet detection come BEFORE any merge or integration. Ambient-RNA cleanup (SoupX/CellBender), mito/gene-count filtering, and doublet calling (expected rate ~0.8% per 1,000 recovered cells) are per-capture operations; running them on a pooled object lets one lane's artifacts contaminate the shared null and lets surviving doublets form fake intermediate clusters. See single-cell/preprocessing and single-cell/doublet-detection.
- Integrate, THEN cluster, THEN annotate - never reorder. For multi-sample designs, batch-correct on a shared embedding (Harmony/scVI/RPCA) and cluster on the corrected graph; clustering before correction makes clusters track samples or lanes instead of biology, and annotation has nothing to label until clusters exist. See single-cell/batch-integration, single-cell/clustering, single-cell/cell-annotation.
- Over-integration erases biology. Batch-mixing metrics (kBET, iLISI) are maximized by destroying structure, so a method that flattens real cell states scores perfectly; pair every batch metric with a bio-conservation metric and re-inspect rare populations after correction. Batch-corrected expression is for embedding and clustering only - never feed it to differential expression. See single-cell/batch-integration.
- Condition-level DE must use pseudobulk, not cells-as-replicates. Testing thousands of cells per donor as independent replicates is pseudoreplication and inflates false positives by orders of magnitude (Squair 2021); aggregate RAW counts per sample x cell type and hand them to DESeq2/edgeR/limma-voom. Step 7 marker p-values are descriptive ranking for labeling, not condition inference. See single-cell/markers-annotation and differential-expression/deseq2-basics.
- Composition shifts are tested separately and masquerade as DE. A cluster-level expression change between conditions can be a pure proportion shift (a mixed cluster's substates re-balance with no gene changing per cell), invisible if only DE is run; always pair condition DE with a differential-abundance test (Milo cluster-free, or scCODA/sccomp/propeller cluster-based). See single-cell/differential-abundance.

The Seurat and Scanpy paths below are written single-sample for clarity. A multiplexed design demultiplexes the pool first (single-cell/hashing-demultiplexing); a multi-sample design then inserts per-sample QC and doublet removal, a merge, and an integration step between normalization (Step 4) and dimensionality reduction (Step 5); the rest of the order is unchanged.

## Workflow Overview

```
10X data (RAW + filtered_feature_bc_matrix)
    |
    v
[0. Ambient removal] ---> SoupX/CellBender on RAW (per sample)   (single-cell/preprocessing)
    |
    v
[1. Load Data] ---------> Read10X / read_10x_h5
    |
    v
[2. QC + Doublets] -----> MAD-adaptive nFeature/percent.mt; scDblFinder/Scrublet (per sample, before merge)
    |
    v
[3. Normalization] -----> SCTransform or LogNormalize
    |
    v
[4. HVG Selection] -----> FindVariableFeatures
    |
    v
[5. Dim Reduction] -----> PCA -> UMAP
    |
    v
[6. Clustering] --------> FindNeighbors -> FindClusters
    |
    v
[7. Markers] -----------> FindAllMarkers
    |
    v
[8. Annotation] --------> Manual or automated
    |
    v
Annotated Seurat/AnnData object
```

## Primary Path: Seurat (R)

The Seurat/Scanpy paths below start from the filtered matrix for clarity. In practice, run ambient-RNA removal FIRST on the RAW droplet matrix per sample (SoupX needs raw+filtered+clusters; CellBender folds calling+denoise and is best for snRNA/high-ambient) — pick ONE (double-correction over-strips). See single-cell/preprocessing.

### Step 1: Load 10X Data

```r
library(Seurat)
library(ggplot2)
library(dplyr)

# Load from Cell Ranger output
data_dir <- 'cellranger_output/filtered_feature_bc_matrix'
counts <- Read10X(data.dir = data_dir)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, project = 'my_project',
                                  min.cells = 3, min.features = 200)
```

### Step 2: Quality Control

```r
# Calculate QC metrics
seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')
seurat_obj[['percent.ribo']] <- PercentageFeatureSet(seurat_obj, pattern = '^RP[SL]')

# Visualize QC metrics
VlnPlot(seurat_obj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

# Filter cells
seurat_obj <- subset(seurat_obj,
                     nFeature_RNA > 200 &
                     nFeature_RNA < 5000 &
                     percent.mt < 20 &
                     nCount_RNA > 500)

cat('Cells after QC:', ncol(seurat_obj), '\n')
```

**QC Checkpoint 1:** Review QC plots
- Remove cells with very low/high gene counts
- Remove cells with high mitochondrial content (dying cells)

### Step 3: Doublet Detection

```r
library(scDblFinder)

# Convert to SCE for scDblFinder
sce <- as.SingleCellExperiment(seurat_obj)
sce <- scDblFinder(sce)

# Add back to Seurat
seurat_obj$doublet_class <- sce$scDblFinder.class
seurat_obj$doublet_score <- sce$scDblFinder.score

# Remove doublets
seurat_obj <- subset(seurat_obj, doublet_class == 'singlet')
cat('Cells after doublet removal:', ncol(seurat_obj), '\n')
```

### Step 4: Normalization with SCTransform

```r
# SCTransform (recommended for most analyses)
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
```

Alternative: Standard normalization
```r
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = 'vst', nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
```

Regressing out `percent.mt` (or `nCount`/cell-cycle) is NOT reflexive: those covariates are confounded with real cell state and regressing them can erase biology, so only pass `vars.to.regress` for a covariate verified not confounded with the signal of interest. See single-cell/preprocessing.

### Step 5: Dimensionality Reduction

```r
# PCA
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

# Determine optimal PCs
ElbowPlot(seurat_obj, ndims = 50)

# UMAP
n_pcs <- 30  # Choose based on elbow plot
seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs, verbose = FALSE)
```

### Step 6: Clustering

```r
# Find neighbors
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs, verbose = FALSE)

# Find clusters (try multiple resolutions)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0), verbose = FALSE)

# Visualize
DimPlot(seurat_obj, reduction = 'umap', group.by = 'SCT_snn_res.0.4', label = TRUE)
```

**QC Checkpoint 2:** Assess clustering
- Clusters should be visually separable on UMAP
- Resolution 0.4-0.8 is often appropriate

### Step 7: Find Marker Genes

```r
# Set identity to chosen resolution
Idents(seurat_obj) <- 'SCT_snn_res.0.4'

# Find markers for all clusters
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top markers per cluster
top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

# Visualize top markers
DoHeatmap(seurat_obj, features = top_markers$gene) + NoLegend()
```

### Step 8: Cell Type Annotation

```r
# Manual annotation based on known markers
# Example for PBMC data:
cluster_annotations <- c(
    '0' = 'CD4 T cells',
    '1' = 'CD14 Monocytes',
    '2' = 'B cells',
    '3' = 'CD8 T cells',
    '4' = 'NK cells',
    '5' = 'CD16 Monocytes',
    '6' = 'Dendritic cells'
)

seurat_obj$cell_type <- cluster_annotations[as.character(Idents(seurat_obj))]

# Final UMAP
DimPlot(seurat_obj, reduction = 'umap', group.by = 'cell_type', label = TRUE)

# Save object
saveRDS(seurat_obj, 'seurat_annotated.rds')
```

## Alternative Path: Scanpy (Python)

```python
import scanpy as sc
import numpy as np

# Load 10X data
adata = sc.read_10x_h5('filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()

# QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter (flat cutoffs are illustrative; prefer MAD-adaptive, tissue-aware thresholds, see single-cell/preprocessing)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

# Doublet detection (rate from recovered cells, not the 0.05 placeholder, see single-cell/doublet-detection)
expected_rate = 0.008 * adata.n_obs / 1000
sc.pp.scrublet(adata, expected_doublet_rate=expected_rate)
adata = adata[~adata.obs['predicted_doublet'], :]

# Normalize and HVGs
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# PCA, neighbors, UMAP -- stash log-normalized values in .raw BEFORE scaling, or rank_genes_groups
# later computes logfoldchanges from z-scores (negative means -> NaN/garbage LFCs)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

# Clustering (pin the backend for reproducibility; default flips to igraph)
sc.tl.leiden(adata, resolution=0.5, flavor='igraph', n_iterations=2, directed=False)

# Markers
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

# Save
adata.write('scanpy_annotated.h5ad')
```

## Parameter Recommendations

| Step | Parameter | Recommendation |
|------|-----------|----------------|
| QC | min.features | 200-500 |
| QC | max.features | 2500-5000 (depends on data) |
| QC | percent.mt | <10-20% (tissue-dependent) |
| SCTransform | vars.to.regress | none by default; only a validated non-confounded covariate |
| PCA | npcs | 30-50 |
| UMAP | dims | 15-30 (check elbow plot) |
| Clustering | resolution | 0.4-0.8 (start with 0.5) |

QC cutoffs above are illustrative starting points, not universal thresholds: set min/max features and mito % per dataset from the data (MAD-adaptive, tissue-aware) rather than porting flat values, since healthy mito fraction varies by tissue. See single-cell/preprocessing.

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| All cells filtered | QC too strict | Relax thresholds |
| Poor UMAP separation | Too few HVGs or PCs | Increase nfeatures, check n_pcs |
| Too many/few clusters | Wrong resolution | Adjust resolution parameter |
| Unknown cell types | Missing markers | Check known marker genes manually |

## Complete R Workflow

```r
library(Seurat)
library(scDblFinder)
library(ggplot2)
library(dplyr)

# Configuration
data_dir <- 'filtered_feature_bc_matrix'
output_dir <- 'results'
dir.create(output_dir, showWarnings = FALSE)

# Load
counts <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
cat('Initial cells:', ncol(seurat_obj), '\n')

# QC
seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')
seurat_obj <- subset(seurat_obj, nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
cat('After QC:', ncol(seurat_obj), '\n')

# Doublets
sce <- as.SingleCellExperiment(seurat_obj)
sce <- scDblFinder(sce)
seurat_obj$doublet <- sce$scDblFinder.class
seurat_obj <- subset(seurat_obj, doublet == 'singlet')
cat('After doublet removal:', ncol(seurat_obj), '\n')

# Normalize (no reflexive vars.to.regress; regress only a validated non-confounded covariate)
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

# Dimension reduction
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

# Cluster
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

# Markers
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(output_dir, 'markers.csv'))

# Save
saveRDS(seurat_obj, file.path(output_dir, 'seurat_object.rds'))

# Plots
pdf(file.path(output_dir, 'umap.pdf'), width = 10, height = 8)
DimPlot(seurat_obj, reduction = 'umap', label = TRUE)
dev.off()

cat('Pipeline complete. Object saved to:', output_dir, '\n')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Ambient markers everywhere (e.g. Hb across all PBMCs) | Ran SoupX/CellBender on the FILTERED matrix (soup already removed) | Run ambient removal on the RAW droplet matrix, before QC |
| Fake intermediate cell states | Doublets removed AFTER clustering/integration | Call doublets per sample before merge/normalize |
| Clusters track sample/lane, not biology | Clustered before integration | Integrate -> cluster -> annotate |
| Inflated DE, thousands of false positives | Tested cells as replicates (on integrated values) | Pseudobulk RAW counts per sample x cell-type -> DESeq2/edgeR/limma (Squair 2021) |
| "DE change" is really a proportion shift | Only ran DE, not abundance | Pair with differential abundance (Milo/scCODA/propeller) |
| Biology collapses to one blob | Reflexively regressed total_counts/cell-cycle | Regress only a validated, non-confounded covariate |
| CITE-seq/hashtag features vanish | `gex_only=True` default | Set False; split by `feature_types` |

## References

- Heumos L, Schaar AC, Lance C, et al (2023) Best practices for single-cell analysis across modalities. *Nature Reviews Genetics* 24:550-572. DOI 10.1038/s41576-023-00586-w. (canonical pipeline order.)
- Squair JW, Gautier M, Kathe C, et al (2021) Confronting false discoveries in single-cell differential expression. *Nature Communications* 12:5692. DOI 10.1038/s41467-021-25960-2. (cells-as-replicates is pseudoreplication; pseudobulk.)
- Fleming SJ, Chaffin MD, Arduini A, et al (2023) Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender. *Nature Methods* 20:1323-1335. DOI 10.1038/s41592-023-01943-7. (ambient removal on the RAW matrix.)

## Related Skills

- database-access/geo-data - Resolve GSE to SRA; detect SuperSeries before processing
- database-access/sra-data - Download 10x records with --include-technical for barcodes/UMIs
- single-cell/data-io - Loading 10X, h5ad, RDS, and h5mu formats
- single-cell/preprocessing - QC thresholds, ambient-RNA removal, normalization choice
- single-cell/doublet-detection - Per-sample doublet calling before integration
- single-cell/hashing-demultiplexing - Assign multiplexed pools to samples and call cross-sample doublets before QC
- single-cell/batch-integration - Multi-sample integration and over-correction diagnosis
- single-cell/clustering - Resolution sweep and cluster validation
- single-cell/markers-annotation - Marker discovery, manual labeling, and pseudobulk condition DE
- single-cell/cell-annotation - Automated reference-based label transfer
- single-cell/differential-abundance - Test whether cell-type proportions shifted between conditions
- single-cell/trajectory-inference - Pseudotime and lineage reconstruction for continuous processes
- single-cell/multimodal-integration - CITE-seq and multiome joint analysis
- differential-expression/deseq2-basics - Pseudobulk condition DE engine for aggregated counts
- differential-expression/de-results - Shrink, filter, and interpret pseudobulk DE results
- pathway-analysis/go-enrichment - Functional interpretation of marker and DE gene lists
- workflows/grn-pipeline - Downstream: infer gene regulatory networks (pySCENIC Path A) from the annotated object
- workflows/spatial-pipeline - Downstream: the annotated reference deconvolves spatial spots
