---
name: bio-single-cell-trajectory-inference
description: Infer developmental trajectories and pseudotime from single-cell RNA-seq data using Monocle3, Slingshot, and scVelo for RNA velocity analysis. Use when inferring developmental trajectories or pseudotime.
tool_type: mixed
primary_tool: Monocle3
---

# Trajectory Inference

## Monocle3 (R)

```r
library(monocle3)

# Create cell_data_set from Seurat
cds <- as.cell_data_set(seurat_obj)

# Preprocess (if not already done)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = 'UMAP')

# Cluster cells
cds <- cluster_cells(cds)

# Learn trajectory graph
cds <- learn_graph(cds)

# Order cells (select root interactively or programmatically)
cds <- order_cells(cds, root_cells = root_cell_ids)

# Plot trajectory with pseudotime
plot_cells(cds, color_cells_by = 'pseudotime', label_branch_points = TRUE, label_leaves = TRUE)

# Get pseudotime values
pseudotime <- pseudotime(cds)
```

## Set Root Programmatically

```r
# Find root by marker gene expression
get_earliest_principal_node <- function(cds, cluster_name) {
    cell_ids <- which(colData(cds)$seurat_clusters == cluster_name)
    closest_vertex <- cds@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[cell_ids, ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[['UMAP']])$name[as.numeric(names(which.max(table(closest_vertex))))]
    root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, 'stem_cluster'))
```

## Slingshot (R)

```r
library(slingshot)
library(SingleCellExperiment)

# From Seurat object
sce <- as.SingleCellExperiment(seurat_obj)
reducedDims(sce)$UMAP <- Embeddings(seurat_obj, 'umap')

# Run slingshot
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')

# Get pseudotime for each lineage
pseudotime_mat <- slingPseudotime(sce)

# Get lineage curves
curves <- slingCurves(sce)

# Plot trajectories
plot(reducedDims(sce)$UMAP, col = sce$seurat_clusters, pch = 16)
lines(SlingshotDataSet(sce), lwd = 2)
```

## Slingshot with Start/End Clusters

```r
# Specify starting cluster
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP', start.clus = 'HSC')

# Specify start and end
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP',
                 start.clus = 'HSC', end.clus = c('Erythroid', 'Myeloid'))
```

## scVelo RNA Velocity (Python)

```python
import scvelo as scv
import scanpy as sc

# Load data with spliced/unspliced counts
adata = scv.read('data.h5ad')

# Or merge loom files from velocyto
ldata = scv.read('velocyto_output.loom')
adata = scv.utils.merge(adata, ldata)

# Preprocess
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Compute velocity (stochastic model)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

# Visualize velocity streams
scv.pl.velocity_embedding_stream(adata, basis='umap', color='clusters')
```

## scVelo Dynamical Model

```python
# More accurate but slower
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

# Latent time (pseudotime)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', cmap='gnuplot')

# Velocity confidence
scv.tl.velocity_confidence(adata)
scv.pl.scatter(adata, color=['velocity_confidence', 'velocity_length'])
```

## Gene Dynamics Along Trajectory

```r
# Monocle3: Find genes varying over pseudotime
graph_test_res <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)
sig_genes <- graph_test_res %>% filter(q_value < 0.05) %>% arrange(desc(morans_I))

# Plot gene expression over pseudotime
plot_genes_in_pseudotime(cds[rownames(cds) %in% top_genes, ], color_cells_by = 'cluster')
```

```python
# scVelo: Top likelihood genes
scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=0.3)
top_genes = adata.uns['rank_velocity_genes']['names']

# Plot phase portraits
scv.pl.velocity(adata, var_names=['gene1', 'gene2'], basis='umap')
```

## Branch Point Analysis

```r
# Monocle3: Genes differentially expressed at branch points
branch_genes <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

# Slingshot + tradeSeq for branch analysis
library(tradeSeq)
sce <- fitGAM(sce, nknots = 6)
branch_res <- earlyDETest(sce, knots = c(3, 4))
```

## Velocyto Preprocessing

```bash
# Generate loom file with spliced/unspliced counts
velocyto run10x -m repeat_mask.gtf /path/to/cellranger_output annotation.gtf

# For SmartSeq2
velocyto run_smartseq2 -o output -m repeat_mask.gtf -e sample bam_files/*.bam annotation.gtf
```

## PAGA Trajectory (Scanpy)

```python
import scanpy as sc

# Compute PAGA
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color='leiden', threshold=0.03)

# PAGA-initialized UMAP
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color='leiden')

# Diffusion pseudotime
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden'] == 'root_cluster')[0]
sc.tl.dpt(adata)
sc.pl.draw_graph(adata, color='dpt_pseudotime')
```

## Related Skills

- single-cell/clustering - Prerequisite clustering
- single-cell/cell-communication - Downstream signaling analysis
- differential-expression/deseq2-basics - DE along trajectory
