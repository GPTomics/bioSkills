---
name: bio-single-cell-clustering
description: Dimensionality reduction and graph-based clustering for single-cell RNA-seq with Scanpy (Python) and Seurat (R). Resolves which algorithm to use (Leiden vs Louvain), how many PCs and neighbors to set, how to sweep and validate resolution, when a split is over-clustering, and why post-clustering marker p-values are not valid inference. Use when clustering cells, choosing a clustering resolution, deciding whether two clusters are one population, building a UMAP/tSNE, or judging whether clusters are real.
tool_type: mixed
primary_tool: Seurat
---

## Version Compatibility

Reference examples tested with: scanpy 1.10+, Seurat 5.0+, anndata 0.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Single-Cell Clustering

**"Cluster my cells"** -> Build a k-nearest-neighbor graph in PCA space, partition it into communities, and embed in 2D for display.
- Python: `sc.pp.neighbors` + `sc.tl.leiden` + `sc.tl.umap` (Scanpy)
- R: `FindNeighbors` + `FindClusters` + `RunUMAP` (Seurat)

## Governing Principle

Resolution is not a truth knob, and clusters are not discoveries. Graph-based clustering partitions a kNN graph built in PCA space; there is no ground truth, no "correct" number of clusters, and resolution selects a scale of description rather than revealing one. Clusters are hypotheses that must be validated by markers, stability, and (where claims are made) a significance test before any partition is named a cell population. Over-clustering is the default failure mode: any homogeneous blob can be bisected, and a higher resolution will always split it further. UMAP/tSNE distances, cluster sizes, and apparent gaps are artifacts of a non-linear neighbor-preserving objective and are not metric data (Chari & Pachter 2023) - cluster on the graph, never on embedding coordinates. The deepest trap: clustering chooses labels to maximize between-group separation, so running a marker test on those same clusters tests a hypothesis built from the data used to test it (double-dipping), and the resulting p-values are not merely inflated, they are invalid inference.

## Leiden vs Louvain

Leiden is the current default for graph community detection because Louvain can return internally disconnected communities (Traag 2019). Seurat still ships Louvain as its default algorithm; Scanpy uses Leiden but is mid-migration between backends. Methodology evolves - verify the current default and backend against the installed package docs before pinning a pipeline.

| Method | Model / assumption | Use when | Fails when |
|--------|-------------------|----------|------------|
| Leiden | Modularity/CPM optimization with a refinement phase guaranteeing connected communities | Default for scRNA-seq; reproducibility matters; large graphs (faster) | Backend/iteration count left unpinned -> silently different labels across Scanpy versions |
| Louvain | Modularity optimization without refinement | Legacy pipelines; Seurat default (`algorithm=1`) | Can yield internally disconnected communities; superseded by Leiden (Traag 2019) |
| SLM | Smart Local Moving refinement | Seurat option (`algorithm=3`) for tighter modularity optima | Slower; rarely needed over Leiden |

Scanpy 1.10 backend migration (pin for reproducibility): `sc.tl.leiden` still defaults to the `leidenalg` backend through 1.10-1.12 and emits a FutureWarning that the default will switch to igraph. The exact flip version is unconfirmed, so pin the backend explicitly: `sc.tl.leiden(adata, flavor='igraph', n_iterations=2, directed=False)`. Switching backend or `n_iterations` changes the labels - a pipeline that pins neither is non-reproducible across versions. Seurat's Leiden (`algorithm=4`) requires the `leidenalg` Python module via reticulate, which is why most Seurat pipelines still run Louvain.

## Parameter Reference

`n_pcs` dominates the result far more than `n_neighbors` and is the largest under-tuned lever - too few collapses real structure, too many reintroduces technical noise.

| Parameter | Typical range | Rationale | Validation |
|-----------|--------------|-----------|------------|
| n_pcs | 30-50 (check elbow) | Captures biological variance while denoising; effect dwarfs n_neighbors | Elbow plot; cluster stability across nearby n_pcs values |
| n_neighbors | 10-30 (15 default) | Higher = smoother, fewer fine clusters; lower = more local, fragmented | Secondary lever; vary only after n_pcs is set |
| resolution | 0.2-2.0 (sweep, do not fix) | Higher = more, smaller clusters; has no biological meaning | clustree across the sweep; marker check; significance test |
| min_dist (UMAP) | 0.1-0.5 | Visualization only; lower = tighter visual clusters | Affects display, never the partition |

Resolution is an unidentifiable nuisance parameter: it cannot be validated internally (no ground truth), so tuning it until clusters "match known cell types" is confirmation bias laundered as analysis. Sweep a range, visualize cell flow with clustree, and pick the coarsest level whose populations are defensible by orthogonal evidence - label finer splits as hypotheses.

## Cluster Cells with Scanpy

**Goal:** Reduce dimensions, build the neighbor graph, partition with Leiden, and embed for display.
**Approach:** PCA -> kNN graph on a chosen `n_pcs` -> Leiden with a pinned backend -> UMAP.

```python
import scanpy as sc

adata = sc.read_h5ad('preprocessed.h5ad')

sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)  # elbow to choose n_pcs

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5, flavor='igraph', n_iterations=2, directed=False)
adata.obs['leiden'].value_counts()

sc.tl.umap(adata, min_dist=0.3)
sc.pl.umap(adata, color=['leiden', 'CD3D', 'MS4A1', 'CD14'])
```

## Sweep and Validate Resolution

**Goal:** Choose a defensible granularity instead of a single tuned-to-taste resolution.
**Approach:** Cluster across a resolution range, inspect cell flow (clustree), and confirm each cluster carries distinct markers.

```python
import scanpy as sc

for res in [0.2, 0.4, 0.6, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_r{res}',
                 flavor='igraph', n_iterations=2, directed=False)
    print(res, adata.obs[f'leiden_r{res}'].nunique(), 'clusters')

sc.pl.umap(adata, color=['leiden_r0.2', 'leiden_r0.6', 'leiden_r1.0'], ncols=3)
# clustree (R) or sc.tl.dendrogram for flow across resolutions; merge clusters
# whose top markers are indistinguishable -> they are one population over-split
```

## Cluster Cells with Seurat

**Goal:** Run PCA, build the SNN graph, partition, and embed in Seurat.
**Approach:** RunPCA -> FindNeighbors on chosen dims -> FindClusters (sweep resolutions) -> RunUMAP.

```r
library(Seurat)

seurat_obj <- readRDS('preprocessed.rds')
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
ElbowPlot(seurat_obj, ndims = 50)  # choose dims

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

library(clustree)
clustree(seurat_obj, prefix = 'RNA_snn_res.')  # cell flow across the sweep
DimPlot(seurat_obj, reduction = 'umap', label = TRUE)
```

`FindClusters` defaults to Louvain (`algorithm=1`); pass `algorithm=4` for Leiden (requires the leidenalg Python module). Resolutions stored as `RNA_snn_res.<r>` columns feed clustree directly.

## Subclustering

**Goal:** Resolve fine states inside a coarse cluster without importing global axes.
**Approach:** Subset the cluster, then recompute HVGs, PCA, and the kNN graph on the subset.

```python
sub = adata[adata.obs['leiden'] == '3'].copy()
sc.pp.highly_variable_genes(sub, n_top_genes=2000)
sub = sub[:, sub.var.highly_variable]
sc.pp.scale(sub, max_value=10)
sc.tl.pca(sub, n_comps=30)
sc.pp.neighbors(sub, n_neighbors=15, n_pcs=20)
sc.tl.leiden(sub, resolution=0.4, flavor='igraph', n_iterations=2, directed=False)
```

Reusing the global PCA imports axes uninformative within a homogeneous subset and manufactures artifactual sub-splits. Subclustering compounds double-dipping (cells selected twice), so stop when splits lose distinct markers or fail a significance test - not when resolution can technically still split (it always can).

## Validating That Clusters Are Real

Stability and significance are separate questions, and both differ from biological reality.

- Stability (necessary, not sufficient): bootstrap cells, re-cluster, measure per-cluster label agreement (Jaccard >= ~0.6-0.7 = stable). A perfectly reproducible split can still be technical - driven by cell-cycle phase, dissociation stress (FOS/JUN/HSPA1A), mitochondrial fraction, ambient RNA, or batch. Stable does not mean real.
- Significance (whether a split is two populations or one): scSHC (Grabski 2023) and CHOIR (2025) test each split under a null with error control; a failed split is over-clustering and the clusters should be merged.
- Biological-vs-technical adjudication: check that a split survives regressing out cycle/mito and carries non-stress markers before claiming a population.

Double-dipping (post-clustering inference is invalid, not just inflated): `rank_genes_groups`/`FindAllMarkers` p-values are conditioned on a clustering chosen to maximize separation, so under one homogeneous population they do not follow their nominal null - type-I error approaches 1 as resolution rises, and BH correction does nothing because the p-values are invalid before correction. Use these marker tests for ranking and labeling only. To make a defensible claim that a cluster is a distinct population, pair a cluster significance test (scSHC/CHOIR) with a double-dipping-robust DE method (ClusterDE's synthetic null, or count splitting where the noise model holds - Poisson thinning on overdispersed counts silently reinstates the bias). See markers-annotation for marker testing and the pseudobulk path for cross-condition DE.

## UMAP and tSNE Are Visualization Only

Cluster on the graph or PCA, never on embedding coordinates. Inter-cluster distances, relative sizes, and apparent gaps in UMAP/tSNE are artifacts of the embedding objective and are not metric (Chari & Pachter 2023) - do not read them as lineage, evolutionary distance, or population separation. "Cells far apart in UMAP are more different" and "UMAP preserves global structure" are folklore; PCA initialization plus high perplexity makes tSNE less misleading but does not make distances trustworthy (Kobak & Berens 2019). Embeddings display graph-derived labels; back every structural claim with the graph, validly-tested markers, or quantitative analysis in PCA space.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| One giant blob, no structure | Too few HVGs or wrong n_pcs (too few), or the sample is genuinely one cell type | Increase HVGs (~2000), raise n_pcs, check the elbow plot; if markers stay uniform across a resolution sweep, the blob may be a single real population, not a parameter bug |
| Far too many clusters | Resolution too high; n_pcs too high (noise) | Lower resolution; sweep with clustree; reduce n_pcs to the elbow |
| Adjacent clusters share all top markers | Over-clustering one population | Merge them; lower resolution; significance-test the split (scSHC/CHOIR) |
| Labels change between runs/versions | Leiden backend or n_iterations unpinned | Pin `flavor='igraph', n_iterations=2, directed=False`; set random_state |
| A cluster maps to one sample/lane only | Batch effect, not biology | Integrate batches first (batch-integration); inspect QC covariates |
| A "stable" cluster of stress/cycle genes | Technical split (dissociation, cycle, mito) | Regress out cycle/mito or score and exclude; require non-stress markers |
| Cluster expresses two lineages' markers | Doublets clustering together | Run doublet detection before clustering (doublet-detection) |
| Marker p-values quoted as proof clusters are real | Double-dipping (selective inference) | Use markers for ranking only; validate with scSHC/CHOIR + ClusterDE |

## Related Skills

- preprocessing - QC, normalization, and HVG selection that must precede clustering
- doublet-detection - Remove doublets before clustering so they do not form fake intermediate clusters
- batch-integration - Integrate batches before clustering when a cluster tracks a single sample
- markers-annotation - Find and test markers per cluster (with the double-dipping caveat)
- cell-annotation - Assign cell-type identities to validated clusters
- single-cell/differential-abundance - Test whether cluster proportions shift across conditions
- data-visualization/dimensionality-reduction-plots - Publication-quality UMAP/tSNE/PCA figures
- pathway-analysis/go-enrichment - Interpret per-cluster marker sets

## References

- Traag, Waltman & van Eck (2019). From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep 9:5233.
- Chari & Pachter (2023). The specious art of single-cell genomics. PLoS Comput Biol 19(8):e1011288.
- Kobak & Berens (2019). The art of using t-SNE for single-cell transcriptomics. Nat Commun 10:5416.
- Grabski, Street & Irizarry (2023). Significance analysis for clustering with single-cell RNA-sequencing data. Nat Methods 20:1196-1202.
- Neufeld, Gao, Popp, Battle & Witten (2024). Inference after latent variable estimation for single-cell RNA-seq (count splitting). Biostatistics 25(1):270-287.
- Zappia & Oshlack (2018). Clustering trees: a visualization for evaluating clusterings at multiple resolutions. GigaScience 7(7):giy083.
