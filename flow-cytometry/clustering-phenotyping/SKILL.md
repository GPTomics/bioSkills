---
name: bio-flow-cytometry-clustering-phenotyping
description: Unsupervised clustering and cell-type identification for high-dimensional flow, spectral, and mass cytometry - FlowSOM, PhenoGraph, FlowSOM-via-CATALYST, with UMAP/tSNE for visualization. Covers the type-vs-state marker distinction (cluster on lineage, test state within clusters), over-provision-then-metacluster, the Weber-Robinson benchmark, seed dependence and metacluster stability, why embeddings are for looking not measuring, and median-heatmap annotation/merging. Use when discovering populations without predefined gates, choosing a clustering algorithm, selecting the number of metaclusters, or annotating clusters into cell types.
tool_type: r
primary_tool: CATALYST
---

## Version Compatibility

Reference examples tested with: CATALYST 1.26+, FlowSOM 2.10+, flowCore 2.14+; Rphenograph (GitHub: JinmiaoChenLab/Rphenograph).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

`Rphenograph` is GitHub-only (`remotes::install_github('JinmiaoChenLab/Rphenograph')`) and returns a list - membership is `igraph::membership(out[[2]])`, not a vector. Adapt rather than retrying.

# Clustering and Phenotyping

**"Cluster my cytometry data to find cell types"** -> Discover populations in high-dimensional data without gates, then annotate them by marker expression.
- R: `CATALYST::cluster()` (wraps FlowSOM + ConsensusClusterPlus) - the field default
- R: `FlowSOM::FlowSOM()` directly, or `Rphenograph()` for graph-based clustering

## The Single Most Important Modern Insight -- Cluster on Type Markers; the Embedding Is for Looking, Not Measuring

Two rules carry most of the correctness here. First, the type-vs-state distinction: LINEAGE/type markers (CD3, CD4, CD8, CD19) DEFINE clusters; functional/STATE markers (phospho-epitopes, cytokines, Ki-67, activation markers) must be WITHHELD from clustering and tested within clusters instead (the DA/DS framework, Nowicka 2017 *F1000Res* 6:748). Clustering on state markers splits "activated CD4" from "resting CD4" and confounds abundance with activation - a classic, silent design error. Second, t-SNE/UMAP embeddings do NOT preserve inter-cluster distances, cluster sizes, or densities (the apparent "UMAP preserves global structure" edge over tSNE is largely an initialization artifact - Kobak & Linderman 2021 *Nat Biotechnol* 39:156). Define populations by clustering in the HIGH-DIMENSIONAL space and COLOR the embedding by cluster; never gate on the embedding or read biology off blob distances.

## Algorithm Taxonomy

| Algorithm | Citation | Mechanism | Speed | Rare-pop | Determinism |
|-----------|----------|-----------|-------|----------|-------------|
| FlowSOM | Van Gassen 2015 *Cytometry A* 87:636 | SOM grid -> MST (viz) -> consensus metaclustering | fastest | good if grid over-provisioned | stochastic; seed-controllable |
| PhenoGraph | Levine 2015 *Cell* 162:184 | kNN graph (Jaccard) + Louvain | moderate | strong (no preset k) | seed-fragile (>40% reassignment reported) |
| X-shift | Samusik 2016 *Nat Methods* 13:493 | weighted kNN density + auto cluster # | slow | excellent | more deterministic |
| flowMeans | Aghaeepour 2011 *Cytometry A* 79:6 | k-means multi-cluster + change-point k | fast | moderate | stochastic |

Benchmark: Weber & Robinson 2016 *Cytometry A* 89:1084 tested 18 methods - FlowSOM (with metaclustering) was a top performer AND by far fastest, hence the field default; but its accuracy depends on supplying the right number of metaclusters.

## Why Over-Provision the Grid, Then Metacluster

Set the SOM grid (e.g. 10x10 = 100 nodes) MUCH larger than the number of populations expected, then metacluster down. The asymmetry: metaclustering can MERGE over-fine nodes into a real population, but can NEVER SPLIT a node that erroneously fused two cell types. Too coarse commits the unrecoverable error; too fine commits only the recoverable one. So over-cluster, then merge by hand off the median heatmap.

## CATALYST Clustering Pipeline

**Goal:** Cluster on type markers and prepare for annotation.

**Approach:** `prepData` builds the SCE (panel `marker_class` flags type vs state); `cluster()` wraps FlowSOM+ConsensusClusterPlus. Defaults `xdim=ydim=10`, `maxK=20` (the metacluster cap people forget); set `seed` on the function.

```r
library(CATALYST)

sce <- prepData(fs, panel, md, transform = TRUE, cofactor = 5)   # cofactor 5 = CyTOF; ~150 for fluorescence
sce <- cluster(sce, features = 'type',                            # type markers only
               xdim = 10, ydim = 10, maxK = 20, seed = 42)        # maxK caps metaclusters at 20 by default
plotExprHeatmap(sce, features = 'type', by = 'cluster_id', k = 'meta20', scale = 'last')
```

## PhenoGraph (graph-based alternative)

**Goal:** Cluster with a kNN graph when a data-driven cluster count is wanted.

**Approach:** `Rphenograph` on the type-marker matrix (cells x markers); extract membership from the list.

```r
library(Rphenograph)
type_expr <- t(assay(sce, 'exprs')[rowData(sce)$marker_class == 'type', ])
out <- Rphenograph(type_expr, k = 30)                  # only knob: k (neighbors)
sce$phenograph <- factor(igraph::membership(out[[2]]))  # list -> membership, not a vector
```

## Dimensionality Reduction (visualization only) and Annotation

**Goal:** Visualize structure and assign cell-type labels.

**Approach:** `runDR` subsamples per sample (`cells=`); color by cluster, never gate on it. Annotate from the median heatmap, then `mergeClusters` with a curated table.

```r
sce <- runDR(sce, dr = 'UMAP', features = 'type', cells = 2000)   # subsampled embedding
plotDR(sce, 'UMAP', color_by = 'meta20')

merging <- data.frame(old_cluster = 1:20,
                      new_cluster = c('CD4 T','CD4 T','CD8 T', '...'))   # curated from the heatmap
sce <- mergeClusters(sce, k = 'meta20', table = merging, id = 'annotated')
```

## Per-Method Failure Modes

### Clustering on state markers
**Trigger:** activation/phospho markers in the clustering feature set. **Mechanism:** state contaminates lineage identity. **Symptom:** "activated" and "resting" versions of a type split as separate clusters. **Fix:** cluster on `type` only; test state markers within clusters (differential-analysis).

### Seed-dependent "novel populations"
**Trigger:** a population that appears at one seed and vanishes at another. **Mechanism:** FlowSOM init / Louvain are stochastic. **Symptom:** non-reproducible clusters. **Fix:** set + report the seed; check multi-seed stability; treat unstable clusters as hypotheses.

### Reading biology off the embedding
**Trigger:** "cluster A is closer to B than C." **Mechanism:** UMAP/tSNE distances are non-metric. **Symptom:** false developmental/relatedness claims. **Fix:** quantify in marker space; embedding for display only.

### Clustering uncompensated/untransformed data
**Trigger:** raw linear input to FlowSOM. **Mechanism:** spillover + scale dominate Euclidean distance. **Symptom:** clusters track intensity, not biology. **Fix:** compensate + transform first.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| over-provision grid (10x10) >> expected pops | Van Gassen 2015 | metacluster can merge, never split |
| maxK = 20 default | CATALYST | metacluster cap; raise if expecting more |
| FlowSOM needs correct K | Weber & Robinson 2016 | accuracy depends on metacluster number |
| use median (not mean) per cluster | Bendall 2011 *Science* 332:687 | robust to doublet/spillover contamination |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| clustering uses scatter/Time/state | `features` not restricted | `features='type'` / `colsToUse=` lineage markers |
| Rphenograph result unusable | it returns a list | `igraph::membership(out[[2]])` |
| `set.seed` doesn't make FlowSOM reproducible | internal reseeding | pass `seed=` to `cluster()` |
| only 20 clusters no matter what | `maxK` default | raise `maxK` |

## References

- Van Gassen 2015 *Cytometry A* 87(7):636-645 — FlowSOM.
- Levine 2015 *Cell* 162(1):184-197 — PhenoGraph.
- Samusik 2016 *Nat Methods* 13(6):493-496 — X-shift.
- Weber & Robinson 2016 *Cytometry A* 89(12):1084-1096 — clustering benchmark (FlowSOM top + fastest).
- Nowicka 2017 *F1000Research* 6:748 — CyTOF workflow; type-vs-state markers.
- Kobak & Linderman 2021 *Nat Biotechnol* 39:156-157 — embedding initialization artifact.
- Bendall 2011 *Science* 332(6030):687-696 — arcsinh-median analysis of CyTOF data.

## Related Skills

- compensation-transformation - Compensate/transform before clustering
- gating-analysis - Supervised alternative; needed for rare populations
- differential-analysis - Test abundance/state of clusters between conditions
- cytometry-qc - Cluster only QC-passed events
- single-cell/clustering - Leiden/Louvain on scRNA-seq (shared graph-clustering ideas)
- imaging-mass-cytometry/phenotyping - Same CATALYST/FlowSOM conventions for imaging
