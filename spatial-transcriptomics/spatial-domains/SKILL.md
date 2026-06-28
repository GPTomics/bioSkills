---
name: bio-spatial-transcriptomics-spatial-domains
description: Identify spatially coherent tissue domains (regions like cortical layers, tumor vs stroma) in Visium, Visium HD, Xenium, MERFISH, Slide-seq, and Stereo-seq data with Squidpy, BANKSY, BayesSpace, STAGATE, and GraphST. Use when distinguishing a domain (a region with many cell types) from a cell type (one cell's identity) and a niche (local cell-type composition); choosing a domain method by tissue geometry (laminar/continuous vs high-resolution imaging vs non-contiguous); tuning the spatial-weight knob (BANKSY lambda, BayesSpace smoothing, SpaGCN histology weight, GNN graph radius) to avoid over-smoothing into blobs or under-smoothing into salt-and-pepper; choosing the number of domains k as a biological decision with k+-1 sensitivity; and reading the Yuan 2024 benchmark with the DLPFC continuous-laminar caveat.
tool_type: python
primary_tool: squidpy
---

## Version Compatibility

Reference examples tested with: squidpy 1.4+, scanpy 1.10+, anndata 0.10+, scikit-learn 1.4+; BANKSY (banksy_py / R Banksy), BayesSpace 1.12+ (R), STAGATE/GraphST optional

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Spatial Domain Detection

**"Identify tissue domains in my section"** -> Partition the tissue into spatially contiguous regions of homogeneous expression/composition, using BOTH transcriptional similarity AND spatial proximity, so the output is a region label per spot/cell that is spatially coherent.
- Python: `squidpy.gr.spatial_neighbors` for the graph, then a neighbor-augmented or graph-based domain method (BANKSY, STAGATE, GraphST)
- R: BayesSpace (`spatialCluster`), BASS, or Banksy

## Governing Principle

A spatial domain is a REGION, not a cell type, and not a niche -- conflating the three is the central conceptual error of this analysis. A CELL TYPE is one cell's transcriptional identity. A NICHE (cellular neighborhood) is the local cell-type COMPOSITION around a cell -- which types co-occur. A SPATIAL DOMAIN is a contiguous tissue region (a cortical layer, tumor core, stromal band) that contains MANY cell types and several niches. If the question is "which cell types co-occur," that is a niche question and belongs to neighborhood-enrichment analysis (spatial-statistics), NOT domain segmentation; running the wrong one answers the wrong question.

Domain methods exist because plain Leiden/Louvain on expression alone ignores coordinates and produces salt-and-pepper, spatially-incoherent labels. Every domain method deliberately adds a spatial term so neighbors tend to share a label. The load-bearing decision is therefore NOT the clustering algorithm -- it is the SPATIAL-WEIGHT knob (BANKSY lambda, BayesSpace smoothing, SpaGCN histology weight, GNN graph radius). Too much spatial weight is the #1 trap, over-smoothing: it erases real boundaries and merges biologically distinct regions into blobs. Too little reverts to salt-and-pepper. The number of domains k is the second decision, and it is a BIOLOGICAL choice (how fine a regionalization the question needs), not a statistic to optimize against a silhouette score -- report k, justify it, and show sensitivity to k+-1.

## Domain vs Niche vs Cell Type

| Concept | What it is | Unit | Right tool |
|---------|-----------|------|-----------|
| Cell type | One cell's transcriptional identity (lineage/state) | a cell | clustering + markers (single-cell/clustering) |
| Niche / cellular neighborhood | Local cell-type composition -- which types co-occur around a cell | a neighborhood | neighborhood enrichment / co-occurrence (spatial-statistics) |
| Spatial domain | Contiguous tissue region of homogeneous expression, containing many types | a region | domain segmentation (this skill) |

A domain can contain several niches; a niche can span domain boundaries. BANKSY and BASS switch between cell-typing and domain detection via a single parameter, which underscores that these are different outputs of related machinery, not the same task.

## Domain Methods by Mechanism

The mechanism for "using the neighbors" is the axis that separates the methods. No method is a universal winner (Yuan 2024 *Nat Methods*); selection is scenario-specific.

| Method | Spatial mechanism | Needs k? | Best when | Fails when |
|--------|-------------------|----------|-----------|------------|
| BayesSpace | MRF (Potts) smoothing prior on a t-mixture in PCA space; also enhances sub-spot resolution | yes (q) | Visium laminar/continuous tissue; want sub-spot enhancement | non-contiguous regions; MCMC is slow; smoothing strength is a fixed prior |
| BASS | Bayesian hierarchical: cell-type AND domain jointly, Potts prior, multi-sample | yes (C, D) | joint cell-type + domain, multi-sample, low-continuity tissue | heavier to run; needs both counts set |
| STAGATE | Graph attention autoencoder; learns per-edge weights, reconstructs from a spatially-smoothed latent | embedding free; downstream k for mclust | continuous tissue; attention down-weights cross-boundary edges (fights over-smoothing); scales to Slide-seq/Stereo-seq | black-box embedding; sensitive to graph radius |
| GraphST | Graph self-supervised contrastive learning | yes (mclust) | low-res Visium; also does integration + deconvolution | sensitive to graph construction |
| SpaGCN | GCN fusing expression + coordinates + histology RGB | searches resolution to hit target count | H&E histology is informative | needs registered histology |
| SEDR | Masked autoencoder + variational graph autoencoder | yes (mclust) | Visium/Slide-seq/Stereo-seq, robust on DLPFC | graph-construction sensitivity |
| BANKSY | Neighbor-AUGMENTED features: own + neighborhood-mean + azimuthal Gabor; then Leiden/k-means | no fixed k | high-res imaging AND sequencing; unifies cell typing (lambda~0.2) and domains (lambda~0.8); very scalable; transparent | lambda mis-set collapses the task it solves |
| UTAG | Message passing: multiply features by normalized adjacency (one-hop average), then Leiden | no fixed k | multiplexed imaging/proteomics (IMC, CODEX, MIBI); fast | one-hop smoothing only |
| stLearn (SME) | Histology-CNN-weighted smoothing of each spot's expression, then clustering | downstream Louvain/k-means | H&E available and well-registered | only as good as image registration |

Mechanism families: MRF/Bayesian smoothing (BayesSpace, BASS, PRECAST) vs graph neural net (STAGATE, GraphST, SpaGCN, SEDR) vs neighbor-augmentation (BANKSY, UTAG) vs histology-guided (stLearn, SpaGCN). Benchmarks evolve fast -- verify the current verdict for the tissue geometry before committing.

## Read the Benchmark With the DLPFC Caveat

The standard ground truth is DLPFC (Maynard 2021 *Nat Neurosci* 24:425-436): 12 Visium sections of human dorsolateral prefrontal cortex, manually annotated into 6 cortical layers + white matter, scored by ARI. The Yuan 2024 (*Nat Methods* 21:712-722) benchmark of 13 methods x 34 datasets concludes there is NO single winner -- methods are complementary across accuracy, spatial continuity, marker detection, scalability, and robustness. On DLPFC, GNN/graph methods (STAGATE, SEDR, DeepST) and BayesSpace are among the most robust, and a technology-stratified benchmark (Chen 2025 *iMeta* 4:e70084) finds STAGATE/GraphST best on low-resolution Visium while BASS/stLearn/BANKSY lead on high-resolution platforms.

The caveat: DLPFC is a CONTINUOUS, LAMINAR tissue, which flatters smoothing-friendly methods. Do not over-generalize these rankings to non-laminar tissue. ALL methods struggle with NON-CONTIGUOUS domains (the same region appearing in separated patches -- scattered tumor nests, immune aggregates) because the spatial prior assumes contiguity; for non-contiguous biology, lower the spatial weight or switch to a niche/neighborhood analysis rather than domain segmentation.

## Build the Spatial Graph

**Goal:** Construct the spatial neighbor graph that every domain method inherits.

**Approach:** Use a hex lattice for Visium (6 neighbors) and a generic kNN/Delaunay graph for imaging point clouds; the graph radius IS a spatial-weight knob (too dense over-smooths, too sparse fragments).

```python
import squidpy as sq
import scanpy as sc

adata = sc.read_h5ad('preprocessed.h5ad')

# Visium hex lattice: 6 immediate neighbors. For imaging point clouds use
# coord_type='generic' with KNN or Delaunay. Build in microns where possible --
# a graph built in pixels has a different radius than one built in microns.
sq.gr.spatial_neighbors(adata, coord_type='grid', n_neighs=6)
```

## Expression-Only Clustering Is the Salt-and-Pepper Baseline

**Goal:** Show why a domain method is needed at all.

**Approach:** Cluster on expression PCA with no spatial term; the result is spatially incoherent and demonstrates the failure domain methods correct.

```python
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5, key_added='expr_leiden',
             flavor='igraph', n_iterations=2, directed=False)

# Plot on tissue: speckled, no contiguous regions -- this is the baseline a
# domain method is built to beat, not a domain result.
sq.pl.spatial_scatter(adata, color='expr_leiden')
```

## Neighbor-Augmented Domains (BANKSY-style) and the Spatial-Weight Knob

**Goal:** Produce spatially coherent domains and make the over-smoothing knob explicit.

**Approach:** BANKSY concatenates each cell's own expression with its neighborhood-mean expression, mixed by lambda; lambda~0.2 yields cell typing, lambda~0.8 yields domains. A transparent neighbor-augmented matrix reproduces the idea with Squidpy when the BANKSY package is unavailable -- the lambda here is the load-bearing decision, not the clustering call.

```python
import numpy as np
from sklearn.preprocessing import normalize

# Mean expression over each spot's spatial neighbors (the neighborhood signal).
W = normalize(adata.obsp['spatial_connectivities'], norm='l1', axis=1)
neighbor_mean = W @ adata.obsm['X_pca']

# lambda is the spatial-weight knob: ~0.8 for domains, ~0.2 for cell typing.
# Too high -> blobs (boundaries erased); too low -> salt-and-pepper.
lam = 0.8
augmented = np.concatenate(
    [np.sqrt(1 - lam) * adata.obsm['X_pca'], np.sqrt(lam) * neighbor_mean], axis=1)
adata.obsm['X_banksy'] = augmented

sc.pp.neighbors(adata, use_rep='X_banksy', key_added='banksy')
sc.tl.leiden(adata, resolution=0.5, key_added='domains', neighbors_key='banksy',
             flavor='igraph', n_iterations=2, directed=False)
sq.pl.spatial_scatter(adata, color='domains')
```

If the BANKSY package is installed, prefer it: `import banksy_py` (Python) or the R `Banksy` Bioconductor package compute the own + neighborhood-mean + azimuthal Gabor (AGF) features and expose lambda directly.

## Sweep the Spatial Weight to Find the Over-Smoothing Edge

**Goal:** Locate the lambda where boundaries are sharp but regions stay contiguous.

**Approach:** Recompute domains across a lambda grid and inspect boundaries against histology; the right lambda is the largest value before distinct regions merge into blobs, judged visually, not by an internal score.

```python
for lam in [0.2, 0.5, 0.8]:
    aug = np.concatenate(
        [np.sqrt(1 - lam) * adata.obsm['X_pca'],
         np.sqrt(lam) * (W @ adata.obsm['X_pca'])], axis=1)
    adata.obsm['X_banksy'] = aug
    sc.pp.neighbors(adata, use_rep='X_banksy', key_added='banksy')
    sc.tl.leiden(adata, resolution=0.5, key_added=f'domains_l{lam}',
                 neighbors_key='banksy', flavor='igraph', n_iterations=2, directed=False)

# Low lambda -> speckled; high lambda -> over-smoothed blobs. Pick by eye
# against known tissue architecture, not by silhouette.
sq.pl.spatial_scatter(adata, color=['domains_l0.2', 'domains_l0.5', 'domains_l0.8'])
```

## Choose k as a Biological Decision, With k+-1 Sensitivity

**Goal:** Set the number of domains to the regionalization the question needs and show the answer is not fragile to it.

**Approach:** When a method takes k directly (BayesSpace q, mclust on a GNN embedding), run k, k-1, k+1 and report all three; do not silently optimize k against a clustering score, which has no biological ground truth.

```python
from sklearn.mixture import GaussianMixture

for k in [6, 7, 8]:  # DLPFC has 6 layers + white matter -> k near 7
    gm = GaussianMixture(n_components=k, covariance_type='full', random_state=0)
    adata.obs[f'domains_k{k}'] = gm.fit_predict(adata.obsm['X_banksy']).astype(str)

# Report k+-1 side by side; a domain that only appears at one k is a weak claim.
sq.pl.spatial_scatter(adata, color=['domains_k6', 'domains_k7', 'domains_k8'])
```

## BayesSpace for Laminar Tissue (R)

**Goal:** Apply an MRF smoothing prior, the robust choice on continuous Visium tissue.

**Approach:** Preprocess, then `spatialCluster` with q domains; the smoothing prior couples neighboring spots. Run in R and import the labels.

```r
library(BayesSpace)

sce <- readRDS('sce.rds')
sce <- spatialPreprocess(sce, platform = 'Visium', n.PCs = 15)
# q is the biological k; nrep is MCMC iterations. spatialEnhance() can split
# spots into subspots for higher resolution after spatialCluster().
sce <- spatialCluster(sce, q = 7, nrep = 10000)
write.csv(data.frame(barcode = colnames(sce), domain = sce$spatial.cluster),
          'bayesspace_domains.csv', row.names = FALSE)
```

## STAGATE for High-Resolution or Large Sections (optional)

**Goal:** Learn a spatially-aware embedding whose attention down-weights cross-boundary edges.

**Approach:** Build the STAGATE radius graph, train the graph-attention autoencoder, then cluster the embedding; set the radius to the over-smoothing knob.

```python
import STAGATE  # optional dependency; pip install STAGATE_pyG or use the TF build

STAGATE.Cal_Spatial_Net(adata, rad_cutoff=150)  # rad_cutoff is the graph radius knob
STAGATE.Stats_Spatial_Net(adata)
adata = STAGATE.train_STAGATE(adata)

sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.leiden(adata, resolution=0.5, key_added='stagate_domains',
             flavor='igraph', n_iterations=2, directed=False)
```

## Name Domains From Markers

**Goal:** Attach anatomical labels to spatially coherent clusters.

**Approach:** Rank per-domain markers, then map cluster IDs to region names; marker p-values from the same data that defined the domains are for ranking and labeling only, not inference (the double-dipping caveat from single-cell/clustering).

```python
sc.tl.rank_genes_groups(adata, groupby='domains', method='wilcoxon')
markers = sc.get.rank_genes_groups_df(adata, group=None)
print(markers.groupby('group').head(5))

domain_names = {'0': 'White matter', '1': 'Layer 1', '2': 'Layer 2/3'}
adata.obs['region'] = adata.obs['domains'].map(domain_names)
sq.pl.spatial_scatter(adata, color='region')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Speckled, spatially incoherent labels | No spatial term (plain Leiden on expression) or spatial weight too low | Use a domain method; raise lambda / smoothing / graph density |
| Regions merged into smooth blobs, boundaries gone | Over-smoothing -- spatial weight too high or graph radius too large | Lower lambda / smoothing strength / rad_cutoff; check boundaries against histology |
| Scattered tumor nests collapse into one domain or vanish | Non-contiguous domain; spatial prior assumes contiguity (Yuan 2024) | Lower the spatial weight, or switch to niche/neighborhood analysis (spatial-statistics) |
| Domain count feels arbitrary / reviewer questions k | k optimized against a silhouette score instead of chosen biologically | Set k from the biology; report k+-1 sensitivity; justify the regionalization |
| Domains track a single sample/section | Batch confounded with biology in multi-sample data | Use a multi-sample method (BASS, PRECAST, GraphST integration) or integrate first |
| Answer changes with no parameter change | Spatial graph rebuilt with different units (pixels vs microns) or different kNN/Delaunay | Pin the graph construction and coordinate units; build once, reuse |
| "Domain" answer to a "which types co-occur" question | Domain segmentation used for a niche question | Use neighborhood enrichment / co-occurrence (spatial-statistics) instead |
| Marker p-values quoted as proof a domain is real | Double-dipping (testing the clustering that defined the groups) | Use markers for ranking/labeling only; validate regions against known architecture |

## Related Skills

- spatial-neighbors - Build and tune the spatial graph every domain method inherits
- spatial-statistics - Niche/neighborhood enrichment and co-occurrence when the question is which cell types co-occur, not which region this is
- spatial-deconvolution - Per-spot cell-type composition; domains are regions, deconvolution is composition within a spot
- single-cell/clustering - Non-spatial clustering, resolution sweeps, and the double-dipping caveat on post-clustering marker tests

## References

- Zhao et al. (2021). Spatial transcriptomics at subspot resolution with BayesSpace. Nat Biotechnol 39:1375-1384.
- Hu et al. (2021). SpaGCN: integrating gene expression, spatial location and histology to identify spatial domains and SVGs by graph convolutional network. Nat Methods 18:1342-1351.
- Dong & Zhang (2022). Deciphering spatial domains from spatially resolved transcriptomics with an adaptive graph attention auto-encoder (STAGATE). Nat Commun 13:1739.
- Long et al. (2023). Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST. Nat Commun 14:1155.
- Singhal et al. (2024). BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis. Nat Genet 56:431-441.
- Kim et al. (2022). Unsupervised discovery of tissue architecture in multiplexed imaging (UTAG). Nat Methods 19:1653-1661.
- Xu et al. (2024). Unsupervised spatially embedded deep representation of spatial transcriptomics (SEDR). Genome Med 16:12.
- Li & Zhou (2022). BASS: multi-scale and multi-sample analysis enables accurate cell type clustering and spatial domain detection in spatial transcriptomic studies. Genome Biol 23:168.
- Yuan et al. (2024). Benchmarking spatial clustering methods with spatially resolved transcriptomics data. Nat Methods 21:712-722.
- Chen et al. (2025). A comprehensive benchmarking for spatially resolved transcriptomics clustering methods across variable technologies, organs, and replicates. iMeta 4:e70084.
- Maynard et al. (2021). Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex. Nat Neurosci 24:425-436.
- Palla et al. (2022). Squidpy: a scalable framework for spatial omics analysis. Nat Methods 19:171-178.
