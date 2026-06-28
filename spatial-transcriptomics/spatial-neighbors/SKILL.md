---
name: bio-spatial-transcriptomics-spatial-neighbors
description: Build the spatial neighbor graph that every downstream spatial statistic (Moran's I, neighborhood enrichment, co-occurrence, spatial domains) inherits, using Squidpy. Use when choosing the graph type (kNN vs Delaunay vs fixed-radius vs Visium hex grid) and understanding why it silently changes every downstream result; handling variable cell density (kNN fixes neighbor COUNT, fixed-radius fixes physical DISTANCE -- each distorts the other); getting coordinate units right (pixels vs microns; Visium array coords are not distance); pruning Delaunay long edges across tissue gaps; running the graph sensitivity analysis almost nobody runs; and knowing when planar section neighbors misrepresent a 3D tissue.
tool_type: python
primary_tool: squidpy
---

## Version Compatibility

Reference examples tested with: squidpy 1.4+, scanpy 1.10+, anndata 0.10+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Spatial Neighbor Graphs

**"Build a spatial neighbor graph for my tissue"** -> Define which cells or spots count as spatial neighbors, encoded as a sparse weights matrix W over `adata.obsm['spatial']`.
- Python: `squidpy.gr.spatial_neighbors()` -> writes `adata.obsp['spatial_connectivities']` and `adata.obsp['spatial_distances']`

The platform-class fork sets the default geometry. Sequencing/capture data on a fixed lattice (Visium hex, Visium HD grid) has a KNOWN adjacency -> use `coord_type='grid'`. Imaging/in-situ point clouds (Xenium, MERFISH, CosMx) have irregular cell positions -> use `coord_type='generic'` with Delaunay or kNN. The first question is always which side of the fork the data is on, because it decides whether "neighbor" is a lattice fact or a modeling choice.

## Governing Principle

The graph is the model that all spatial statistics inherit. Moran's I, Geary's C, neighborhood enrichment, co-occurrence, and every graph-neural-network spatial-domain method are functions `f(expression, W)` of the spatial weights matrix W -- the graph is not preprocessing, it is literally an argument to the statistic. Change k from 6 to 30, switch Delaunay to kNN, or row-standardize W instead of leaving it binary, and the Moran's I value, the enrichment z-scores, the SVG ranking, and the domain boundaries all move. The analyst is never measuring "spatial structure"; they are measuring spatial structure as seen through one particular definition of adjacency.

This is the most under-reported researcher degree of freedom in the field, and it has a name geographers settled decades ago: the modifiable areal unit problem -- aggregate or re-adjacency the units and the answer changes. There is no canonical W. The honest workflow therefore does what almost no paper does: build the graph under at least two definitions, rerun the downstream statistic, and report which genes/pairs/domains are graph-robust versus graph-fragile. A result that survives only one graph choice is a result about that graph, not about the tissue.

Two failure directions bound the choice. Too dense a graph (large k, large radius, many rings) over-smooths -- it inflates apparent autocorrelation, washes out local detail, and merges distinct domains into blobs. Too sparse a graph fragments the tissue into disconnected components, flags spurious local outliers, and misses real medium-scale structure. The right density is the one whose downstream conclusion is stable; the specific k is not the deliverable, the stability across k is.

## The Graph-Construction Decision

Each family silently assumes something different about the tissue, and that assumption -- not the algorithm -- is what fails.

| Graph type | Degree behavior | Density bias | Best when | Fails when |
|---|---|---|---|---|
| Visium hex grid (`coord_type='grid'`, `n_rings`) | Fixed 6 per ring; known lattice | None (regular lattice) | Visium / capture grids where geometry is fixed and exact | A spot is treated as a cell -- it is a 1-10-cell mixture, so spot adjacency mixes deconvolution error with real contact |
| kNN (`coord_type='generic'`, `n_neighs`) | Constant COUNT k | Radius implicitly stretches in sparse regions, connecting distant cells | Single-cell platforms (Xenium/MERFISH/CosMx) where fixed degree is wanted | Density varies sharply; asymmetric by default (A is B's neighbor but not vice versa) |
| Fixed-radius (`radius=r`) | Variable -- more neighbors where dense | STRONG: dense regions get more neighbors of EVERYTHING -> inflated enrichment that is pure density artifact | A real physical interaction range exists (ligand diffusion ~tens of microns) AND density is ~uniform | Density gradients; radius set in the wrong coordinate unit |
| Delaunay (`delaunay=True`) | Variable; parameter-free "who touches whom" | Mild | Single-cell data wanting a parameter-free contact graph | Tissue has gaps/holes/folds -> long spurious edges leap across empty space; needs distance pruning |

Three points the naive analyst misses. `squidpy.gr.nhood_enrichment` builds NO graph of its own -- it consumes whatever graph `spatial_neighbors` already stored in `obsp` (and errors if none exists); the Squidpy non-grid default is kNN with `n_neighs=6` (`delaunay=False`), so a published z-score is specific to whichever graph produced it and would change under a different graph -- always know which graph produced the number. kNN is asymmetric; "mutual kNN" (edge only if both cells are in each other's k-set) is more conservative and stops hub cells in dense regions from dominating. Unpruned Delaunay over tissue with necrotic holes or folds connects cells micrometers apart on the slide but biologically unrelated -- pruning by a max edge length is the standard fix.

## Variable Cell Density: There Is No Free Lunch

A fixed-radius neighborhood gives dense regions more neighbors and sparse regions fewer. Because neighborhood enrichment, co-occurrence, and local statistics all depend on neighbor COUNTS, a pure density gradient masquerades as biological signal: a dense lymphoid follicle gets inflated "enrichment" of everything simply because every cell there has more neighbors. kNN fixes the count (it adapts the radius to local density) but then distorts physical distance -- a "neighbor" in sparse stroma may sit far away. The density structure of the tissue dictates which distortion is tolerable: use kNN/Delaunay when density varies (the common case in real tissue); reserve fixed-radius for roughly uniform tissue where an absolute physical interaction range is the actual biological question.

## The Coordinate-Unit Trap

A radius, a co-occurrence interval, and a Delaunay pruning cutoff are all in PHYSICAL distance units. If `adata.obsm['spatial']` holds pixels, array row/col indices, or arbitrary units, a "50-unit radius" is silently meaningless. Visium array row/col is a lattice index, not microns; full-resolution Visium pixel coordinates need the Space Ranger scale factor (`spot_diameter_fullres`, `tissue_hires_scalef`) to convert to physical distance. Imaging platforms store microns or pixels depending on the reader. Confirm the unit BEFORE setting any distance parameter -- the single cheapest check is to measure nearest-neighbor spacing and compare it to the known platform pitch (Visium ~100 microns center-to-center).

**Goal:** Confirm the coordinate unit so that any radius is physically meaningful.

**Approach:** Measure median nearest-neighbor distance from a temporary kNN graph and compare it to the known platform geometry; a Visium grid in microns should read ~100, in pixels it reads hundreds-to-thousands.

```python
import squidpy as sq
import scanpy as sc
import numpy as np

sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=1)   # nearest neighbor only, just to read spacing
nn = adata.obsp['spatial_distances'].data
print(f'median nearest-neighbor spacing: {np.median(nn):.1f} units')
# Visium pitch is ~100 microns; a value of hundreds-to-thousands means coords are in PIXELS -> rescale or use grid mode
```

## Build the Graph by Platform Class

**Goal:** Construct the adjacency that matches the platform geometry rather than a one-size default.

**Approach:** Use grid mode for Visium hex (the lattice is exact and known); use generic Delaunay or kNN for imaging point clouds; store under named keys so multiple graphs coexist for the sensitivity check below.

```python
# Visium hex lattice: 6 immediate neighbors per ring; n_rings=2 widens the neighborhood deliberately
sq.gr.spatial_neighbors(adata, coord_type='grid', n_neighs=6, n_rings=1, key_added='visium_hex')

# Imaging point cloud, fixed-degree: constant k, radius adapts to local density
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=10, key_added='knn10')

# Imaging point cloud, parameter-free contact graph (delaunay=True is opt-in; the
# generic default is kNN with n_neighs=6)
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True, key_added='delaunay')
```

## Prune Delaunay Long Edges Across Tissue Gaps

**Goal:** Stop Delaunay from inventing long-range "neighbors" that leap across necrotic holes, folds, or slide background.

**Approach:** Build Delaunay, then prune to a physically sensible maximum edge length using `radius` as a `(min, max)` interval -- edges longer than `max` (in microns) are dropped.

```python
# radius as a (min, max) tuple prunes the graph to edges within that physical-distance interval;
# choose max from the tissue: a few cell diameters (e.g. 50 microns) kills cross-gap edges, keeps true contacts
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True, radius=(0.0, 50.0), key_added='delaunay_pruned')

pruned = adata.obsp['delaunay_pruned_connectivities']
print(f'edges after pruning: {pruned.nnz}; mean degree: {pruned.nnz / adata.n_obs:.1f}')
```

## Run the Graph Sensitivity Analysis (the one almost nobody runs)

**Goal:** Decide whether a downstream conclusion is a property of the tissue or an artifact of the graph choice.

**Approach:** Build the graph under several adjacency definitions, store each under its own key, then recompute the downstream statistic on each and flag results that are not stable across graphs.

```python
graphs = {}
for k in (6, 15, 30):
    sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=k, key_added=f'knn{k}')
    graphs[f'knn{k}'] = adata.obsp[f'knn{k}_connectivities']
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True, key_added='delaunay')
graphs['delaunay'] = adata.obsp['delaunay_connectivities']

# Recompute the downstream statistic per graph (Moran's I shown); compare rankings, not single values.
# A gene/pair/domain that only appears under one graph is graph-fragile -- report it as such.
import scanpy as sc
for name, W in graphs.items():
    adata.obsp['spatial_connectivities'] = W            # spatial_autocorr reads the active 'spatial' graph
    adata.obsp['spatial_distances'] = adata.obsp[f'{name}_distances'] if f'{name}_distances' in adata.obsp else adata.obsp['spatial_distances']
    sq.gr.spatial_autocorr(adata, mode='moran', genes=adata.var_names[:50].tolist())
    adata.uns[f'moranI_{name}'] = adata.uns['moranI'].copy()
```

## Inspect the Graph Before Trusting It

**Goal:** Catch fragmentation (too sparse) and over-connection (too dense) before they corrupt every downstream number.

**Approach:** Summarize degree distribution and connected components; a healthy graph is one connected component with a tight degree distribution, not many islands or a few hub cells.

```python
import numpy as np
conn = adata.obsp['spatial_connectivities']
degree = np.asarray((conn > 0).sum(axis=1)).ravel()
print(f'mean degree {degree.mean():.1f}; min {degree.min()}; max {degree.max()}')
# isolated cells (degree 0) signal fragmentation; a heavy max-degree tail signals density-driven hubs
print(f'isolated cells: {(degree == 0).sum()}')
```

## The 2D-Section vs 3D-Tissue Caveat

A tissue section is one ~5-10 micron optical/physical plane of a three-dimensional organ. Two cells that are planar neighbors in the section may be far apart in the intact tissue, and two true 3D neighbors may sit in different sections and never appear adjacent in the graph. Cells truncated at the section's top or bottom surface carry partial transcript profiles (only the captured fraction of the cell), which depresses their counts and distorts their degree. Any neighbor graph built from a single section is a planar slice of the real 3D adjacency -- adequate for in-plane analysis, but it does not license 3D-contact claims. Reconstructing true 3D neighbors from serial sections (registration, z-stacking, alignment across planes) is a different problem that this graph does not solve. Layered, ducted, or crypted tissue is also anisotropic (covariance is direction-dependent), so an isotropic graph that ignores orientation underpowers directional structure -- a caveat to keep when neighbor counts feed directional or layer-aware statistics.

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| Every cell has wildly different neighbor counts; dense regions show "enrichment" of everything | Fixed-radius graph on density-varying tissue -- pure density artifact | Use kNN or Delaunay (constant or contact-based degree); reserve `radius` for ~uniform tissue with a real physical range |
| A radius of 50 captures all cells or none | `adata.obsm['spatial']` is in pixels or array units, not microns | Confirm the unit (measure nearest-neighbor spacing vs platform pitch); apply the Visium scale factor or use `coord_type='grid'` |
| Long edges cross empty space / necrotic holes; spurious long-range neighbors | Unpruned Delaunay over tissue with gaps or folds | Prune with `radius=(0, max)` at a few cell diameters; inspect the overlaid graph |
| Visium neighbors look irregular instead of a clean hex lattice | `coord_type='generic'` used on a Visium grid | Use `coord_type='grid'` with `n_rings`; the lattice adjacency is exact and known |
| Downstream Moran's I / enrichment z-scores change when k is changed | Expected -- the statistic is `f(expression, W)`; the graph is the model | Run the sensitivity analysis across k and Delaunay; report only graph-robust results |
| Graph splits into many connected components | Graph too sparse (k too small, radius too short) | Increase k or radius, or switch to Delaunay; check isolated-cell count |
| Spot-level neighborhood enrichment over-interpreted as cell-cell contact | A Visium spot is a 1-10-cell mixture, not a cell | Treat spot adjacency as spot-level; deconvolve (spatial-transcriptomics/spatial-deconvolution) before cell-level claims |
| 3D-contact conclusion drawn from one section | Planar neighbors are a slice of 3D adjacency; truncated cells have partial profiles | Restrict claims to in-plane; use serial-section reconstruction for true 3D neighbors |

## Related Skills

- spatial-transcriptomics/spatial-statistics - the neighbor graph is the W fed to Moran's I, neighborhood enrichment, and co-occurrence
- spatial-transcriptomics/spatial-domains - graph-based domain methods inherit this adjacency and its over-smoothing/fragmentation tradeoff
- spatial-transcriptomics/spatial-communication - ligand-receptor proximity tests run on this graph and inherit its density bias
- spatial-transcriptomics/spatial-data-io - load coordinates and confirm their unit before building any graph
- single-cell/clustering - expression-space kNN graphs, the non-spatial counterpart

## References

- Palla G, Spitzer H, Klein M, et al. (2022) Squidpy: a scalable framework for spatial omics analysis. Nature Methods 19(2):171-178. DOI 10.1038/s41592-021-01358-2
- Moran PAP (1950) Notes on continuous stochastic phenomena. Biometrika 37(1/2):17-23. DOI 10.1093/biomet/37.1-2.17
- Geary RC (1954) The contiguity ratio and statistical mapping. The Incorporated Statistician 5(3):115-145. DOI 10.2307/2986645
- Getis A, Ord JK (1992) The analysis of spatial association by use of distance statistics. Geographical Analysis 24(3):189-206. DOI 10.1111/j.1538-4632.1992.tb00261.x
- Ripley BD (1977) Modelling spatial patterns. Journal of the Royal Statistical Society Series B 39(2):172-212. DOI 10.1111/j.2517-6161.1977.tb01615.x
- Dos Santos Peixoto R, Miller BF, Brusko MA, et al. (2025) Characterizing cell-type spatial relationships across length scales in spatially resolved omics data. Nature Communications 16:350. DOI 10.1038/s41467-024-55700-1
