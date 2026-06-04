---
name: bio-imaging-mass-cytometry-spatial-analysis
description: Analyze spatial cell-cell interactions, neighborhoods, and niches in IMC/MIBI data with squidpy and imcRtools, covering neighborhood-enrichment permutation nulls, the abundance-vs-density confound, inhomogeneous Ripley's K, cellular-neighborhood discovery, graph-construction (contact vs proximity), and edge effects. Use when testing whether cell types co-locate, choosing a spatial null, building a neighbor graph, discovering tissue niches, or deciding whether a spatial pattern is real or a density/segmentation artifact.
tool_type: python
primary_tool: squidpy
---

## Version Compatibility

Reference examples tested with: squidpy 1.3+, scanpy 1.10+, anndata 0.10+, numpy 1.26+, imcRtools 1.8+ (R), spatstat 3.0+ (R)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: squidpy's `nhood_enrichment` z-score is unbounded and scales with graph degree, so it is NOT comparable across images of different cell counts. squidpy does not implement Ripley's K (its `ripley` and `co_occurrence` are loose analogs); the inhomogeneous K that conditions on tissue density lives in R `spatstat::Kinhom`/`Kcross.inhom`. Build the graph with `gr.spatial_neighbors(coord_type='generic')` before any test.

# Spatial Analysis for IMC

**"Analyze spatial cell interactions in my IMC data"** -> Build a neighbor graph, then test cell-type co-location against an explicitly chosen null, or discover recurrent niches.
- Python: `squidpy.gr.spatial_neighbors`, `squidpy.gr.nhood_enrichment`, `squidpy.gr.co_occurrence`
- R: `imcRtools::buildSpatialGraph`, `spatstat::Kcross.inhom`

## The Single Most Important Modern Insight -- a spatial interaction is a hypothesis test, and the null silently decides what is measured

A "spatial interaction" or "niche" is not an observation; it is a test against a null model, and which null is chosen (label-shuffle vs CSR vs density-conditioned vs patient-level) decides whether the result is real biology, a density gradient, a segmentation artifact, or upstream clustering. The label-permutation null (histoCAT, squidpy `nhood_enrichment`) holds cell positions fixed and shuffles identities: it DOES control global abundance (a rare type cannot show spurious enrichment just from being rare) but it does NOT control local density or tissue architecture -- so two cell types that merely share an anatomical compartment (both enriched in a follicle or at an invasive margin) score as strongly "interacting" with no direct affinity. This density confound is the dominant false-positive engine, and a finding significant under CSR but null under inhomogeneous Ripley's K is a density artifact, not an interaction. Two further facts compound it: the `nhood_enrichment` z-score is unbounded and graph-degree-dependent, so a z of 30 in a 50k-cell image and a z of 8 in a 5k-cell image are not on the same scale (never threshold a fixed z across heterogeneous images); and the replicate is the patient, not the cell or the image, so a per-cell test over tens of thousands of cells reports p~0 for trivial effects (pseudoreplication). The skill forces two questions onto every analysis: against which null, and at which unit.

## Spatial Methods Taxonomy

| Method | Measures | Null / inference | Confound it does NOT control |
|--------|----------|------------------|------------------------------|
| histoCAT / squidpy `nhood_enrichment` | A neighbors B more/less than chance | within-image label shuffle; z-score | local density, architecture, edges |
| squidpy `interaction_matrix` | raw cluster-cluster edge counts | none (descriptive) | everything -- just counts |
| squidpy `co_occurrence` | P(B \| A at distance d) / P(B) | none (descriptive curve) | needs its own CSR null |
| Ripley's K/L, cross-K, K_inhom | clustering vs dispersion over radii | CSR; K_inhom conditions on intensity | tissue inhomogeneity (unless K_inhom); ROI shape |
| Cellular Neighborhoods (CN) | recurrent local compositions (niches) | none -- clustering, not a test | window size = scale; doubly-derived |
| Spatial-LDA / UTAG | microenvironment topics / tissue domains | generative / unsupervised | topic/domain count arbitrary |
| Moran's I / Geary's C | spatial autocorrelation of a continuous mark | permutation / analytic | graph definition; global stat masks local |

## Decision Tree by Scenario

| Question | Method | Why |
|----------|--------|-----|
| Do two named types co-locate more than chance? (fast screen, one image) | `nhood_enrichment` / histoCAT | abundance-aware permutation z |
| Same, but worried about density/architecture | inhomogeneous cross-K (`Kcross.inhom`) | conditions on the tissue's own intensity surface |
| At what scale is a type clustered/dispersed? | Ripley's K/L over radii (edge-corrected) | second-order structure across distance |
| What recurrent multicellular niches exist? (discovery) | Cellular Neighborhoods (sweep window k) | recurrent compositions; no built-in test |
| Is a continuous marker/score spatially structured? | Moran's I (global) / Geary's C (local) | autocorrelation |
| Does an interaction/niche differ between conditions? | hand off to differential-analysis | per-image summary -> patient unit -> mixed model + FDR |

## Build the Spatial Graph (contact vs proximity)

**Goal:** Construct the graph whose definition matches the biological claim.

**Approach:** Delaunay approximates physical adjacency (contact/juxtacrine) but invents long edges across lumen/necrosis, so prune by a max distance; fixed radius gives true proximity (paracrine) at a stated micron scale; kNN silently mixes contact and proximity because fixed k spans microns in dense regions and hundreds of microns in sparse ones. Build the graph per image.

```python
import squidpy as sq
import numpy as np

# contact graph: Delaunay, then prune edges longer than a biological max distance (um)
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
dist = adata.obsp['spatial_distances']
keep = dist.copy(); keep.data[keep.data > 30] = 0; keep.eliminate_zeros()   # cap at ~30 um
adata.obsp['spatial_connectivities'] = (keep > 0).astype(float)

# OR proximity graph: fixed radius at a justified micron scale (paracrine range)
# sq.gr.spatial_neighbors(adata, coord_type='generic', radius=30.0)
```

## Neighborhood Enrichment with a Named Null

**Goal:** Test co-location per image with the abundance-aware permutation null, knowing its blind spot.

**Approach:** Run `nhood_enrichment` per image (so the shuffle is within-image), keep the z as a per-image summary, and never threshold a fixed z across images of different size. Cross-check density-driven hits against inhomogeneous K.

```python
per_image_z = {}
for img_id, idx in adata.obs.groupby('image_id').groups.items():
    sub = adata[idx].copy()
    sq.gr.spatial_neighbors(sub, coord_type='generic', delaunay=True)
    sq.gr.nhood_enrichment(sub, cluster_key='cell_type', seed=0)
    per_image_z[img_id] = sub.uns['cell_type_nhood_enrichment']['zscore']
# aggregate these per-image summaries to the PATIENT unit in differential-analysis, not here
```

## Cellular Neighborhoods (niche discovery)

**Goal:** Find recurrent local cell-type compositions, treating them as exploratory.

**Approach:** Per cell, summarize the composition of its window of neighbors, then cluster the windows. The window size IS the spatial scale and is almost always unjustified, so sweep it and report that the biological conclusion survives k in {10, 20, 30}. A CN has no built-in significance test; significance enters only as a cross-condition comparison (differential-analysis).

```python
from sklearn.cluster import KMeans
import pandas as pd

def cellular_neighborhoods(adata, k_window=20, n_cn=10):
    sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=k_window)   # window = k neighbors
    conn = adata.obsp['spatial_connectivities']
    onehot = pd.get_dummies(adata.obs['cell_type']).values
    comp = (conn @ onehot)                       # neighbor composition per cell
    comp = comp / comp.sum(axis=1, keepdims=True).clip(min=1)
    return KMeans(n_clusters=n_cn, random_state=0).fit_predict(comp)

adata.obs['CN'] = cellular_neighborhoods(adata, k_window=20)   # sweep k_window to check stability
```

## Per-Method Failure Modes

### nhood_enrichment -- density false positive
**Trigger:** reporting a z-score as a "significant interaction." **Mechanism:** the label-shuffle null absorbs abundance but not local density, so co-compartmentalized types score as interacting. **Symptom:** every type pair sharing a region looks attracted. **Fix:** name the null; cross-check with inhomogeneous cross-K; or shuffle within a compartment, not the whole image.

### Fixed z threshold across images
**Trigger:** calling z>2 "significant" across images of different cell counts. **Mechanism:** the z is unbounded and scales with graph degree. **Symptom:** large images dominate; small images never reach threshold. **Fix:** rank within image or convert to an effect size; aggregate per-image summaries to the patient unit.

### Unpruned Delaunay / kNN mislabeled as contact
**Trigger:** kNN graph results called "contact," or Delaunay across tissue gaps. **Mechanism:** fixed k mixes contact and proximity by density; Delaunay invents long edges across voids. **Symptom:** phantom interactions across lumen/necrosis. **Fix:** prune Delaunay by a max distance; use fixed radius for paracrine claims; state the micron scale.

### CSR Ripley's K on inhomogeneous tissue
**Trigger:** homogeneous-Poisson K on structured tissue. **Mechanism:** CSR assumes constant intensity, so everything tests as clustered. **Symptom:** universal "clustering". **Fix:** inhomogeneous K (`Kinhom`/`Kcross.inhom`) that divides by the estimated intensity surface.

### Niche taken as a discovered fact
**Trigger:** treating a CN as established biology. **Mechanism:** it is doubly-derived (cluster cells -> cluster windows) with an arbitrary window and k, and no significance test. **Symptom:** different labs report different niches on the same tissue. **Fix:** sweep window k; report stability (NMI/ARI); validate niche-defining markers against raw images for spillover.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| CN window = 10 nearest neighbors, 9 CNs retained | Schurch 2020 *Cell* 182:1341 | the canonical CN convention -- sweep it, do not adopt blindly |
| CODEX i-niches: Delaunay 1st-tier, k-means = 100 | Goltsev 2018 *Cell* 174:968 | window/scale is a choice, not a default |
| n_perm >= 10,000 for small corrected p | Schapiro 2017 *Nat Methods* 14:873 | n_perm=1000 floors p at ~1/1001, too coarse after FDR |
| BH-FDR across ~C(C+1)/2 type pairs x radii | multiplicity | ~200 pairs at p<0.05 guarantees false positives |
| Prune Delaunay at a biological max distance (e.g. ~30 um) | graph hygiene | removes edges across acellular voids |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "Significant interaction" that is two types in one compartment | density confound under label-shuffle | inhomogeneous cross-K; shuffle within compartment |
| p~0 across 50k cells | cell-level pseudoreplication | per-image summary -> patient unit (differential-analysis) |
| Niche changes with the window/k | window IS the scale | sweep k_window and n_cn; report stability |
| Boundary cells dominate a small ROI | ignored edge effects | edge-corrected K; erode/buffer the ROI interior |
| Many "significant" pairs | no multiple-testing correction | BH-FDR across all pairs and radii |

## References

- Schapiro D, Jackson HW, Raghuraman G, et al. 2017. histoCAT: analysis of cell phenotypes and interactions in multiplex image cytometry data. *Nat Methods* 14(9):873-876. — neighborhood permutation test.
- Goltsev Y, Samusik N, Kennedy-Darling J, et al. 2018. Deep Profiling of Mouse Splenic Architecture with CODEX Multiplexed Imaging. *Cell* 174(4):968-981.e15. — cellular neighborhoods (i-niches).
- Schurch CM, Bhate SS, Barlow GL, et al. 2020. Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front. *Cell* 182(5):1341-1359.e19. — canonical CN convention.
- Palla G, Spitzer H, Klein M, et al. 2022. Squidpy: a scalable framework for spatial omics analysis. *Nat Methods* 19(2):171-178. — squidpy spatial functions.
- Chen Z, Soifer I, Hilton H, Keren L, Jojic V. 2020. Modeling Multiplexed Images with Spatial-LDA Reveals Novel Tissue Microenvironments. *J Comput Biol* 27(8):1204-1218. — Spatial-LDA.
- Kim J, Rustam S, Mosquera JM, et al. 2022. Unsupervised discovery of tissue architecture in multiplexed imaging. *Nat Methods* 19(12):1653-1661. — UTAG domains.
- Ripley BD. 1977. Modelling Spatial Patterns. *J R Stat Soc Series B* 39(2):172-212. — the K-function.
- Moran PAP. 1950. Notes on Continuous Stochastic Phenomena. *Biometrika* 37(1/2):17-23. — Moran's I.
- Geary RC. 1954. The Contiguity Ratio and Statistical Mapping. *Incorporated Statistician* 5(3):115-145. — Geary's C.

## Related Skills

- phenotyping - cell-type labels are the input to every spatial test
- differential-analysis - testing whether interactions/niches differ between conditions at the patient level
- cell-segmentation - over-segmentation and lateral spillover create fake niches
- data-preprocessing - uncompensated spillover manufactures false cell-cell interactions
- spatial-transcriptomics/spatial-statistics - shared squidpy neighborhood and autocorrelation methods
