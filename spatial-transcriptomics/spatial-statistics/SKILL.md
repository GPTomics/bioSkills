---
name: bio-spatial-transcriptomics-spatial-statistics
description: Detects spatially variable genes, spatial autocorrelation, and cell-type colocalization for spatial transcriptomics using Squidpy with PySAL/esda for local statistics. Use when choosing an SVG method by its null and scaling (SpatialDE/SPARK GP variance-component vs SPARK-X/nnSVG linear vs Moran/Geary graph autocorrelation); separating genes that are spatially variable because of cell-type composition from genes regulated within a cell type; choosing the right autocorrelation statistic (global Moran/Geary vs Getis-Ord hot/cold spots vs local LISA and its FDR trap); and choosing a colocalization null strong enough to defeat the abundance/compartment confound (conditional or toroidal vs the weak Squidpy default permutation).
tool_type: python
primary_tool: squidpy
---

## Version Compatibility

Reference examples tested with: squidpy 1.4+, scanpy 1.10+, anndata 0.10+, esda 2.5+, libpysal 4.9+ (SPARK, SPARK-X, and nnSVG are R/Bioconductor packages)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Spatial Statistics

**"Find spatially variable genes / run Moran's I / test which cell types colocalize"** -> Quantify spatial structure in expression or in cell-type arrangement against an explicit null model.
- Python: `squidpy.gr.spatial_autocorr` (Moran/Geary), `squidpy.gr.nhood_enrichment`, `squidpy.gr.co_occurrence`, `esda.Moran_Local`/`esda.getisord.G_Local` (LISA, Getis-Ord)
- R: SPARK / SPARK-X / nnSVG (SVG), `spdep` (Moran/Geary/Getis-Ord/LISA)

## Governing Principle

Two reframes decide whether a spatial-statistics result means anything. Both are silent failures: the code runs, the numbers look clean, and the interpretation is wrong.

A "spatially variable gene" is not necessarily spatially REGULATED. Moran's I, Geary's C, SPARK-X, and SpatialDE all detect that a gene's expression is spatially autocorrelated -- but a gene that is simply a marker of a spatially clustered cell type scores as "spatially variable" with zero cell-intrinsic spatial regulation, purely because the CELL TYPE is spatially organized. A hepatocyte gene in zonated liver, MAG in white matter, KRT17 in epithelium: all top SVGs, none of them regulated in space. This cell-type-driven signal swamps the genuinely interesting cell-type-INDEPENDENT signal (a gene graded across a niche WITHIN one cell type). Sample-wide SVG lists therefore largely re-derive marker genes and overlap heavily with HVGs -- if the SVG list is roughly the HVG list, the spatial test added almost nothing. The interesting question (within-type spatial regulation) needs cell-type-aware methods (C-SIDE, CTSV, CELINA/Celina), which are themselves unsettled and carry their own false positives. Decision: if the question is "where is tissue organized," sample-wide SVG is correct (the cell-type structure IS the answer); if the question is "which genes are regulated beyond cell identity," sample-wide SVG is the WRONG tool -- test within cell type or regress out composition first.

The null and the graph define the result. Every spatial statistic is computed on a neighbor graph (a weights matrix W) against a null distribution, and both are researcher choices, not properties of the tissue. Change kNN k from 6 to 30 and Moran's I, the enrichment z-scores, and the SVG ranking all move. Change the colocalization null from global label-shuffle to within-compartment and most "A is near B" claims evaporate. SVG methods disagree heavily across the literature precisely because each tests a DIFFERENT null (variance-component-zero vs covariance-independence vs graph-autocorrelation-zero) -- that disagreement is expected, not a bug. The honest workflow names the null, names the graph (cross-ref spatial-neighbors), and reports whether a hit survives a second graph or a stronger null.

## Choosing an SVG Method

**Goal:** Pick a spatially-variable-gene test whose null hypothesis and computational scaling match the platform and the biological question.

**Approach:** Match GP variance-component methods to small Gaussian/count data, linear-scaling methods to single-cell-resolution data, and treat the SVG list as method-conditional -- cross-method intersection is more trustworthy than any single ranking, though it is small.

| Method | Null it tests | Scaling | Best when | Fails when |
|--------|---------------|---------|-----------|------------|
| SpatialDE (GP) | spatial variance component = 0 at the tested length scale | O(n^3); infeasible past ~1e4 locations | Small Visium-scale, Gaussian on log-normalized | Sparse/low counts violate Gaussian; fixed length-scale grid misses other scales |
| SPARK (count GLSM) | no pattern matching any of 10 fixed kernels | O(n^3)-ish (PQL); slow at large n | Small count data; want Poisson model, not Gaussian | Pattern unlike its 10 kernels; still cell-type-confounded |
| SPARK-X | expression covariance independent of location covariance | LINEAR in n and genes | 1e4-1e6 cells (MERFISH/Xenium/CosMx) needing scalability | Fixed location kernels miss unusual length scales; low power on small focal hotspots |
| nnSVG (NNGP) | spatial variance = 0, with a per-gene length scale | LINEAR in n | Length scales genuinely differ across genes; large single-cell data | Still cell-type-confounded; more compute per gene than SPARK-X; needs adequate counts |
| Moran's I / Geary's C | no autocorrelation on graph W | Fast (sparse W) | Quick screen on an existing neighbor graph | Single fixed scale (the graph); misses multi-focal/small hotspots |

There is no uniformly best SVG method; power is pattern-specific (SPARK-X, nnSVG, and Moran's I all have LOW power for genes high in small focal areas). Methods evolve fast -- verify the current benchmark before committing. Threshold on EFFECT SIZE (fraction of spatial variance), not p alone: with thousands of cells, trivial autocorrelation reaches tiny p-values.

## Computing Spatial Autocorrelation with Squidpy

**Goal:** Rank genes by graph-based spatial autocorrelation as a fast SVG screen.

**Approach:** Build a neighbor graph, run Moran's I (or Geary's C) per gene with a permutation/analytic p-value and FDR, then read effect size before significance.

```python
import squidpy as sq
import scanpy as sc

# The graph IS the model: k, coord_type, and units all change the result (see spatial-neighbors)
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)   # 6 mimics the Visium hex lattice

# genes=None defaults to highly_variable if present; n_perms adds a permutation null, corr_method applies FDR
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100, corr_method='fdr_bh')   # statsmodels name, not 'benjamini-hochberg'
moran = adata.uns['moranI']                  # columns: I, pval_norm, pval_norm_fdr_bh, ...

# Threshold on effect size (I) AND FDR, not p alone -- large n makes trivial autocorrelation 'significant'
svg = moran[(moran['I'] > 0.1) & (moran['pval_norm_fdr_bh'] < 0.05)].sort_values('I', ascending=False)
```

A top-ranked SVG here is a hypothesis about spatial structure, NOT evidence of spatial regulation. Before interpreting, compare the SVG list to the HVG list: the overlap is cell-type marker genes; the SVG-not-HVG subset (modest-amplitude gradients) is where spatial information actually lives.

## Choosing an Autocorrelation Statistic

**Goal:** Match the statistic to the spatial question -- "is this gene clustered" is a different question from "where is it HIGH" and from "is THIS region a cluster."

**Approach:** Use a global statistic for one tissue-wide number, Getis-Ord when the sign (hot vs cold) matters, and local LISA for non-stationary tissue -- but pay the local multiple-testing tax correctly.

| Statistic | Global / local | Hot vs cold? | Use when |
|-----------|----------------|--------------|----------|
| Moran's I | global | NO (clustering of like values only) | One number: is this gene spatially structured across the whole section |
| Geary's C | global | NO | Same as Moran but more sensitive to LOCAL/short-range differences; disagreement with Moran is informative |
| Getis-Ord Gi* | local | YES -- separates high-clusters from low-clusters | "Where is this gene HIGH" -- hot/cold spot mapping |
| Local Moran / LISA | local | partial (HH/LL/HL/LH quadrants) | Non-stationary tissue: per-location clusters and spatial outliers |

Moran's I and Geary's C cannot tell a hot spot from a cold spot -- both flag "similar values cluster" regardless of high or low. Choosing Moran when the question is "where is this gene HIGH" is a category error; use Getis-Ord Gi*. A non-significant GLOBAL Moran's I does NOT mean "no spatial structure": over heterogeneous tissue, positive autocorrelation in one region cancels negative in another, so use local statistics for non-stationary sections.

```python
from esda.getisord import G_Local
from esda.moran import Moran_Local
from libpysal.weights import KNN

coords = adata.obsm['spatial']
w = KNN.from_array(coords, k=6)
w.transform = 'r'                            # row-standardized; changes the value AND its variance vs binary W

gene = adata[:, 'GENE1'].X.toarray().ravel()
gi = G_Local(gene, w, transform='B', star=True, permutations=999)   # star=True -> Gi* (includes self); binary weights for Getis-Ord
lisa = Moran_Local(gene, w, permutations=999)            # conditional-permutation local null
adata.obs['GENE1_hotspot'] = gi.Zs                       # positive Z = hot spot, negative = cold spot
adata.obs['GENE1_lisa_q'] = lisa.q                       # 1=HH, 2=LH, 3=LL, 4=HL
```

Local statistics carry a DOUBLE trap. There are n tests (one per location), so uncorrected LISA/Gi* maps are mostly false positives -- FDR is mandatory, and Anselin recommends stricter base cutoffs (0.01/0.005/0.001), not 0.05. Worse, the local statistics are themselves spatially autocorrelated (adjacent locations share neighbors, so adjacent I_i values are correlated), which violates the independence assumption of standard BH-FDR; the effective number of tests is far below n. Conditional permutation gives the correct local null but does not fix the cross-location dependence. Treat the cluster map as exploratory, not a set of independent discoveries.

## Testing Cell-Type Colocalization

**Goal:** Decide whether two cell types are SPECIFICALLY associated in space, not merely both abundant or both in the same compartment.

**Approach:** Choose a permutation null strong enough to defeat the abundance/compartment confound; the Squidpy default answers only the weak question, and a co-occurrence distance profile is more informative than a single z-score.

| Null model | What it permutes | Controls for | Misses |
|------------|------------------|--------------|--------|
| Global label permutation (Squidpy `nhood_enrichment` default) | all labels over all positions | graph topology, marginal counts | tissue compartmentalization -- two abundant co-compartment types pass trivially |
| Conditional / within-compartment permutation | labels within a region only | shared-compartment forcing | cross-compartment questions; the region choice is itself a degree of freedom |
| Toroidal shift | whole label field translated (wrapped) | each type's first-order density pattern | anisotropy; boundary realism (wrapping a bounded tissue is artificial) |
| Grid-tile shuffle across scales (CRAWDAD) | labels within tiles of size s | structure above scale s -- isolates colocalization AT scale s | within-tile structure below s |

The single most common error in spatial-omics papers is reading a positive `nhood_enrichment` z-score as a specific A-B interaction. Under the weak global-permutation null it usually reflects co-compartmentalization plus abundance: two stromal populations, or tumor plus tumor-associated macrophages both in the tumor bed, pass trivially. A specific-affinity claim must SURVIVE a stronger null (conditional/within-compartment, toroidal shift preserving each type's density, or CRAWDAD's scale-explicit tiles). Rare-type enrichment is the least trustworthy: few edges give high-variance, often spuriously large |z| -- be most skeptical exactly where the biology is most exciting (a rare type near the tumor).

For clustering as a function of distance rather than a single graph z, Ripley's K/L (`squidpy.gr.ripley`, mode `'L'`) counts within-cluster neighbors within radius r against a complete-spatial-randomness expectation, so it reads as clustering vs dispersion ACROSS scale per cell type (squidpy computes the univariate L per cluster; a true bivariate cross-K answering "are A and B closer than chance, at what radius" needs a dedicated point-process tool such as spatstat). Its assumptions are the geostatistics ones tissue violates: a bounded, holey, non-stationary window. EDGE CORRECTION is mandatory and usually omitted -- without it, counts near the tissue boundary or a necrotic hole are biased DOWN and read as false depletion, so an ROI with gaps needs an edge-corrected estimator (or restrict analysis to the interior).

```python
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)   # nhood_enrichment reuses this stored graph

# Global label-permutation null -- the WEAK question: more adjacent than complete spatial randomness?
sq.gr.nhood_enrichment(adata, cluster_key='cell_type', n_perms=1000)
z = adata.uns['cell_type_nhood_enrichment']['zscore']

# co_occurrence gives a DISTANCE PROFILE (at what scale colocalization appears/vanishes), unlike a single z
sq.gr.co_occurrence(adata, cluster_key='cell_type')
```

Cellular Neighborhoods (Schurch/Nolan: cluster per-cell windows of neighbor composition) have NO inferential null at all -- they are descriptive k-means clusters whose number and window size are user knobs. They are useful summaries but routinely over-read as tested findings; the number of neighborhoods is chosen, not discovered, and identity shifts with window size. Test them downstream (neighborhood composition vs outcome), do not report them as significant in themselves.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Top SVGs are all known cell-type markers; SVG list ~= HVG list | Sample-wide SVG re-derives markers of spatially clustered cell types (cell-type-driven, not regulated) | Ask within-type: regress out cell-type composition or use ctSVG (C-SIDE/CTSV/CELINA); report the SVG-not-HVG subset |
| Moran's I finds "clustering" but cannot locate where the gene is HIGH | Moran/Geary detect clustering of like values, blind to high vs low | Use Getis-Ord Gi* for hot/cold spots |
| LISA/Gi* map is mostly "significant"; thousands of hits | n tests, AND local statistics are spatially autocorrelated so naive BH-FDR is invalid | FDR with stricter cutoffs (0.001); conditional permutation; treat map as exploratory, not independent discoveries |
| Confident "cell type A interacts with B" that vanishes on a second look | Default `nhood_enrichment` global-permutation null only beats complete randomness; abundant co-compartment types pass trivially | Demand survival under a conditional/within-compartment or toroidal null; report abundances |
| Global Moran's I near zero but tissue is clearly structured | Non-stationarity: opposite-sign local regions cancel in one global number | Use local statistics (LISA/Gi*); stratify by region |
| SVG ranking changes completely between two runs/tools | Different graph (k, Delaunay vs kNN) or different null -- methods test different hypotheses | Name the graph and null; report graph-robust hits; expect cross-method intersection to be small |
| Radius/length-scale statistic gives nonsense | Coordinates in pixels/array units, parameter in microns; Visium array coords are not microns | Convert to microns via scale factors before any distance parameter (see spatial-neighbors) |
| Rare cell type shows a huge enrichment z-score | Few edges -> high-variance estimate -> large |z| by chance | Report cell-type abundances; discount enrichment involving rare types |

## Related Skills

- spatial-neighbors - builds the graph W that every statistic here inherits; the choice propagates
- spatial-domains - region-level structure; a domain is not a colocalization result
- spatial-communication - ligand-receptor in space; the spillover/abundance confounds recur there
- single-cell/markers-annotation - cell-type labels feeding colocalization, and the marker overlap that confounds SVG

## References

- Svensson V, Teichmann SA, Stegle O (2018) SpatialDE: identification of spatially variable genes. Nature Methods 15(5):343-346. DOI 10.1038/nmeth.4636
- Sun S, Zhu J, Zhou X (2020) Statistical analysis of spatial expression patterns for spatially resolved transcriptomic studies (SPARK). Nature Methods 17(2):193-200. DOI 10.1038/s41592-019-0701-7
- Zhu J, Sun S, Zhou X (2021) SPARK-X: non-parametric modeling enables scalable and robust detection of spatial expression patterns for large spatial transcriptomic studies. Genome Biology 22:184. DOI 10.1186/s13059-021-02404-0
- Weber LM, Saha A, Datta A, Hansen KD, Hicks SC (2023) nnSVG for the scalable identification of spatially variable genes using nearest-neighbor Gaussian processes. Nature Communications 14:4059. DOI 10.1038/s41467-023-39748-z
- Palla G, Spitzer H, Klein M, et al. (2022) Squidpy: a scalable framework for spatial omics analysis. Nature Methods 19(2):171-178. DOI 10.1038/s41592-021-01358-2
- Schurch CM, Bhate SS, Barlow GL, et al. (2020) Coordinated cellular neighborhoods orchestrate antitumoral immunity at the colorectal cancer invasive front. Cell 182(5):1341-1359. DOI 10.1016/j.cell.2020.07.005
- Dos Santos Peixoto R, Miller BF, Brusko MA, et al. (2025) Characterizing cell-type spatial relationships across length scales in spatially resolved omics data (CRAWDAD). Nature Communications 16:350. DOI 10.1038/s41467-024-55700-1
- Moran PAP (1950) Notes on continuous stochastic phenomena. Biometrika 37(1/2):17-23. DOI 10.2307/2332142
- Geary RC (1954) The contiguity ratio and statistical mapping. The Incorporated Statistician 5(3):115-145. DOI 10.2307/2986645
- Getis A, Ord JK (1992) The analysis of spatial association by use of distance statistics. Geographical Analysis 24(3):189-206. DOI 10.1111/j.1538-4632.1992.tb00261.x
- Anselin L (1995) Local indicators of spatial association -- LISA. Geographical Analysis 27(2):93-115. DOI 10.1111/j.1538-4632.1995.tb00338.x
