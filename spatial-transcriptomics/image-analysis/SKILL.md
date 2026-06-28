---
name: bio-spatial-transcriptomics-image-analysis
description: Segments cells/nuclei and extracts image features from imaging spatial transcriptomics (Xenium, MERFISH/MERSCOPE, CosMx) and H&E/IF tissue images using Cellpose, StarDist, Baysor, and Squidpy. Use when choosing a segmentation strategy (DAPI nucleus + expansion vs membrane-stain whole-cell vs transcript-aware Baysor/proseg vs segmentation-free SSAM) given the available stain; judging whether transcript spillover is fabricating false co-expression and short-range cell-cell signal; and deciding whether the derived cell-by-gene matrix is trustworthy before downstream typing, DE, or ligand-receptor analysis.
tool_type: python
primary_tool: squidpy
---

## Version Compatibility

Reference examples tested with: squidpy 1.7+, scanpy 1.10+, scikit-image 0.22+, numpy 1.26+, pandas 2.2+, cellpose 4.0+ (CLI tools: Baysor 0.6+, proseg 1.0+)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Image Analysis for Spatial Transcriptomics

**"Segment cells from my imaging data"** -> Draw cell/nucleus boundaries on an image (or on the transcript point cloud) and assign each molecule to one cell, producing a cell-by-gene matrix.
- Python image-based: `cellpose.models.CellposeModel().eval()`, StarDist for nuclei, `squidpy.im.segment()` (watershed baseline)
- CLI transcript-based: Baysor, proseg run on the molecule table (x, y, gene)

**"Extract image features for my spots"** -> Summarize pixel intensity/texture under each Visium spot (a DIFFERENT operation from segmentation -- no cell boundaries are drawn).
- Python: `squidpy.im.calculate_image_features()` (summary, histogram, texture/GLCM)

Keep these two operations distinct. `squidpy.im` feature extraction (texture/summary on H&E) describes the image patch under a spot for domain detection. Segmentation manufactures the cells themselves. This skill leads with segmentation because that is the dominant error source; feature extraction is a downstream convenience.

## Governing Principle

In imaging spatial data the cell is a segmentation HYPOTHESIS, not an observation. The raw data are images plus a table of decoded molecules with x, y (and sometimes z); there is no native cell. A cell-by-gene matrix exists only AFTER an algorithm draws boundaries and assigns each molecule to one cell or to background. Every row of that matrix is produced by the segmentation step, which is the dominant, irreducible upstream error source for every imaging platform -- it confounds typing, DE, and communication downstream (Mitchel et al. 2026 *Nat Genet* 58:434, who find segmentation errors "dominate the results" for context-dependent DE and ligand-receptor inference).

Segmentation fails three ways, and each fabricates a specific downstream lie:
- Over-segmentation splits one cell into many -> inflated cell count; fragments look low-quality and get filtered or mis-typed.
- Under-segmentation merges neighbors into a spatial DOUBLET whose profile is a mixture of two types. Unlike droplet doublets these are spatially structured -- adjacent types merge preferentially -- so standard doublet detectors (Scrublet, DoubletFinder) MISS them, because the synthetic doublets those tools simulate are random pairs, not neighbors.
- Transcript mis-assignment / spillover: a molecule of cell A is assigned to neighbor B (diffusion, optical PSF bleed, z-collapse, boundary error). Spillover is DISTANCE-DEPENDENT and strongest between adjacent heterotypic cells, so it manufactures spatially-structured contamination, not uniform noise.

The cascade is worst for exactly the analyses people prize. A T cell abutting epithelium picks up keratin spillover -> false marker co-expression -> a spurious "transitional/hybrid" state that is pure artifact; rare cells in dense parenchyma (TILs, neutrophils) are swamped by neighbor spillover and lost. Crucially, distance-dependent spillover FABRICATES the short-range co-localization that ligand-receptor and cell-cell communication tools detect -- so an L-R "hit" between two touching types can be a pure segmentation artifact (a circularity; see spatial-communication). Treat the derived matrix as PROVISIONAL and run a contamination QC step before trusting any single-cell-resolution claim.

## The decision: choose by available signal, not by reflex

The first question is not "which tool" but "what boundary signal do I have?" A nuclear stain says where the nucleus is; a membrane/boundary stain says where the cell ENDS. Whole-cell segmentation is boundary-finding, and DAPI carries no boundary information -- so DAPI-only whole-cell is always inference (expansion). Adding a membrane/boundary stain converts that inference into measurement, and is the single highest-leverage change available -- it beats swapping algorithms on DAPI-only data.

| Tool | Class | Input signal | Best when | Fails / weak when |
|------|-------|--------------|-----------|-------------------|
| StarDist | image, star-convex polygons | DAPI, 1 channel | round, crowded NUCLEI; fast | non-convex shapes (whole cells, neurons) cannot be represented; not a whole-cell tool |
| Cellpose | image, DL generalist | 1-2 ch (nucleus +/- membrane) | generalist; 2-channel nucleus+membrane = true whole-cell; retrainable | no transcript-only mode; over/under-segments on DAPI-only without a membrane channel |
| Watershed (`squidpy.im.segment`) | classic flooding | DAPI seeds + intensity | fast baseline; seeded splitting of touching nuclei | over-segments textured nuclei; seed/threshold-sensitive; no shape prior |
| Mesmer / DeepCell | image, DL whole-cell (TissueNet) | 2 ch: nuclear + membrane | any platform WITH a membrane stain (also CODEX/MIBI/IMC); human-level whole-cell | needs a membrane channel; DAPI-only -> nuclear only |
| Baysor | transcript MRF+EM, optional prior | molecule table (+ optional DAPI) | refining/replacing image segmentation by transcriptional composition; recovers cells images miss; runs with or without a prior | sparse/low-plex panel + no prior -> unstable; compute-heavy |
| proseg | transcript, cell-simulation membership | molecule table | transcript-only whole-cell WITHOUT a membrane stain; fast; recovers hard immune cells | sparse-panel limits of all transcript-only methods; newer/less battle-tested |
| SSAM / ClusterMap | segmentation-FREE molecule density | molecule table | cell-type/domain MAPPING when boundaries are hopeless; recovers low-density types | produce NO cell objects -> no per-cell composition, counts, or neighbor graph |

The ladder of trust runs: DAPI-only StarDist/Cellpose-nuclei (accept nuclear sensitivity loss or expansion bias) < membrane-stain whole-cell Mesmer or 2-channel Cellpose < transcript-aware Baysor/proseg (molecules, not a fixed radius, set boundaries) < segmentation-free SSAM/ClusterMap (best mapping, but the cell unit is lost -- no per-cell matrix, neighbor graph, or QC). Methods evolve fast here; verify current best practice against the latest benchmarks before committing.

### Nucleus-only and the expansion trap

Nuclear (DAPI) segmentation is robust because nuclei are round, separated, and high-contrast -- but the nucleus holds only a minority of mRNA, so cytoplasmic transcripts fall OUTSIDE the mask and are discarded (large sensitivity loss) or must be reassigned. The cheap substitute is nucleus expansion: dilate each nuclear mask by a fixed radius until it hits a neighbor. This assumes round, equal-sized, isotropically-arranged cells -- false for almost all tissue. In dense tissue expanded disks collide and partition intercellular space by a Voronoi-like rule unrelated to true membranes (the worst region for spillover); elongated or large-cytoplasm cells (neurons, muscle, glia, macrophages) are badly served -- a fixed disk captures none of their projections and steals neighbors' transcripts. No single radius is correct for a heterogeneous tissue. Expansion is a baseline, not a solution.

Xenium makes this concrete: XOA v1.0-1.9 used DAPI + 15 um nucleus expansion; v2.0+ cut the default to 5 um "for improved accuracy" -- an admission that 15 um over-assigned in dense tissue. The vendor changed the answer, so do not treat any expansion radius as ground truth. Note that Cellpose on Xenium is a community path via Xenium Ranger `import-segmentation`, NOT the built-in XOA default -- do not conflate them.

## Segment nuclei from a DAPI/IF image

**Goal:** Produce instance masks (one integer label per cell) from a nuclear-stain image as the starting cell hypotheses.

**Approach:** Run Cellpose's generalist model; in v4 (Cellpose-SAM) there is one model, channels are no longer an input, and `diameter` is optional because the model is size-invariant. With a membrane channel available, pass it as a second channel for true whole-cell masks instead of nuclei + expansion.

```python
from cellpose import models

model = models.CellposeModel(gpu=False)            # v4 Cellpose-SAM single generalist model; v3 used models.Cellpose(model_type='nuclei')
masks, flows, styles = model.eval(dapi_image, diameter=None)   # v4 returns 3 values + drops channels=; v3 returned masks, flows, styles, diams and took channels=[0,0]
# masks: integer label image; 0 = background, 1..N = cells. This is a HYPOTHESIS, not ground truth.
n_cells = int(masks.max())
```

DAPI-only masks are nuclei. Approximating whole cells without a membrane stain requires expansion (round-cell bias above) or a transcript-aware method. With a membrane/boundary channel, stacking `[nucleus, membrane]` and passing both lets Cellpose-SAM use the first channels in any order -- that converts boundary inference into measurement.

## Re-segment from the transcript table (membrane-free whole-cell)

**Goal:** Recover cells that image-based nuclear segmentation drops (small, irregular, immune) by letting transcript composition and density define boundaries.

**Approach:** Run Baysor or proseg on the per-molecule table (x, y, gene). These are CLI tools; the molecule table is the source of truth and the only object that permits re-segmentation. Optionally seed Baysor with the vendor nuclear masks as a prior.

```bash
# Baysor: molecule-table segmentation; -m = min transcripts/cell, -s = expected cell scale (um), :gene names the gene column
baysor run -x x_location -y y_location -g feature_name -m 30 -s 10 \
  --prior-segmentation-confidence 0.5 transcripts.csv nucleus_id

# proseg: transcript-only whole-cell, reads Xenium/CosMx/MERSCOPE molecule tables directly
proseg --xenium transcripts.csv.gz --output-counts counts.csv.gz --output-cell-polygons cells.geojson
```

Baysor and proseg output per-molecule cell assignments -> rebuild a cell-by-gene matrix from those. Where the panel is sparse and no membrane stain exists, all transcript-only methods become unstable -- check cell-yield and size distributions against the image before trusting them.

## QC the segmentation before trusting the matrix

**Goal:** Detect the segmentation failure modes (over/under-segmentation, spillover) BEFORE they propagate into typing and communication results.

**Approach:** Treat the matrix as provisional. Inspect cell-size and transcripts-per-cell distributions (bimodality flags merged doublets or fragments), check for impossible co-expression of mutually exclusive lineage markers (a spillover signature), and where possible run a dedicated contamination tool.

```python
import numpy as np

counts = np.asarray(adata.X.sum(axis=1)).ravel()       # transcripts per cell
area = adata.obs['cell_area'].to_numpy()               # from the segmentation polygons
# Over-segmentation: a spike of tiny, low-count fragments. Under-segmentation: a tail of huge, high-count "cells".
print('transcripts/cell pct [5,50,95]:', np.percentile(counts, [5, 50, 95]))
print('cell area pct [5,50,95]:', np.percentile(area, [5, 50, 95]))

# Spillover signature: cells co-expressing markers of two mutually exclusive lineages (e.g. epithelial KRT + T-cell CD3).
# Distance-dependent -> these false double-positives concentrate at heterotypic boundaries.
epi = np.asarray(adata[:, 'EPCAM'].X).ravel() > 0
tcell = np.asarray(adata[:, 'CD3E'].X).ravel() > 0
print('suspicious EPCAM+CD3E+ cells:', int((epi & tcell).sum()))
```

Dedicated correction/QC tools target the distance-dependent contamination directly: SPLIT and neighborhood factorization (Mitchel et al. 2026) for contamination, FastReseg for transcript-based re-segmentation (CosMx), and ovrlpy for vertical/z-collapse doublets. SOPA runs Cellpose and Baysor on the same data with patch-based conflict resolution. Run a contamination step before any rare-state, hybrid-state, or short-range L-R claim.

## Extract image features per spot (NOT segmentation)

**Goal:** Summarize the tissue image under each Visium spot (intensity, texture) to augment expression-based spatial-domain detection.

**Approach:** Wrap the image in a Squidpy `ImageContainer` and call `calculate_image_features`. This draws no cell boundaries -- it describes the pixel patch under each spot. Pass `layer='image'` explicitly if a segmentation layer already exists on the container.

```python
import squidpy as sq

img = sq.datasets.visium_hne_image_crop()              # ImageContainer; pair with the matching adata
adata = sq.datasets.visium_hne_adata_crop()
sq.im.calculate_image_features(adata, img, layer='image', features=['summary', 'texture'],
                               key_added='img_features', n_jobs=1, show_progress_bar=False)
# texture = GLCM (contrast, homogeneity, correlation, ASM); summary = per-channel intensity stats
feats = adata.obsm['img_features']                     # rows = spots, columns = features
```

Watershed via `sq.im.segment(img, layer='image', method='watershed')` is a fast classical baseline that over-segments H&E; use it for a quick look, not for production cell calling. Morphology per mask comes from `skimage.measure.regionprops_table` (area, eccentricity, solidity).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Many tiny, low-count "cells" | Over-segmentation split single cells | Lower model sensitivity / raise `min_size`; check cell-area histogram for a fragment spike; prefer a learned model over watershed |
| Cluster of cells co-expressing exclusive lineage markers (KRT + CD3) | Transcript spillover fabricating co-expression at heterotypic boundaries | Run a contamination QC (SPLIT, neighborhood factorization); re-segment with a membrane stain or Baysor/proseg; do not interpret as a "hybrid state" |
| A short-range ligand-receptor hit between two touching types | Distance-dependent spillover manufactures the short-range co-localization (circularity) | Validate against segmentation quality; constrain L-R by distance; treat as hypothesis (see spatial-communication) |
| Standard doublet detector finds almost nothing, yet merged cells exist | Spatial doublets are neighbor merges, not random pairs the detector simulates | Inspect transcripts/cell and area tails; re-segment; do not rely on Scrublet/DoubletFinder for imaging merges |
| Cytoplasmic markers nearly absent from every cell | Nucleus-only mask discarded cytoplasmic mRNA | Expand the mask, add a membrane stain, or use a transcript-aware method |
| Sharp drop in transcripts/cell after a vendor software update | Xenium expansion default cut 15 um -> 5 um (v2.0) | Expected; the smaller radius assigns fewer (and fewer mis-assigned) transcripts -- re-run downstream, do not "fix" |
| `model.eval` returns 3 values but code unpacks 4 | Cellpose v4 dropped `diams` and the `channels=` argument | Unpack `masks, flows, styles`; remove `channels=`; `diameter` is optional in v4 |
| `Unable to determine which layer to use` | A segmentation layer was added, so the container has >1 layer | Pass `layer='image'` to `calculate_image_features` / `segment` |
| Transcript-only segmentation yields implausible cell shapes/yield | Sparse/low-plex panel with no prior -> Baysor/proseg unstable | Add a nuclear prior; compare yield + size to the image; fall back to image segmentation |

## Related Skills

- spatial-transcriptomics/spatial-preprocessing - QC floors and non-gene-count normalization for the post-segmentation matrix
- spatial-transcriptomics/spatial-communication - ligand-receptor inference, where segmentation spillover fabricates short-range signal (the circularity)
- spatial-transcriptomics/spatial-proteomics - whole-cell segmentation on membrane markers for CODEX/IMC/MIBI (Mesmer/DeepCell)
- imaging-mass-cytometry/cell-segmentation - segmentation for multiplexed-imaging proteomics
- spatial-transcriptomics/spatial-data-io - load the molecule table (the only object that permits re-segmentation) and the derived matrix

## References

- Stringer C, Wang T, Michaelos M, Pachitariu M (2021) Cellpose: a generalist algorithm for cellular segmentation. Nature Methods 18:100-106. DOI 10.1038/s41592-020-01018-x
- Pachitariu M, Stringer C (2022) Cellpose 2.0: how to train your own model. Nature Methods 19:1634-1641. DOI 10.1038/s41592-022-01663-4
- Schmidt U, Weigert M, Broaddus C, Myers G (2018) Cell Detection with Star-Convex Polygons (StarDist). MICCAI 2018, LNCS 11071:265-273. DOI 10.1007/978-3-030-00934-2_30
- Greenwald NF, Miller G, Moen E, et al. (2022) Whole-cell segmentation of tissue images with human-level performance using large-scale data annotation and deep learning (Mesmer/DeepCell). Nature Biotechnology 40:555-565. DOI 10.1038/s41587-021-01094-0
- Petukhov V, Xu RJ, Soldatov RA, et al. (2022) Cell segmentation in imaging-based spatial transcriptomics (Baysor). Nature Biotechnology 40:345-354. DOI 10.1038/s41587-021-01044-w
- Jones DC, Elz AE, Hadadianpour A, et al. (2025) Cell simulation as cell segmentation (proseg). Nature Methods 22:1331-1342. DOI 10.1038/s41592-025-02697-0
- Park J, Choi W, Tiesmeyer S, et al. (2021) Cell segmentation-free inference of cell types from in situ transcriptomics data (SSAM). Nature Communications 12:3545. DOI 10.1038/s41467-021-23807-4
- Mitchel J, Gao T, Petukhov V, et al. (2026) Impact and correction of segmentation errors in spatial transcriptomics. Nature Genetics 58:434-444. DOI 10.1038/s41588-025-02497-4
- Palla G, Spitzer H, Klein M, et al. (2022) Squidpy: a scalable framework for spatial omics analysis. Nature Methods 19:171-178. DOI 10.1038/s41592-021-01358-2
- Janesick A, Shelansky R, Gottscho AD, et al. (2023) High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis (Xenium). Nature Communications 14:8353. DOI 10.1038/s41467-023-43458-x
