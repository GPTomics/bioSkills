---
name: bio-spatial-transcriptomics-spatial-data-io
description: Loads spatial transcriptomics data from Visium, Visium HD, Xenium, MERFISH/MERSCOPE, CosMx, Slide-seq/Curio, and Stereo-seq into AnnData or SpatialData using spatialdata-io and Squidpy. Use when deciding which platform class is in hand (imaging/in-situ vs sequencing/capture), which reader matches the platform (spatialdata_io.xenium/merscope/cosmx vs squidpy.read.visium/vizgen/nanostring), whether to work from the per-transcript molecule table (the re-segmentable source of truth) or the segmentation-derived per-cell matrix (quality-filtered, inherits all segmentation error), whether a molecule table even exists (spot platforms have none), and how to keep coordinate frames and units (pixel vs micron) registered to histology.
tool_type: python
primary_tool: spatialdata
---

## Version Compatibility

Reference examples tested with: spatialdata 0.2+, spatialdata-io 0.1.5+, squidpy 1.4+, scanpy 1.10+, anndata 0.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Spatial Data I/O

**"Load my spatial data"** -> Parse a platform's output bundle into one coordinate frame holding the expression matrix, coordinates, images, and (for imaging) the molecule table and segmentation shapes.
- Imaging/in-situ (Xenium, MERSCOPE/MERFISH, CosMx, seqFISH): `spatialdata_io.{xenium, merscope, cosmx}` -> SpatialData with a per-transcript `points` table AND a derived per-cell `tables` matrix.
- Sequencing/capture (Visium, Visium HD, Slide-seq/Curio, Stereo-seq): `squidpy.read.visium` or `spatialdata_io.{visium, visium_hd, curio, stereoseq}` -> spot/bin matrix + coordinates; NO molecule table.

## Governing Principle

The single most consequential I/O fact is that the two platform classes emit different primary objects, and one class emits two of them that are easy to confuse.

Imaging/in-situ platforms emit TWO physically distinct primary objects. The first is a per-TRANSCRIPT molecule table -- one row per decoded molecule with x, y (and often z), gene, a decoding-quality value, and a cell-assignment-or-unassigned. The second is a per-CELL expression matrix, genes-by-cells, DERIVED by overlaying a segmentation and counting the molecules that fall inside each boundary. The matrix looks exactly like scRNA-seq and is therefore wrongly trusted as ground truth, but it is a downstream product: it inherits every segmentation error and is usually quality-filtered (Xenium keeps Q>=20 in the matrix while the transcript table keeps everything). The molecule table is the source of truth and the ONLY object that lets the analyst re-segment, recover unassigned molecules, or do subcellular work. A loader that returns only the cell matrix has silently discarded the re-segmentable layer.

Sequencing/capture platforms (Visium, Slide-seq, Stereo-seq) have NO molecule table -- a spot/bead/bin is mini-bulk over the cells beneath it, captured as a single barcoded profile. Do not go looking for a transcript table that does not exist; the only objects are the barcode-by-gene matrix, the coordinates, and the tissue image.

A second trap is coordinate frames. Images, spot/cell coordinates, and molecule points each live in their own intrinsic pixel or array axes; overlaying transcripts on histology, or building a neighbor graph with a micron radius, requires the right transform (Visium scalefactors; imaging micron-to-pixel matrices). SpatialData makes the frames explicit (intrinsic vs a shared extrinsic "global" system); the legacy AnnData layout hides them in `uns['spatial'][library_id]['scalefactors']`. Mixing frames silently places points off the image or builds a graph at the wrong scale.

## The Platform-Class Fork

The first question of any spatial dataset is which side of the fork it sits on, because it decides what objects exist and what the rest of the pipeline must do.

| Class | Platforms | Primary objects | Molecule table? | Cell unit | Downstream |
|---|---|---|---|---|---|
| Imaging / in-situ | Xenium, MERSCOPE/MERFISH, CosMx, seqFISH | molecule table + segmentation-derived cell matrix + images + shapes | YES (source of truth) | from segmentation (a hypothesis) | segment, then label-transfer typing |
| Sequencing / capture | Visium, Visium HD, Slide-seq/Curio, Stereo-seq, GeoMx | barcode/bin-by-gene matrix + coordinates + image | NO | spot/bin = 1-10-cell MIXTURE (GeoMx ROI = many-cell region; sub-cell bins = fraction of a cell) | deconvolution (or bin/segment-up for sub-cell bins) |

## Object-Model Landscape

Each toolkit is a strategy for co-storing four things in one frame: an expression matrix, geometry (spot circles, cell/nucleus polygons, centroids), raster images (H&E, DAPI, multiplex IF -- often gigapixel), and (imaging only) the molecule point cloud. They differ in how separately they keep these and which is the forward path.

| Framework (language) | Core object | Molecule table | Per-cell matrix | Segmentation geometry | Images | Best when |
|---|---|---|---|---|---|---|
| SpatialData / scverse (Python) | `SpatialData` of elements | `points` (dask -> Parquet) | `tables` (AnnData) | `labels` (masks) + `shapes` (geopandas polygons) | `images` (xarray, OME-NGFF/Zarr, lazy) | Multimodal, larger-than-memory, multiple platforms in one store; re-segmentation |
| Squidpy + AnnData (Python) | `AnnData` | via SpatialData/readers | `adata.X` | in `obs`/external | `uns['spatial']` (legacy) | Standard AnnData spatial graph stats on spot or cell matrices |
| Seurat v5 (R) | `Seurat`, geometry in `@images` | `FOV@molecules` | `Assay5` counts | `FOV@boundaries` (Centroids + Segmentation) | platform `SpatialImage` | R single-cell users; v5 integration |
| SpatialExperiment / SFE-Voyager (Bioc, R) | `SpatialExperiment` / `SpatialFeatureExperiment` | `rowGeometries` (sf points, SFE only) | SCE assay | `colGeometries` (cellSeg/nucSeg/spotPoly) | `imgData()` | Bioconductor scran/scater; geospatial ESDA (Moran's I) |
| Giotto Suite (R) | `giotto` (multi-scale) | subcellular molecule layer | aggregated cell layer | polygon/cell layers | image layers | One technology-agnostic object spanning molecule -> cell -> region |

SpatialData is the forward-path Python standard because it is the only framework that natively keeps the molecule table, multiscale OME-NGFF images, and segmentation shapes as first-class, frame-aware elements. Because a SpatialData `table` IS an AnnData, all of `squidpy.gr.*` runs on it unchanged.

## Platform-to-Reader Map

| Platform | spatialdata-io reader | squidpy.read | Key I/O fact |
|---|---|---|---|
| Visium | `visium` | `visium` | spot, no molecule table; `tissue_positions.csv` gained a header at Space Ranger v2.0 (readers handle both) |
| Visium HD | `visium_hd` | -- | bins (2/8/16um); `tissue_positions.parquet` (PARQUET, not CSV); 8um bin still spans ~2 cells |
| Xenium | `xenium` | -- (none) | molecule table `transcripts.parquet` (all Q) + Q>=20 cell matrix; needs `experiment.xenium` manifest |
| MERSCOPE / MERFISH | `merscope` | `vizgen` | there is NO `merfish` reader -- `merscope` handles both; boundaries went hdf5-folder -> single `cell_boundaries.parquet` at instrument SW v232 |
| CosMx | `cosmx` | `nanostring` | flat CSVs with a run/slide prefix; `tx_file` (molecule table) absent for protein-only panels |
| Slide-seq / Curio | `curio` | -- (none) | bead, no molecule table |
| Stereo-seq | `stereoseq` | -- | DNB sub-cellular; binned up to cells |

`squidpy.read` provides ONLY `visium`, `vizgen`, and `nanostring` -- it has no `xenium` or `slideseq` reader. For Xenium, Slide-seq/Curio, Stereo-seq, and Visium HD, use the `spatialdata_io` reader. `scanpy.read_visium` is deprecated as of scanpy 1.11 -- prefer `squidpy.read.visium` (identical obsm/uns layout) or `spatialdata_io.visium`.

## Load Visium (Spot/Capture)

**Goal:** Read a Space Ranger bundle into an AnnData with coordinates, image, and scalefactors, without reaching for the deprecated scanpy reader.

**Approach:** Use `squidpy.read.visium`; coordinates land in `obsm['spatial']` (pixels of the full-res image), image and scalefactors in `uns['spatial'][library_id]`.

```python
import squidpy as sq

adata = sq.read.visium('spaceranger_out/')         # filtered matrix + spatial/; NOT a cell -- each spot is a 1-10-cell mixture
library_id = list(adata.uns['spatial'].keys())[0]
scalef = adata.uns['spatial'][library_id]['scalefactors']
# obsm['spatial'] is in FULL-RES pixels; multiply by tissue_hires_scalef to index the hires image
print(adata.n_obs, 'spots', adata.n_vars, 'genes', '| spot diameter (px):', scalef['spot_diameter_fullres'])
```

## Load Imaging Data and Keep the Molecule Table

**Goal:** Load Xenium (or MERSCOPE/CosMx) so BOTH the per-transcript molecule table and the derived cell matrix are available, not just the matrix.

**Approach:** Use the `spatialdata_io` reader, which returns a SpatialData object; the molecule table lives in `sdata.points`, the cell matrix in `sdata.tables`, segmentation polygons in `sdata.shapes`, images in `sdata.images`. Inspect element names with `print(sdata)` -- they vary by platform and reader version.

```python
import spatialdata_io as sdio

sdata = sdio.xenium('xenium_out/')                 # needs experiment.xenium manifest
print(sdata)                                       # lists points/tables/shapes/images element names

transcripts = sdata.points['transcripts']          # dask DataFrame: x, y, z, feature_name, qv, cell_id -- ALL Q-scores
adata = sdata.tables['table']                      # AnnData cell matrix -- Q>=20 filtered, inherits segmentation error
# the matrix is a DERIVED product; the molecule table is the re-segmentable source of truth
print('molecules:', transcripts.shape[0].compute(), '| cells in matrix:', adata.n_obs)
```

For MERSCOPE substitute `sdio.merscope('merscope_out/')`; for CosMx `sdio.cosmx('cosmx_out/')`. As an AnnData-only alternative for MERSCOPE, `sq.read.vizgen(path, counts_file='cell_by_gene.csv', meta_file='cell_metadata.csv')` returns the cell matrix but discards the molecule table.

## Load Other Capture Platforms

**Goal:** Load Visium HD bins, Slide-seq/Curio beads, or Stereo-seq into a SpatialData object.

**Approach:** Use the matching `spatialdata_io` reader; none of these has a molecule table, and Visium HD / Stereo-seq bins are smaller than a cell (the inverse-of-deconvolution regime -- see spatial-deconvolution).

```python
import spatialdata_io as sdio

sdata_hd = sdio.visium_hd('visium_hd_out/')        # tissue_positions are PARQUET; pick a bin (8um default still ~2 cells)
sdata_ss = sdio.stereoseq('stereoseq_out/')        # DNB sub-cellular spots, binned up to cells
sdata_bead = sdio.curio('slideseq_out/')           # Slide-seq/Curio beads; ~1 cell but ~1/3 carry >=2 types
```

## Inspect and Register Coordinate Frames

**Goal:** Confirm whether coordinates are in pixels or microns before building a graph or overlaying on histology, so a neighbor radius or a plotted point lands at the right scale.

**Approach:** In SpatialData read the element transformations (intrinsic vs the shared "global" extrinsic frame); in the AnnData layout read the scalefactors. Never assume `obsm['spatial']` units -- Visium is full-res pixels, most imaging readers place a micron "global" frame.

```python
from spatialdata.transformations import get_transformation

# SpatialData: every element carries transforms into shared coordinate systems
print(sdata.coordinate_systems)                    # e.g. ['global']
print(get_transformation(sdata['transcripts'], get_all=True))   # intrinsic -> global (often micron scaling)
```

## Convert SpatialData to AnnData

**Goal:** Extract the cell/spot matrix as a plain AnnData for tools that expect one, while keeping coordinates.

**Approach:** Copy the `table`, set `obsm['spatial']` from the matching shapes/centroids. Persist via `sdata.write(...)` ONLY to a scratch path -- a `.zarr` store is a DIRECTORY, not a file, so it must never be committed.

```python
adata = sdata.tables['table'].copy()
# write a zarr STORE (a directory) to scratch, never the repo; rm -rf when done
# sdata.write('/tmp/scratch/store.zarr')
```

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| `AttributeError: module 'squidpy.read' has no attribute 'xenium'` (or `slideseq`) | `squidpy.read` only has `visium`, `vizgen`, `nanostring` | Use `spatialdata_io.xenium` / `spatialdata_io.curio` for those platforms |
| `spatialdata_io` has no `merfish` reader | The reader is named for the instrument, not the chemistry | Use `spatialdata_io.merscope` (handles MERFISH and MERSCOPE) |
| `DeprecationWarning` / future removal on `scanpy.read_visium` | Deprecated as of scanpy 1.11 | Use `squidpy.read.visium` or `spatialdata_io.visium` (same layout) |
| Cell matrix has far fewer transcripts than the molecule table | Imaging cell matrix is Q>=20 filtered and segmentation-derived | Treat the matrix as provisional; use `sdata.points` (all Q) to re-segment or audit |
| Trusting the cell matrix as ground truth; weird co-expression | The matrix inherits all segmentation/spillover error | Validate against the molecule table; re-segment (see image-analysis) |
| Looking for a transcript table in Visium/Slide-seq and finding none | Capture platforms have no molecule table | Stop -- a spot is mini-bulk; there is nothing to re-segment |
| Visium HD reader fails reading `tissue_positions.csv` | Visium HD positions are PARQUET (`tissue_positions.parquet`) | Use `spatialdata_io.visium_hd`, which expects the parquet bundle |
| Older Visium positions parse with a shifted header | `tissue_positions.csv` gained a header at Space Ranger v2.0 | Current readers handle both; upgrade `spatialdata-io`/`squidpy` if parsing legacy files |
| MERSCOPE boundaries not found | hdf5-folder boundaries became single `cell_boundaries.parquet` at SW v232 | Match reader version to instrument software; point at the parquet if present |
| Seurat `@coordinates` slot missing (R interop) | Seurat 5.1 `VisiumV2` has no `coordinates` slot (V1 did) | Use `GetTissueCoordinates()`; there is no V1->V2 converter |
| Transcripts plot off the image | Points and image in different frames/units (pixel vs micron) | Apply the reader's transform / scalefactor before overlaying |
| A stray `.zarr` directory left in the repo after writing | `sdata.write` makes a directory store; `git status` hides untracked dirs | Write to scratch; `rm -rf` the store; run `git status --porcelain | grep '/$'` |

## Related Skills

- spatial-preprocessing - QC floors and normalization that differ by platform class after loading
- image-analysis - re-segment the molecule table; the cell matrix is a segmentation hypothesis
- spatial-deconvolution - recover cell-type proportions from spot mixtures that have no molecule table
- high-resolution-binning - bin/segment-up sub-cellular Visium HD and Stereo-seq captures
- spatial-visualization - plot spots vs imaging FOVs with the correct coordinate frame
- single-cell/data-io - non-spatial scRNA-seq loading for the deconvolution/label-transfer reference

## References

- Marconato L, Palla G, Yamauchi KA, et al. (2025) SpatialData: an open and universal data framework for spatial omics. Nature Methods 22(1):58-62. DOI 10.1038/s41592-024-02212-x
- Palla G, Spitzer H, Klein M, et al. (2022) Squidpy: a scalable framework for spatial omics analysis. Nature Methods 19(2):171-178. DOI 10.1038/s41592-021-01358-2
- Wolf FA, Angerer P, Theis FJ (2018) SCANPY: large-scale single-cell gene expression data analysis. Genome Biology 19:15. DOI 10.1186/s13059-017-1382-0
- Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA (2024) anndata: Access and store annotated data matrices. Journal of Open Source Software 9(101):4371. DOI 10.21105/joss.04371
- Hao Y, Stuart T, Kowalski MH, et al. (2024) Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology 42(2):293-304. DOI 10.1038/s41587-023-01767-y
- Righelli D, Weber LM, Crowell HL, et al. (2022) SpatialExperiment: infrastructure for spatially-resolved transcriptomics data in R using Bioconductor. Bioinformatics 38(11):3128-3131. DOI 10.1093/bioinformatics/btac299
- Moore J, Allan C, Besson S, et al. (2021) OME-NGFF: a next-generation file format for expanding bioimaging data-access strategies. Nature Methods 18:1496-1498. DOI 10.1038/s41592-021-01326-w
- Janesick A, Shelansky R, Gottscho AD, et al. (2023) High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis (Xenium). Nature Communications 14:8353. DOI 10.1038/s41467-023-43458-x
