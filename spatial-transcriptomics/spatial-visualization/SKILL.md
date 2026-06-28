---
name: bio-spatial-transcriptomics-spatial-visualization
description: Plots spatial transcriptomics expression, clusters, and annotations on tissue using Squidpy and Scanpy. Use when choosing the plotter and spot size by platform fork (sc.pl.spatial / sq.pl.spatial_scatter with real scalefactors and capture diameter for spot/capture data like Visium and Slide-seq, versus molecule/segmentation overlays for imaging/FOV data like Xenium, MERFISH, and CosMx); getting the histology coordinate-frame transform right (micron<->pixel, scalefactors) so points land on the image; and avoiding the honest-visualization traps where interpolation/KDE manufactures spatial pattern not in the data, oversized markers fake tissue coverage, jet and other non-uniform colormaps distort structure, and non-metric UMAP/tSNE distances are misread as spatial conclusions.
tool_type: python
primary_tool: squidpy
---

## Version Compatibility

Reference examples tested with: squidpy 1.4+, scanpy 1.10+, anndata 0.10+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Spatial Visualization

**"Plot expression / clusters / a score on my tissue section"** -> Render per-feature values at their real spatial coordinates, optionally over the histology image, without inventing structure the assay did not measure.
- Spot/capture fork (Visium, Visium HD, Slide-seq): `scanpy.pl.spatial` or `squidpy.pl.spatial_scatter` WITH the dataset's `scalefactors` and a `spot_size`/`size` set to the real capture diameter.
- Imaging/FOV fork (Xenium, MERFISH/MERSCOPE, CosMx): `squidpy.pl.spatial_scatter` with `shape=None` for the cell/molecule point cloud, or polygon shapes for segmentation overlays, or a platform viewer (Xenium Explorer, napari, Vitessce, TissUUmaps).

## Governing Principle

Plotting differs by the platform-class fork, and rendering choices can manufacture pattern that is not in the data.

The first decision is which side of the fork the data sits on, because it selects the plotter and the meaning of marker size. Spot/capture data carries a histology image and a `scalefactors` block that maps array coordinates to image pixels; the plotted spot stands for a real capture spot (a 55 um Visium spot is a 1-10-cell mixture, not a cell) and the marker size should reflect that capture diameter. Imaging/FOV data is a point cloud of segmented cells or individual transcript molecules with no Visium-style hex lattice; forcing it through a spot plotter or oversizing markers paints continuous tissue coverage over what is actually sparse, discrete detections. Using the wrong plotter or an arbitrary spot size silently misrepresents the tissue.

The deeper trap is that several common rendering choices invent structure. Smoothing, kernel-density, kriging, or contouring a sparse spatial field produces a continuous surface that looks like high-resolution biology but is interpolated -- the apparent gradients and domains can be artifacts of the kernel, and any spatial statistic (Moran's I, domain calls) computed on the smoothed field is partly circular. Oversized markers are rhetorical: in `scanpy.pl.spatial` `size` is a scaling factor on the spot diameter, so inflating it merges neighbors and fakes contiguity the assay never resolved. A perceptually non-uniform colormap (jet/rainbow) invents banding and edges in a smooth gradient and is unreadable under color-vision deficiency, and silent `vmin`/`vmax` clipping can erase or exaggerate differences. Finally, UMAP/tSNE distances are NOT metric (the same caveat as single-cell/clustering) -- gaps and cluster spacing in an embedding carry no spatial meaning and must not be read as tissue conclusions. Honest spatial visualization shows the raw points, names the transform and any clipping, and never lets a plotting parameter assert biology the measurement did not contain.

## Plot genre by platform fork

| Plot genre | Spot/capture fork (Visium, Slide-seq) | Imaging/FOV fork (Xenium, MERFISH, CosMx) | Honesty pitfall to avoid |
|------------|----------------------------------------|--------------------------------------------|--------------------------|
| Expression / cluster on tissue | `sc.pl.spatial` (uses `scalefactors`) or `sq.pl.spatial_scatter` | `sq.pl.spatial_scatter(shape=None)` point cloud | Oversized `spot_size`/`size` faking coverage |
| Histology overlay | `sc.pl.spatial(img_key='hires')`; scalefactor maps coords->pixels | `sq.pl.spatial_scatter(img=True, img_res_key=...)` with the registered image | Wrong coordinate frame (micron vs pixel) -> points off image |
| Single-molecule / transcript map | not applicable (no molecule table) | scatter the transcript x,y table, or Xenium Explorer / napari | Treating segmented matrix as raw molecules |
| Segmentation / boundary overlay | not applicable | `sq.pl.spatial_scatter` polygon shapes, or napari/TissUUmaps | Hiding segmentation error behind tidy cell polygons |
| Continuous field / heatmap | per-spot color, NO interpolation | per-cell color, NO interpolation | KDE/kriging/contour manufacturing gradients |
| Embedding (UMAP/tSNE) | `sc.pl.umap` for QC only | `sc.pl.umap` for QC only | Reading non-metric embedding distance as spatial |

When competing rendering options exist (point cloud vs polygon overlay, sequential vs diverging colormap), verify the current platform viewer and Squidpy plotting docs before committing -- spatial tooling and platform exports change quickly.

## Spot/Capture Plot with Real Scalefactors and Spot Size

**Goal:** Show expression or cluster labels at true spot positions on a spot/capture section with a marker size that reflects the capture geometry, not a guess.

**Approach:** Let `sc.pl.spatial` read the `uns['spatial']` `scalefactors` so spot coordinates align to the histology image; size markers from the recorded spot diameter rather than an arbitrary constant.

```python
import scanpy as sc
import squidpy as sq

# scalefactors live in adata.uns['spatial'][library_id]; sc.pl.spatial reads them automatically.
sc.pl.spatial(adata, color=['leiden', 'total_counts'], img_key='hires', alpha_img=0.6, ncols=2)

# A spot is a 1-10-cell MIXTURE, not a cell -- do not relabel spot clusters as cell types.
# squidpy resolves the scalefactor from library_id; size here is relative to the spot diameter.
sq.pl.spatial_scatter(adata, color='leiden', library_id='V1_Human_Lymph_Node', size=1.0)
```

## Imaging/FOV Overlay (Point Cloud and Segmentation)

**Goal:** Render imaging-platform cells or molecules in their real micron coordinates without imposing a spot lattice they do not have.

**Approach:** Use `sq.pl.spatial_scatter` with `shape=None` for the segmented-cell point cloud (or polygon shapes when boundaries are stored), and overlay the registered image only when its transform is known.

```python
# Imaging data is a point cloud, not a hex grid: shape=None plots cells as points in micron space.
# With no image, `size` is the ACTUAL dot size, not a scaling factor -- keep it small so sparse
# detections do not visually merge into fake continuous tissue.
sq.pl.spatial_scatter(adata, color='cell_type', shape=None, size=8, img=False)

# Overlay the registered morphology image only when the coordinate frame is trusted.
sq.pl.spatial_scatter(adata, color='EPCAM', shape=None, size=8, img=True, img_alpha=0.5)
```

## Histology Coordinate-Frame Overlay

**Goal:** Place transcripts/spots on the H&E or DAPI image so each point lands on the histological structure it came from.

**Approach:** Map array/micron coordinates into image-pixel space with the correct scalefactor (or platform affine); never plot raw micron coordinates onto a pixel image. Inspect the alignment before trusting any structure read off the overlay.

```python
# Spot/capture: hires-image pixel coords = spatial coords * tissue_hires_scalef.
library_id = list(adata.uns['spatial'].keys())[0]
scalef = adata.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef']
img = adata.uns['spatial'][library_id]['images']['hires']

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(img)                                   # image is in pixel space
coords_px = adata.obsm['spatial'] * scalef       # transform microns/array units -> pixels
ax.scatter(coords_px[:, 0], coords_px[:, 1], s=6, c='red')
ax.set_axis_off()                                # a small misalignment puts expression in the wrong structure
```

## Honest Continuous Field (Colormap and No Interpolation)

**Goal:** Show a continuous score across the section truthfully -- visible raw points, a perceptually uniform colormap, and disclosed clipping.

**Approach:** Color each measured spot/cell directly (never interpolate between them), pick a perceptually uniform map, and state any `vmin`/`vmax` clip rather than letting it silently reshape the gradient.

```python
# Color the MEASURED points only. Do NOT KDE/kriging/contour a sparse field -- that manufactures
# gradients and any Moran's I / domain call computed on the smoothed surface is partly circular.
sc.pl.spatial(adata, color='CD3D', cmap='viridis', vmin=0, vmax='p99')   # 'p99' clip is disclosed, not hidden

# Avoid jet/rainbow: perceptually non-uniform maps invent banding and fail color-vision-deficiency
# readers. Scientific colour maps (Crameri) are perceptually uniform; install cmcrameri to use them.
# import cmcrameri.cm as cmc; sc.pl.spatial(adata, color='CD3D', cmap=cmc.batlow)
```

## Interactive Exploration

Large imaging sections and multi-resolution images are better explored interactively than in static panels. napari (image + points + shapes layers), Vitessce (web, multimodal), TissUUmaps (large image-plus-marker viewing), and the vendor Xenium Explorer / Xenium Panel viewer all pan-and-zoom over the full-resolution data. The same coordinate-frame discipline applies: points must be transformed into the viewer's pixel space (for spot/capture, multiply spatial coordinates by the relevant `tissue_*_scalef`).

```python
import napari
library_id = list(adata.uns['spatial'].keys())[0]
img = adata.uns['spatial'][library_id]['images']['hires']
scalef = adata.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef']
viewer = napari.Viewer()
viewer.add_image(img, name='tissue')
viewer.add_points(adata.obsm['spatial'] * scalef, size=10, name='spots')   # transform into pixel space
napari.run()
```

## Embedding Caveat

`sc.pl.umap`/`sc.pl.tsne` are legitimate for QC and cluster sanity-checks, but UMAP/tSNE distances are not metric: the size of gaps between clusters and the apparent spacing of points carry no quantitative meaning, and nothing spatial can be concluded from them. Read tissue structure off the spatial plot, never off the embedding (see single-cell/clustering for the full caveat).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Spots overlap into a solid sheet; sparse signal looks continuous | `spot_size`/`size` set far above the real capture diameter | Size markers from the spot diameter; for imaging keep `size` small (it is the actual dot size when no image) |
| Points land off the image or in the wrong tissue region | Plotting micron/array coordinates onto a pixel image without the scalefactor/affine | Transform coords -> pixels (`* tissue_hires_scalef`, or the platform affine) before overlay |
| Smooth gradients/domains that vanish on the raw points | Field was KDE/kriged/contoured/imputed; pattern is the kernel, not the tissue | Plot measured points only; show raw alongside any smoothed view and disclose the kernel |
| Banding/edges appear in a smooth field; figure unreadable in grayscale | jet/rainbow or other perceptually non-uniform colormap | Use a perceptually uniform map (viridis, or Crameri scientific colour maps via cmcrameri) |
| Two conditions look very different for the same expression | Inconsistent or silent `vmin`/`vmax` between panels | Fix and disclose the color scale across panels (shared `vmin`/`vmax`) |
| Imaging cells plotted on a hex/grid lattice or with empty image background | Spot plotter (`sc.pl.spatial`) or default `shape` used on imaging point-cloud data | Use `sq.pl.spatial_scatter(shape=None)`; pass the registered image only with a known transform |
| Conclusions drawn from gaps between UMAP clusters | Treating non-metric embedding distance as spatial/quantitative | Restrict spatial claims to the spatial plot; use UMAP for QC only |
| Spot clusters labeled as cell types | A capture spot is a 1-10-cell mixture, not a cell | Label spot clusters as regions/niches; deconvolve for composition (spatial-deconvolution) |
| Per-spot proportion/scatterpie map read as measured composition | Deconvolution output is a model estimate carrying reference and fit uncertainty | Present proportion maps as estimates; rare-type fractions are least reliable, so corroborate before reading them off the map (spatial-deconvolution) |

## Related Skills

- spatial-data-io - load the platform data and the histology image plus scalefactors that plotting depends on
- spatial-domains - produce the region labels rendered on the section
- spatial-statistics - compute Moran's I / neighborhood enrichment whose results are plotted here
- data-visualization/heatmaps-clustering - general perceptually-uniform colormap and figure conventions
- single-cell/clustering - the non-metric UMAP/tSNE distance caveat that applies to embeddings

## References

- Palla G, Spitzer H, Klein M, et al. (2022) Squidpy: a scalable framework for spatial omics analysis. Nature Methods 19(2):171-178. DOI 10.1038/s41592-021-01358-2
- Wolf FA, Angerer P, Theis FJ (2018) SCANPY: large-scale single-cell gene expression data analysis. Genome Biology 19:15. DOI 10.1186/s13059-017-1382-0
- Marconato L, Palla G, Yamauchi KA, et al. (2025) SpatialData: an open and universal data framework for spatial omics. Nature Methods 22(1):58-62. DOI 10.1038/s41592-024-02212-x
- Crameri F, Shephard GE, Heron PJ (2020) The misuse of colour in science communication. Nature Communications 11:5444. DOI 10.1038/s41467-020-19160-7
- Chari T, Pachter L (2023) The specious art of single-cell genomics. PLoS Computational Biology 19(8):e1011288. DOI 10.1371/journal.pcbi.1011288
