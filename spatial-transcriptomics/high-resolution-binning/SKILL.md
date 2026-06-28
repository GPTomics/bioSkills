---
name: bio-spatial-transcriptomics-high-resolution-binning
description: Reconstructs single cells from sub-cellular spatial capture units (Visium HD 2um bins, Stereo-seq DNB spots, Slide-seqV2 beads) by aggregating bins UP into cells rather than deconvolving a mixture DOWN. Use when choosing a bin size and recognizing the sparsity-vs-mixture dilemma (2um bins are too sparse to cluster, but binning to 8/16um re-creates the multi-cell mixture deconvolution was meant to escape); deciding between morphology-driven cell reconstruction (Bin2cell -- StarDist/Cellpose nuclei on a registered H&E/DAPI image, then assign 2um bins to nuclei) and fixed-bin aggregation by whether a co-registered cell image exists; recognizing this as the INVERSE of deconvolution (bin UP, not mix DOWN -- this is the AMBIGUOUS regime of the resolution fork); and handling each platform (Visium HD has an image so reconstruct, Slide-seqV2 has no per-bead image so aggregate or deconvolve, Stereo-seq depends on a registered stain).
tool_type: python
primary_tool: bin2cell
---

## Version Compatibility

Reference examples tested with: bin2cell 0.3+, scanpy 1.10+, anndata 0.10+, spatialdata 0.1+, squidpy 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# High-Resolution Binning

**"Turn my Visium HD 2um bins into cells"** -> Aggregate sub-cellular capture features UP into single-cell profiles, using a registered nucleus image to decide which bins belong to which cell when one exists.
- Python: Bin2cell (`b2c.read_visium` -> `b2c.stardist` -> `b2c.insert_labels` -> `b2c.bin_to_cell`) for image-guided reconstruction; `scanpy`/`squidpy` for fixed-bin aggregation when no image exists

## Governing Principle

Binning is the INVERSE of deconvolution. Deconvolution takes a capture unit that is LARGER than a cell (a 55um Visium spot holding 1-10 cells) and mixes it DOWN into the cell-type fractions inside it. High-resolution platforms have the opposite geometry: a Visium HD 2um bin, a Stereo-seq ~220nm DNB spot, and a Slide-seqV2 10um bead are SMALLER than or comparable to a single cell, so each unit is a fragment of one cell, not a mixture of several. The task is to aggregate fragments UP into whole cells, never to deconvolve a mixture that does not exist. Running deconvolution on 2um bins invents fractional cell-type mixtures inside features that hold only part of one cell.

The trap that defeats the naive fix is coarse binning. The 2um bins are far too sparse to cluster directly -- most bins capture a handful of transcripts or none, so a per-bin expression vector carries no cell-type signal. The reflex is to bin up to a coarser grid (Visium HD ships 8um and 16um bins for exactly this reason). But an 8um bin still spans roughly two cells, so coarse binning trades resolution for the precise multi-cell-mixture problem the high resolution was meant to escape -- it lands back in the DECONVOLVE regime, now needing a reference and a deconvolution method. This is a genuine dilemma, not a tunable knob: too fine is unclusterably sparse, too coarse is a mixture.

The escape is to define the cell from morphology instead of from a fixed grid. When a high-quality registered image exists (Visium HD ships an H&E or DAPI image co-registered to the bin coordinates), segment nuclei on the IMAGE, then assign each 2um bin to the nucleus whose territory contains it, and sum the bins per nucleus into a real single-cell profile. The cell boundary comes from morphology, not from an arbitrary square. Without a per-feature registered cell image (Slide-seqV2 beads have no co-registered cell morphology), morphology reconstruction is impossible and fixed-bin aggregation or bead-level deconvolution (RCTD doublet-mode is common for Slide-seqV2) remains the standard. Platform plus image availability decides the approach -- not the tool.

## The reconstruction decision

This skill IS the AMBIGUOUS regime of the resolution fork named in spatial-deconvolution: the near-single-cell middle where a unit holds part of, or roughly, one cell. The fork there sorts platforms into DECONVOLVE (spot >> cell), SEGMENT (imaging, already single cells), and AMBIGUOUS; everything below is the AMBIGUOUS branch.

| Platform | Native unit | Co-registered cell image? | Recommended approach | Pitfall |
|----------|-------------|---------------------------|----------------------|---------|
| Visium HD | 2um square bins (gapless lawn) | YES -- H&E or DAPI from CytAssist, registered to bins | Morphology-driven reconstruction (Bin2cell: StarDist/Cellpose nuclei -> assign 2um bins -> per-cell sum) | Treating 8um bins as the unit; an 8um bin still mixes ~2 cells |
| Stereo-seq | ~220nm DNB spots, binned (bin20 ~10-14um, bin50 ~25-36um) | Sometimes -- ssDNA/nuclei stain if acquired and registered | Reconstruct from the stain if registered (StereoCell/Cellpose); else fixed-bin aggregation | Default bin50 spans several cells -> a mixture, not a cell |
| Slide-seqV2 | 10um beads (random close-pack) | NO -- beads carry no co-registered cell morphology | Fixed-bin/bead aggregation, or bead deconvolution (RCTD doublet-mode) | Reconstructing cells from morphology -- there is no image to segment |

The discriminating axis is the registered cell image, not the platform name. A Visium HD run without a usable image collapses to the Slide-seqV2 row; a Stereo-seq run with a clean registered ssDNA stain behaves like the Visium HD row. Confirm the image is registered to the bin coordinate frame before trusting any morphology reconstruction; a misregistered image assigns bins to the wrong nuclei silently. Methods here evolve quickly -- verify the current best practice and the tool's registration assumptions against its latest documentation before committing.

## Loading Visium HD bins

**Goal:** Read the 2um bin matrix together with the registered morphology image into one object whose bin coordinates and image pixels share a frame.

**Approach:** Use the Bin2cell reader, which wraps the Space Ranger 2um output and attaches the full-resolution source image; Visium HD tissue positions are PARQUET, not CSV, and the reader handles that. Inspect the bin sparsity before deciding fine-reconstruct vs coarse-aggregate.

```python
import bin2cell as b2c
import numpy as np

# square_002um is the 2um bin output; source_image_path is the full-res H&E/DAPI registered to the bins
adata = b2c.read_visium('visium_hd_outs/binned_outputs/square_002um/',
                        source_image_path='Visium_HD_tissue_image.tif',
                        spaceranger_image_path='visium_hd_outs/spatial/')

median_counts = np.median(np.asarray(adata.X.sum(axis=1)).ravel())   # 2um bins are sparse: often single-digit median UMIs
print(f'bins: {adata.n_obs}, median UMI/bin: {median_counts:.1f}')   # too sparse to cluster -> reconstruct, do not cluster bins
```

## Morphology-driven cell reconstruction (Bin2cell)

**Goal:** Build true single-cell profiles by segmenting nuclei on the registered image and summing the 2um bins that fall inside each nucleus territory.

**Approach:** Scale the H&E to the segmentation resolution, destripe the Visium HD per-row/per-column count artifact, run StarDist for nuclei, insert the labels onto the bin coordinates, expand each nucleus to capture cytoplasmic bins, then collapse bins per label into a cell-level AnnData. Each cell records how many bins it absorbed.

```python
import bin2cell as b2c

mpp = 0.5                                                            # microns-per-pixel for the scaled image; sets StarDist's effective resolution
b2c.scaled_he_image(adata, mpp=mpp, save_path='stardist/he.tiff')
b2c.destripe(adata)                                                 # corrects Visium HD per-row/per-column total-count striping before it biases segmentation

b2c.stardist(image_path='stardist/he.tiff', labels_npz_path='stardist/he.npz',
             stardist_model='2D_versatile_he', prob_thresh=0.01)    # H&E nuclei; '2D_versatile_fluo' for DAPI
b2c.insert_labels(adata, labels_npz_path='stardist/he.npz', basis='spatial',
                  spatial_key='spatial_cropped_150_buffer', mpp=mpp, labels_key='labels_he')
b2c.expand_labels(adata, labels_key='labels_he', expanded_labels_key='labels_he_expanded')   # nucleus -> cell territory for cytoplasmic bins

cdata = b2c.bin_to_cell(adata, labels_key='labels_he_expanded',
                        spatial_keys=['spatial', 'spatial_cropped_150_buffer'])
# cdata is cell-level: bins summed per label; cdata.obs['bin_count'] = bins absorbed per cell -> a QC handle
```

Bins assigned to no nucleus (label 0) are dropped -- they are inter-cellular space or unsegmented territory, and forcing them into a cell fabricates expression. A cell built from very few bins is a low-confidence reconstruction; filter on `bin_count` the way single-cell QC filters on UMIs. When the H&E nuclei miss sparse regions, a second StarDist pass on a gene-expression-derived image (`b2c.grid_image` -> `2D_versatile_fluo`) plus `b2c.salvage_secondary_labels` rescues cells the H&E alone missed.

## Fixed-bin aggregation when no image exists

**Goal:** Produce a workable cell-scale matrix from Slide-seqV2 beads or an imageless Stereo-seq run, accepting that each unit is approximate rather than a morphology-defined cell.

**Approach:** Aggregate to a cell-scale grid (or treat beads as the unit) and pass the result downstream as APPROXIMATE cells; if the bins clearly mix types, hand them to bead-level deconvolution instead of pretending they are pure. Choose the grid in microns, not in bins, so the physical scale is explicit.

```python
import scanpy as sc
import numpy as np

# coords are in microns; choose a grid near one cell diameter (~10um) -- coarser re-creates the multi-cell mixture
bin_um = 10
coords = adata.obsm['spatial']
gx = np.floor(coords[:, 0] / bin_um).astype(int)
gy = np.floor(coords[:, 1] / bin_um).astype(int)
adata.obs['grid'] = [f'{x}_{y}' for x, y in zip(gx, gy)]            # aggregate bins/beads sharing a grid cell

agg = sc.get.aggregate(adata, by_key='grid', func='sum')           # sum counts per grid cell -> approximate cell-scale matrix
agg.X = agg.layers['sum']
# a grid cell spanning two real cells is a MIXTURE -> if so, deconvolve it (see spatial-deconvolution) rather than typing it
```

The honest caveat: a fixed grid is a compromise, and the coarser it is the more it is a deconvolution problem wearing a cell label. If the downstream question is cell typing and the beads visibly mix types, route to spatial-deconvolution (RCTD doublet-mode for Slide-seqV2) instead of clustering the grid.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Clustering on 2um bins yields noise / empty clusters | 2um bins are far too sparse (single-digit UMIs) to carry cell-type signal | Do not cluster bins; reconstruct cells (Bin2cell) or aggregate to a cell-scale grid first |
| "Cell types" from 8um/16um bins look like blends | An 8um bin still spans ~2 cells -- it is a mixture, not a cell | Reconstruct from morphology, or treat the bin as a mixture and deconvolve (spatial-deconvolution) |
| Deconvolution "runs" on 2um bins but fractions are nonsense | Deconvolved a sub-cellular fragment as if it were a multi-cell mixture (inverted the geometry) | Aggregate UP into cells; deconvolution applies to spot >> cell, not bin << cell |
| Bin2cell assigns bins to the wrong nuclei | Source image not registered to the bin coordinate frame | Verify image-to-bin registration before reconstruction; a misregistered image fails silently |
| Reconstruction wanted but there is no image to segment | Slide-seqV2 (and imageless Stereo-seq) have no per-bead cell morphology | Use fixed-bin aggregation or bead deconvolution; morphology reconstruction needs a registered image |
| Reconstructed cells have tiny `bin_count` and erratic profiles | Cells built from too few bins are low-confidence | Filter on `bin_count` as single-cell QC filters on UMIs; consider salvage_secondary_labels |
| Striping artifacts bias nuclei or counts | Visium HD per-row/per-column total-count striping left uncorrected | Run `b2c.destripe` before segmentation and before downstream normalization |

## Related Skills

- spatial-deconvolution - the resolution fork that sends the AMBIGUOUS regime here; deconvolution is the opposite (mix DOWN) geometry to this skill's bin UP
- image-analysis - nucleus/cell segmentation (StarDist, Cellpose) that morphology-driven reconstruction depends on
- spatial-data-io - load Visium HD PARQUET bin positions and the registered image before reconstruction
- spatial-preprocessing - QC and normalize the reconstructed cells once they exist (cell-scale, not bin-scale, thresholds)
- single-cell/cell-annotation - annotate the reconstructed cells with markers or label transfer
- single-cell/clustering - cluster reconstructed cells, which now carry cell-scale signal that raw bins lacked

## References

- Polanski K, Bartolome-Casado R, Sarropoulos I, et al. (2024) Bin2cell reconstructs cells from high resolution visium HD data. Bioinformatics 40(9):btae546. DOI 10.1093/bioinformatics/btae546
- Chen A, Liao S, Cheng M, et al. (2022) Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball-patterned arrays (Stereo-seq). Cell 185(10):1777-1792. DOI 10.1016/j.cell.2022.04.003
- Stickels RR, Murray E, Kumar P, et al. (2021) Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2. Nature Biotechnology 39:313-319. DOI 10.1038/s41587-020-0739-1
- Schmidt U, Weigert M, Broaddus C, Myers G (2018) Cell detection with star-convex polygons (StarDist). MICCAI, Lecture Notes in Computer Science 11071:265-273. DOI 10.1007/978-3-030-00934-2_30
- Stringer C, Wang T, Michaelos M, Pachitariu M (2021) Cellpose: a generalist algorithm for cellular segmentation. Nature Methods 18:100-106. DOI 10.1038/s41592-020-01018-x
- Cable DM, Murray E, Zou LS, et al. (2022) Robust decomposition of cell type mixtures in spatial transcriptomics (RCTD). Nature Biotechnology 40:517-526. DOI 10.1038/s41587-021-00830-w

Visium HD (2um bins, registered CytAssist H&E/DAPI image, PARQUET tissue positions) is a 10x Genomics product; 10x provides the Space Ranger output specification and onboard image registration but no primary peer-reviewed platform paper, so it is attributed to 10x Genomics rather than a citation.
