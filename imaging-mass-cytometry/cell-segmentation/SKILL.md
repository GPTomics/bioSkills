---
name: bio-imaging-mass-cytometry-cell-segmentation
description: Segment single cells from multiplexed IMC/MIBI tissue images using Mesmer/DeepCell, Cellpose, or ilastik+CellProfiler, covering whole-cell vs nuclear segmentation, the summed-membrane-channel decision, nuclear-expansion bias, lateral spillover, resolution-floor parameters, and downstream-proxy evaluation. Use when delineating cells after preprocessing, choosing a segmentation model, building a cell mask for quantification, diagnosing impossible double-positive populations, or troubleshooting over/under-segmentation.
tool_type: python
primary_tool: deepcell
---

## Version Compatibility

Reference examples tested with: steinbock 0.16+, DeepCell 0.12+ (Mesmer), Cellpose 3.0+, numpy 1.26+, scikit-image 0.22+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: Mesmer expects `(batch, y, x, 2)` with channel 0 = nuclear, channel 1 = membrane, and `image_mpp` set to the TRUE acquisition resolution (~1.0 for IMC) -- it was trained at `model_mpp ~= 0.5` and rescales the input, so a wrong mpp degrades everything. Cellpose `cyto` trains at 30-px and `nuclei` at 17-px diameter; steinbock feeds nuclear-first (reversed vs native Cellpose). Recent steinbock cellpose containers default to the `cpsam` (Cellpose-SAM) model -- pin the version.

# Cell Segmentation for IMC

**"Segment cells from my IMC images"** -> Draw a per-cell boundary mask so that averaging the channels inside each mask yields single-cell expression.
- Python: `deepcell.applications.Mesmer().predict(...)`, `cellpose.models`
- CLI: `steinbock segment deepcell`, `steinbock segment cellpose`

## The Single Most Important Modern Insight -- segmentation is the largest irreversible error source, not a preprocessing step

A single-cell table is literally `for each mask_id: mean(pixels_in_mask, every_channel)`, so the mask defines the support of every measurement and no downstream step -- clustering, batch correction, differential abundance -- can recover a cell the mask merged or split. Two failure modes, and their asymmetry dictates how to tune. Under-segmentation (two cells in one mask) produces LOUD, catchable artifacts: a mask spanning a T cell and a macrophage reports CD3+CD68+, so biologically-impossible co-expression is a segmentation diagnosis until proven otherwise, not a discovery. Over-segmentation (one cell fragmented) is the QUIET, dangerous error: each fragment still looks like a plausible cell, but counts inflate and spatial-neighborhood statistics corrupt without obvious tells. Tuning a watershed or threshold until masks "look clean" usually trades the loud error for the quiet one, which is worse for spatial work. A second, independent problem rides on top: lateral (spatial) spillover -- real signal from a neighbor's membrane bleeding across the shared boundary at ~1 um resolution -- produces the same impossible-co-expression signature even with flawless masks and perfect channel compensation, so the two are confounded and must be addressed separately (REDSEA after segmentation; channel compensation before aggregation).

## Methods Landscape

| Tool / model | Class | Input it consumes | Strength | Fails when |
|--------------|-------|-------------------|----------|------------|
| Mesmer / DeepCell (Greenwald 2022) | deep, trained on TissueNet (incl. IMC/MIBI) | 2-ch: nuclear + summed membrane | purpose-built for multiplexed tissue; the IMC default | summed membrane channel is weak/patchy; wrong `image_mpp` |
| Cellpose / `cpsam` (Stringer 2021; Pachitariu 2025) | deep, flow-field / SAM backbone | 1-2 ch (cyto +- nuclear) | generalist; `cpsam` needs no diameter | default models carry a non-IMC size prior; wrong `diameter` |
| StarDist (Schmidt 2018) | star-convex polygon regression | single nuclear ch | excellent for crowded round nuclei | nuclear-only; breaks on irregular/elongated cells |
| ilastik + CellProfiler (Berg 2019; McQuin 2018) | random-forest pixels -> watershed | painted nucleus/cyto/background | transparent, tunable, no GPU; original IMC pipeline | semantic not instance; seed-threshold-sensitive; manual tuning |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Whole-cell phenotyping with a good broadly-expressed membrane marker set | Mesmer whole-cell, `image_mpp` = true resolution | TissueNet includes this modality; first choice for IMC/MIBI |
| Membrane staining weak/patchy/cell-type-specific | Nuclei (StarDist/Mesmer-nuclear) + small constrained expansion | a poor membrane sum systematically under-segments types lacking a marker |
| Only nuclear/intracellular markers needed (TFs, Ki-67) | Nuclear segmentation, quantify directly | nuclear markers barely suffer lateral spillover -- sidesteps the boundary problem |
| Mesmer struggles on the tissue | Cellpose / `cpsam`, optionally retrain (Cellpose 2.0) | a panel-specific learned prior beats a wrong generalist prior |
| Legacy / no-GPU / need full transparency | ilastik -> CellProfiler watershed | the original Bodenmiller pipeline; fully tunable |
| Impossible co-expression appears after any path | re-tune and/or REDSEA before clustering | the rate is the headline under-segmentation/spillover metric |

## Whole-Cell Segmentation with Mesmer

**Goal:** Produce whole-cell instance masks for surface-marker phenotyping.

**Approach:** Stack nuclear and summed-membrane channels as `(batch, y, x, 2)` and pass the true acquisition resolution as `image_mpp`. Mesmer internally rescales to its training resolution, so the mpp is load-bearing, not cosmetic.

```python
import numpy as np
from deepcell.applications import Mesmer

nuclear = img[dna_idx]                       # DNA/Ir channel
membrane = build_membrane(img, membrane_idx) # broadly-expressed membrane sum (see below)
stack = np.stack([nuclear, membrane], axis=-1)[np.newaxis, ...]  # (1, y, x, 2)

app = Mesmer()
masks = app.predict(stack, image_mpp=1.0, compartment='whole-cell')[0, ..., 0]  # ~1.0 for IMC
```

## Build the Summed Membrane Channel

**Goal:** Construct channel 2 so whole-cell masks are not biased against cell types lacking a marker.

**Approach:** Sum BROADLY-expressed membrane markers chosen to cover every cell type present (not only the types of interest), because a cell-type-specific sum is bright on some types and dark on others, systematically under-segmenting the dark ones.

```python
def build_membrane(img, membrane_idx):
    # sum pan-membrane markers covering ALL populations (e.g. pan-cytokeratin for
    # epithelium, CD45 for immune, E-cadherin, Na/K-ATPase) -- inspect the result
    # before trusting whole-cell masks; a patchy sum collapses into nuclear-like masks
    return img[membrane_idx].sum(axis=0)
```

## Orchestrate via steinbock

```bash
# Mesmer/DeepCell (nuclear-first); membrane channels are aggregated per the panel column
steinbock segment deepcell --minmax -o masks

# Cellpose container (current default model is cpsam; channel order is reversed vs native)
steinbock segment cellpose --minmax -o masks

# aggregate per-cell mean intensities (mean is the default and the right phenotyping choice)
steinbock measure intensities -o intensities
```

## Nuclear Segmentation with Constrained Expansion (fallback)

**Goal:** Approximate whole cells when membrane staining is absent, without the bias of free dilation.

**Approach:** Segment nuclei, then expand with a small radius under a competitive/watershed constraint so pixels are owned by exactly one cell. Fixed isotropic dilation is a cell-type-correlated bias (under-captures macrophages, over-captures small cells) and free dilation double-counts boundary pixels into two masks.

```python
from skimage.segmentation import expand_labels, watershed

# expand_labels grows each label into background but stops at the midline between
# labels (no overlap), so each pixel is assigned once -- a partition, unlike free dilation
expanded = expand_labels(nuclear_masks, distance=3)   # ~3 px at 1 um; report the radius
assert expanded.max() == nuclear_masks.max()          # no cells created/destroyed
```

## Per-Tool Failure Modes

### Mesmer -- wrong image_mpp
**Trigger:** leaving the default mpp on 1 um IMC. **Mechanism:** Mesmer rescales the image to its ~0.5 um training resolution; a wrong mpp rescales cells to the wrong learned size. **Symptom:** systematic over/under-segmentation across the whole image. **Fix:** pass `image_mpp` = true acquisition resolution (~1.0 IMC, ~0.5 or finer MIBI).

### Cellpose -- auto-diameter on tiny cells
**Trigger:** auto-diameter on ~5-px IMC nuclei. **Mechanism:** `cyto` rescales to a 30-px target; at single-digit diameters the rescale factor is large and unstable. **Symptom:** merged or fragmented masks. **Fix:** set `diameter` from known cell size in pixels, or use `cpsam` (no diameter dependence).

### Generalist model -- wrong size prior
**Trigger:** default Cellpose `cyto`/`cyto3` on IMC, trusted blindly. **Mechanism:** at ~5 px of evidence the learned PRIOR, not the image, draws the boundary, and a non-IMC prior is wrong. **Symptom:** plausible-looking but systematically biased masks. **Fix:** prefer Mesmer (TissueNet includes IMC/MIBI) or fine-tune Cellpose on the panel; evaluate on downstream proxies.

### Channel compensation in the wrong order
**Trigger:** REDSEA before segmentation, or channel compensation after aggregation. **Mechanism:** channel spillover is pixel-level (must be corrected before the per-cell average); lateral spillover is defined on segmented neighbors (must be corrected after). **Symptom:** residual impossible co-expression. **Fix:** pixel-compensate -> segment -> aggregate -> REDSEA.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `image_mpp` ~= 1.0 (IMC) | Greenwald 2022; Mesmer `model_mpp` ~0.5 | match acquisition resolution to the rescaler |
| Cellpose `cyto` 30 px / `nuclei` 17 px | Stringer 2021 | the trained-diameter targets the model rescales to |
| Nuclear expansion ~3 px @ 1 um | ImcSegmentationPipeline convention | approximates a thin cytoplasm without crossing into neighbors |
| Lymphocyte ~5-7 px @ 1 um/px | Giesen 2014 | the resolution floor that makes size priors load-bearing |
| Impossible-co-expression rate | Bai 2021 | per-slide under-segmentation/lateral-spillover monitor; not an F1 substitute |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| CD3+CD68+ "hybrid" cluster | under-segmentation or lateral spillover | treat as QC failure; re-tune + REDSEA before clustering |
| Whole-cell masks collapse to nuclei | weak/patchy summed membrane channel | broaden the membrane sum or fall back to nuclei + expansion |
| Same pixel counted in two cells | free dilation expansion | use `expand_labels`/watershed (exclusive ownership); assert label count unchanged |
| Macrophages under-captured | fixed isotropic nuclear dilation | constrained expansion; accept and report the bias; don't cross-compare with whole-cell data |
| Native Cellpose channel args do nothing in steinbock | steinbock reverses channel order | configure channels via the steinbock panel column, not native `--chan` semantics |
| High IoU but wrong biology | optimizing a pixel-overlap metric | accept on downstream proxies (impossible-co-expression rate, count/density sanity, positive-fraction stability); audit dense regions |

## References

- Giesen C, Wang HAO, Schapiro D, et al. 2014. Highly multiplexed imaging of tumor tissues with subcellular resolution by mass cytometry. *Nat Methods* 11(4):417-422. — IMC ~1 um resolution floor.
- Berg S, Kutra D, Kroeger T, et al. 2019. ilastik: interactive machine learning for (bio)image analysis. *Nat Methods* 16(12):1226-1232. — pixel classification stage.
- McQuin C, Goodman A, Chernyshev V, et al. 2018. CellProfiler 3.0: Next-generation image processing for biology. *PLoS Biol* 16(7):e2005970. — watershed instance segmentation.
- Schmidt U, Weigert M, Broaddus C, Myers G. 2018. Cell Detection with Star-Convex Polygons. *MICCAI 2018*, LNCS 11071:265-273. — StarDist nuclear baseline.
- Stringer C, Wang T, Michaelos M, Pachitariu M. 2021. Cellpose: a generalist algorithm for cellular segmentation. *Nat Methods* 18(1):100-106. — Cellpose diameters.
- Pachitariu M, Rariden M, Stringer C. 2025. Cellpose-SAM: superhuman generalization for cellular segmentation. *bioRxiv* doi:10.1101/2025.04.28.651001. — the cpsam SAM-backbone model (preprint).
- Greenwald NF, Miller G, Moen E, et al. 2022. Whole-cell segmentation of tissue images with human-level performance using large-scale data annotation and deep learning. *Nat Biotechnol* 40(4):555-565. — Mesmer/DeepCell, TissueNet, image_mpp.
- Bai Y, Zhu B, Rovira-Clave X, et al. 2021. Adjacent Cell Marker Lateral Spillover Compensation and Reinforcement for Multiplexed Images. *Front Immunol* 12:652631. — REDSEA boundary compensation; lateral spillover signature.
- Windhager J, Zanotelli VRT, Schulz D, et al. 2023. An end-to-end workflow for multiplexed image processing and analysis. *Nat Protoc* 18(11):3565-3613. — steinbock segmentation/measurement.

## Related Skills

- data-preprocessing - channel spillover compensation precedes segmentation
- phenotyping - consumes the single-cell mask and intensities; double-positives diagnose segmentation
- spatial-analysis - over-segmentation corrupts neighborhood statistics
- quality-metrics - segmentation QC metrics and the impossible-co-expression monitor
- interactive-annotation - overlay masks on channels to audit boundaries
