---
name: bio-imaging-mass-cytometry-data-preprocessing
description: Load and preprocess imaging mass cytometry (IMC) and MIBI data from raw MCD/TXT through hot-pixel removal, spillover compensation, and variance-stabilizing transformation, covering readimc/steinbock ingestion, NNLS spillover compensation (CATALYST), IMC-Denoise, and the IMC arcsinh-cofactor question. Use when starting analysis from raw MCD files, building per-channel TIFF stacks, compensating channel spillover, choosing an arcsinh cofactor, or preparing single-cell intensities for phenotyping.
tool_type: mixed
primary_tool: steinbock
---

## Version Compatibility

Reference examples tested with: steinbock 0.16+, readimc 0.7+, numpy 1.26+, CATALYST 1.28+ (R/Bioconductor), cytomapper 1.16+ (R)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: IMC pixels are integer ion COUNTS, not fluorescence intensities. `CATALYST::compCytof` defaults to `method='nnls'` (non-negativity preserved); pass `cofactor=1` explicitly for IMC single-cell means (the compCytof default is `NULL`/5, the suspension value). The spillover-matrix and SCE channel names must both be `(metal)(mass)Di` (e.g. `Sm152Di`) or compensation silently no-ops. steinbock's hot-pixel filter is `steinbock preprocess imc images --hpf 50` (a signed 8-neighbor difference, not a median filter).

# IMC Data Preprocessing

**"Preprocess my imaging mass cytometry data"** -> Ingest raw acquisitions, suppress acquisition noise, compensate channel spillover, and variance-stabilize counts so each cell's measured intensity reflects real antigen abundance.
- Python/CLI: `readimc.MCDFile`, `steinbock preprocess imc images --hpf 50`
- R: `CATALYST::compCytof`, `cytomapper::compImage` for spillover compensation

## The Single Most Important Modern Insight -- IMC pixels are ion counts, and non-negativity is physics, not a preference

Every IMC/MIBI pixel is an integer number of detected metal ions from one ~1 um laser shot, drawn from a low-count, zero-inflated, near-Poisson regime where most pixels read 0-2 counts and the limit of detection is ~6 counts. Four consequences govern every preprocessing decision, and importing fluorescence-microscopy habits violates all four. (1) There is no continuous Gaussian background to subtract -- the floor is a count floor, and a true-negative pixel still reads 1-2 counts by Poisson chance. (2) Negative values are physically meaningless, so flow-style spillover compensation (exact matrix inverse) is wrong because it manufactures negatives; non-negative least squares (NNLS) is mandatory (Chevrier 2018 *Cell Syst* 6:612). (3) Per-pixel "expression" is mostly shot noise -- signal emerges only after segmentation sums a cell's pixels, so pixel maps are for localization and QC, never quantification. (4) Spillover is SPATIAL: a bright cell bleeding into a neighboring mass channel contaminates adjacent pixels, fabricating co-localization and false marker positivity at cell borders -- which means uncompensated spillover manufactures the exact cell-cell interactions IMC exists to measure. The corollary that trips up suspension-CyTOF veterans: the arcsinh cofactor is NOT 5 (cofactor 1 is the modern IMC default, Hunter 2024 *Cytometry A* 105:36), and there are no in-stream calibration beads in ablated tissue.

## Preprocessing Pipeline Order (load-bearing)

```
read .mcd (readimc) -> panel filter+sort (keep column) -> hot-pixel removal (DIMR or --hpf 50)
  -> [optional] DeepSNiF on low-SNR channels only -> spillover compensation (NNLS)
  -> segmentation (on compensated nuclear/membrane channels) -> per-cell aggregation
  -> arcsinh(cofactor 1) -> z-score / cohort-anchored normalization
```

Order is not cosmetic: denoise operates on RAW counts (the Poisson noise model is defined there), compensate BEFORE segmentation when spatial fidelity matters (corrupted membrane channels yield wrong boundaries that no later step recovers), and transform/normalize LAST.

## Denoising Taxonomy

| Method | Targets | Risk | When to use | Fails when |
|--------|---------|------|-------------|------------|
| steinbock `--hpf 50` | hot pixels | low | default fast hot-pixel pass | absolute-count threshold over-clips bright channels, under-cleans dim ones |
| IMC-Denoise DIMR (Lu 2023) | hot pixels | low | self-calibrating hot-pixel removal | misclassifies large multi-pixel hot-pixel clusters as signal |
| IMC-Denoise DeepSNiF (Lu 2023) | shot noise | HIGH | only channels with mean positive intensity < ~7 | over-smooths sparse/punctate markers and sub-1-2 um structure; biases extreme-low-count regions |
| 3x3 median filter | (do not use) | severe | never | returns 0 for isolated real single-positive pixels -- erases sparse biology |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Single-stain QC shows >~2% off-target spillover at relevant masses | Compensate (NNLS) before phenotyping | leak corrupts type calls and (spatially) neighborhood stats |
| Spatial neighborhood / interaction analysis is the endpoint | Pixel-level compensation (`compImage`) before segmentation | spillover is spatial and fakes interactions; cell-mean compensation comes too late |
| Means-only phenotyping, segmentation channels uncontaminated | Cell-level `compCytof` on SCE means | cheaper, less per-pixel Poisson noise |
| Channel driven above ~5,000 dual counts | Re-titrate at acquisition; do not compensate | linearity (and thus the matrix) breaks above saturation |
| Panel pre-designed to avoid bright/dim mass adjacencies, QC near-clean | Compensation ~ identity; skipping is defensible | avoids NNLS noise on near-zero pixels |
| Channel mean positive intensity < ~7, cannot phenotype on it | DIMR + DeepSNiF | shot noise dominates; denoising is the only way to use it |
| Channel clean, or punctate, or the one segmentation runs on | DIMR / `--hpf` only; skip DeepSNiF | DeepSNiF over-smooths real isolated signal and blurs boundaries |

## Read Raw Acquisitions

**Goal:** Ingest the multi-ROI MCD (the canonical source) and keep the metal-to-target panel mapping intact.

**Approach:** Open the MCD with readimc, iterate slides and acquisitions, and carry both `channel_names` (metals) and `channel_labels` (targets) -- conflating them is the most common ingestion bug. Prefer MCD over TXT (one ROI per file, and absent on Hyperion XTi).

```python
from readimc import MCDFile

with MCDFile('slide.mcd') as f:
    slide = f.slides[0]
    for acq in slide.acquisitions:
        img = f.read_acquisition(acq)              # (channels, y, x) float32 ion counts
        metals = acq.channel_names                 # e.g. 'Sm152' -- the mass channel
        targets = acq.channel_labels               # e.g. 'CD3'  -- the antibody target
        print(acq.id, img.shape, dict(zip(metals, targets)))
```

## Extract and Hot-Pixel Filter with steinbock

**Goal:** Build per-channel TIFF stacks, filtered to the analysis panel and de-spiked.

**Approach:** The panel CSV is both a filter and a sort key -- only `keep==1` rows are written and their row order defines channel order in the stack, so it must be pinned as a versioned artifact. The `--hpf` filter compares each pixel to its 8 neighbors with a signed difference and replaces spikes with the neighbor maximum (a conservative, valley-preserving operation), not a median.

```bash
# generate the panel template (edit the keep column before extracting)
steinbock preprocess imc panel

# extract TIFFs (keep-filtered, panel-ordered) with hot-pixel removal; 50 is a count
# difference, not a universal constant -- raise it for high-dynamic-range markers
steinbock preprocess imc images --hpf 50
```

## Spillover Compensation (NNLS)

**Goal:** Remove channel crosstalk (oxide M+16, abundance-sensitivity M+-1, isotopic impurity) without introducing negative counts.

**Approach:** Estimate the spillover matrix from single-stain controls per positive EVENT then take the median (population-summary estimation overcompensates because IMC's zero background biases ratios), then apply NNLS. Compensate pixels (`compImage`) before segmentation for spatial work, or cell means (`compCytof`) after segmentation otherwise. Channel names must be `(metal)(mass)Di`.

```r
library(CATALYST)
library(imcRtools)

# estimate the matrix from spotted single-stain TXTs (filenames carry the metal)
sce <- readSCEfromTXT('spillover/')
sce <- prepData(sce, transform = TRUE, cofactor = 5)   # cofactor 5 for PIXEL spot data
sce <- assignPrelim(sce); sce <- estCutoffs(sce); sce <- applyCutoffs(sce)
sm  <- computeSpillmat(sce)                            # the spillover matrix

# cell-level compensation on segmented single-cell means (NNLS is the default)
sce_cells <- compCytof(sce_cells, sm, method = 'nnls', cofactor = 1, overwrite = FALSE)
```

```r
# OR pixel-level compensation on the image stack, before segmentation (spatial fidelity)
library(cytomapper)
images <- compImage(images, adaptSpillmat(sm, channelNames(images)))
```

## Transform and Normalize

**Goal:** Variance-stabilize counts and remove batch offset without destroying cross-sample comparability.

**Approach:** Apply arcsinh with cofactor 1 on single-cell means (the OPTIMAL-derived IMC default, not the suspension-CyTOF 5), then z-score per channel against cohort-wide statistics. Per-image percentile or min-max scaling is a one-way door that makes equal biology look unequal across samples -- reserve it for visualization and always retain raw/compensated counts.

```python
import numpy as np

def arcsinh_cofactor1(cell_means):
    # cofactor 1 for IMC single-cell means (Hunter 2024); state the cofactor explicitly --
    # no field-wide standard exists, so reproducibility requires reporting it
    return np.arcsinh(cell_means / 1.0)

def zscore_per_channel(expr, mean, std):
    # mean/std computed COHORT-WIDE (not per-image) so scales stay comparable across samples
    return (expr - mean) / std
```

## Per-Method Failure Modes

### Flow-style compensation -- negative counts
**Trigger:** exact matrix inversion (`method='flow'` or generic linear unmixing). **Mechanism:** the inverse violates non-negativity and produces negative ion counts. **Symptom:** negative compensated values, downstream stats corrupted. **Fix:** `compCytof(..., method='nnls')` (the default); negatives are a solver artifact, not evidence compensation is wrong.

### Channel-name mismatch -- silent no-op
**Trigger:** SCE channel names are `Sm152` but the spillover matrix uses `Sm152Di` (or vice versa). **Mechanism:** name-based mapping finds no match. **Symptom:** compensation runs without error but changes nothing. **Fix:** enforce `(metal)(mass)Di`; reconcile with `adaptSpillmat()`.

### DeepSNiF on every channel -- invented structure
**Trigger:** blanket denoising "to clean things up." **Mechanism:** the Hessian continuity prior imposes spatial smoothness on genuinely sparse/punctate markers. **Symptom:** rare-population or punctate signal blurred into neighbors; biased low-count regions. **Fix:** denoise only channels with mean positive intensity < ~7; validate against the un-denoised image; never denoise the segmentation channel without checking boundary integrity.

### Per-image normalization before cross-sample comparison
**Trigger:** independent per-image 99th-percentile or min-max scaling, then comparing samples. **Mechanism:** image A's 99th percentile (40 counts) and image B's (400 counts) map to the same [0,1]. **Symptom:** a dim positive in A reads like a bright positive in B; differential abundance is spurious. **Fix:** derive normalization from cohort-wide statistics or a shared anchor; keep raw counts to re-derive.

### Median filtering as a denoiser
**Trigger:** `ndimage.median_filter` on the count image. **Mechanism:** a 3x3 median over sparse single-positive pixels returns 0. **Symptom:** real isolated membrane/punctate signal erased. **Fix:** use the neighbor-spike filter (`--hpf`) or DIMR; never blanket-median IMC.

### MIBI data run through the IMC pipeline unchanged
**Trigger:** applying this steinbock/CATALYST flow to MIBI-TOF data as if it were IMC. **Mechanism:** MIBI is SIMS on Au/Ta conductive slides, so it carries a 197Au slide background and crosstalk classes the CATALYST single-stain-bead model does not target. **Symptom:** gold-background contamination and uncorrected MIBI crosstalk. **Fix:** remove the Au/native-background channels and run MAUI (Baranski 2021) before this pipeline; treat the CATALYST spillover step as IMC-specific.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| arcsinh cofactor 1 (single-cell means) | Hunter 2024 *Cytometry A* 105:36 | maximizes positive/negative separation (Fisher ratio) for IMC counts; 5 over-compresses |
| arcsinh cofactor 5 (pixel/spot data) | CATALYST IMC workflow | pixel spot counts are higher-scale than cell means |
| Linearity ceiling ~5,000 dual counts | Chevrier 2018 *Cell Syst* 6:612 | above it count->abundance bends and the spillover matrix is invalid |
| `--hpf 50` (count difference) | steinbock convention | a per-experiment heuristic, not a universal constant -- tune to dynamic range |
| DeepSNiF only if mean positive intensity < ~7 | Lu 2023 *Nat Commun* 14:1601 | above ~7 the channel is effectively noise-immune; denoising is pure risk |
| Detection limit ~6 ion counts | Lu 2023 *Nat Commun* 14:1601 | below this, signal and shot noise are indistinguishable per pixel |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `compCytof` runs but values unchanged | channel-name format mismatch | enforce `Sm152Di`; `adaptSpillmat()` |
| Negative compensated counts | flow-style inversion | use `method='nnls'` |
| "Channel 12 is a different antibody for collaborators" | panel keep/order changed after segmentation | pin `panel.csv` as a versioned artifact; never reorder post-segmentation |
| Rare population vanished after denoising | DeepSNiF on a sparse channel | restrict DeepSNiF to low-SNR non-punctate channels |
| TXT loads only one ROI | analyzing TXT instead of MCD | ingest the multi-ROI `.mcd` with `readimc` |
| Cross-sample differences disappear or explode | per-image normalization | cohort-anchored normalization; keep raw counts |

## References

- Giesen C, Wang HAO, Schapiro D, et al. 2014. Highly multiplexed imaging of tumor tissues with subcellular resolution by mass cytometry. *Nat Methods* 11(4):417-422. — IMC ~1 um ion-count origin.
- Chevrier S, Crowell HL, Zanotelli VRT, Engler S, Robinson MD, Bodenmiller B. 2018. Compensation of Signal Spillover in Suspension and Imaging Mass Cytometry. *Cell Syst* 6(5):612-620.e5. — three spillover sources, single-stain beads, NNLS, ~5,000 dual-count linearity, CATALYST.
- Lu P, Oetjen KA, Bender DE, et al. 2023. IMC-Denoise: a content aware denoising pipeline to enhance Imaging Mass Cytometry. *Nat Commun* 14:1601. — DIMR hot-pixel and DeepSNiF shot-noise removal; mean>7 noise-immune guide.
- Windhager J, Zanotelli VRT, Schulz D, et al. 2023. An end-to-end workflow for multiplexed image processing and analysis. *Nat Protoc* 18(11):3565-3613. — the steinbock workflow and readimc/imcRtools ingestion.
- Hunter B, Nicorescu I, Foster E, et al. 2024. OPTIMAL: An OPTimized Imaging Mass cytometry AnaLysis framework for benchmarking segmentation and data exploration. *Cytometry A* 105(1):36-53. — arcsinh cofactor 1 and z-score-after-arcsinh for IMC.
- Baranski A, Milo I, Greenbaum S, et al. 2021. MAUI (MBI Analysis User Interface): An image processing pipeline for Multiplexed Mass Based Imaging. *PLoS Comput Biol* 17(4):e1008887. — MIBI-specific crosstalk, aggregate, and gold-background removal.

## Related Skills

- quality-metrics - reading the spillover matrix and gating channels before compensation
- cell-segmentation - segmentation runs on compensated nuclear/membrane channels
- phenotyping - consumes the arcsinh-transformed single-cell matrix
- flow-cytometry/compensation-transformation - suspension spillover and arcsinh background
- single-cell/preprocessing - AnnData conventions for the single-cell matrix
