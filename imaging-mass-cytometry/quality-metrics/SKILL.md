---
name: bio-imaging-mass-cytometry-quality-metrics
description: Quality control for IMC/MIBI data across pixel, channel, image, slide, and batch levels, covering Poisson-count SNR (cell-level Gaussian-mixture and empty-channel comparison), spillover-matrix QC (the three physical sources), drift and the missing EQ-bead analog, acquisition artifacts, and sample-of-origin batch effects. Use when deciding whether to keep or drop a channel, ROI, or slide, distinguishing a dim antibody from a failed one, reading a spillover matrix, or diagnosing batch-driven clustering before analysis.
tool_type: mixed
primary_tool: CATALYST
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, scikit-learn 1.4+, scanpy 1.10+, CATALYST 1.28+ (R), spillR 1.0+ (R)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: IMC values are integer ion (dual) counts, so SNR must respect Poisson statistics, not fluorescence intuition; biology lives in 1-2 count differences. There is no EQ-bead in-line drift normalizer for ablated tissue. CATALYST `normCytof` is for SUSPENSION bead normalization, not IMC images -- do not apply it to image data. Compensate raw pixels before transformation.

# IMC Quality Metrics

**"Assess the quality of my IMC acquisition"** -> Gate the data at the level each failure lives at -- pixel, channel, image, slide, batch -- before analysis, not by normalizing after.
- Python: `numpy`/`scikit-learn` for SNR, artifacts, batch diagnosis
- R: `CATALYST::plotSpillmat`, `spillR` for spillover QC

## The Single Most Important Modern Insight -- QC is multi-level, and every metric is blind at some level

IMC/MIBI QC is not one number. The data live in a Poisson ion-count regime where "noise" has a defined statistical meaning, and the failures that actually destroy an experiment -- a dead antibody, an unbalanced batch, cells that cluster by which slide they came from rather than by phenotype -- are panel/staining/batch problems that are invisible to per-image SNR. So a single metric is always blind at some level, and the discipline is to gate (drop a channel, ROI, or slide) before analysis rather than normalize after, because correction moves a distribution but never creates the positive/negative separation that staining never produced. Three corollaries a postdoc internalizes. (1) Counts are Poisson: a real floor-abundance epitope genuinely yields a few counts, so "2 counts" is signal or noise depending on dwell, area, and the aggregation level -- SNR grows as sqrt(N) when pixels are pooled into a cell. (2) Dim is not failed: a correctly-titrated antibody to a sparse antigen is supposed to be dim; failure is INSEPARABILITY of positive and negative populations, judged against a known-empty channel, not low absolute intensity. (3) IMC has no EQ-bead drift normalizer -- ablated fixed tissue cannot be spiked with calibration beads, so the only honest cross-batch yardstick is an anchor reference sample included in every run (Casanova 2025), and the absence of an in-line standard is itself expert knowledge.

## Multi-Level QC Framework

| Level | What to measure | Characteristic failure | Blind to |
|-------|-----------------|------------------------|----------|
| Pixel | hot pixels, shot noise, dynamic range | detector spikes; Poisson noise on dim signal | whether the channel is biologically real |
| Channel (marker) | cell-level SNR, spillover in/out, vs empty channel | dead antibody, crosstalk, oxide/+-1 leak | spatial artifacts, batch |
| Image / ROI | mean intensity, cell coverage, ablation completeness | failed ablation, folding, off-target ROI | cross-sample comparability |
| Slide / acquisition | detector drift over time, tune (Lu duals) | within-run sensitivity decay, mis-tune | between-slide offset |
| Batch / cohort | sample-of-origin clustering, lot effects | the cohort clusters by batch not biology | nothing -- the top level |

## Decision Tree by Scenario

| Observation | Decision | Basis |
|-------------|----------|-------|
| Cell-level positive/negative mixture won't separate; signal ~ empty/80ArAr channel | DROP (failed antibody) | inseparability, not intensity |
| Low absolute counts but clean separation, pattern matches biology + control tissue | KEEP (dim-but-real); use at aggregated levels | dim != failed |
| High signal, low SNR (everything "positive") | DROP or re-titrate | non-specific binding |
| Heavy +16 oxide or impurity from a co-expressed partner | DROP or re-mass the panel | unrescuable by compensation |
| Striping / incomplete-ablation banding | DROP the ROI | physical failure, not correctable noise |
| DNA/Ir dropout over a region | MASK the region, keep the rest | non-ablation/tissue loss |
| Tune fails (Lu duals below panel criterion) | RE-TUNE / re-acquire | instrument not in spec |
| Cells cluster by slide/patient not phenotype | batch-correct; if it won't mix, the contrast is confounded | sample-of-origin effect |

## Cell-Level SNR (the decision-relevant number)

**Goal:** Judge marker adequacy at the unit of analysis (the cell), in a count-aware way.

**Approach:** Fit a two-component Gaussian mixture to per-cell mean counts (on non-transformed counts) and take mean(positive)/mean(negative); a failed antibody is one whose components do not separate, regardless of brightness. Compare the distribution to a known-empty channel as the operational "did this antibody work" test.

```python
import numpy as np
from sklearn.mixture import GaussianMixture

def cell_snr(per_cell_counts):
    # two-component mixture on raw per-cell means: separation, not brightness, defines success
    gm = GaussianMixture(n_components=2, random_state=0).fit(per_cell_counts.reshape(-1, 1))
    pos, neg = np.sort(gm.means_.ravel())[::-1]
    return pos / neg if neg > 0 else np.inf

def matches_empty(channel_counts, empty_channel_counts, q=95, tol=2.0):
    # compare the POSITIVE tail (q-th percentile), not the median: a real-but-sparse marker
    # carries its signal in the tail while a failed channel's tail sits at the empty floor.
    # True -> indistinguishable from 80ArAr / an unconjugated lanthanide -> the honest "failed" test
    return np.percentile(channel_counts, q) - np.percentile(empty_channel_counts, q) <= tol
```

## Spillover Matrix QC

**Goal:** Decide whether a panel's crosstalk is acceptable before compensating.

**Approach:** Generate the matrix from single-stain controls and read whole rows, not just neighbors -- spillover has three physically distinct sources with different mass signatures and different fixes. Acceptability is co-expression-dependent: the same percentage is fine between unrelated markers and fatal between co-expressed ones.

```r
library(CATALYST)

sce <- readSCEfromTXT('spillover/')         # single-metal-spotted slides; filenames carry the metal
sce <- prepData(sce, transform = TRUE, cofactor = 5)
sce <- assignPrelim(sce); sce <- applyCutoffs(estCutoffs(sce))
sm  <- computeSpillmat(sce)
plotSpillmat(sce, sm)                        # inspect M+-1 (abundance), M+16 (oxide), and
                                             # any bright off-diagonal at a NON-adjacent mass (impurity)
```

## Batch / Sample-of-Origin QC

**Goal:** Catch the dominant real-world failure -- cells clustering by slide/patient rather than phenotype -- which no per-image metric reports.

**Approach:** Embed cells and color by patient, slide, day, and antibody lot; if cells separate by sample, there is a batch problem. Diagnose before correcting, and correct at the batch layer with an anchor reference sample.

```python
import scanpy as sc

sc.pp.pca(adata); sc.pp.neighbors(adata); sc.tl.umap(adata)
sc.pl.umap(adata, color=['patient', 'slide', 'acquisition_day', 'antibody_lot'])
# separation by these = batch, not biology; a dead/unbalanced channel cannot be normalized into life
```

## Per-Source Failure Modes

### "2 counts is noise"
**Trigger:** discarding low-count channels by fluorescence intuition. **Mechanism:** at 1 um^2/~1 ms dwell a floor-abundance epitope yields a few counts; biology lives in 1-2 count differences. **Symptom:** real dim markers dropped. **Fix:** judge adequacy at the aggregation level analyzed; SNR scales sqrt(N) with pooled pixels; mean expression > ~7 is effectively noise-immune.

### Pixel correlation read as spillover
**Trigger:** flagging high pixel-channel Pearson correlation as spillover. **Mechanism:** co-expressed real markers correlate too; spillover is a directional, mass-structured leak. **Symptom:** false spillover calls, missed real ones. **Fix:** read the single-stain spillover matrix; diagnose by mass signature (+-1, +16, named impurity mass), not correlation.

### Compensating a saturated or co-expressed channel
**Trigger:** trusting compensation on very bright donors or co-expressed pairs. **Mechanism:** the matrix is linear only in the linear range and cannot separate real co-expression from leak. **Symptom:** over/under-shoot; subtracted real biology. **Fix:** keep total per-pair spillover low by panel design; use NNLS (CATALYST) or flag-and-replace (spillR); the real fix is upstream mass assignment.

### Per-image QC declared sufficient
**Trigger:** passing per-image SNR and skipping cross-sample QC. **Mechanism:** FFPE/ischemia/lot variation shifts baselines by batch. **Symptom:** unsupervised analysis groups by sample-of-origin. **Fix:** UMAP/ridgeline by patient/slide/day/lot; include anchor samples; correct at the batch layer.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Mean expression > ~7 ~ noise-immune | Lu 2023 *Nat Commun* 14:1601 | below it shot noise perturbs per-cell values |
| Abundance sensitivity (M+-1) < 0.3% (Tb) | Han 2018 *Nat Protoc* 13:2121 | TOF peak-tail spec; tuning target, not guarantee |
| Oxide (M+16) < 3% (La) | Han 2018 *Nat Protoc* 13:2121 | plasma-oxide spec; worst for abundant structural markers |
| Isotopic impurity up to ~4% at a named mass | Han 2018 *Nat Protoc* 13:2121 | not predictable from mass proximity -- read the lot |
| Tune ~ Lu >= 1500 dual counts | panel-specific convention | a pass criterion, stated in dual counts (unit matters) |
| Pixel foreground (Otsu) signal < 2/image -> flag | steinbock/IMCDataAnalysis | image-level marker filter |

Where the field has NO accepted threshold (itself expert knowledge): a universal SNR cutoff for a "good marker"; a single spillover percentage defining an "acceptable panel" (acceptability is co-expression-dependent); an in-line pixel-level drift-normalization standard equivalent to EQ beads; a fixed hot-pixel count threshold (DIMR/KNN are adaptive precisely because a fixed cutoff fails across brightnesses).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Dropped a real sparse marker | judged by absolute intensity | drop on inseparability + match-to-empty + wrong spatial pattern |
| Spillover "fixed" but double-positives persist | compensated co-expressed/saturated channel | re-mass the panel; NNLS/spillR; compensate raw pixels |
| Cohort clusters by patient | unaddressed batch | anchor reference sample; diagnose before correcting |
| Threshold "1500" or "2" ambiguous | unit omitted | always state dual counts; thresholds are panel/instrument-specific |
| Striped ROI "denoised" | physical ablation failure treated as noise | drop the ROI |
| Indium nuclear signal taken as a marker (MIBI) | In localizes to nuclei | treat as artifact unless validated; use 197Au + background masking |

## References

- Giesen C, Wang HAO, Schapiro D, et al. 2014. Highly multiplexed imaging of tumor tissues with subcellular resolution by mass cytometry. *Nat Methods* 11(4):417-422. — IMC origin; ~50 copies/um^2 detection floor.
- Chevrier S, Crowell HL, Zanotelli VRT, et al. 2018. Compensation of Signal Spillover in Suspension and Imaging Mass Cytometry. *Cell Syst* 6(5):612-620.e5. — spillover sources, single-stain beads, less accurate at high ion load, CATALYST.
- Han G, Spitzer MH, Bendall SC, Fantl WJ, Nolan GP. 2018. Metal-isotope-tagged monoclonal antibodies for high-dimensional mass cytometry. *Nat Protoc* 13(10):2121-2148. — M+-1/M+16/impurity specs and panel design.
- Finck R, Simonds EF, Jager A, et al. 2013. Normalization of mass cytometry data with bead standards. *Cytometry A* 83A(5):483-494. — EQ four-element bead normalization (suspension).
- Ijsselsteijn ME, Somarakis A, Lelieveldt BPF, Hollt T, de Miranda NFCC. 2021. Semi-automated background removal limits data loss and normalizes imaging mass cytometry data. *Cytometry A* 99(12):1187-1197. — sample-of-origin clustering, FFPE/ischemia variation.
- Baranski A, Milo I, Greenbaum S, et al. 2021. MAUI: An image processing pipeline for Multiplexed Mass Based Imaging. *PLoS Comput Biol* 17(4):e1008887. — MIBI artifacts, gold/indium, 1-2 count biology.
- Lu P, Oetjen KA, Bender DE, et al. 2023. IMC-Denoise: a content aware denoising pipeline to enhance Imaging Mass Cytometry. *Nat Commun* 14:1601. — Poisson noise model, mean>7 noise-immune.
- Guazzini M, Reisach AG, Weichwald S, Seiler C. 2024. spillR: spillover compensation in mass cytometry data. *Bioinformatics* 40(6):btae337. — flag-and-replace compensation, preserves correlations.
- Windhager J, Zanotelli VRT, Schulz D, et al. 2023. An end-to-end workflow for multiplexed image processing and analysis. *Nat Protoc* 18(11):3565-3613. — image/cell-level QC and SNR conventions.
- Casanova C, et al. 2025. Standardization of Suspension and Imaging Mass Cytometry Single-Cell Readouts for Clinical Decision Making. *Cytometry A* 107(6):390-403. — anchor/reference samples for batch-level drift correction.

## Related Skills

- data-preprocessing - hot-pixel removal, denoising, and NNLS spillover compensation
- cell-segmentation - segmentation QC and the impossible-co-expression monitor
- phenotyping - failed channels and batch corrupt cell-type calls
- differential-analysis - batch as a covariate when comparing across conditions
- flow-cytometry/cytometry-qc - suspension bead normalization and channel QC background
