# XCMS Untargeted LC-MS Preprocessing Usage Guide

## Overview

This skill drives programmatic untargeted LC-MS feature extraction in R with the modern xcms 4.x API: raw centroided mzML in, a features-by-samples intensity table out, via CentWave peak detection, retention-time alignment, peak-density correspondence, gap-filling, CAMERA redundancy collapse, and built-in QC feature filtering. It guards against the traps that make untargeted tables untrustworthy: parameters that silently delete real peaks, gap-filling that fabricates intensities for absent compounds, alignment that looks perfect in QCs while mis-registering rare features, and the one-compound-many-features redundancy that inflates apparent dimensionality.

## Prerequisites

```r
if (!require('BiocManager')) install.packages('BiocManager')
BiocManager::install(c('xcms', 'MsExperiment', 'Spectra', 'CAMERA'))
```

Conceptual prerequisites: data must be centroided (centroid in software with a documented algorithm, e.g. ProteoWizard msconvert vendor peakPicking, rather than relying on irreversible on-instrument centroiding); the instrument class and chromatography must be known to set `ppm` and `peakwidth`; the sample design (biological groups, pooled QCs, blanks, injection order) must be recorded as a per-file table before processing; and the analyst should understand that the resulting feature table is a parameterized result that must be reported with its full processing specification.

## Quick Start

Tell your AI agent what you want to do:
- "Process my centroided mzML files with xcms into a feature table"
- "Detect peaks with CentWave using parameters appropriate for my Orbitrap UHPLC data"
- "Align retention times with obiwarp and group peaks across samples"
- "Gap-fill, then report the fraction of filled values per feature"
- "Collapse adduct and isotope redundancy with CAMERA before annotation"
- "Filter features by QC CV and D-ratio with xcms filterFeatures"

## Example Prompts

### Peak Detection
> "Set up CentWaveParam for qTOF UHPLC data and explain why ppm should not be the spec mass accuracy."
> "My trace metabolites are missing from the table - which of prefilter, noise, and snthresh should I lower first?"

### Alignment and Grouping
> "Align retention times to a pooled QC reference rather than the first file, then regroup."
> "Choose between obiwarp and peakGroups for a heterogeneous case/control cohort with few shared peaks."
> "Set the grouping bw relative to my residual post-alignment RT scatter on UHPLC."

### Gap-Filling and QC
> "Gap-fill with ChromPeakAreaParam but flag the filled values and report per-feature filled fraction."
> "Filter features with RsdFilter and DratioFilter using my QC and study sample indices."

### Redundancy and Export
> "Run CAMERA in the correct order to collapse adducts and isotopes before annotation."
> "Export the feature table with mz, rt, and intensities for downstream statistics."

## What the Agent Will Do

1. Load centroided mzML into an `MsExperiment` via `readMsExperiment` with a per-file sample table.
2. Run `findChromPeaks` with a `CentWaveParam` (or `MatchedFilterParam` for low-res/profile data) tuned to the instrument and chromatography.
3. Align retention times with `ObiwarpParam` or `PeakGroupsParam`, then regroup.
4. Group peaks into features with `PeakDensityParam`.
5. Gap-fill with `ChromPeakAreaParam`, tracking `is_filled`.
6. Collapse redundancy with CAMERA and optionally filter features with `filterFeatures`.
7. Export the features-by-samples matrix with mz/rt annotations and the full parameter set.

## Tips

- Use CentWave for centroided high-res data (most modern instruments); MatchedFilter only for low-res/quadrupole/profile data.
- Measure `peakwidth` from real EIC base-widths; the default c(20,50) silently drops sharp UHPLC peaks.
- Set `ppm` from empirical across-scan centroid scatter (~2-3x), not the datasheet - too tight fragments mass traces and deletes features.
- Align to a pooled QC or representative mid-batch run; aligning to file #1 propagates an outlier's idiosyncrasy into every warp.
- Inspect RT-deviation plots and a few EICs of significant features; alignment can look perfect in QCs while corrupting rare features.
- Treat filled values as imputations: track `is_filled`, report filled fractions, and prefer MNAR-aware imputation for inferential statistics.
- A finding that survives only one software/parameter set is a candidate, not a result; report the full processing specification.

## Related Skills

- metabolomics/normalization-qc - Drift correction, CV/D-ratio filtering, and feature-table normalization
- metabolomics/metabolite-annotation - Identification of features into named metabolites
- metabolomics/msdial-preprocessing - GUI/MS-DIAL alternative with MS2Dec deconvolution and GC-EI support
- metabolomics/statistical-analysis - Differential and multivariate statistics on the feature table
