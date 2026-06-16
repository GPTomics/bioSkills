---
name: bio-proteomics-proteomics-qc
description: Quality control for bottom-up proteomics across three levels -- instrument/raw-signal (mass accuracy, RT/iRT fit, FWHM, TIC vs injection time, % MS2 identified), identification/run (missed cleavages, charge states, PTM handling artifacts, contaminants), and experiment/quantitative (replicate correlation on log2, CV on the linear scale, completeness, MNAR-vs-MCAR missingness, PCA/batch, TMT channel balance, DIA q-values). Frames QC as a control chart against a per-instrument rolling baseline, not fixed cutoffs, and mandates inspecting raw boxplots, per-sample ID counts, total signal, and contaminant removal BEFORE normalizing -- because median normalization erases loading failures. Use when assessing proteomics data quality, diagnosing outlier samples, or deciding which samples to exclude before differential testing. The statistical test itself is differential-abundance; normalization mechanics are quantification; DIA q-value internals are dia-analysis.
tool_type: mixed
primary_tool: pandas
---

## Version Compatibility

Reference examples tested with: pandas 2.2+, numpy 1.26+, scipy 1.12+, matplotlib 3.8+, scikit-learn 1.4+, limma 3.58+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Proteomics Quality Control -- A Three-Level Funnel Where the Matrix Sees the Failure Last

**"Check the quality of my proteomics data"** -> Read instrument, identification, and quantitative metrics as a descending funnel of silent failures, and inspect raw signal BEFORE normalizing -- because by the time a fault reaches the deliverable matrix, normalization has usually erased the evidence.
- Python: `pandas` for matrix QC; `matplotlib`/`seaborn` for raw boxplots, correlation heatmaps, PCA
- R: PTXQC `createReport()` for MaxQuant search-table QC; `limma::plotMDS()`/`plotDensities()`; MSstatsTMT `dataProcessPlotsTMT()` for TMT channel balance

Scope: This skill OWNS QC diagnosis across all three levels -- which metric localizes which fault, what threshold means trouble, and the mandatory inspect-before-normalize ordering. Normalization mechanics route to quantification; the differential test routes to differential-abundance; DIA q-value computation routes to dia-analysis. OUT OF SCOPE: running the statistical test, the normalization algorithms themselves, and DIA q-value/FDR internals.

## The Single Most Important Modern Insight -- QC Is a Three-Level Funnel and Normalization Hides the Evidence

1. **QC is a three-level funnel of silent failures, and the deliverable matrix is the LAST place a problem becomes visible.** Faults originate at the instrument (spray, calibration, column) or in identification (digestion, contamination, PTM artifacts), but a protein matrix only shows the downstream symptom -- a low correlation or an outlier sample. A matrix-only QC pass is one-third of the job and blind to where faults actually start. Localize by descending: read every metric together with its co-readouts, never alone.

2. **Almost no metric has a universal pass/fail cutoff; the defensible practice is a per-instrument, per-method control chart.** Deviation from a lab's own rolling baseline (Levey-Jennings, +/-2 SD warn, +/-3 SD action) detects faults that a constant threshold misses or false-flags (Neely and Palmblad 2024). The numbers below seed a control chart, they are not standards-body limits.

3. **Median normalization HIDES loading problems that must be SEEN first.** Median/quantile normalization works by forcing a chosen summary statistic of every sample equal. A sample that genuinely loaded 3x low sits visibly shifted down in a RAW boxplot -- an obvious, diagnosable defect. The instant the matrix is median-normalized, the algorithm shifts that sample up by a constant to match everyone's median; boxplots line up perfectly; the evidence is mathematically erased. Worse, the low-loaded sample's noisy low-abundance signal gets stretched up to mid-range and injected into the differential test while QC plots look pristine. MANDATE: inspect raw/un-normalized boxplots plus per-sample ID counts, total signal, and missing fraction BEFORE normalizing; remove loading/injection failures and contaminants; THEN normalize and re-plot on the survivors. The identical principle governs TMT channel-loading balance.

## The Three Levels in One Table

| Level | Question | Inputs | Faults localized |
|-------|----------|--------|------------------|
| 1. Instrument / raw-signal | Is the LC-MS hardware performing? | Vendor `.raw`/`.d`; RawTools/RawBeans/rawrr/rawDiag/QuaMeter; Panorama AutoQC | Column, spray/emitter, mass analyzer/calibration |
| 2. Identification / run | Did this run identify peptides correctly? | Search tables (MaxQuant `txt/`, FragPipe `*.tsv`, DIA-NN report); PTXQC | Digestion, sample-handling PTM artifacts, contamination, FDR efficiency |
| 3. Experiment / quantitative | Are the numbers reproducible and comparable? | Protein/peptide intensity matrix; MSstatsTMT, pandas/limma | Loading/pipetting, batch, outliers, missingness, sample swaps |

The most integrative metric (% MS2 identified / ID count) is the first alarm but the LEAST specific -- it moves whenever anything upstream degrades. The same protein-count drop means spray (erratic TIC + maxed injection time), column (lost RT + broad peaks + rising backpressure), or sample (high contaminant fraction) depending on what co-moves.

## Tool Taxonomy

| Tool / method | Level | Citation | Mechanism / role | When |
|---------------|-------|----------|------------------|------|
| PTXQC (R, CRAN) | 2+3 | Bielow 2016 | `createReport()` over MaxQuant `txt/` or mzTab; per-metric scores in [0,1], QC heatmap PDF | MaxQuant output, fast multi-metric report |
| RawTools / RawBeans | 1 | Kovalchik 2019; Morgenstern 2021 | Parse Thermo `.raw` for IT, TIC, FWHM, scan timing | Diagnose instrument faults from raw files |
| rawrr / rawDiag | 1 | Kockmann 2021; Trachsel 2018 | R access to Orbitrap scan metadata | Custom Level-1 plots / method optimization |
| QuaMeter | 1+2 | Ma 2012 | Vendor-independent ID-free and ID-based metrics | Cross-vendor Level-1 QC |
| Skyline + Panorama AutoQC | 1 (longitudinal) | Bereman 2016 | Levey-Jennings + CUSUM/Moving-Range, SD-band flagging | System-suitability trending over time |
| MSstatsTMT | 3 (TMT) | Huang 2020 | `proteinSummarization()`, `dataProcessPlotsTMT()`; filters isolation interference on import | TMT channel-balance and QC plots |
| pandas / limma matrix QC | 3 | (this skill) | Correlation, CV, completeness, PCA on the matrix | Experiment-level QC (the code below) |
| differential-abundance | (route OUT) | -- | The moderated test itself | Hit calling after QC passes |
| quantification | (route OUT) | -- | Normalization and imputation mechanics | The how of normalizing |
| dia-analysis | (route OUT) | -- | DIA q-value/FDR internals | DIA-NN/Spectronaut report computation |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| MaxQuant `txt/` folder, want fast multi-metric report | PTXQC `createReport(txt_folder=...)` | Scores Level-2/3 metrics vs a representative file; one PDF |
| Protein-count drop, cause unknown | Descend to Level 1: read TIC + injection time + RT/FWHM together | Co-readouts localize spray vs column vs sample |
| Replicate correlation low for one sample | Check if it correlates better with a DIFFERENT group | Distinguishes sample swap from prep failure |
| Boxplots flat but a sample feels wrong | Re-plot the RAW (un-normalized) matrix | Normalization erased the loading evidence |
| Deciding how to impute | Diagnose MNAR (left tail) vs MCAR (all-abundance) from the histogram FIRST | Wrong imputer corrupts present/absent calls |
| TMT data, channel looks off | MSstatsTMT QC plots on RAW reporter intensities | See the imbalance before global median rescales it |
| DIA matrix, how many proteins are real | Filter Global.Q.Value and Global.PG.Q.Value, route q internals to dia-analysis | Precursor q != protein q; both needed |
| Long sample queue, drift suspected | Interspersed QC every 4th-5th injection + Levey-Jennings | Turns one check into a time series |

Default when uncertain: plot the RAW per-sample boxplots, ID counts, total signal, and missing fraction first; remove loading/injection failures and contaminants; only then normalize, re-plot, and proceed to correlation/CV/PCA on the survivors.

## Inspect Raw Signal and Remove Contaminants Before Normalizing

**Goal:** Catch loading/injection failures and strip contaminant/decoy rows while they are still visible -- before normalization erases them.

**Approach:** Load the un-normalized matrix, plot per-sample boxplots plus ID counts and total signal, filter MaxQuant `Potential contaminant`/`Reverse`/`Only identified by site` rows, THEN log-transform and normalize on the survivors.

```python
import pandas as pd
import numpy as np

contaminant_flags = ['Potential contaminant', 'Reverse', 'Only identified by site']

def strip_contaminant_rows(protein_groups):
    keep = pd.Series(True, index=protein_groups.index)
    for col in contaminant_flags:
        match = next((c for c in protein_groups.columns if c.lower() == col.lower()), None)  # MaxQuant casing varies by version -- match case-insensitively
        if match is not None:
            keep &= protein_groups[match].fillna('') != '+'  # MaxQuant marks flagged rows with a literal '+'
    return protein_groups[keep]

def raw_sample_qc(raw_intensities):
    return pd.DataFrame({
        'n_quantified': raw_intensities.notna().sum(),
        'total_signal': raw_intensities.sum(),
        'median_intensity': raw_intensities.median(),
        'missing_pct': 100 * raw_intensities.isna().sum() / len(raw_intensities)})
```

Read the boxplots before normalizing: a sample shifted >=2-3x below its group median is a loading/injection failure to exclude, not to rescale. The contaminant fraction of summed intensity should be small (PTXQC default flags >1%); keratin and trypsin autolysis dominate LOW-INPUT samples (single-cell, IPs, gel bands) because they are a roughly fixed absolute amount whose fractional share explodes as load shrinks.

## Replicate Correlation on log2

**Goal:** Quantify reproducibility without letting a few abundant proteins fake agreement.

**Approach:** Correlate on log2 intensities (variance-stabilized, high-abundance tail compressed), report within-group pairs, and flag a sample correlating better with another group as a possible swap.

```python
from itertools import combinations

def replicate_correlation(log2_intensities, sample_groups):
    corr = log2_intensities.corr(method='pearson')  # log2 first: Pearson on raw is a high-abundance artifact
    rows = []
    for group in sample_groups.unique():
        members = sample_groups[sample_groups == group].index
        for s1, s2 in combinations(members, 2):
            rows.append({'group': group, 's1': s1, 's2': s2, 'r': corr.loc[s1, s2]})
    return pd.DataFrame(rows)
```

Technical replicates r > 0.98 (instrument noise only); biological r ~ 0.90-0.98 (genuine variance, lower is expected and correct); soft floor r > 0.8 to retain a biological replicate. A Spearman check is a robustness aid only -- ranks discard the magnitude that quant QC cares about.

## Coefficient of Variation on the Linear Scale

**Goal:** Summarize per-condition precision with a number that means what it says.

**Approach:** Compute CV = SD/mean on LINEAR (non-log) intensities; if only logged values exist use the geometric-CV formula. Report the median CV per condition (the per-protein distribution is right-skewed).

```python
def median_cv_linear(linear_intensities, sample_groups):
    rows = []
    for group in sample_groups.unique():
        block = linear_intensities[sample_groups[sample_groups == group].index]
        per_protein_cv = block.std(axis=1) / block.mean(axis=1)  # base CV formula REQUIRES linear scale
        rows.append({'group': group, 'median_cv_pct': 100 * per_protein_cv.median()})
    return pd.DataFrame(rows)

def geometric_cv_from_log(log_intensities):
    sigma = log_intensities.std(axis=1) * np.log(2)  # convert log2 SD to natural-log SD
    return 100 * np.sqrt(np.expm1(sigma ** 2))  # gCV = sqrt(exp(sigma^2) - 1)
```

Applying the base formula to log-transformed data compresses CV ~14x (most proteins appear to have CV < 1%) -- meaningless (Brenes 2024). State normalization state, transform, and software params or the CV is uninterpretable: DIA-NN "High precision" mode silently median-normalizes, halving median CV vs "High accuracy". Technical median CV < ~10-20%, biological ~20-40%; a LOWER CV is not automatically better (loose FDR or faulty MS1 extraction produce artificially low CVs).

## Missingness Mechanism and Completeness

**Goal:** Decide how to impute by first deciding why values are missing.

**Approach:** Diagnose the missingness profile -- left-tail concentration means MNAR (left-censored, abundance-dependent), all-abundance scatter means MCAR -- and filter on completeness before imputing only the shallow remainder.

```python
def missingness_profile(log2_intensities, n_bins=20):
    observed = log2_intensities.stack()
    abundance_bins = pd.qcut(observed, n_bins, duplicates='drop')
    present_per_protein = log2_intensities.notna().mean(axis=1)
    mean_abundance = log2_intensities.mean(axis=1)
    return mean_abundance, present_per_protein  # plot present-fraction vs abundance: rising-with-abundance = MNAR

def completeness_filter(log2_intensities, sample_groups, min_valid_frac=0.7):
    keep = pd.Series(False, index=log2_intensities.index)
    for group in sample_groups.unique():
        block = log2_intensities[sample_groups[sample_groups == group].index]
        keep |= block.notna().mean(axis=1) >= min_valid_frac  # valid in >=70% of >=1 condition
    return log2_intensities[keep]
```

kNN-imputing a genuinely-absent (MNAR) value invents mid-range abundance and KILLS a real present/absent difference; a left-shifted draw (Perseus down-shifted normal, downshift=1.8 SD below the observed mean, width=0.3 of observed SD) on an MCAR gap FABRICATES a false low and inflates a difference. Match the imputer to the mechanism. The imputation mechanics themselves are quantification.

## PCA and Batch Detection

**Goal:** See whether the dominant variance is biology or batch, and flag outlier samples.

**Approach:** On the normalized survivors, run PCA, color by condition and by batch, and test whether top PCs associate with batch.

```python
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import f_oneway

def pca_batch_check(normalized_log2, sample_info, batch_col='batch'):
    imputed = normalized_log2.apply(lambda r: r.fillna(r.median()), axis=1)  # temporary, for PCA only
    pcs = PCA(n_components=5).fit(StandardScaler().fit_transform(imputed.T))
    coords = pd.DataFrame(pcs.transform(StandardScaler().fit_transform(imputed.T)),
                          columns=[f'PC{i+1}' for i in range(5)], index=normalized_log2.columns).join(sample_info)
    for pc in ['PC1', 'PC2', 'PC3']:
        groups = [coords[coords[batch_col] == b][pc] for b in coords[batch_col].unique()]
        _, p = f_oneway(*groups)
        print(f'{pc} ~ {batch_col}: p={p:.4f}')
    return coords, pcs.explained_variance_ratio_
```

A sample isolated from its group is a removal/re-run candidate. If batch is PC1, correct it explicitly (ComBat, or include batch in the design matrix downstream) and re-inspect; never let batch be the dominant axis going into differential testing. Visualization of the projection routes to data-visualization/dimensionality-reduction-plots.

## Per-Method Failure Modes

### Median normalization hides loading failures
**Trigger:** Normalizing the matrix before inspecting raw per-sample signal. **Mechanism:** median-centering shifts each sample by a constant to equalize the very statistic that was the symptom of a low load. **Symptom:** flat, clean boxplots that hide a 3x-low sample now stretched into mid-range. **Fix:** plot RAW boxplots + ID counts + total signal first; exclude failures; then normalize.

### Normalizing with contaminants still in the matrix
**Trigger:** Contaminant/decoy rows left in before log + normalize. **Mechanism:** keratin/trypsin/albumin inflate the denominator and shift the median; when their load differs across groups the differential gets normalized into the real proteins. **Symptom:** spurious fold changes; a contaminant fraction that varies by group. **Fix:** filter `Potential contaminant` + `Reverse` + `Only identified by site` BEFORE log + normalize.

### Wrong imputer for the missingness mechanism
**Trigger:** kNN on MNAR, or left-shift on MCAR. **Mechanism:** kNN borrows mid-range neighbors for a value that is low because it is absent; left-shift draws a deep low for a value missing at random. **Symptom:** killed present/absent calls (kNN-on-MNAR) or inflated false lows (left-shift-on-MCAR). **Fix:** diagnose left-tail vs all-abundance from the histogram first; mechanics route to quantification.

### CV computed on log-transformed data
**Trigger:** Base CV formula applied after log2. **Mechanism:** SD/mean is defined for linear intensity; logging compresses it ~14x. **Symptom:** most proteins appear to have CV < 1%. **Fix:** compute on linear intensity, or use the geometric-CV formula on logged values; always state transform + normalization + software.

### Pearson r on raw intensity
**Trigger:** Correlating un-logged intensities. **Mechanism:** a few high-abundance proteins dominate the covariance. **Symptom:** r = 0.99 while the bulk disagrees. **Fix:** log2 before correlating; Spearman as a robustness check only.

### Bimodal mass-error histogram read as calibration
**Trigger:** A two-peaked ppm-error distribution. **Mechanism:** almost always monoisotopic mis-assignment or co-isolation (a search/sample problem), NOT calibration drift; a generous tolerance still IDs the mis-assigned precursors so ID rate looks fine. **Symptom:** bimodal histogram, normal ID rate. **Fix:** correct isotope-error tolerance/deisotoping, not recalibration.

### Charge-state distribution off the platform baseline
**Trigger:** The fully-tryptic 2+ fraction drifts from the rolling baseline. **Mechanism:** a tryptic peptide carries two basic sites (C-terminal K/R + N-terminus) so 2+ dominates; excess 3+/4+ comes from internal basic residues left by missed cleavages, excess 1+ from poor ionization, short peptides, or contaminants. **Symptom:** raised high-charge fraction (a digestion/chemistry signal) or raised 1+ fraction (an ionization/spray signal). **Fix:** read the charge distribution together with the missed-cleavage rate (high charge co-moving with missed cleavages = under-digestion) to separate a chemistry problem from a spray problem a raw protein-count drop cannot resolve alone.

### TMT ratio compression from co-isolation
**Trigger:** Isobaric quant with a wide isolation window. **Mechanism:** near-isobaric co-eluting precursors are co-isolated and add their own reporters across all channels, a uniform pedestal. **Symptom:** every fold-change compressed toward 1 (a real 10:1 reads ~5:1). **Fix:** filter isolation interference < 50%; prefer SPS-MS3 (McAlister 2014) / FAIMS / narrow windows; run a TKO or empty-channel control to measure the floor.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Mass accuracy (internal/lock) median \|err\| < 1-3 ppm, single-mode centered 0 | CONVENTION | width approaching MS1 tolerance loses real IDs |
| iRT RT-fit R^2 > 0.99 (warn below) | CONVENTION (mechanism firm) | residual is just LC noise on a stable gradient |
| FWHM alarm on > 20-30% rise vs baseline; >= 8-10 points across FWHM (floor ~5) | Kocher 2011 | peak capacity tracks peptide IDs; points needed for accurate AUC |
| % MS2 identified: < 20% bad, 20-35% ok, >= 35% great | PTXQC `createYaml.R` | generic pass marks, not a biological ceiling |
| Missed cleavages: >= 75-85% at 0 MC; flag > 25-30% with >= 1 MC | OPERATIONAL (PTXQC-style) | porcine trypsin ~78% efficient even ideally |
| Replicate Pearson r (log2): technical > 0.98, biological 0.90-0.98, floor 0.8 | CONVENTION (mechanism firm) | log2 variance-stabilizes; biological variance is real |
| Median CV (linear): technical < 10-20%, biological 20-40% | CONVENTION | DIA < DDA (no stochastic sampling); lower not always better |
| Completeness: valid in >= 50-70% of replicates in >= 1 condition | CONVENTION | filter before imputing |
| Perseus MNAR imputation: downshift = 1.8 SD, width = 0.3 | Tyanova 2016 | deep left tail simulates below-LOD, narrowed so not mistaken for real |
| TMT channel deviation: investigate > ~2x, flag > ~3-4x | CONVENTION | pipetting/labeling vs biology |
| Isolation interference < 50% PSM filter (< 30% stricter) | CONVENTION (PD practice) | above ~50% contaminant dominates, ratios uninterpretable |
| DIA precursor + protein q both <= 0.01, GLOBAL q, PICKED protein estimator | STANDARD | precursor != protein FDR; route internals to dia-analysis |
| Levey-Jennings: +/-2 SD warn, +/-3 SD action; QC every 4th-5th injection | Bereman 2016 | ~95/99.7% of points under stable normal; catch drift early |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `KeyError: 'Mass Error [ppm]'` | MaxQuant column casing varies by version | match case-insensitively (`Mass error` vs `Mass Error`); never hard-code |
| Contaminant rows have `True`/`False`, filter keeps all | MaxQuant flags with a literal `'+'`, not a boolean | filter `col != '+'` |
| CV unexpectedly tiny (< 1%) | base CV formula applied to log2 data | compute on linear intensity or use geometric CV |
| `r = 0.99` but samples clearly differ | Pearson on raw (un-logged) intensity | log2 transform before correlating |
| PCA dominated by injection day | batch effect, not biology | correct (ComBat / batch in design) and re-inspect; do not proceed |
| PTXQC "not found" via `BiocManager` | PTXQC is on CRAN, not Bioconductor | `install.packages('PTXQC')` |
| `createReport()` errors on a dataframe arg | it takes a txt-folder path / mzTab / YAML, not dataframes | pass `txt_folder=` (the MaxQuant `txt/` directory) |

## References

- Bielow C, Mastrobuoni G, Kempa S. Proteomics Quality Control: Quality Control Software for MaxQuant Results. *J Proteome Res* 2016;15(3):777-787.
- Kovalchik KA, Colborne S, Spencer SEP, et al. RawTools: Rapid and Dynamic Interrogation of Orbitrap Data Files for Mass Spectrometer System Management. *J Proteome Res* 2019;18(2):700-708.
- Morgenstern D, Barzilay R, Levin Y. RawBeans: A Simple, Vendor-Independent, Raw-Data Quality-Control Tool. *J Proteome Res* 2021;20(4):2098-2104.
- Kockmann T, Panse C. The rawrr R Package: Direct Access to Orbitrap Data and Beyond. *J Proteome Res* 2021;20(4):2028-2034.
- Trachsel C, Panse C, Kockmann T, et al. rawDiag: An R Package Supporting Rational LC-MS Method Optimization for Bottom-up Proteomics. *J Proteome Res* 2018;17(8):2908-2914.
- Ma ZQ, Polzin KO, Dasari S, et al. QuaMeter: Multivendor Performance Metrics for LC-MS/MS Proteomics Instrumentation. *Anal Chem* 2012;84(14):5845-5850.
- Bereman MS, Beri J, Sharma V, et al. An Automated Pipeline to Monitor System Performance in Liquid Chromatography-Tandem Mass Spectrometry Proteomic Experiments. *J Proteome Res* 2016;15(12):4763-4769.
- Kocher T, Swart R, Mechtler K. Ultra-High-Pressure RPLC Hyphenated to an LTQ-Orbitrap Velos Reveals a Linear Relation between Peak Capacity and Number of Identified Peptides. *Anal Chem* 2011;83(7):2699-2704.
- Tyanova S, Temu T, Sinitcyn P, et al. The Perseus computational platform for comprehensive analysis of (prote)omics data. *Nat Methods* 2016;13(9):731-740.
- Brenes AJ. Calculating and Reporting Coefficients of Variation for DIA-Based Proteomics. *J Proteome Res* 2024;23(12):5274-5278.
- McAlister GC, Nusinow DP, Jedrychowski MP, et al. MultiNotch MS3 Enables Accurate, Sensitive, and Multiplexed Detection of Differential Expression across Cancer Cell Line Proteomes. *Anal Chem* 2014;86(14):7150-7158.
- Huang T, Choi M, Tzouros M, et al. MSstatsTMT: Statistical Detection of Differentially Abundant Proteins in Experiments with Isobaric Labeling and Multiple Mixtures. *Mol Cell Proteomics* 2020;19(10):1706-1723.
- Neely BA, Palmblad M, et al. Quality Control in the Mass Spectrometry Proteomics Core: A Practical Primer. *J Biomol Tech* 2024;35(3).

## Related Skills

- data-import - Load search-engine output and intensity matrices before QC
- quantification - Normalization and imputation mechanics that QC mandates running AFTER inspection
- differential-abundance - The moderated statistical test QC gates
- dia-analysis - DIA q-value/FDR internals behind the protein-count QC
- data-visualization/dimensionality-reduction-plots - PCA/MDS projection plotting
- workflows/proteomics-pipeline - End-to-end pipeline placing QC before differential testing
