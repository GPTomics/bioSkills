---
name: bio-longitudinal-monitoring
description: Tracks ctDNA across serial liquid-biopsy timepoints for molecular residual disease (MRD) and treatment-response monitoring, treating MRD as a binary integrated detection call across the patient's full variant set (with a defined LoD95 and per-sample specificity) rather than a per-timepoint VAF threshold, and handling undetectable samples as left-censored at the per-sample limit of detection rather than true zeros. Covers tumor-informed bespoke vs tumor-naive design, landmark vs surveillance sampling, molecular-response definitions and their non-standardization, censoring-aware clearance kinetics, and the multiple-testing structure of repeated surveillance. Use when monitoring ctDNA during therapy, calling molecular relapse before imaging, or estimating clearance half-life from serial samples.
tool_type: python
primary_tool: pandas
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, pandas 2.2+, scipy 1.12+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: this skill is statistical, not tool-bound. The hard parts are interpretive (left-censoring, multiple testing, lead-time bias), not API calls. `scipy.stats.linregress` returns a named tuple whose `.slope`/`.pvalue` attributes are stable across recent versions; the censoring-aware fit below uses only `linregress` on the uncensored decay phase plus a manual interval check, so version drift is low-risk.

# Longitudinal Monitoring

**"Track ctDNA over this course of treatment"** -> Integrate serial plasma measurements into a binary detected/not-detected trajectory plus censoring-aware burden kinetics for MRD and response monitoring.
- Python: `pandas` for the per-timepoint table, `scipy.stats` for censoring-aware decay/trend, `matplotlib` for log-scale trajectory plots

## The Single Most Important Modern Insight -- MRD is a binary integrated detection call, and "undetectable" is left-censored, not zero

A tumor-informed MRD assay does not ask "is the VAF at locus X above a threshold?" It integrates signal across the patient's entire personal variant set (16 to 500+ loci) into ONE detected/not-detected call with a defined LoD95 (the tumor fraction detected 95% of the time at a given input) and a per-sample specificity. Signal invisible at any single 0.001%-VAF locus becomes significant when summed across hundreds of loci against a modeled error background; this is why bespoke assays reach 10^-4 to 10^-6 tumor fraction. The detected/not-detected call is the unit of analysis -- per-locus VAF is plumbing, not the readout. Re-deriving a per-timepoint "VAF < X" cutoff throws away the multi-locus integration that makes MRD work and inflates false positives from a single noisy locus.

The second half of the insight: an "undetectable" result is conditional on how many genome-equivalents were interrogated. VAF=0 is LEFT-CENSORED at the per-sample LoD, not a true zero. A 10 mL tube yields roughly 50 ng cfDNA, around 15,000 haploid genome-equivalents; at 0.01% tumor fraction that is roughly 1.5 expected tumor molecules, squarely in the Poisson-limited regime (lambda < 3) where detection is stochastic. "Undetectable" at a low-input draw may simply mean the assay could not have seen the burden it saw at a higher-input draw. Every undetectable must carry its per-sample LoD; plugging 0 into a log-fit or fold-change biases everything and log(0) breaks the fit outright.

## Design Decision: tumor-informed vs tumor-naive, landmark vs surveillance

| Axis | Tumor-informed bespoke | Tumor-naive (fixed panel) |
|------|------------------------|---------------------------|
| Variant set | Patient-specific, designed from tumor/normal WES/WGS | Fixed gene panel, identical across patients |
| Examples | Signatera (16 SNVs, Reinert 2019), RaDaR (up to ~48 amplicons), INVAR (hundreds-thousands of loci, Wan 2020) | Broad cfDNA panels, sWGS |
| MRD sensitivity | Very high (10^-4 to 10^-6 TF); LoD scales with #loci x input | Lower for MRD; few loci per region |
| Needs tumor tissue | Yes (design step, weeks of turnaround) | No (tissue-free, faster) |
| CHIP confounding | Low (tracks known tumor somatic variants) | High (de novo calls include clonal hematopoiesis) |
| Best use | Defined-burden MRD/surveillance after curative intent | No tissue available, or broad genotyping in metastatic disease |

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Post-curative-intent MRD / recurrence surveillance | Tumor-informed bespoke, binary call | Reaches ppm LoD by integrating across the personal variant set; CHIP-resistant |
| No tissue available, metastatic response monitoring | Tumor-naive panel, track aggregated burden | Tissue-free and immediate; accept higher LoD and mandatory CHIP control |
| Post-surgical landmark (single decisive timepoint) | One draw at ~2-10 weeks post-op | Avoids the surgical cfDNA surge; conventional Week 4 default |
| Serial surveillance over months-years | Trend over >=2 consecutive draws, confirm before acting | Each draw is another false-positive opportunity (multiple testing) |
| Defining "molecular response" | Use the assay's own validated cutoff; do not import "2-log" or "90%" blindly | Cutoffs are non-harmonized across assays (see below) |

Methodology evolves: verify the current best practice and the assay's validated definitions against the latest tool/vendor documentation before fixing any threshold in code.

## ctDNA Kinetics Biology -- why timing and shedding gate interpretation

ctDNA has a plasma half-life of roughly 2 h (114 min, Diehl 2008); broader literature spans ~16 min to 2.5 h. This fast turnover is the entire reason serial monitoring works: plasma concentration tracks CURRENT tumor flux, not a weeks-old average, while imaging tumor volume lags. The same fast clearance makes landmark timing fragile. Surgery dumps a transient cfDNA surge into plasma (tissue trauma, wound healing, neutrophil extracellular traps) that dilutes tumor fraction and can transiently raise total cfDNA. Drawing at post-op day 1-3 reads this surge, not residual disease; the conventional landmark window is ~2-10 weeks (Week 4 a frequent default). A clearance fit that includes a post-op surge point mis-estimates the half-life.

Shedding is not uniform. "ctDNA-negative" does NOT equal "disease-free": some early lung adenocarcinomas and indolent/low-volume tumors shed below detectable thresholds, and brain metastases behind the blood-brain barrier shed poorly into plasma (CSF is the better CNS compartment). A patient can have radiographic progression with clean plasma. Negativity has high negative predictive value for relapse in shedding tumors but is never a guarantee -- imaging stays mandatory for low/non-shedders and sanctuary sites.

## Molecular-Response Definitions -- and the non-standardization caveat

| Term | Representative operationalization | Caveat |
|------|-----------------------------------|--------|
| Molecular response (MR) | >= 90% drop (ctMoniTR-style) or >= 2-log/100x (immuno/heme heritage) from baseline | 2-log and 90% are different magnitudes; cutoff is study/assay-specific |
| Molecular complete response (mCR) | ctDNA becomes undetectable | "Undetectable" is LoD-conditional, not zero |
| ctDNA clearance | Sustained detectable -> undetectable | Depends on input/depth of the clearing draw; confirm with re-draw |
| Molecular progression / relapse | Confirmed re-detection or rise-from-nadir | Require trend over >=2 draws (multiple testing) |

These definitions are NOT harmonized. "2-log reduction," "90% reduction," and "molecular complete response" are assay- and study-dependent, not interchangeable. The Friends of Cancer Research ctMoniTR project is the field's standardization attempt (pooling ctDNA-change data across NSCLC immunotherapy studies to validate ctDNA change as an intermediate endpoint), not a settled standard. The FDA ctDNA guidance for curative-intent solid-tumor drug development was issued as a draft in May 2022 and finalized in November 2024; it endorses ctDNA for patient selection, MRD-based enrichment, and as a measure of response, but does NOT yet endorse ctDNA change as a validated surrogate endpoint for DFS/EFS/OS. Code should accept the assay's own validated cutoff rather than baking one in.

## Clinical Evidence and Lead Time

ctDNA MRD predicts relapse months before imaging across tumor types: breast median ~8 mo (Garcia-Murillas 2015), NSCLC median ~5.2 mo (Chaudhuri 2017), CRC mean ~8.7 mo (Reinert 2019); TRACERx phylogenetic ctDNA tracks clonal evolution and metastatic seeding (Abbosh 2017, 2023). The interventional landmark is DYNAMIC (Tie 2022): a ctDNA-guided strategy in stage II colon cancer reduced adjuvant chemotherapy use (15% vs 28%) without compromising 2-year recurrence-free survival, proving an MRD-negative call can justify de-escalation. Caveat -- lead-time bias: "ctDNA detects relapse N months before imaging" is a real analytic-sensitivity advantage, but measuring survival from molecular detection vs clinical detection merely moves the clock back and inflates apparent survival. Demonstrating clinical utility (that acting on the earlier signal improves outcomes) requires an interventional design like DYNAMIC, not earlier detection alone. Flag lead-time bias wherever lead time is reported.

## Tumor-Fraction Trend with Baseline and Nadir

**Goal:** Summarize a serial trajectory into baseline, nadir, and baseline-referenced change, with below-LoD points marked as censored, not zero.

**Approach:** Sort by time, carry a per-sample LoD column, flag any point at-or-below its LoD as left-censored, and compute log-fold change from baseline only on the uncensored estimates (substituting the LoD bound, never 0, for censored points).

```python
import numpy as np
import pandas as pd

def summarize_trajectory(df):
    '''df columns: timepoint, tumor_fraction, per_sample_lod (genome-equivalent-aware).'''
    df = df.sort_values('timepoint').copy()
    df['censored'] = df['tumor_fraction'] <= df['per_sample_lod']
    df['tf_for_log'] = np.where(df['censored'], df['per_sample_lod'], df['tumor_fraction'])
    baseline = df.iloc[0]['tf_for_log']
    df['log2_fc_baseline'] = np.log2(df['tf_for_log'] / baseline)
    detected = df[~df['censored']]
    nadir = detected['tumor_fraction'].min() if len(detected) else np.nan
    return df, {'baseline_tf': baseline, 'nadir_tf': nadir, 'n_censored': int(df['censored'].sum())}
```

## Mutation Tracking and Censoring-Aware Clearance Kinetics

**Goal:** Pivot per-mutation VAF over time and estimate a clearance half-life only over the genuine decay phase, treating below-LoD timepoints as censored.

**Approach:** Build a timepoint-by-mutation pivot, mark cleared loci as below-LoD (not missing-equals-zero), then fit ln(VAF) ~ time by OLS over the monotonic-decay phase only, excluding the surgical-surge point, any post-nadir rebound, and all censored points; half-life = ln(2)/(-slope).

```python
from scipy import stats

def clearance_half_life(df, lod):
    '''df columns: timepoint, vaf for one mutation. lod = per-sample detection bound.
       Fits the uncensored decay phase up to the nadir (drops post-nadir rebound and
       every below-LoD point); never feeds log(0) into the OLS.'''
    df = df.sort_values('timepoint')
    uncensored = df[df['vaf'] > lod].reset_index(drop=True)
    if len(uncensored) < 3:
        return None
    decay = uncensored.iloc[:uncensored['vaf'].idxmin() + 1]
    if len(decay) < 3:
        return None
    fit = stats.linregress(decay['timepoint'].values, np.log(decay['vaf'].values))
    half_life = np.log(2) / -fit.slope if fit.slope < 0 else np.inf
    return {'half_life_days': half_life, 'slope': fit.slope, 'r_squared': fit.rvalue ** 2,
            'n_points': len(decay), 'n_censored_excluded': int((df['vaf'] <= lod).sum())}
```

## Molecular-Relapse Calling

**Goal:** Call molecular relapse from confirmed re-detection or sustained rise-from-nadir, not a single excursion.

**Approach:** Find the nadir, then require detection (above per-sample LoD) on >=2 consecutive post-nadir draws, or a rise above an assay-defined margin above nadir confirmed on a re-draw; annotate every call with the draw's per-sample LoD so a low-LoD draw is not mistaken for new biology.

```python
def call_molecular_relapse(df, rise_factor=2.0, min_consecutive=2):
    '''df columns: timepoint, tumor_fraction, per_sample_lod. Requires a confirmed trend.'''
    df = df.sort_values('timepoint').copy()
    df['detected'] = df['tumor_fraction'] > df['per_sample_lod']
    nadir_time = df.loc[df['tumor_fraction'].idxmin(), 'timepoint']
    nadir_tf = max(df['tumor_fraction'].min(), df['per_sample_lod'].min())  # floor at LoD, not a censored value
    post = df[df['timepoint'] > nadir_time]
    consec = (post['detected'] & (post['tumor_fraction'] > nadir_tf * rise_factor)).astype(int)
    run = consec.groupby((consec == 0).cumsum()).cumsum().max() if len(consec) else 0
    relapse = bool(run >= min_consecutive)
    return {'relapse': relapse, 'nadir_tf': nadir_tf, 'confirmed_consecutive': int(run)}
```

## Per-Method Failure Modes

### Naive per-timepoint VAF thresholding
Trigger: applying "VAF < X" per timepoint to a multi-locus assay. Mechanism: discards the integration that makes MRD work; one noisy locus calls positive. Symptom: inflated false positives, jumpy trajectory. Fix: use the assay's integrated binary detected/not-detected call across the full variant set.

### Treating undetectable as a true zero
Trigger: plugging 0.0 (or VAF/2) into a log-fit or fold-change. Mechanism: log(0) is undefined; substituting a small number biases slope and fold-change. Symptom: NaN/inf fits or implausibly fast clearance. Fix: treat below-LoD as left-censored at the per-sample LoD; report "below LoD = X," never "0%."

### Ignoring per-timepoint LoD changes with input mass
Trigger: comparing "undetectable" across draws of different cfDNA input. Mechanism: LoD is input-conditional; a low-input draw could not have seen the prior burden. Symptom: spurious "clearance" or "relapse" at draws with anomalous input. Fix: carry per-sample LoD/genome-equivalents as a covariate; down-weight low-input negatives.

### CHIP rising over time read as relapse
Trigger: tumor-naive longitudinal panel without matched WBC sequencing. Mechanism: clonal hematopoiesis clones expand over time and under chemotherapy (Razavi 2019: majority of plasma variants are CHIP-derived), producing a rising non-tumor "ctDNA" signal. Symptom: false molecular progression in DNMT3A/TET2/ASXL1 hotspots. Fix: tumor-informed tracking, or matched serial WBC sequencing / known-CHIP-gene blacklisting.

### Lead-time bias in outcome claims
Trigger: reporting survival from molecular detection vs clinical detection. Mechanism: moving the detection clock back inflates apparent survival without changing outcome. Symptom: a "benefit" that is an artifact of earlier detection. Fix: claim clinical utility only from interventional designs (DYNAMIC); label lead time as analytic sensitivity, not benefit.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| ctDNA plasma half-life ~114 min (~2 h) | Diehl 2008, *Nat Med* 14:985 | Single-patient post-op estimate; broader range 16 min-2.5 h. Sets why monitoring works and why post-op timing matters |
| Post-op landmark window ~2-10 weeks (Week 4 common) | Convention across DYNAMIC/Signatera/RaDaR | Late enough for the surgical cfDNA surge to clear before reading residual disease |
| ~15,000 haploid genome-equivalents per 10 mL tube (~50 ng cfDNA, ~300 GE/ng) | Standard biophysical constants | Hard sampling floor: at 0.01% TF that is ~1.5 expected tumor molecules |
| Poisson-limited regime: lambda < 3 expected tumor molecules | Poisson detection theory (1-e^-3=0.95) | Below this, detection is stochastic; small input/recovery shifts flip a result across the LoD |
| Bespoke MRD sensitivity 10^-4 to 10^-6 tumor fraction | INVAR (Wan 2020), Signatera (Reinert 2019), RaDaR | Achieved only by integrating across the personal variant set, not per-locus |
| Molecular response ~2-log (100x) or ~90% reduction | ctMoniTR (Vega 2022); immuno/heme heritage | Non-harmonized convention -- surface the assay's own validated cutoff, do not hard-code |
| Per-course false-positive risk = 1 - s^n (s = per-sample specificity, n = draws) | Multiple-testing arithmetic | s=0.995,n=12 -> ~5.8%; s=0.99,n=12 -> ~11%. Report specificity per monitoring course, confirm positives by re-draw |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| log(0) / inf in clearance fit | Censored (below-LoD) point fed into ln(VAF) | Fit only the uncensored decay phase; carry LoD as the censoring bound |
| Half-life implausibly short or fits the surge | Post-op surge or post-nadir rebound point included | Restrict the OLS to the monotonic on-treatment decay phase |
| "Relapse" at one noisy draw | Acting on a single positive | Require >=2 consecutive confirmed draws; re-draw before acting |
| Rising signal in a tissue-free panel mistaken for tumor | CHIP drift over time | Matched WBC sequencing or tumor-informed tracking |
| "Negative = cured" overcall | Low-shedder / sanctuary site / low-input draw | Frame negative as below-LoD; keep imaging for low/non-shedders |

## References

- Diehl F, et al. 2008. Circulating mutant DNA to assess tumor dynamics. *Nature Medicine* 14:985-990. -- ctDNA half-life ~114 min; post-resection clearance kinetics.
- Tie J, et al. 2022. Circulating Tumor DNA Analysis Guiding Adjuvant Therapy in Stage II Colon Cancer (DYNAMIC). *New England Journal of Medicine* 386:2261-2272. -- Interventional MRD-guided de-escalation.
- Abbosh C, et al. 2017. Phylogenetic ctDNA analysis depicts early-stage lung cancer evolution. *Nature* 545:446-451. -- TRACERx tumor-informed clonal tracking.
- Abbosh C, et al. 2023. Tracking early lung cancer metastatic dissemination in TRACERx using ctDNA. *Nature* 616:553-562. -- Deep tumor-informed surveillance (~200 mutations).
- Garcia-Murillas I, et al. 2015. Mutation tracking in circulating tumor DNA predicts relapse in early breast cancer. *Science Translational Medicine* 7:302ra133. -- Lead time ~8 mo in breast cancer.
- Wan JCM, et al. 2020. ctDNA monitoring using patient-specific sequencing and integration of variant reads (INVAR). *Science Translational Medicine* 12:eaaz8084. -- Per-sample LoD derived from informative reads; the formal binary-integration framing.
- Reinert T, et al. 2019. Analysis of Plasma Cell-Free DNA by Ultradeep Sequencing in Patients With Stages I to III Colorectal Cancer. *JAMA Oncology* 5:1124-1131. -- Signatera; serial HR ~43.5; mean lead ~8.7 mo.
- Chaudhuri AA, et al. 2017. Early Detection of Molecular Residual Disease in Localized Lung Cancer by Circulating Tumor DNA Profiling. *Cancer Discovery* 7:1394-1403. -- CAPP-Seq; median lead ~5.2 mo.
- Razavi P, et al. 2019. High-intensity sequencing reveals the sources of plasma circulating cell-free DNA variants. *Nature Medicine* 25:1928-1937. -- Majority of plasma cfDNA variants are CHIP-derived; matched WBC sequencing essential.
- Merino Vega D, et al. 2022. Changes in Circulating Tumor DNA Reflect Clinical Benefit Across Multiple Studies of Patients With Non-Small-Cell Lung Cancer Treated With Immune Checkpoint Inhibitors. *JCO Precision Oncology* 6:e2100372. -- Friends of Cancer Research ctMoniTR Step 1; the molecular-response standardization effort, not a settled threshold.

The RaDaR/LUCID early-NSCLC residual-ctDNA study (*Annals of Oncology* 2022, 33:500-510) is referenced generically above; verify first-author attribution and the exact article identifier against the journal record before citing it formally.

## Related Skills

- ctdna-mutation-detection - detect the variant set that is then tracked
- tumor-fraction-estimation - per-timepoint tumor burden
- analytical-validation - per-timepoint LoD and left-censoring of undetectable samples
- fragment-analysis - fragmentomic trends as a complementary monitoring signal
- clinical-biostatistics/survival-analysis - relapse, lead-time, and endpoint analysis
