# Longitudinal Monitoring - Usage Guide

## Overview
Tracks ctDNA across serial liquid-biopsy timepoints for molecular residual disease (MRD) and treatment-response monitoring. Treats MRD as a binary integrated detected/not-detected call across the patient's full variant set rather than a per-timepoint VAF threshold, handles undetectable samples as left-censored at the per-sample limit of detection, and estimates clearance kinetics and molecular relapse with that censoring respected.

## Prerequisites
```bash
pip install pandas numpy scipy matplotlib
```

## Quick Start
Tell your AI agent what you want to do:
- "Build a detected/not-detected MRD trajectory from my serial samples with per-sample LoD carried through"
- "Estimate ctDNA clearance half-life, treating undetectable draws as censored not zero"
- "Call molecular relapse only on a confirmed rising trend over consecutive draws"
- "Generate a monitoring report with baseline, nadir, and below-LoD points flagged"

## Example Prompts

### MRD Trajectory and Detection
> "Integrate these serial draws into a binary MRD call per timepoint, keeping the per-sample LoD as a covariate."

> "Mark every undetectable draw as left-censored at its limit of detection rather than plotting it at zero."

### Clearance Kinetics
> "Estimate the clearance half-life from the on-treatment decay phase, excluding the post-op surge and any rebound."

> "Fit log-linear decay only over the uncensored points and report the slope and R-squared."

### Relapse and Response
> "Call molecular relapse only if ctDNA rises above nadir on at least two consecutive confirmed draws."

> "Classify molecular response using my assay's validated cutoff, not a generic 2-log rule."

### Interpretation
> "Flag where lead time is being reported so I do not confuse it with a survival benefit."

> "Check whether any rising signal could be CHIP drift rather than tumor in this tumor-naive panel."

## What the Agent Will Do
1. Assemble serial measurements with timepoint, tumor fraction/VAF, and per-sample LoD (genome-equivalents)
2. Flag below-LoD draws as left-censored and compute baseline/nadir-referenced change without using zero
3. Estimate clearance half-life over the genuine monotonic-decay phase only
4. Call molecular relapse from a confirmed rising trend, annotated with each draw's LoD
5. Surface caveats: non-harmonized response definitions, lead-time bias, CHIP drift, low-shedder negatives

## Tips
- MRD is a binary integrated call across the full variant set, not a per-locus VAF threshold
- Undetectable is left-censored at the per-sample LoD, never a true zero; log(0) breaks decay fits
- Per-sample LoD depends on cfDNA input; a low-input negative is weak evidence and not comparable to a high-input one
- Use the same assay and tube type across timepoints; mixing them confounds the trajectory
- Plot tumor fraction on a log y-axis with a shaded below-LoD band; draw censored points at the LoD, not at zero
- Molecular-response cutoffs ("2-log" vs "90%" vs molecular CR) are non-harmonized; use the assay's own validated definition
- Require a confirmed trend over >=2 consecutive draws before calling relapse; repeated surveillance is a multiple-testing problem
- For tumor-naive panels, control for CHIP drift with matched WBC sequencing before reading a rise as relapse
- Lead time over imaging is analytic sensitivity, not proven benefit; clinical utility needs an interventional design

## Related Skills
- ctdna-mutation-detection - detect the variant set that is then tracked
- tumor-fraction-estimation - per-timepoint tumor burden
- analytical-validation - per-timepoint LoD and left-censoring of undetectable samples
- fragment-analysis - fragmentomic trends as a complementary monitoring signal
- clinical-biostatistics/survival-analysis - relapse, lead-time, and endpoint analysis
