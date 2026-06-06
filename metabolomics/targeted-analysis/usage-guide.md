# Targeted Metabolomics Analysis Usage Guide

## Overview

Targeted metabolomics quantifies a closed, pre-defined panel of known metabolites and reports absolute concentrations with units, using MRM/SRM on a triple quadrupole or PRM on a high-resolution instrument. The central enemy is the matrix effect: co-eluting matrix suppresses analyte ionization, and only a co-eluting stable-isotope-labeled internal standard truly corrects it. This guide covers building a calibration assay, choosing the internal-standard and weighting strategy, confirming identity by ion ratio, and validating to a depth that matches the decision the number supports.

## Prerequisites

```bash
# Skyline (free, vendor-neutral) for transition lists, integration, export
# Download from: https://skyline.ms/
# R for post-export curve fitting and validation metrics
Rscript -e 'install.packages(c("ggplot2", "dplyr"))'
# Optional Python path
pip install pandas numpy scipy
```

Conceptual prerequisites: an authentic reference standard for each analyte (its purity scales every reported number), ideally one stable-isotope-labeled internal standard per analyte, a defined transition list (quantifier + qualifier per analyte with collision energies), and a decision about how much validation the application demands.

## Quick Start

Tell your AI agent what you want to do:
- "Build a 1/x^2-weighted calibration curve and accept it by back-calculated %RE, not R-squared"
- "Normalize each analyte to its co-eluting internal standard and quantify the samples"
- "Confirm identity with the quantifier/qualifier ion ratio and flag interferences"
- "Compute LOD and LLOQ and mark samples below the validated range as not reportable"
- "Help me decide an internal-standard and validation strategy for a clinical assay"

## Example Prompts

### Calibration and Weighting
> "Fit unweighted, 1/x, and 1/x^2 calibration curves and pick the weighting that minimizes low-end back-calculated relative error."
> "Show the per-level %RE for each calibrator and tell me whether the curve passes ICH M10."
> "Set the LLOQ to the lowest calibrator within +/-20% accuracy."

### Quantification
> "Normalize analyte areas to the SIL-IS and back-calculate concentrations from the response-ratio curve."
> "Apply the dilution factor and report concentrations in micromolar, flagging anything below the LLOQ."

### Identity and Quality
> "Compute the qualifier/quantifier ion ratio per sample and flag any outside +/-30% of the calibrator ratio."
> "Check carryover in the blank injected after the top calibrator."
> "Estimate the matrix factor and the IS-normalized matrix factor across matrix lots."

### Strategy and Validation
> "My panel has 40 chemically diverse metabolites and one global IS -- where is my accuracy at risk?"
> "Lay out the ICH M10 parameters I need for a regulated PK assay versus an exploratory study."
> "Should I use a deuterated or 13C internal standard, and what do I verify before trusting it?"

## What the Agent Will Do

1. Establish the panel, transitions, and internal-standard strategy, and pick a validation tier by application.
2. Fit weighted calibration on the analyte/IS response ratio and accept it by per-level back-calculated %RE.
3. Set the LLOQ from accuracy and noise, not from an extrapolated curve.
4. Normalize samples to the internal standard and back-calculate concentrations within the validated range.
5. Confirm identity by ion ratio and flag isobaric interferences.
6. Compute validation metrics (accuracy, precision, matrix factor, recovery, carryover) and export concentrations with quality flags.

## Tips

- Judge a calibration by back-calculated %RE at the low end, never by R-squared alone.
- One SIL-IS per analyte is the gold standard; the further an analyte is from its IS in retention time and chemistry, the larger the uncorrected matrix error.
- Prefer 13C/15N internal standards over deuterium; if forced to deuterium, verify co-elution by overlaying analyte and IS chromatograms.
- Low CV is not evidence of a correct number -- precision and accuracy decouple in shared-IS panels.
- Pre-analytics is upstream of every assay safeguard and invisible to all of them; quench fast and measure stability for labile metabolites.
- A single transition has no defense against isobaric interference -- always carry a qualifier where sensitivity allows.

## Related Skills

- metabolomics/xcms-preprocessing - Upstream feature detection for untargeted discovery before targeted validation
- metabolomics/statistical-analysis - Group comparison and multivariate analysis of quantified concentrations
- metabolomics/isotope-tracing - Stable-isotope tracing and flux (MID), the adjacent discipline this skill hands off to
- metabolomics/normalization-qc - QC-sample-driven drift correction and RSD filtering
- clinical-biostatistics/cdisc-data-handling - Regulated-trial bioanalysis data handling when targeted numbers feed a clinical study
