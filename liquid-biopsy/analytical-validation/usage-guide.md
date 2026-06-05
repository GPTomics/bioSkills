# Analytical Validation and Detection Limits - Usage Guide

## Overview
Establishes and audits the sensitivity claim of a liquid-biopsy assay using measurement science: the Poisson genome-equivalent sampling ceiling, the error-suppression floor, and the CLSI EP17 LoB/LoD/LoD95/LoQ framework, including the per-locus-vs-panel-integrated distinction and honest LoD reporting conditioned on input mass.

## Prerequisites
```bash
pip install numpy scipy statsmodels
```

## Quick Start
Tell your AI agent what you want to do:
- "How many genome equivalents do I need to detect a 0.1% variant 95% of the time?"
- "Estimate the LoD95 from my dilution series of detection calls"
- "Compute the panel-integrated LoD for a 30-variant bespoke MRD assay"
- "Audit this 'detects 0.1% VAF' sensitivity claim for missing input mass"

## Example Prompts
### Detection Limits
> "Given 25 ng of cfDNA input, what is the probability of detecting a true 0.05% VAF variant at a single locus?"
> "What input mass puts lambda at 3 for a 1e-4 variant, and why is that the 95% sampling threshold?"

### Validation Study
> "Fit a probit detection curve to my dilution series and report LoD95 with the LoB from my blank replicates"
> "Distinguish the LoD from the LoQ for this assay and tell me which VAFs I can report as quantitative"

### Panel Design
> "Compare the per-locus LoD to the panel-integrated LoD for tracking 16 vs 48 clonal variants"
> "Explain why a bespoke panel reaches ppm sensitivity when each locus alone floors at 1e-4"

### Reporting and Audit
> "Rewrite this LoD claim to condition on input genome equivalents, consensus depth, and replicate detection rate"
> "Flag whether this sensitivity number is sampling-limited or error-limited at the stated VAF"

## What the Agent Will Do
1. Convert input mass to haploid genome equivalents (~330/ng) and compute lambda = input_GE x VAF.
2. Return the Poisson detection probability and the minimum GE needed for lambda >= 3 (95% sampling-detection).
3. Identify whether the achieved LoD is sampling-limited or error-limited at the target VAF.
4. Estimate LoB from blank replicates and fit a probit/logistic detection curve for LoD95 from a dilution series.
5. Combine per-locus detection probabilities into the panel-integrated LoD under a >=k-of-N positivity rule.
6. Reframe any bare-VAF sensitivity claim into an LoD conditioned on input mass, consensus depth, and replicate detection rate.

## Tips
- **Anchor to molecules** - Always pair a VAF with input genome equivalents; lambda = GE x VAF, not the VAF, determines detection.
- **Stop buying depth** - Past the deduplication plateau the assay is sampling-saturated; add plasma volume or conversion efficiency instead.
- **Separate detect from quantify** - LoD is binary detection; LoQ (CV<=20%) is quantitation. Report near-floor VAFs as detected/not-detected, not as measured values.
- **Mind the integration level** - A per-locus LoD and a panel-integrated LoD differ by orders of magnitude; never quote one for the other.
- **Duplex for the deepest claims** - Single-strand UMI consensus does not remove template-resident deamination/oxidation damage; sub-1e-5 claims need duplex strand-concordance.
- **Use commutable standards** - Validate against contrived SEQC2 Sample A / HCC1395 admixtures fragmented to ~160 bp, and treat commutability with real plasma as the open caveat.

## Related Skills
- ctdna-mutation-detection - applies these limits to low-VAF somatic calls
- longitudinal-monitoring - per-timepoint LoD and left-censoring of undetectable samples
- tumor-fraction-estimation - the ~3% CNA-based detection floor as an LoD
- experimental-design/multiple-testing - repeated-surveillance specificity and FDR
- clinical-biostatistics/power-and-sample-size - validation-study design
