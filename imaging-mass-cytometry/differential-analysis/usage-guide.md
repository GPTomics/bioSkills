# IMC Differential Analysis - Usage Guide

## Overview

Compares cell-type composition and spatial features across conditions in an IMC/MIBI cohort. The one rule that shapes every test here: the replicate is the patient, not the cell or the image -- testing at the cell level over millions of correlated cells manufactures significance, so every differential question follows one spine (per-image summary, aggregate to patient, test across patients), with cell-type proportions treated as compositional and acquisition batch as a covariate.

## Prerequisites

```bash
# Python
pip install statsmodels scanpy sccoda pandas numpy

# R / Bioconductor (diffcyt; mixed models; SpaceANOVA)
# BiocManager::install('diffcyt'); install.packages(c('lme4', 'SpaceANOVA'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Test whether Tregs differ between responders and non-responders, at the patient level"
- "Compare cell-type composition with a method that handles the sum-to-one constraint"
- "Test whether CD8-tumor contact differs between groups, aggregated to patients"
- "Set up a mixed model with patient as a random effect and batch as a covariate"
- "Tell me my real sample size given 4 patients and 2 million cells"

## Example Prompts

### Differential abundance
> "I have per-cell type labels for 30 patients in two groups. Compare cell-type proportions, accounting for multiple ROIs per patient and the fact that proportions sum to one."

### Differential spatial feature
> "I computed per-image neighborhood enrichment z-scores. Test whether the CD8-tumor interaction differs between conditions, aggregating to the patient and correcting across all cell-type pairs."

### Power and unit
> "My collaborator ran a Wilcoxon over 500,000 cells and got p < 1e-50. Why is that wrong, and what is the correct test?"

### Compositional
> "Several of my cell types shift between groups. Use scCODA so I don't get spurious opposite-direction changes."

## What the Agent Will Do

1. Aggregate the cell table to per-image summaries, then to one value per patient.
2. Choose the test from the decision tree (mixed model for nested ROIs; scCODA for compositional abundance; diffcyt-DA for cluster counts; pseudobulk for within-type marker expression).
3. Enter batch as a covariate and check it is not confounded with condition.
4. For spatial differences, treat the per-image spatial statistic as the summary and test it at the patient unit with FDR across pairs.
5. Report the patient sample size and honest power, never rescuing significance with cell count.

## Tips

- Cell count is not replication; the effective n is the number of patients.
- Proportions are compositional -- a real rise in one type mechanically depresses others; use scCODA or a CLR transform with a reference type.
- Multiple ROIs per patient are nested, not independent; use a mixed model with patient random effect.
- Over-correcting batch can erase the disease signal; correct minimally and keep integration out of the inference path.
- Randomize acquisition order against condition so run drift does not confound the comparison.
- Apply BH-FDR across all cell-type pairs and radii tested.

## Related Skills

- phenotyping - supplies the cell-type labels whose proportions are compared
- spatial-analysis - supplies the per-image spatial statistics that become patient-level summaries
- quality-metrics - batch must be diagnosed and entered as a covariate
- experimental-design/randomization-blocking - the experimental-unit and pseudoreplication foundation
- clinical-biostatistics/subgroup-analysis - multiplicity and effect estimation in clinical cohorts
- flow-cytometry/differential-analysis - diffcyt-DA/DS for suspension cytometry
