# Epigenetic Clocks (DNAm Age) - Usage Guide

## Overview

This skill computes DNA methylation age (DNAm age) and pace of aging by applying frozen, pre-trained epigenetic clocks to a clean beta or M-value matrix. It covers choosing a clock by the question being asked, computing age acceleration (the real endpoint), handling the per-CpG reliability crisis with principal-component (PC) clocks, and dealing with clock-CpG dropout on EPICv2.

The central idea: a clock is a frozen elastic-net predictor, not a mechanism. Its CpGs are prediction features (do not GO-enrich them), the endpoint is age ACCELERATION rather than raw age, and on noisy first-gen clocks the same sample can age several years between technical replicates. This skill applies clocks; it does not train them.

## Prerequisites

Install the clock packages:

```r
BiocManager::install('methylclock')
# methylclockData is pulled in as a dependency for the coefficient sets
remotes::install_github('yiluyucheng/dnaMethyAge')        # author-year clock IDs, EAA in one call
remotes::install_github('MorganLevineLab/methylCIPHER')   # PC clocks for longitudinal/trial use
remotes::install_github('danbelsky/DunedinPACE')          # pace of aging
```

Conceptual prerequisites:
- A clean beta or M-value matrix with CpGs as rows and samples as columns, already normalized and QC'd (this skill assumes the matrix is ready).
- For EPICv2 data, replicate probes (multiple beads per CpG, suffixed names) must be collapsed to one value per CpG before any clock call.
- A phenotype table with chronological age (and sex for GrimAge variants) is needed to compute age acceleration.
- Clock coefficient sets are fixed and version-pinned; the package version mostly controls which clocks ship and what the clock-name strings are.

## Quick Start

Tell your AI agent what you want to do:
- "Compute Horvath and Hannum DNAm age for my beta matrix and report age acceleration"
- "Estimate DunedinPACE pace of aging for these blood samples"
- "Which clock should I use for a mortality endpoint?"
- "Check how many clock CpGs survive on my EPICv2 data before estimating ages"
- "Use a PC clock so my longitudinal trial estimates are reliable"

## Example Prompts

### Choosing and applying a clock
> "I have an EPIC beta matrix for 200 blood samples with chronological ages. Compute Horvath, Hannum, and PhenoAge DNAm age, then give me the age acceleration (residual on chronological age) for each clock."

### Pace of aging
> "Run DunedinPACE on my blood methylation data and report the pace value per sample. Keep it on its own scale, do not convert it to an age."

### Reliability for a trial
> "This is a longitudinal intervention study with two timepoints per person. Use PC clocks so the measurement noise does not swamp the intervention effect, and document the reliability."

### EPICv2 coverage check
> "Before estimating any clock on my EPICv2 data, tell me what fraction of each clock's CpGs is actually present, and warn me about clocks that lose more than 10% of their CpGs."

### Choosing by question
> "I want a clock that predicts mortality and healthspan, not just chronological age. Recommend one and explain why it differs from Horvath."

## What the Agent Will Do

1. Confirm the input is a clean beta/M-value matrix and, for EPICv2, that replicate probes have been collapsed.
2. Select clocks by the question: first-gen for chronological age, second-gen (PhenoAge/GrimAge) for health/mortality, DunedinPACE for pace, tissue-matched clocks for pediatric/gestational, mitotic clocks for cancer-risk questions.
3. Check clock-CpG coverage with `checkClocks` and report the fraction present per clock, refusing or flagging high-missingness samples.
4. Apply the clocks with `DNAmAge` / `methyAge` / `PACEProjector`, supplying chronological age so age-acceleration columns are produced.
5. Return age acceleration (the residual endpoint), not raw DNAm age, and keep DunedinPACE on its own rate scale.
6. Recommend PC clocks (methylCIPHER) for any longitudinal or interventional design and hand off cell-count adjustment (IEAA), survival modeling, and predictor training to the related skills.

## Tips

- Report age ACCELERATION, not raw DNAm age. Raw age is dominated by chronological age (r often > 0.9); the residual is the signal. DunedinPACE is the exception (it is already a rate).
- Never GO-enrich clock CpGs. They are penalty-selected prediction features, one arbitrary representative per correlated cluster, not an aging pathway.
- Two clocks for the same outcome can share almost no CpGs. That is expected from elastic-net selection, not a contradiction.
- Always report how many of each clock's CpGs were present. A large imputed fraction biases age acceleration toward zero. On EPICv2, GrimAge, Hannum, and DunedinPACE lose more than 10% of their CpGs, and Hannum can return negative ages.
- For repeated-measures designs, use PC clocks (methylCIPHER) or at minimum document an ICC/reliability assessment; first-gen per-CpG noise can exceed a small intervention effect.
- Match the clock to tissue and age: skin&blood for fibroblasts/skin, PedBE for pediatric buccal, Knight/Bohlin for newborn cord blood. A blood clock on another tissue needs a known-age calibration check.
- For IEAA, deconvolve first (cell-type-deconvolution) and residualize the clock on the estimated cell counts.

## Related Skills

- array-preprocessing - Provides the clean beta matrix clocks consume
- cell-type-deconvolution - IEAA: adjust age acceleration for cell composition
- machine-learning/model-validation - Predictor training, cross-validation, leakage (clocks apply frozen models)
- clinical-biostatistics/survival-analysis - Survival/mortality modeling of age acceleration
- ewas-design - Age-acceleration association study design
- workflows/methylation-pipeline - End-to-end methylation pipeline
