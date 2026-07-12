# Differential Splicing - Usage Guide

## Overview
Detect alternative splicing changes between conditions using rMATS-turbo (binomial LRT), leafcutter (Dirichlet-multinomial GLM), MAJIQ V3 deltapsi/HET (Bayesian posterior), SUPPA2 (empirical null), or Shiba (2025 SOTA at low coverage). Tools differ in statistical model, annotation dependence, and calibration regime - these determine which is right for a given experimental design.

## Prerequisites
```bash
# Python tools (rMATS-turbo is not on PyPI; install via bioconda below)
pip install suppa pandas numpy

# R / Bioconductor
BiocManager::install(c('DESeq2'))
# leafcutter (GitHub, not Bioconductor)
# devtools::install_github('davidaknowles/leafcutter/leafcutter')

# CLI tools
conda install -c bioconda rmats regtools
# MAJIQ V3 from majiq.biociphers.org (academic license)
```

## Quick Start
Tell your AI agent what you want to do:
- "Find differential splicing between tumor and normal samples"
- "Run rMATS and leafcutter in parallel and report concordant hits"
- "Use MAJIQ-HET for differential splicing across a heterogeneous patient cohort"
- "Compare splicing between treatment and control with low replicate count"
- "Prioritize differential splicing events by combined statistical and biological significance"

## Example Prompts

### Standard Replicate Designs
> "I have n=3 vs n=3 RNA-seq BAMs; run rMATS-turbo with FDR<0.05 and |dPSI|>0.10, then filter for >=10 junction reads per replicate."

> "Use leafcutter Dirichlet-multinomial GLM on intron clusters from regtools junctions for annotation-free differential splicing."

### Heterogeneous Cohorts
> "I have 30 tumor and 30 normal samples from heterogeneous patients; use MAJIQ V3 HET module with posterior threshold P(|dPSI|>0.2)>0.95."

### Low Replicate Count
> "n=2 vs n=2 design - avoid SUPPA2; use leafcutter or Shiba with junction-imbalance correction."

### Result Prioritization
> "Compute combined score (-log10(FDR) * |dPSI|) and pull the top 50 events with NMD-status and protein-domain annotation."

> "Cross-reference top hits with eCLIP RBP target databases to identify candidate trans-regulators."

## What the Agent Will Do
1. Choose tools based on replicate count, cohort heterogeneity, and annotation availability
2. Configure stat tests (LRT, Dirichlet-multinomial GLM, Bayesian posterior, or empirical null)
3. Apply replicate-aware filtering (per-replicate minimum coverage)
4. Adjust for batch effects and confounders if covariates available
5. Produce ranked results with FDR, ΔPSI, and biological annotation
6. Recommend orthogonal validation (e.g. SpliceAI for variant-driven hits, FRASER2 for outliers)

## Tips
- Run two complementary tools (rMATS + leafcutter) and require concordance for high-confidence calls
- For n=2 vs n=2, avoid SUPPA2 (poorly calibrated null); prefer leafcutter or Shiba
- Document `--b1`/`--b2` order - IncLevelDifference sign matches this
- Check strandedness with RSeQC `infer_experiment.py` before quantification
- Increased PSI of a poison exon decreases functional protein due to NMD - always check ORF
- Report SF3B1 cancer signal as cryptic 3'ss ~10-30nt upstream of canonical
- TDP-43 ALS/FTD signature: cryptic exons in UNC13A, STMN2, ATG4B
- MAJIQ posterior P>0.95 ~ FDR<0.05 in many regimes but not exactly equivalent

## Related Skills

- splicing-quantification - Per-event PSI estimation
- isoform-switching - DTU framework for functional consequences
- sashimi-plots - Visualize significant events
- outlier-splicing-detection - Single-sample-vs-cohort testing for rare disease
- splice-variant-prediction - SpliceAI for variant-driven differential splicing
- long-read-splicing - Differential analysis from full-isoforms
- read-alignment/star-alignment - STAR 2-pass alignment required
