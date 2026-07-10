# GRN Inference and TF Activity - Usage Guide

## Overview

Infer gene regulatory networks from bulk (or general) expression data and read transcription-factor protein activity from the resulting regulons. Edge inference uses mutual information (ARACNe-AP) or tree ensembles (GENIE3, GRNBoost2); activity inference uses VIPER/msVIPER. The central claim is the Califano-lab "activity, not edges" paradigm: an expression-derived GRN is an undirected association graph whose individual edges are unreliable, but a regulon read as a multiplexed reporter gives a robust estimate of a TF's protein activity -- which can be high even when the TF's own mRNA is unchanged. DREAM5 showed no single inference method dominates (ensembles win) and that synthetic-benchmark accuracy does not transfer to real eukaryotic data. This is the bulk counterpart to the single-cell scenic-regulons skill.

## Prerequisites

```r
BiocManager::install(c('viper', 'GENIE3'))
```

```bash
# ARACNe-AP is a Java tool built from source (https://github.com/califano-lab/ARACNe-AP)
# Python tree-ensemble alternative:
pip install arboreto
```

A candidate TF list matched to the expression matrix's gene namespace is required for edge inference.

## Quick Start

Tell your AI agent what you want to do:
- "Infer a TF-target network from my bulk RNA-seq with GENIE3"
- "Build an ARACNe network and find master regulators with VIPER"
- "Score per-sample TF activity for patient stratification"
- "Which transcription factors are the master regulators of my tumor vs normal contrast?"

## Example Prompts

### Edge Inference
> "I have a bulk RNA-seq matrix and a TF list. Reverse-engineer a regulatory network with GENIE3."

> "Run ARACNe-AP with bootstrapping and consolidation to build a mutual-information network."

### Master Regulators and Activity
> "Build a regulon from my ARACNe network and run msVIPER to rank master regulators of disease."

> "Compute a per-sample VIPER activity matrix and cluster patients by TF activity."

### Robustness
> "Infer the network with multiple methods and combine them into a consensus (wisdom of crowds)."

## What the Agent Will Do

1. Orient the expression matrix correctly (genes in rows for GENIE3) and load the TF list
2. Run edge inference (GENIE3 random forests, or ARACNe-AP bootstraps + consolidation)
3. Assemble a regulon, assigning each target a Mode-of-Regulation sign
4. Build a null model and run msVIPER (master regulators) or VIPER (per-sample activity)
5. Report activity by NES, not node degree, and evaluate edges with AUPRC where a gold standard exists

## Tips

- Activity beats edges - individual inferred edges are noisy; the regulon-as-reporter VIPER activity is the robust, interpretable output.
- Direction is assumed, not inferred - MI/correlation are symmetric; treat edges as associations unless perturbation/time/sequence supports direction.
- Mode of Regulation is required - build regulons with `aracne2regulon` so target signs are assigned, or VIPER degrades to plain enrichment.
- ARACNe needs bootstrapping + consolidation - a single run is unstable; ~100 bootstraps then `--consolidate`.
- GENIE3 orientation - genes in rows, samples in columns (the transpose of WGCNA); set a seed for reproducibility.
- Evaluate with AUPRC - AUROC near 1 can hide near-random AUPRC under sparse true edges; gold standards are incomplete, so absent edges are not true negatives.
- No single best method - ensemble multiple inferences for robustness (DREAM5); be skeptical of methods validated only on synthetic data.

## Related Skills

- scenic-regulons - single-cell regulons with motif-pruning directness
- coexpression-networks - undirected co-expression modules (no TF privileging)
- differential-networks - VIPER differential activity / rewiring between conditions
- multiomics-grn - enhancer-driven directed GRNs from accessibility
- differential-expression/de-results - signatures that feed msVIPER
