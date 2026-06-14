# Rare-Variant Association - Usage Guide

## Overview

Rare variants have too few carriers for single-variant GWAS to detect, so the unit of analysis becomes a gene or region and the question becomes which aggregation test to use. That choice is a bet about the unobserved effect architecture: a burden test assumes every variant in the mask acts in the same direction and is most powerful when that holds but cancels to a null under mixed directions, SKAT is a variance-component test robust to mixed directions but weaker when the truth is unidirectional, and SKAT-O optimally blends the two. The mask - which variants enter, by functional class and MAF cutoff - is the hypothesis, not a preprocessing detail, and the null model (SPA/Firth for imbalance, an LMM for relatedness) is what keeps the aggregate calibrated.

## Prerequisites

- regenie (`conda install -c bioconda regenie`) and/or SAIGE (`conda install -c bioconda saige`) for biobank-scale set tests
- The SKAT and STAAR R packages (`install.packages('SKAT')`; STAAR from `xihaoli/STAAR` via remotes) for direct/scriptable analysis
- QC'd genotypes (sequenced calls preferred for rare-variant masks) and per-variant functional annotations (LoF, missense, CADD, regulatory)
- Covariates (age, sex, principal components) and, for related samples, the inputs to fit an LMM null
- Conceptual prerequisites and big notes:
  - The mask (functional class x MAF cutoff) IS the hypothesis; report it alongside any p-value and test a small grid of masks rather than one default.
  - Burden assumes one effect direction; SKAT does not; SKAT-O picks between them. Pick from the architecture prior, not habit.
  - Case/control imbalance and low MAC make naive set tests anti-conservative; use SPA (SAIGE-GENE+) or Firth (regenie).
  - Population structure and relatedness still need an LMM null; PCs alone cannot remove relatedness.
  - The multiple-testing burden is per-gene-per-mask (exome-wide ~2.5e-6, tighter with multiple masks); combine masks via ACAT-O to limit it.
  - Imputed or low-callrate variants corrupt masks; filter by INFO/R2 and genotype quality first.

## Quick Start

Tell your AI agent what you want to do:
- "Run gene-based burden, SKAT, and SKAT-O tests on my exome data"
- "Build LoF-only and LoF+missense masks at MAF 0.001 and 0.01 and test each gene"
- "Use SAIGE-GENE+ for an imbalanced binary trait in a related cohort"
- "Run regenie gene tests with the SKAT-O and ACAT-O omnibus and Firth fallback"
- "Decide whether burden or SKAT fits my gene's effect architecture"
- "Set an exome-wide significance threshold for my gene tests"

## Example Prompts

### Test selection
> "I have a gene with both predicted gain-of-function and loss-of-function rare variants - which aggregation test should I run, and why not a plain burden test?"

### Biobank pipeline
> "Run regenie step 2 gene-based tests on my whole-exome data with LoF and LoF+missense masks at AAF 0.001 and 0.01, request SKAT-O and ACAT-O, and turn on Firth for the imbalanced binary trait."

### Imbalance and relatedness
> "My case:control ratio is 1:50 and the sample has cryptic relatedness - set up SAIGE-GENE+ to test each gene across MAF cutoffs 0.0001, 0.001, and 0.01."

### Direct R analysis
> "For a small case-control cohort, fit a SKAT null model on covariates and report burden, SKAT, and SKAT-O p-values for one gene with Beta(1,25) MAF weighting."

### Masks and thresholds
> "Explain how to build the annotation, set-list, and mask-definition files for regenie and what exome-wide significance threshold to use across multiple masks."

## What the Agent Will Do

1. Confirm the scope is gene/region aggregation, not single-variant GWAS (which routes to association-testing), and identify the trait type and sample structure.
2. Choose the engine: regenie or SAIGE-GENE+ for biobank scale (SAIGE-GENE+ when imbalance/relatedness dominate), or the SKAT R package for small scriptable cohorts.
3. Build the variant masks from functional annotations and MAF cutoffs, treating the mask as the hypothesis and testing a small grid rather than one default.
4. Select the test from the architecture prior: burden for a single-direction mask, SKAT for mixed directions, SKAT-O when unknown, ACAT-O/STAAR for omnibus or annotation-weighted analysis.
5. Calibrate the null: SPA/Firth for imbalance, an LMM null (LOCO predictor or sparse GRM) for relatedness, and QC out imputed/low-quality variants first.
6. Set the per-gene-per-mask significance threshold (exome-wide ~2.5e-6, tighter with multiple masks or combined via ACAT-O) and report the mask with each result.

## Tips

- Validate mask inputs before running: regenie `--check-burden-files` flags set-list variants missing from the annotation file, a silent cause of empty masks.
- Pass nested MAF cutoffs in one run (regenie `--aaf-bins 0.0001,0.001,0.01`, SAIGE `--maxMAF_in_groupTest 0.0001,0.001,0.01`) rather than rerunning per cutoff.
- A burden hit with a null SKAT (or vice versa) is diagnostic of the architecture, not a contradiction - report both and let SKAT-O/ACAT-O adjudicate.
- ACAT-O is dominated by its smallest input p-value by design, so a single artifactual variant can drive a gene; QC the inputs (INFO/R2, genotype quality) before trusting a hit.
- Reuse the step-1 null across single-variant and gene tests (regenie, SAIGE) so the relatedness/structure model is identical for both.
- For genome-wide gene scans in R, use SSD files (`Generate_SSD_SetID`, `Open_SSD`, `SKAT.SSD.All`) instead of holding every gene's matrix in memory.

## Related Skills

- association-testing - single-variant GWAS (linear/logistic/LMM/SPA per marker) that this skill aggregates beyond
- plink-basics - genotype QC and format conversion before masking
- population-structure - PCs and relatedness for the null model that calibrates the set test
- variant-calling/variant-annotation - functional annotations (LoF, missense, CADD, regulatory) that define masks
- clinical-databases/variant-prioritization - clinical interpretation of variants flagged by gene tests
- causal-genomics/fine-mapping - resolving which variants in a significant gene/region carry the signal
