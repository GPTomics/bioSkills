# Genomic SEM Usage Guide

## Overview

Genomic SEM fits structural equation models to GWAS summary statistics, treating each GWAS as an indicator of one or more latent genetic factors. The framework (Grotzinger 2019 Nat Hum Behav 3:513) uses an LDSC-derived genetic covariance matrix (S) and its sampling covariance (V) as inputs to a `lavaan`-based SEM, supporting:

- Common-factor confirmatory models across correlated traits
- User-specified factor structures with arbitrary lavaan syntax
- Exploratory structural equation modeling (ESEM)
- Multivariate GWAS in which a SNP is regressed on a latent factor (`commonfactorGWAS`) with the Q_SNP heterogeneity test
- Stratified GenomicSEM partitioning factor heritability across functional annotations
- Cross-checks against MTAG multi-trait analysis

The skill emphasizes when GenomicSEM is the right multivariate framework, when MTAG is a better fit, and how to read fit indices and Q_SNP to avoid over-claiming "common-factor SNPs" that are in fact trait-specific.

## Prerequisites

- R 4.4+ with `GenomicSEM` (GitHub `GenomicSEM/GenomicSEM`), `lavaan`, `Matrix`, `gdata`
- Munged GWAS summary statistics (HapMap3 SNPs, allele-aligned) -- typically produced by `GenomicSEM::munge()` or LDSC's `munge_sumstats.py`
- LDSC reference panel matching the GWAS ancestry (e.g., `eur_w_ld_chr/` from alkesgroup)
- For multivariate GWAS: SNP-level reference panel (1000G with MAF filter) for SNP harmonization
- For MTAG comparison: Python LDSC fork + MTAG (`JonJala/mtag`). Prefer `abdenlab/ldsc-python3` (v2.0.0) over `belowlab/ldsc` v3.0.1 (whose CLI is broken per the README; Docker `jtb114/ldsc:latest` is the belowlab fallback).
- For stratified GenomicSEM: `baselineLD_v2.2` annotations + per-ancestry frequency files
- Minimum 3 traits to identify a single common factor; more traits stabilize loadings

## Quick Start

Tell an AI agent what to model:

- "Fit a common-factor model to my 5 psychiatric GWAS using GenomicSEM and report CFI and RMSEA"
- "Run a common-factor GWAS across MDD, anxiety, and PTSD; flag Q_SNP-significant SNPs"
- "Compare GenomicSEM common-factor results against MTAG on the same input traits"
- "Fit a two-factor user model: F1 loads on lipid traits, F2 loads on glycemic traits, with correlated factors"
- "Diagnose this Heywood case (negative residual variance) in my GenomicSEM output"
- "Partition the heritability of the common factor across baseline-LD annotations"
- "Check whether the sampling covariance V is positive definite before fitting"

## Example Prompts

### Common-Factor Model

> "Take these three educational-attainment GWAS sumstats, run `GenomicSEM::ldsc()` against the EUR reference, fit a common-factor CFA, and report CFI, RMSEA, SRMR, and standardized loadings. Diagnose any Heywood case if present."

### Common-Factor GWAS with Q_SNP

> "Run a common-factor GWAS across MDD, BIP, and SCZ using `commonfactorGWAS` with `DWLS` estimation. Report SNPs with factor p < 5e-8 AND Q_SNP p > Bonferroni-corrected threshold. Exclude Q_SNP-significant SNPs from the common-factor SNP list."

### Confirmatory Two-Factor Model

> "Fit a two-factor user model where F1 = LDL + HDL + triglycerides and F2 = fasting glucose + HbA1c + 2hr glucose, with F1 ~~ F2 free to estimate factor correlation. Use `usermodel` with DWLS estimation. Report fit indices and factor correlation."

### MTAG Comparison

> "Run MTAG via CLI on the same three lipid GWAS, then compare per-trait top hits with GenomicSEM common-factor SNPs. Check that MTAG `maxFDR` is < 5% for each trait. Where they disagree, classify the SNP as factor-mediated vs trait-specific."

### Stratified GenomicSEM

> "Use `s_ldsc()` to partition the heritability of the latent externalizing factor across baseline-LD annotations + ENCODE/Roadmap cell-type marks. Report per-annotation factor enrichment and the cell type with the largest factor-tau coefficient."

### Sample Overlap Diagnosis

> "These four GWAS are from UK Biobank. Verify the bivariate LDSC intercepts and confirm the V matrix off-diagonals reflect overlap. Fit the common-factor model and check that SEs differ from a naive analysis that ignores V."

### ESEM Exploratory

> "Run an ESEM with two free factors on these 8 cognition + personality GWAS to discover whether a single g factor captures the shared variance or whether g + p_factor emerges. Report rotation: geomin oblique."

## What the Agent Will Do

1. Verify each input GWAS sumstats file is HapMap3-aligned and has columns SNP, A1, A2, BETA/OR, SE, P, N (and EAF when available); call `munge()` if not already munged.
2. Run `ldsc()` to produce the genetic covariance matrix S and sampling covariance V across all input traits; verify mean chi-square per trait > 1.02.
3. Inspect S for pairwise rg near +/-1 (multicollinearity warning); inspect V eigenvalues for positive-definiteness.
4. Fit the requested model:
   - `commonfactor()` for single-factor CFA
   - `usermodel()` with lavaan syntax for arbitrary structure
   - `commonfactorGWAS()` for SNP-level multivariate GWAS with Q_SNP
5. Report model fit indices: CFI, TLI, RMSEA, SRMR, chi-square + df + p-value, AIC, BIC.
6. Inspect standardized loadings; flag any Heywood cases (residual variance < 0 OR standardized loading > 1).
7. For common-factor GWAS, parse factor effect (Est, SE, P) AND Q_SNP (chi-square, df, P) per SNP; apply Bonferroni for heterogeneity threshold.
8. Classify discovered SNPs into "common-factor SNPs" (factor sig + Q_SNP non-sig) vs "heterogeneous SNPs" (factor sig + Q_SNP sig).
9. When MTAG is also requested, build per-trait sumstats input, run `mtag.py`, verify MaxFDR < 5%, and compare top hits with GenomicSEM factor hits.
10. For stratified analysis, run `s_ldsc()` with baseline-LD and target annotations, then `enrich()` for factor enrichment per category.
11. Reconcile any disagreement between GenomicSEM, MTAG, and per-trait univariate GWAS using the reconciliation table in SKILL.md.
12. Produce a manuscript-ready summary: fit indices, factor diagram, SNP table with Q_SNP, and a methodological statement covering sample overlap handling via V.

## Tips

- **Q_SNP is mandatory.** Reporting common-factor SNPs without Q_SNP heterogeneity testing is the single most common reviewer-flagged error. Always include the Q_SNP p-value column.
- **Sample overlap.** Whenever any two input GWAS share even partial sample overlap (e.g., UK Biobank phenotypes), do not bypass V. The full `ldsc()` output -- both S and V -- is the input GenomicSEM expects.
- **Heywood cases require diagnosis.** Negative residual variance is a model-data conflict, not a fitting nuisance. Identify the offending trait (usually highly correlated with another), constrain the variance to be non-negative, or drop the trait. Document the choice.
- **MaxFDR is MTAG's Q_SNP.** When running MTAG alongside GenomicSEM, treat MaxFDR > 5% the same way as significant Q_SNP -- a flag that the assumption is violated for that trait.
- **DWLS vs ML.** Use `estimation = 'DWLS'` unless the V matrix is essentially diagonal. ML is faster but produces wrong SE when sample overlap exists (i.e., almost always in practice).
- **Identification.** A single factor needs >= 3 traits to be identified. A two-factor model needs >= 6 traits or 3+ per factor without cross-loadings, plus a constraint (typically fixing factor variance to 1 OR fixing one loading per factor to 1).
- **Cross-ancestry.** As of 2026, GenomicSEM has no validated cross-ancestry pipeline. Run per-ancestry, compare loadings qualitatively, and avoid pooled cross-ancestry V matrices.
- **CFI vs chi-square.** At large N the chi-square test almost always rejects exact fit. Use CFI >= 0.95 and RMSEA <= 0.05 as the primary fit anchor.
- **Power for ESEM.** Exploratory analysis with rotation needs >= 6 traits and stable loadings across rotation choices (geomin, target, quartimin); report sensitivity.
- **Downstream MR.** The factor sumstats from `commonfactorGWAS` (excluding Q_SNP-significant SNPs) is a valid MR exposure for testing the causal effect of the latent factor on a separate outcome. See causal-genomics/mendelian-randomization.

## Related Skills

causal-genomics/mendelian-randomization - Run MR with the common-factor sumstats as exposure
causal-genomics/heritability-partitioning - sLDSC foundations + LDAK comparison; required upstream for stratified GenomicSEM
causal-genomics/colocalization-analysis - Resolve overlap of factor-significant SNPs with eQTL/pQTL signals
causal-genomics/pleiotropy-detection - Q_SNP is the per-SNP pleiotropy diagnostic in GenomicSEM
causal-genomics/fine-mapping - Construct credible sets at factor-significant loci
causal-genomics/mediation-analysis - SEM mediation overlaps with `usermodel` indirect path coefficients
causal-genomics/transcriptome-wide-association - TWAS on factor GWAS output
population-genetics/association-testing - GWAS workflow upstream of GenomicSEM
