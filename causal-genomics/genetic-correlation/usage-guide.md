# Genetic Correlation Usage Guide

## Overview

Genetic correlation (rg) is the central cross-trait statistic in causal genomics: a single number summarizing how much two traits share genetic architecture. It is computed from GWAS summary statistics (cross-trait LDSC, HDL, LAVA, rho-HESS, Popcorn) or individual-level genotypes (GREML-bivariate), and underpins multi-trait analyses (GenomicSEM, MTAG), latent-causal inference (LCV gcp), and Mendelian randomization validity assessment (CHP-aware sensitivity is required whenever rg between exposure and outcome is large).

This skill covers global rg, local rg per locus, trans-ancestry rg, and the operational rules for choosing among methods. The two postdoc-level traps are (a) misreading the cross-trait LDSC intercept as a bias on rg (it absorbs sample overlap without contaminating the rg estimate) and (b) running HDL on overlapping cohorts (HDL is biased under any non-trivial overlap).

## Prerequisites

- GWAS summary statistics with `SNP`, `A1`, `A2`, `N`, `BETA` (or `Z`), `SE`, `P`, ideally `EAF`
- Ancestry-matched LD-score reference (EUR: `eur_w_ld_chr/` from alkesgroup; non-EUR: matched references at the same site)
- For HDL: UKB-array SVD eigen99 reference; for LAVA: 1KG-based LD reference + LDetect partition file
- For Popcorn: cross-population LD scores (one-time per ancestry pair)
- Effective N per trait of at least ~50k (LDSC mean chi-square > 1.02)
- For LAVA / HESS: GRCh37 vs GRCh38 build must match between sumstats, LD reference, and partition file

## Quick Start

Tell the AI agent what to estimate:

- "Compute genetic correlation between trait A and trait B from these two sumstats files"
- "Run LDSC rg with the EUR LD-score reference and report the intercept as a sample-overlap diagnostic"
- "Use HDL for maximum precision but only if both GWAS are independent (overlap < 5%)"
- "Run LAVA local genetic correlation across all LDetect loci and report Bonferroni-significant hits"
- "Compute trans-ancestry rg between EUR and EAS T2D with Popcorn"
- "Check whether genetic correlation between exposure and outcome is large enough to require CHP-aware MR sensitivity"
- "Estimate LCV gcp to distinguish partial causation from pure correlation"
- "Reconcile a near-zero global LDSC rg against significant LAVA local hits"

## Example Prompts

### Standard global rg

> "Compute genetic correlation between LDL-cholesterol and coronary artery disease using cross-trait LDSC with the EUR reference panel. Report rg, SE, p-value, and the gcov intercept. Explain whether the intercept indicates sample overlap and whether that affects the rg estimate."

### High-precision rg under verified independence

> "Both GWAS come from non-overlapping cohorts (verified by participant IDs). Use HDL for maximum precision and report rg, SE, and the variance reduction vs LDSC. Confirm overlap is < 5% before running."

### Local rg when global rg is near zero

> "Global LDSC rg between major depressive disorder and rheumatoid arthritis is ~0. Run LAVA across the ~2495 EUR LDetect loci, filter loci with significant univariate local h2 in both traits, and report Bonferroni-significant bivariate local rg. Annotate any hit loci with GWAS catalog and pathway enrichment."

### Trans-ancestry rg

> "Compute genetic correlation between EUR-ancestry and EAS-ancestry type-2 diabetes GWAS using Popcorn. Report rho_ge (effect correlation) and rho_gi (impact correlation). Explain whether rho_ge < 1 reflects population-specific causal architecture or methodological artifact."

### rg as MR validity check

> "Before running Mendelian randomization of BMI on cardiovascular disease, compute LDSC rg. If |rg| > 0.3, add CAUSE (or LHC-MR) to the sensitivity battery because IVW, Egger, and MR-PRESSO are all blind to correlated horizontal pleiotropy. Report rg + intercept alongside the MR results."

### LCV gcp to distinguish causation from correlation

> "Genetic correlation between two traits is 0.6. Run LCV to estimate gcp; interpret gcp = 0 as pure correlation (no partial causation in tested direction) and gcp near 1 as full causation. Use the LCV scripts on LDSC-merged sumstats."

### Reconciliation of LDSC vs HDL

> "I have rg estimates from both LDSC and HDL that disagree. Diagnose by inspecting LDSC intercept and verifying sample overlap. If overlap is non-trivial, report LDSC as primary and exclude HDL; if overlap is verified < 5%, prefer HDL as primary (lower variance)."

### Conditional / partial local rg

> "Run LAVA's run.pcor() to estimate the local partial rg between schizophrenia and bipolar disorder (as the `target` pair) conditional on major depressive disorder (`phenos`), using up to 4 conditioning traits before identification fails."

## What the Agent Will Do

1. Inspect both GWAS sumstats for column completeness (SNP, A1, A2, N, BETA, SE, P, ideally EAF) and effective sample size
2. Verify ancestry matches the LD-score reference panel; warn if non-EUR sumstats are paired with EUR LD scores
3. Munge each sumstats file with `munge_sumstats.py` (LDSC) to a harmonized format restricted to HapMap3 SNPs
4. Compute global rg with cross-trait LDSC; report rg, SE, p-value, intercept, gcov_int, and ratio
5. Diagnose sample overlap from the LDSC intercept; decide whether HDL is appropriate (overlap < 5%) or whether to skip it
6. If HDL is selected, format sumstats as HDL data frames and run `HDL.rg` with `N0 = 0` (or the overlap count)
7. If local rg is requested, run LAVA: load LDetect partitioning, filter to loci with significant univariate local h2 in both traits, run bivariate per locus, apply Bonferroni for ~2495 loci
8. If trans-ancestry rg is requested, run Popcorn with cross-population LD scores; report rho_ge and rho_gi separately
9. If MR is downstream, automatically flag whether |rg| > 0.3 should trigger CHP-aware sensitivity (CAUSE / LHC-MR per pleiotropy-detection skill)
10. Reconcile any disagreement between methods using the operational rules in the SKILL.md reconciliation table
11. Annotate hit local-rg loci with GWAS catalog and tissue / pathway enrichment when requested

## Tips

- **LDSC intercept is not a bias.** A non-zero cross-trait LDSC intercept absorbs sample overlap; the rg estimate (slope) remains unbiased. Do NOT switch to HDL because the LDSC intercept is non-zero; HDL is the method that breaks under overlap, not LDSC.
- **HDL is biased under sample overlap > 5%.** Verify cohort independence by participant identifier, not by file source. Two GWAS from UK Biobank are almost always overlapping.
- **Use ancestry-matched LD scores.** EUR `eur_w_ld_chr/` for EUR GWAS, EAS / AFR for non-EUR; using EUR LD scores on non-EUR sumstats produces inflated intercepts and unreliable rg.
- **Global rg can hide local rg.** Two traits with biologically plausible shared etiology can show global rg ~ 0 due to locus-level cancellation. Run LAVA over the standard ~2495 LDetect loci before concluding "no shared genetic architecture".
- **LAVA univariate filter is mandatory.** Run `run.univ()` first on every locus and filter at `p < 0.05 / N_loci` in BOTH traits before running `run.bivar()`. Skipping this filter produces spurious boundary-case rg values.
- **Popcorn for trans-ancestry, not within-pop LDSC.** Within-population LDSC is invalid for cross-population rg; the LD scores and MAF distributions differ.
- **LCV gcp is genome-wide, not MR.** LCV uses LDSC-style moments on ALL genome-wide SNPs (not the instrument set) and complements but does not replace instrument-based MR.
- **|rg| > 0.3 + significant MR estimate REQUIRES CHP-aware sensitivity.** IVW, Egger, and MR-PRESSO are all blind to correlated horizontal pleiotropy. Add CAUSE (if sig SNPs >= 100) or LHC-MR per the pleiotropy-detection skill.
- **Mean chi-square > 1.02 is the underpower floor.** Below this, no method can rescue rg precision. Meta-analyze first; do not switch methods.
- **Build matching matters for LAVA and HESS.** Sumstats, LD reference, and partition file must all be on the same build (GRCh37 or GRCh38). Cross-build mixing fails silently with `Insufficient SNPs at locus` errors.
- **Original LDSC is Python 2.7.** Prefer the `abdenlab/ldsc-python3` fork (v2.0.0) for a working Python 3 CLI; `belowlab/ldsc` v3.0.1 README states the `--h2/--rg/--h2-cts` CLI is currently broken (use Docker `jtb114/ldsc:latest` for the belowlab fallback). The original `bulik/ldsc` repository has not been updated since 2019.
- **rho_ge and rho_gi can diverge.** Popcorn's effect correlation (causal-effect-size scale) and impact correlation (MAF-weighted) capture different biology; report both.

## Computational Footprint

| Method | Per-pair runtime | Hardware |
|--------|-----------------|----------|
| LDSC rg | ~1 min | laptop |
| HDL.rg | ~10-30 sec | laptop |
| LAVA bivariate (full ~2495 loci) | ~hours | server (parallelizable by chromosome) |
| SUPERGNOVA | ~hours | server |
| Popcorn / Popcorn-2 | ~minutes | laptop |
| GenomicSEM ldsc() (multi-trait) | minutes to ~hour | depends on trait count |
| rho-HESS | ~30-60 min per chromosome | server |
| GREML-bivariate (GCTA) | hours to days | server with large RAM |

## Related Skills

causal-genomics/mendelian-randomization - Primary causal estimation; |rg| > 0.3 motivates CHP-aware sensitivity
causal-genomics/pleiotropy-detection - CAUSE, LHC-MR, LCV; CHP-aware MR battery triggered by high rg
causal-genomics/heritability-partitioning - Partner method; LDSC stack for univariate h2 and partitioned enrichment
causal-genomics/genomic-sem - GenomicSEM ldsc() is the multivariate extension of bivariate LDSC
causal-genomics/colocalization-analysis - Locus-level shared causal variant; complements LAVA hits
causal-genomics/fine-mapping - Credible-set construction at LAVA-significant loci
causal-genomics/mediation-analysis - MVMR for X -> M -> Y after rg motivates causal hypothesis
population-genetics/association-testing - GWAS summary statistics underlying all rg methods
clinical-biostatistics/effect-measures - Translate genetic-architecture findings to clinical effect measures
