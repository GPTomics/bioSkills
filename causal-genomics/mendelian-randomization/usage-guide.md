# Mendelian Randomization - Usage Guide

## Overview

Mendelian randomization (MR) estimates the causal effect of an exposure on an outcome using genetic variants as instrumental variables (IVs). Because alleles are randomized at conception and fixed across the lifecourse, well-chosen IVs are independent of confounders that bias observational associations, allowing causal inference from GWAS summary statistics rather than randomized trials.

This skill covers the full sensitivity battery (IVW, MR-Egger, weighted median/mode, MR-RAPS, MR-PRESSO, CAUSE, MVMR, MR-Clust, LCV, LHC-MR), the bias structure of one-sample vs two-sample designs, drug-target / cis-MR for pQTL/eQTL exposures, and STROBE-MR-compliant reporting. The agent selects the primary method based on the experimental regime, runs a concordant sensitivity battery, and flags weak-IV / pleiotropy / sample-overlap issues automatically.

## Prerequisites

```r
install.packages(c('remotes', 'MendelianRandomization', 'MVMR', 'coloc', 'simex'))
remotes::install_github('MRCIEU/TwoSampleMR')        # 0.6.0+ for JWT auth
remotes::install_github('MRCIEU/ieugwasr')           # 1.0+ for JWT auth
remotes::install_github('rondolab/MR-PRESSO')
remotes::install_github('qingyuanzhao/mr.raps')      # CRAN-archived 2025-03-01
remotes::install_github('jean997/cause')
remotes::install_github('cnfoley/mrclust')
remotes::install_github('LizaDarrous/lhcMR')
remotes::install_github('n-mounier/MRlap')           # overlap + winner's curse + weak-IV joint correction
remotes::install_github('HDTian/DRMR')               # non-linear MR (with Hamilton 2023 caveat)
```

OpenGWAS API: OAuth was deprecated May 2024; JWT token required. Generate at api.opengwas.io, set `OPENGWAS_JWT=<token>` in `~/.Renviron`, restart R, verify with `ieugwasr::get_opengwas_jwt()`. For production, skip the API entirely: use `ld_clump_local()` with local 1KG bfile.

For local LD clumping (preferred since OpenGWAS rate limits tightened in 2024):

```bash
# plink2 binary + 1000 Genomes EUR reference panel
# https://www.cog-genomics.org/plink/2.0/
# 1KG bfile prebuilt at https://mrcieu.github.io/ieugwasr/
```

## Quick Start

Tell the AI agent what to do:
- "Run a two-sample MR of LDL cholesterol on coronary heart disease from local GWAS files"
- "Test whether circulating IL-6R protein causally affects rheumatoid arthritis using cis-pQTL instruments"
- "Run multivariable MR of BMI on T2D adjusting for waist-hip-ratio"
- "Estimate the indirect effect of BMI on CHD mediated through LDL cholesterol"
- "Distinguish a causal BMI to CHD effect from correlated horizontal pleiotropy using CAUSE"
- "Produce a STROBE-MR-compliant sensitivity battery for my exposure-outcome pair"

## Example Prompts

### Standard Two-Sample MR

> "I have GWAS summary statistics for systolic blood pressure (exposure) and stroke (outcome) from independent cohorts. Run TwoSampleMR with IVW primary, plus MR-Egger, weighted median, weighted mode, MR-PRESSO, and Steiger directionality. Report effect sizes with 95% CIs and flag any pleiotropy concerns."

> "Test the causal effect of educational attainment on Alzheimer's disease using GWAS summary stats. Use local plink clumping and apply Steiger filtering."

### Drug-Target / cis-MR

> "Run a cis-MR analysis of PCSK9 inhibition on coronary artery disease using cis-pQTL instruments from UKB-PPP within 500 kb of the PCSK9 gene. Cross-validate with coloc PP.H4."

> "I want to predict the effect of an IL-6 receptor antagonist on rheumatoid arthritis. Use cis-pQTLs for IL6R from deCODE plus colocalization."

### MVMR and Mediation

> "Run multivariable MR of LDL, HDL, and triglycerides on CHD jointly. Report conditional F for each lipid trait. Use MVMR's `strength_mvmr()`."

> "Test whether BMI's effect on T2D is mediated through fasting glucose. Implement the two-step MR and MVMR difference methods."

### Non-Linear MR

> "Test for a J-shaped relationship between alcohol consumption and all-cause mortality using doubly-ranked stratification (Tian 2023)."

### Polygenic Exposure with CHP Concern

> "Run CAUSE on BMI -> coronary heart disease to separate a causal effect from correlated horizontal pleiotropy. Compare against IVW and MR-PRESSO."

### One-Sample / Sample-Overlap Concern

> "Both my exposure and outcome GWAS are from UK Biobank. Apply MRlap for joint correction of sample overlap, winner's curse, and weak-IV bias (Mounier 2023); fall back to MR-RAPS plus Burgess 2016 correction if LDSC h^2 < 0.05."

### Binary Outcomes

> "Run MR of LDL on T2D and explicitly handle the non-collapsibility of the logistic OR; report the per-SD genetically-predicted exposure effect on the marginal log-OR scale (Burgess & Labrecque 2018)."

## What the Agent Will Do

1. Diagnose the experimental regime: one-sample, two-sample independent, partial overlap, drug-target, polygenic, or non-linear
2. Read exposure GWAS, select genome-wide significant SNPs (P < 5e-8), compute F-statistics from exposure (not outcome)
3. Clump instruments locally (plink + 1KG bfile preferred); cis-MR uses r2 < 0.1 within window, polygenic MR uses r2 < 0.001 across 10 Mb
4. Read outcome GWAS, extract matched SNPs
5. Harmonise alleles with `action = 2` (infer from EAF; drops MAF~0.5 palindromes) and document drops
6. Run IVW (random effects) as primary, plus MR-Egger (with `I^2_GX` check + SIMEX if needed), weighted median, weighted mode, MR-RAPS
7. Run MR-PRESSO with NbDistribution >= 10000 for global, outlier, and distortion tests
8. Run CAUSE if exposure is polygenic and >= 100 sig SNPs; report ELPD difference
9. For multivariable designs, run MVMR with `strength_mvmr()` per-exposure conditional F
10. Steiger directionality test (with the Hemani Tilling 2022 confounder caveat noted)
11. Generate scatter, forest, leave-one-out, and funnel plots
12. Produce STROBE-MR-compliant report with all 20 items and explicit causal-claim language

## Tips

- Compute F-statistic from EXPOSURE GWAS, not outcome; using outcome F is a common error that produces meaningless instrument-strength claims
- MR-RAPS CRAN-archived 2025-03-01; install the GitHub source `qingyuanzhao/mr.raps`. TwoSampleMR's `mr_raps()` wrapper is a thin caller into that package; the `MendelianRandomization` package does NOT export `mr_raps()`
- OpenGWAS authentication has tightened repeatedly since 2024; use local plink + 1KG bfile via `ieugwasr::ld_clump()` for production work
- When both ends use UK Biobank, treat as one-sample-equivalent: weak-IV bias is TOWARD the observational estimate, not the null
- Palindromic SNP handling: `action = 2` drops MAF~0.5 palindromes (default); `action = 3` drops all palindromes (most conservative); never `action = 1` unless strand convention is guaranteed
- MR-Egger requires `I^2_GX >= 0.9` for NOME; below that apply SIMEX correction or downgrade Egger to exploratory
- CAUSE needs >= 100 significant SNPs; for sparse exposures use IVW + Egger + MR-PRESSO + LCV instead
- MVMR: report CONDITIONAL F per exposure, not total F; total F can be high while one exposure's conditional F is below 10
- CAUSE explicitly models correlated horizontal pleiotropy; use it when prior on CHP is high (BMI / lipids / smoking on cardiometabolic outcomes)
- Bonferroni-correct when running pheWAS-MR across many outcomes; FDR only for explicitly exploratory screens
- Report all 20 STROBE-MR items (Skrivankova 2021); missing items trigger first-pass desk rejection at most epi/genetics journals since 2022
- The Steiger filter has a known failure mode under unmeasured confounding (Hemani Tilling 2022 Wellcome Open Res 7:14); cross-validate direction with bidirectional MR or LHC-MR, not Steiger alone
- MRlap (Mounier 2023) jointly corrects sample overlap, winner's curse, and weak-IV bias from sumstats; prefer when (a) any sample overlap suspected, (b) only sumstats available, (c) UKB-on-UKB or FinnGen-on-FinnGen designs
- For non-linear MR (J-curves, U-shapes): Hamilton 2023 medRxiv 23293658 shows BOTH doubly-ranked and residual stratification produce stratum-specific bias from age/sex effects; pre-specify the non-linear hypothesis, run both methods side-by-side, report negative-control outcomes (genotype vs sex within strata)
- One-sample F-statistic floor: F >= 20 (not 10) per Bowden 2019 IJE 48:728; jackknife SE preferred over analytic; never run exposure and outcome GWAS on the same individuals and claim two-sample (Hartwig 2021 collider bias)
- Binary outcomes: MR returns the population-averaged log-OR, NOT the conditional log-OR; for common diseases OR diverges from RR/HR (Burgess & Labrecque 2018); phrase as "per 1-SD increase in genetically-predicted X, OR for Y = ..."
- MVMR with conditional F < 10: switch from `ivw_mvmr()` to `qhet_mvmr(r_input, pcor, CI = TRUE, iterations = 1000)` (Sanderson 2021 Stat Med 40:5434) which Q-minimizes rather than weighting by inverse variance

## Related Skills

causal-genomics/colocalization-analysis - Required for cis-MR drug-target work to confirm a shared causal variant
causal-genomics/pleiotropy-detection - Deep MR-PRESSO / Egger / contamination-mixture diagnostics
causal-genomics/fine-mapping - Credible-set construction at instrument loci
causal-genomics/mediation-analysis - Two-step MR and MVMR difference method for X -> M -> Y
causal-genomics/transcriptome-wide-association - TWAS as MR-adjacent framework for gene-level inference
causal-genomics/proteome-mr-drug-target - Cis-pQTL drug-target MR workflow with UKB-PPP/deCODE/Fenland and coloc triangulation
population-genetics/association-testing - GWAS source for instrument selection
clinical-databases/clinvar-lookup - Annotate instrument SNPs for downstream interpretation
