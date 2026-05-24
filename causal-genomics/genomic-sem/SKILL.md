---
name: bio-causal-genomics-genomic-sem
description: Fits structural equation models to GWAS summary statistics using GenomicSEM (Grotzinger 2019), including common-factor models, confirmatory factor models, ESEM, common-factor GWAS with Q_SNP heterogeneity, multivariate Wald tests, and stratified GenomicSEM partitioned heritability. Reconciles results against MTAG multi-trait analysis. Handles sample overlap via the LDSC sampling-covariance matrix, identifies and resolves Heywood cases, and verifies model fit with CFI / RMSEA. Use when modeling latent genetic architecture across correlated traits, running multivariate GWAS on a shared factor, distinguishing factor-mediated from trait-specific SNP effects, or comparing GenomicSEM common-factor results against MTAG when both depend on accurate sampling covariance.
tool_type: r
primary_tool: GenomicSEM
---

## Version Compatibility

Reference examples tested with: GenomicSEM 0.0.5+ (GitHub `GenomicSEM/GenomicSEM`), lavaan 0.6-17+, LDSC v1.0.1+ (Python 3; prefer `abdenlab/ldsc-python3` v2.0.0 -- `belowlab/ldsc` v3.0.1 README states the CLI is broken; Docker `jtb114/ldsc:latest` is the belowlab fallback), baselineLD_v2.2 annotations (alkesgroup.broadinstitute.org/LDSCORE), MTAG 1.0.8+ (Python; `JonJala/mtag`), R 4.4+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('GenomicSEM')` then `?ldsc`, `?commonfactor`, `?usermodel`, `?commonfactorGWAS`
- Python (LDSC, MTAG): `<tool>.py -h` and inspect the source under `ldsc/` or `mtag/`

GenomicSEM is GitHub-only (never on CRAN). If `ldsc()` or `usermodel()` throws an error about lavaan syntax or non-positive-definite covariance, introspect the installed API (`getMethod('ldsc')`) and adapt rather than retrying.

# Genomic SEM

**"Model the latent genetic architecture across several correlated GWAS"** -> Treat each GWAS as a measured indicator of one or more latent genetic factors and fit a structural equation model to the LDSC-derived genetic covariance matrix S and its sampling covariance V (Grotzinger 2019 Nat Hum Behav 3:513). The framework extends naturally to a multivariate GWAS in which a SNP is regressed on a latent factor (common-factor GWAS), with Q_SNP testing whether the SNP effect is homogeneous across factor loadings. Sample overlap between input GWAS is absorbed by the off-diagonals of V; ignoring V inflates Type-I.

- R: `GenomicSEM::ldsc()` produces the (S, V) covariance pair from munged sumstats
- R: `GenomicSEM::commonfactor()` fits a single-factor CFA across all traits in S
- R: `GenomicSEM::usermodel()` fits an arbitrary lavaan-syntax model
- R: `GenomicSEM::commonfactorGWAS()` runs SNP -> factor multivariate GWAS with Q_SNP
- R: `GenomicSEM::userGWAS()` runs arbitrary multivariate SNP regression with per-path Q_SNP
- Python (alternative): `mtag.py --sumstats t1,t2,t3 --out mtag_out` (multi-trait power boost on individual traits)

## Statistical Model Taxonomy

| Method | Latent structure | Min traits | SNP-level test | Strength | Fails when |
|--------|------------------|-----------|----------------|----------|------------|
| Common-factor CFA (Grotzinger 2019) | Single F loading all traits | 3 | None (model-fit only) | Tests whether shared variance is unidimensional | Heterogeneous architecture; CFI < 0.9; near-zero loadings |
| User-specified CFA (`usermodel`) | Pre-specified lavaan syntax | 3 | None | Confirmatory; arbitrary structure | Misspecified model; identification under-determined |
| ESEM | Exploratory rotation; cross-loadings allowed | 6+ | None | When factor count and structure unknown | Few traits; collinear traits; rotation arbitrary |
| Common-factor GWAS (`commonfactorGWAS`) | SNP -> F -> trait1..k | 3 | Wald on F + Q_SNP heterogeneity | Discovers SNPs acting via the common factor; flags Q_SNP outliers | Q_SNP-significant SNPs not interpretable as factor SNPs |
| User GWAS (`userGWAS`) | Arbitrary SNP-path lavaan | 3 | Wald per path + Q_SNP | Tests SNP on any specified path | Highly parameterized models lose power |
| Multivariate Wald test | Joint test across SNP -> trait paths | 2+ | Joint chi-square | Boost power when SNP affects multiple traits | Heterogeneous SNP effects collapse joint test |
| Stratified GenomicSEM (Grotzinger 2022 Behav Genet 52:180) | Factor model with sLDSC-partitioned annotations | 3 | Per-annotation factor tau | Localizes heritability of the factor to functional categories | Same sLDSC failure modes (small annotation, collinearity) |
| MTAG (Turley 2018 Nat Genet 50:229) | Empirical-Bayes shrinkage across correlated traits | 2 | Per-trait shrunk z-score | Boosts marginal power for any input trait | MaxFDR > 5% indicates heterogeneity violates MTAG assumption |

Methodology evolves; verify the current Grotzinger 2023+ tutorials at `github.com/GenomicSEM/GenomicSEM/wiki` before locking a method. ESEM rotation choice (geomin vs target rotation) is an active area; report sensitivity to rotation.

## MTAG vs GenomicSEM Common-Factor GWAS

Both methods exploit genetic correlation among input GWAS, but their goals and outputs differ.

| Property | MTAG | GenomicSEM commonfactorGWAS |
|----------|------|------------------------------|
| Output | Per-trait shrunk z-scores | SNP effect on latent factor |
| Sample-overlap handling | Bivariate LDSC intercept | Full LDSC sampling-covariance matrix V |
| Heterogeneity diagnostic | MaxFDR (Turley 2018) | Q_SNP (Grotzinger 2019) |
| Interpretation | "Boosted power for trait k" | "Effect on what the traits share" |
| Min traits | 2 | 3 (otherwise factor not identified) |
| Best when | Power-boost an individual trait | Common factor hypothesized |

Both depend on accurate sampling covariance. MTAG fails (MaxFDR > 5%) under the same heterogeneity that produces large Q_SNP in GenomicSEM. The two methods should be reported together when the prior on a common factor is non-trivial; agreement increases confidence, disagreement points to architecture-specific SNPs.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Multi-trait GWAS power boost for one focal trait | MTAG | Optimized for per-trait marginal power |
| Common-factor architecture hypothesized | GenomicSEM `commonfactorGWAS` | Tests SNP -> factor; reports Q_SNP heterogeneity |
| Heterogeneous architecture (>1 latent factor) | ESEM, then confirmatory `usermodel` | Exploratory first, then confirm |
| Confirming a pre-specified factor structure | `usermodel` with lavaan syntax | Confirmatory factor analysis |
| Partition heritability of factor across annotations | Stratified GenomicSEM | Combines sLDSC + factor model |
| Mediation in a SEM framework | `usermodel` with indirect path | Path coefficients + delta-method SE |
| Sample overlap unknown or any-overlap suspected | Always use `ldsc()` output as input | V matrix off-diagonals absorb overlap |
| Cross-ancestry common-factor analysis | Run per-ancestry, compare loadings; no published cross-ancestry SEM as of 2026 | Method not yet validated for mixed-ancestry V |
| Single biobank for all traits (e.g., UKB only) | GenomicSEM with `ldsc()`; the V matrix will reflect overlap | Equivalent to one-sample MR -- the V matrix is the correction |
| Comparing GenomicSEM and MTAG on the same traits | Run both; compare top hits + heterogeneity | Concordance increases confidence; divergence flags heterogeneity |

## Per-Method Failure Modes

### Heywood case (negative residual variance)

**Trigger:** A residual variance estimate is < 0, or a standardized loading exceeds 1.

**Mechanism:** Empirical underidentification; the genetic covariance matrix S is near-singular OR a trait has near-zero specific variance under the model. The maximum-likelihood / DWLS estimator runs past the boundary of the parameter space.

**Symptom:** `lavaan` warning "some estimated lv variances are negative" or "covariance matrix is not positive definite"; standardized loading > 1; non-convergence.

**Fix:** First, inspect the LDSC S matrix for genetic correlations near 1 (multicollinearity). Drop or merge near-identical traits. Second, constrain the offending residual variance to be non-negative in the lavaan syntax (`trait1 ~~ a*trait1; a > 0`). Third, verify the V matrix is positive definite via `chol(V_LD)`; if not, the bivariate LDSC inputs disagree on intercept signs and need re-munging. Never re-fit without diagnosing the cause.

### Sample overlap mis-specified

**Trigger:** Using LDSC intercept manually or supplying covariance from non-`ldsc()` source.

**Mechanism:** GenomicSEM's `ldsc()` function returns a list with `S` (genetic covariance) AND `V` (sampling covariance of the lower-triangle of S). The V off-diagonals capture sample overlap via cross-trait LDSC intercept. Skipping V and supplying only S treats all inputs as independent samples; Type-I error inflates because the sampling distribution under H0 is wrong.

**Symptom:** SE on factor loadings far too small; many SNPs significant in common-factor GWAS that don't replicate; comparison to MTAG shows disagreement consistent with overlap.

**Fix:** Always pass the full output of `ldsc()` -- both S and V -- to `commonfactor()`, `usermodel()`, and `commonfactorGWAS()`. Never construct S manually from rg estimates.

### Q_SNP not reported in commonfactorGWAS

**Trigger:** Running `commonfactorGWAS()` and reporting only the factor p-value per SNP.

**Mechanism:** Q_SNP tests heterogeneity of the SNP's effect across factor loadings (Grotzinger 2019 supplement). A SNP with significant Q_SNP violates the common-factor assumption: its effect is NOT mediated by the factor, and the factor estimate is meaningless for that SNP.

**Symptom:** Top "common-factor SNPs" are dominated by trait-specific effects; replication in independent cohorts is poor for SNPs with high Q_SNP.

**Fix:** Always report Q_SNP p-value alongside the factor p-value. Flag SNPs with Q_SNP p < 5e-8 / N_factor_SNPs (Bonferroni for the discovered set) as architecture-violating and exclude from "common-factor SNP" claims. Re-fit those SNPs in `userGWAS()` with separate paths to each trait.

### Trait inclusion under heterogeneous factor structure

**Trigger:** Forcing a common-factor model on traits that don't share a single latent factor.

**Mechanism:** When two or more traits load on a different factor than the rest, the single-factor model misfits. lavaan still returns parameter estimates but model fit is poor.

**Symptom:** CFI < 0.9; RMSEA > 0.08; some standardized loadings near 0 while others near 1; chi-square highly significant even after accounting for N.

**Fix:** Run ESEM first (`commonfactor` then `usermodel` with cross-loadings allowed) to discover structure. If two factors emerge, fit a two-factor `usermodel`. Drop traits with near-zero loadings on all factors. Document the model search.

### MTAG MaxFDR > 5%

**Trigger:** Running MTAG on traits with low pairwise genetic correlation or with one trait that has a very different architecture.

**Mechanism:** MTAG assumes a homogeneous variance-covariance structure across SNPs. When heterogeneity dominates, the empirical-Bayes shrinkage can over-claim SNPs in the focal trait. Turley 2018 defines MaxFDR as the maximum estimated false discovery rate under worst-case heterogeneity; > 5% invalidates the published trait-specific summary statistics.

**Symptom:** MTAG output file reports `maxFDR` > 0.05; per-trait MTAG hits don't replicate in independent cohorts.

**Fix:** Check pairwise rg via LDSC; if any pair is < 0.7, MTAG is risky. Drop the most heterogeneous trait and re-run. Alternatively, switch to GenomicSEM's `commonfactorGWAS` which models heterogeneity explicitly via Q_SNP.

### Non-positive-definite V_LD matrix

**Trigger:** LDSC inputs from different ancestry GWAS, or one trait with very low mean chi-square (< 1.02).

**Mechanism:** V is the sampling covariance of vech(S); when individual entries of S have huge SE relative to off-diagonal covariance, the resulting V is not positive definite (negative eigenvalues).

**Symptom:** `commonfactor()` errors with "matrix is not positive definite" before fitting; `eigen(LDSCoutput$V)$values` shows negative values.

**Fix:** Verify per-trait mean chi-square via `ldsc()` log; below 1.02, exclude that trait. Verify all GWAS are EUR ancestry (or match ancestry of LD scores). Apply nearest-PD smoothing via `Matrix::nearPD(V)$mat` ONLY as a last resort and document the approximation in methods.

## Model Fit Diagnostics

| Index | Acceptable | Good | Source |
|-------|-----------|------|--------|
| CFI | >= 0.90 | >= 0.95 | Hu & Bentler 1999 Struct Equ Model 6:1 |
| TLI / NNFI | >= 0.90 | >= 0.95 | Hu & Bentler 1999 |
| RMSEA | <= 0.08 | <= 0.05 | Hu & Bentler 1999 |
| SRMR | <= 0.08 | <= 0.05 | Hu & Bentler 1999 |
| chi-square p-value | (less informative at large N) | n/a | Penalize for N inflation |

AIC / BIC are used for nested-model comparison (lower is better); only compare nested models fit on the same S.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| CFI >= 0.95 | Hu & Bentler 1999 | Conventional good fit; SEM literature default |
| RMSEA <= 0.05 | Hu & Bentler 1999 | Conventional good fit |
| RMSEA 0.05 - 0.08 | Hu & Bentler 1999 | Acceptable; flag as "adequate" not "good" |
| Q_SNP p < 0.05 / N_SNP_factor | Grotzinger 2019 Nat Hum Behav 3:513 | Bonferroni for heterogeneity among factor-significant SNPs |
| MTAG MaxFDR < 5% | Turley 2018 Nat Genet 50:229 | Above this, MTAG marginal trait results invalid |
| Per-trait LDSC mean chi-square > 1.02 | LDSC documentation | Below this, V entries too noisy; factor SE inflated |
| Standardized loading 0.3 - 0.9 typical | SEM conventions | < 0.3 trait loads weakly; > 0.95 may indicate over-fit / collinearity |
| Min 3 traits for common factor | SEM identification | Single factor with k traits has k(k+1)/2 moments; needs k>=3 to identify |

## Standard Workflow

**Goal:** Fit a common-factor model across correlated GWAS and run a multivariate GWAS on the factor with Q_SNP.

**Approach:** Munge sumstats -> LDSC for (S, V) -> common-factor CFA -> inspect fit -> prepare SNPs -> common-factor GWAS -> report factor effects with Q_SNP flags.

```r
library(GenomicSEM)

# Step 1: Munge sumstats (one-time; produces .sumstats.gz files)
files <- c('raw/trait1.txt', 'raw/trait2.txt', 'raw/trait3.txt')
hm3 <- 'w_hm3.snplist'  # HapMap3 SNP list
trait_names <- c('trait1', 'trait2', 'trait3')
N <- c(150000, 200000, 175000)
munge(files = files, hm3 = hm3, trait.names = trait_names, N = N)

# Step 2: LDSC produces both S (genetic covariance) and V (sampling covariance)
traits <- c('trait1.sumstats.gz', 'trait2.sumstats.gz', 'trait3.sumstats.gz')
ldsc_results <- ldsc(
    traits = traits,
    sample.prev = c(0.5, 0.5, NA),    # case prevalence; NA for continuous
    population.prev = c(0.05, 0.05, NA),
    ld = 'eur_w_ld_chr/',
    wld = 'eur_w_ld_chr/',
    trait.names = trait_names
)
# ldsc_results$S = genetic covariance; ldsc_results$V = sampling covariance

# Step 3a: Common-factor CFA via DWLS
cf_fit <- commonfactor(covstruc = ldsc_results, estimation = 'DWLS')
print(cf_fit$modelfit)  # CFI, RMSEA, SRMR, chi-square
print(cf_fit$results)   # loadings + SEs

# Step 3b: Alternative -- user-specified two-factor model (>=3 indicators per factor)
# Identification rule: each factor needs >= 3 indicators OR one anchor loading fixed
# to 1 plus factor variance free. A factor with a single indicator is NOT identified.
model_syntax <- '
    F1 =~ NA*trait1 + trait2 + trait3
    F2 =~ NA*trait4 + trait5 + trait6
    F1 ~~ 1*F1
    F2 ~~ 1*F2
    F1 ~~ F2
'
user_fit <- usermodel(covstruc = ldsc_results, model = model_syntax, estimation = 'DWLS')
```

`estimation = 'DWLS'` (diagonally weighted least squares) is the default and is required when V is large; `'ML'` is faster but assumes a known V and can produce wrong SE under sample overlap.

### ESEM (Exploratory Factor Structure)

When the factor structure is unknown, fit `usermodel()` with all loadings free across all factors, then apply a rotation post-fit. Rotation choices: **geomin oblique** (default; allows factor correlation), **target rotation** (Browne 2001 Multivariate Behav Res 36:111; uses a hypothesized loading template), **quartimin** (orthogonal; assumes factors are uncorrelated).

```r
model_esem <- '
    F1 =~ NA*trait1 + trait2 + trait3 + trait4
    F2 =~ NA*trait1 + trait2 + trait3 + trait4
    F1 ~~ 1*F1
    F2 ~~ 1*F2
    F1 ~~ F2
'
esem_fit <- usermodel(covstruc = ldsc_results, model = model_esem, estimation = 'DWLS')
# Rotate post-fit via GPArotation::GPForth/GPFoblq or lavaan::rotate()
```

**Decision:** ESEM for K-factor exploration when structure is unknown; CFA via `usermodel()` once a structure is confirmed. Report rotation sensitivity (geomin vs target vs quartimin) and treat as exploratory.

**Cross-loadings.** Cross-loadings (one trait loads on > 1 factor) are common in psychiatric and behavioral GWAS. Brown 2015 *Confirmatory Factor Analysis for Applied Research* recommends allowing cross-loadings first and using modification indices to guide simplification. Allow a cross-loading when constraining residual variance otherwise forces a Heywood case. Constrain when CFI < 0.9 and modification indices instead suggest a correlated residual between two indicators (which is the more parsimonious fix).

### userGWAS for Custom Path Models

`userGWAS()` fits arbitrary lavaan-syntax SNP regressions and is the right tool when the SNP needs to be tested on multiple paths simultaneously (e.g. factor-mediated effect AND a direct effect on one indicator).

```r
# Test SNP -> F path + SNP -> trait1 direct path simultaneously
model <- '
    F =~ NA*trait1 + trait2 + trait3
    F ~~ 1*F
    F ~ SNP
    trait1 ~ SNP    # direct effect on trait1, partialed out of F
'
user_results <- userGWAS(covstruc = ldsc_results,
                         SNPs = ss,
                         estimation = 'DWLS',
                         model = model,
                         sub = c('F~SNP', 'trait1~SNP'),
                         parallel = TRUE,
                         cores = 8)
# Output columns include: lhs, op, rhs, est, SE, Z, Pvalue, Q_pval, Q_df.
# Q_pval per SNP measures heterogeneity across loadings AFTER conditioning on
# the direct path; remaining Q_SNP signal indicates a third path is needed.
```

### Higher-Order / Bifactor / p-Factor Models

Use case: psychiatric genetics p-factor (Caspi 2014 Clin Psychol Sci 2:119; Grotzinger 2022 Nat Genet 54:548 cross-disorder), cognitive g-factor (de la Fuente 2021 Nat Hum Behav 5:49).

Hierarchical template -- first-order factors load on a single second-order p-factor:

```r
model_pfactor <- '
    # First-order factors
    INT =~ NA*trait_anx + trait_dep + trait_neuro       # internalizing
    EXT =~ NA*trait_adhd + trait_alc + trait_subst       # externalizing
    THT =~ NA*trait_scz + trait_bp                       # thought-disorder
    # Second-order p-factor
    p =~ NA*INT + EXT + THT
    INT ~~ 1*INT
    EXT ~~ 1*EXT
    THT ~~ 1*THT
    p ~~ 1*p
'
```

Bifactor alternative: `p =~` all traits directly, with `INT`/`EXT`/`THT` as orthogonal residual factors. Bifactor typically gives tighter CFI/RMSEA but the substantive interpretation of the residual factors is harder; bifactor is also prone to over-fitting at modest trait counts (Bonifay 2017 Struct Equ Model 24:378). Cite Grotzinger 2022 Nat Genet 54:548 and Mallard 2022 Am J Psychiatry 179:528 for the canonical psychiatric implementations.

## Common-Factor GWAS with Q_SNP

**Goal:** Identify SNPs that affect the latent factor and flag SNPs whose effect is heterogeneous across loadings.

**Approach:** Build SNP-by-trait effect-and-SE matrix via `sumstats()`, then fit the SNP-augmented model genome-wide via `commonfactorGWAS()`. Report both factor p-value and Q_SNP p-value per SNP.

`sumstats()` flips effect signs based on the reference panel A1/A2 to enforce consistent allele coding across input GWAS. Required input columns (case-sensitive) are `SNP, A1, A2, BETA/Z/OR, SE, P, N, MAF`; some are conditional on `se.logit` / `OLS`. Silent failures are almost always column-name mismatches (e.g. `EA`/`NEA` instead of `A1`/`A2` -> 0% SNPs retained) or a missing `N` column -> SNPs dropped. Always run `head(read.table(file, header=TRUE, nrow=2))` per input before the `sumstats()` call.

```r
# Prepare per-SNP betas and SEs across all input GWAS
ss <- sumstats(
    files = c('raw/trait1.txt', 'raw/trait2.txt', 'raw/trait3.txt'),
    ref = 'reference.1000G.maf.0.005.txt',
    trait.names = trait_names,
    se.logit = c(TRUE, TRUE, FALSE),     # TRUE if trait is logistic-scale (case-control); FALSE if continuous
    OLS = c(FALSE, FALSE, TRUE),
    linprob = c(FALSE, FALSE, FALSE),
    N = N,
    info.filter = 0.9,                   # standard imputation INFO threshold
    maf.filter = 0.01
)

# GenomicSEM internally detects the OS via Sys.info()[['sysname']] and chooses
# PSOCK (Windows) vs FORK (Linux/Mac) clusters automatically; there is no user
# `Operating=` argument. MPI=TRUE switches to an mpirun-based strategy for
# cluster job submission. parallel=TRUE is the default.
cfgwas <- commonfactorGWAS(
    covstruc = ldsc_results,
    SNPs = ss,
    estimation = 'DWLS',
    parallel = TRUE,
    cores = 8,
    MPI = FALSE
)

# cfgwas columns include rsID/chr/BP/MAF/A1/A2/est/se_c/Z_Estimate/Pval_Estimate (factor effect)
# plus heterogeneity columns Q / Q_df / Q_pval. Q_pval IS the per-SNP heterogeneity
# test (often referred to as Q_SNP in the literature; the column name in the data.frame is Q_pval).
cfgwas$factor_sig <- cfgwas$Pval_Estimate < 5e-08
cfgwas$qsnp_sig <- cfgwas$Q_pval < (0.05 / sum(cfgwas$factor_sig))
cfgwas$factor_only <- cfgwas$factor_sig & !cfgwas$qsnp_sig
```

The "factor-only" subset (factor-significant AND Q_SNP non-significant) is the publication-grade set of common-factor SNPs.

## MTAG Comparison

**Goal:** Cross-check GenomicSEM common-factor results against MTAG per-trait shrunk z-scores.

**Approach:** Run MTAG CLI on the same input sumstats; compare top hits with GenomicSEM factor hits. Report MaxFDR.

```bash
# MTAG CLI (Python)
python mtag.py \
    --sumstats trait1.txt,trait2.txt,trait3.txt \
    --n_min 0 \
    --use_beta_se \
    --out mtag_results

# Check MaxFDR per trait
grep -i maxfdr mtag_results.log
# Each per-trait MTAG file: mtag_results_trait<k>.txt
```

If MaxFDR > 0.05 for any trait, MTAG results for that trait are unreliable; GenomicSEM with Q_SNP filtering is the more defensible report.

## Stratified GenomicSEM (Partitioned Heritability of Factor)

For partitioning the heritability of the latent factor across functional annotations, use `s_ldsc()` (stratified LDSC inside GenomicSEM) and pass the multi-annotation output to a stratified model fit.

```r
# Stratified LDSC across baseline + custom annotations
s_results <- s_ldsc(
    traits = traits,
    sample.prev = c(0.5, 0.5, NA),
    population.prev = c(0.05, 0.05, NA),
    ld = 'baselineLD_v2.2.',
    wld = 'weights.hm3_noMHC.',
    frq = '1000G.EUR.QC.',
    trait.names = trait_names
)

# enrich() inventory:
#   params: lavaan syntax of the parameter under enrichment (loading, residual var, or F~~F latent var)
#   fix='regressions': hold regression paths fixed at the genome-wide estimate during stratified fit
#   std.lv=FALSE: do not standardize the latent variance
#   rm_flank=TRUE: drop flanking-window contributions (default)
#   tau=FALSE: do not return per-annotation tau coefficients
#   base=TRUE: include baseline annotation contributions in the partition
#   toler=NULL: matrix-inversion tolerance (let GenomicSEM choose; supply a small value when S is near-singular)
strat_factor <- enrich(s_covstruc = s_results,
                       model = '',
                       params = 'F =~ trait1',
                       fix = 'regressions',
                       std.lv = FALSE,
                       rm_flank = TRUE,
                       tau = FALSE,
                       base = TRUE,
                       toler = NULL)
```

The output gives per-annotation enrichment of the factor h2 -- the analog of cell-type S-LDSC for the latent factor (Grotzinger 2022 Behav Genet 52:180).

## Computational Footprint

| Step | Runtime | Hardware |
|------|---------|----------|
| `ldsc()` multi-trait sampling covariance | minutes | laptop |
| `commonfactor()` / `usermodel()` (no SNP loop) | seconds | laptop |
| `commonfactorGWAS()` over 6-8M SNPs | 4-24h depending on cores | cluster recommended |
| `userGWAS()` over 6-8M SNPs with complex path model | 8-48h | cluster recommended |
| Stratified GenomicSEM with 50+ annotations | 1-3 days | cluster |
| MTAG over 6-8M SNPs | 1-2h | laptop or cluster |

Cluster runs of `commonfactorGWAS()` / `userGWAS()` should use `MPI=TRUE` when submitting via mpirun; GenomicSEM detects the OS internally (via `Sys.info()[['sysname']]`) and selects FORK (Linux/Mac) vs PSOCK (Windows) cluster types automatically -- there is no `Operating=` user argument. On a Mac/Windows workstation, reduce `cores` to the physical-core count to avoid PSOCK fork failures.

## Reconciliation: When GenomicSEM and MTAG Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| GenomicSEM factor SNP sig, MTAG sig for all traits | Genuine common-factor SNP | Report; high confidence |
| GenomicSEM factor SNP sig, MTAG sig in only 1 trait | Q_SNP heterogeneity likely; one-trait-dominant | Check Q_SNP; if sig, this is NOT a factor SNP |
| MTAG sig, GenomicSEM factor null, Q_SNP sig | Trait-specific SNP captured by MTAG shrinkage | Report as trait-specific, not common-factor |
| Both null but per-trait univariate sig | Power loss from multivariate parameterization | Re-check sample overlap V matrix |
| GenomicSEM and MTAG both sig but opposite direction | Sample-overlap mis-specification OR sign error in munging | Re-munge with same allele convention; re-run `ldsc()` |
| MTAG MaxFDR > 5%, GenomicSEM with Q_SNP works | MTAG assumption violated | Prefer GenomicSEM as primary |
| One-trait GWAS sig but common-factor not | Trait-specific architecture | Don't force into common-factor frame |

**Operational rule for publication:** A common-factor SNP claim requires (1) factor p < 5e-8, (2) Q_SNP p > 0.05 / N_factor_SNPs (non-heterogeneous), and (3) replication in an independent set of traits or cohorts. Trait-specific SNPs from MTAG require MaxFDR < 5% for the trait. Reporting only the factor effect without Q_SNP is the most common reviewer-flagged error.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `commonfactor()` complains "S not positive definite" | Genetic correlations near +/-1 among inputs | Drop redundant traits; verify rg < 0.95 pairwise |
| Standardized loading > 1 | Heywood case; under-identification | Constrain residual variance >= 0; inspect S for collinearity |
| Factor p-value reported, Q_SNP not reported | Default focus is on factor effect | Always report Q_SNP from `commonfactorGWAS` output |
| `ldsc()` fails with "category not found" | Wrong LD score column names (legacy format) | Use Python 3 LDSC fork; download `eur_w_ld_chr/` from alkesgroup |
| `lavaan` says "model not identified" | Too few traits for too many parameters | Need >= 3 traits per factor; constrain factor variance to 1 |
| MTAG `MaxFDR` not in log | Older MTAG version (< 1.0.7) | Update MTAG; MaxFDR reporting added late 2019 |
| Singular V matrix on smooth `nearPD` | One trait has near-zero h2 or mean chi-square < 1.02 | Drop the trait; do not smooth as a fix |
| `usermodel()` slow or fails | Complex syntax + many traits | Simplify model; estimate with `estimation = 'DWLS'`, not `'ML'`, when V is informative |
| Sumstats output has zero overlap with reference | Allele coding mismatch in `sumstats()` | Check `se.logit` and `OLS` settings per trait; align A1/A2 |
| GenomicSEM and TwoSampleMR give different rg | TwoSampleMR uses bivariate LDSC; GenomicSEM uses the same S | Match the underlying LDSC reference panel and weights |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Q_SNP reported?" | `Q_pval` column reported alongside SNP->factor effect; Bonferroni threshold 5e-8 / N_factor_SNPs applied |
| "MaxFDR > 5%?" | MTAG `maxFDR` reported per trait; > 5% invalidates MTAG for that trait -> GenomicSEM common-factor used instead |
| "Sample overlap absorbed?" | Full V matrix from `ldsc()` is the input to all model fits; S is never constructed manually from pairwise rg estimates |
| "Model fit?" | CFI >= 0.95, RMSEA <= 0.08, SRMR <= 0.08 reported; AIC / BIC for nested comparison; chi-square reported but treated as inflated at large N |
| "Why a factor model and not MTAG?" | Factor structure tested first; if traits load on a common factor with CFI > 0.95, common-factor GWAS preferred -- explicitly models heterogeneity via Q_SNP |
| "Heywood case?" | Negative residual variance constrained >= 0; OR the offending indicator dropped and the model re-specified; the choice is documented in methods |
| "Cross-ancestry?" | Run per-ancestry; no validated cross-ancestry V matrix as of 2026; loadings compared qualitatively |
| "Why DWLS and not ML?" | ML assumes V is known; DWLS uses the empirical V from `ldsc()` and is the appropriate estimator under sample overlap |

## Tool Installation

```r
# GenomicSEM is GitHub-only
remotes::install_github('GenomicSEM/GenomicSEM')

# Dependencies
install.packages(c('lavaan', 'Matrix', 'gdata'))

# Optional companions
remotes::install_github('MRCIEU/TwoSampleMR')  # for downstream MR using factor GWAS as exposure
```

For Python tools:

```bash
# LDSC python3 fork (GenomicSEM input format). belowlab/ldsc v3.0.1 broke the
# --h2/--rg/--h2-cts CLI per its README; use abdenlab/ldsc-python3 (v2.0.0)
# for a working CLI. Docker jtb114/ldsc:latest is the belowlab fallback.
git clone https://github.com/abdenlab/ldsc-python3.git
cd ldsc-python3 && pip install -r requirements.txt

# MTAG
git clone https://github.com/JonJala/mtag.git
cd mtag && pip install -r requirements.txt
```

Pre-downloaded reference files: `eur_w_ld_chr/`, `baselineLD_v2.2.*`, `w_hm3.snplist`, and 1000G allele-frequency files are hosted at `alkesgroup.broadinstitute.org/LDSCORE/`.

## References

- Grotzinger AD et al 2019 Nat Hum Behav 3:513 (GenomicSEM, common-factor GWAS, Q_SNP)
- Grotzinger AD et al 2022 Behav Genet 52:180 (Stratified GenomicSEM)
- Turley P et al 2018 Nat Genet 50:229 (MTAG; MaxFDR)
- Bulik-Sullivan B et al 2015 Nat Genet 47:291 (LDSC for genetic covariance)
- Bulik-Sullivan B et al 2015 Nat Genet 47:1236 (bivariate LDSC, sample overlap)
- Rosseel Y 2012 J Stat Softw 48:1-36 (lavaan package)
- Hu LT & Bentler PM 1999 Struct Equ Model 6:1 (CFI / RMSEA cutoffs)
- Asparouhov T & Muthen B 2009 Struct Equ Model 16:397 (ESEM framework)
- Finucane HK et al 2015 Nat Genet 47:1228 (S-LDSC, foundation for stratified GenomicSEM)
- Gazal S et al 2017 Nat Genet 49:1421 (baseline-LD annotations)
- Demange PA et al 2021 Nat Genet 53:35 (GenomicSEM applied to cognitive traits; Q_SNP in practice)
- Mallard TT et al 2022 Am J Psychiatry 179:528 (multivariate GWAS of externalizing via GenomicSEM)
- de la Fuente J et al 2021 Nat Hum Behav 5:49 (GenomicSEM for cognitive g factor)
- Skrivankova VW et al 2021 JAMA 326:1614 (STROBE-MR; relevant when downstream MR uses factor GWAS)

## Related Skills

- causal-genomics/mendelian-randomization - Use factor-GWAS effect sizes as MR exposure
- causal-genomics/genetic-correlation - Bivariate LDSC produces the off-diagonals of the S matrix; GenomicSEM is the multi-trait extension
- causal-genomics/heritability-partitioning - LDSC and S-LDSC foundations for stratified GenomicSEM
- causal-genomics/colocalization-analysis - Cross-trait colocalization at common-factor loci
- causal-genomics/pleiotropy-detection - Q_SNP is a per-SNP pleiotropy diagnostic; sibling concept
- causal-genomics/fine-mapping - Resolve factor-significant loci to credible sets
- causal-genomics/mediation-analysis - SEM mediation paths overlap with `usermodel` indirect effects
- causal-genomics/transcriptome-wide-association - TWAS on factor sumstats from `commonfactorGWAS`
- population-genetics/association-testing - GWAS sumstats are the input format
