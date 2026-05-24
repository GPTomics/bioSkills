---
name: bio-causal-genomics-mendelian-randomization
description: Estimate causal effects of an exposure on an outcome from GWAS summary statistics using genetic instruments. Implements IVW (fixed/random), MR-Egger, weighted median/mode, MR-RAPS, CAUSE, GSMR-HEIDI, MR-PRESSO, MVMR, MR-Clust, LCV, and LHC-MR via TwoSampleMR, MendelianRandomization, MR-PRESSO, cause, and lhcMR. Use when testing causal direction between traits, evaluating drug-target effects via cis-pQTL/cis-eQTL, performing multivariable mediation MR, distinguishing causation from correlated horizontal pleiotropy, or producing STROBE-MR-compliant sensitivity batteries.
tool_type: mixed
primary_tool: TwoSampleMR
---

## Version Compatibility

Reference examples tested with: TwoSampleMR 0.6.0+, MendelianRandomization 0.10+, MR-PRESSO 1.0+, cause 1.2+, MVMR 0.4+, ieugwasr 1.0+, MRlap 0.0.3.2+, coloc 5.2+, mrclust 0.1+, lhcMR 0.0.1+, R 4.4+. Both TwoSampleMR 0.6.0 and ieugwasr 1.0 are the JWT-transition versions; older versions still expect deprecated OAuth.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI (plink, GCTA-GSMR): `<tool> --version` then `<tool> --help`

If code throws an error referencing a function that has moved (e.g. `ieugwasr::ld_clump` vs `TwoSampleMR::clump_data`) or an OAuth token failure, introspect the installed API and adapt the example rather than retrying.

# Mendelian Randomization

**"Test whether trait X causally affects trait Y from GWAS summary statistics"** -> Use genetic variants as instrumental variables (IVs) that satisfy three assumptions (relevance, independence, exclusion restriction) to estimate `beta_causal = beta_outcome / beta_exposure` under the IV framework (Davey Smith & Ebrahim 2003 IJE 32:1; Burgess Thompson 2017 SAGE textbook). Tool choice is a decision about the **regime** (one-sample vs two-sample, sparse vs polygenic, drug-target vs polygenic exposure) and the **pleiotropy model** (balanced, directional InSIDE, correlated horizontal). Wrong tool inflates Type-I error or attenuates true effects in a direction predictable from the bias structure.

- R: `TwoSampleMR::mr()` orchestrates IVW + Egger + weighted median + weighted mode in one call
- R: `MendelianRandomization::mr_ivw / mr_egger / mr_median / mr_mbe / mr_conmix` per-method API (S4 objects; MR-RAPS is NOT in this package -- use `TwoSampleMR::mr_raps()` which wraps the GitHub `mr.raps`)
- R: `MRPRESSO::mr_presso()` global / outlier / distortion tests
- R: `cause::cause()` correlated horizontal pleiotropy mixture
- R: `MVMR::strength_mvmr() + MVMR::ivw_mvmr()` multivariable conditional-F + IVW

## Statistical Model Taxonomy

| Method | Pleiotropy assumption | Min instruments | Strength | Fails when |
|--------|------------------------|-----------------|----------|------------|
| IVW (fixed) | All IVs valid | 2 | Most efficient under no pleiotropy | Any directional or balanced pleiotropy inflates Type-I |
| IVW (random effects) | Balanced + InSIDE | 3 | Standard primary; absorbs heterogeneity into wider SE | Directional pleiotropy biases the point estimate |
| MR-Egger | Directional pleiotropy + InSIDE | 10+ for power | Detects + corrects directional pleiotropy via intercept (Bowden 2015 IJE 44:512) | NOME violated (`I^2_GX < 0.9`); SIMEX correction required; underpowered <10 SNPs |
| Weighted median | Up to 50% invalid IVs | 3 | Robust to a minority of bad instruments (Bowden 2016 Genet Epidemiol 40:304) | >50% invalid IVs |
| Weighted mode | Zero modal pleiotropy (ZEMPA) | 3 | Robust if the modal estimate is unbiased (Hartwig 2017 IJE 46:1985) | Bimodal pleiotropy; small numbers |
| MR-RAPS | Balanced pleiotropy + weak instruments | 10+ | Profile-score robust to weak-IV + balanced horizontal pleiotropy (Zhao 2020 Ann Stat 48:1742) | Strong directional pleiotropy; CRAN-archived 2025-03-01 |
| CAUSE | Correlated horizontal pleiotropy (CHP) | 100+ sig SNPs | Explicit shared-factor mixture; protects against CHP-driven false positives (Morrison 2020 Nat Genet 52:740) | Sparse polygenic exposures; <100 sig SNPs |
| GSMR + HEIDI-outlier | Outlier removal under InSIDE | 10+ | Alternative outlier detection; integrates with LD reference (Zhu 2018 Nat Commun 9:224) | Requires individual-level LD; HEIDI conservative |
| MR-PRESSO | Outlier-driven horizontal pleiotropy | 4+ | Global / outlier / distortion three-step (Verbanck 2018 Nat Genet 50:693) | Blind to CHP; computationally heavy at large NbDistribution |
| MVMR (IVW) | Conditional independence after measured pleiotropy | 1+ per exposure | Accounts for measured horizontal pleiotropy via multivariable regression (Sanderson 2019 IJE 48:713) | Conditional F < 10 on any exposure |
| MR-Clust | Heterogeneous causal effects (multiple mechanisms) | 30+ | Clusters SNPs by their causal-effect estimate (Foley 2020 Bioinformatics 37:531) | Single causal mechanism; small instrument sets |
| Contamination mixture | Mixture of valid + invalid IVs | 10+ | Profile-likelihood mixture (Burgess 2020 Nat Commun 11:376) | Sparse signal |
| LCV | Genome-wide; distinguishes causation vs genetic correlation | All SNPs | Tests `gcp` parameter using LDSC-style block jackknife (O'Connor & Price 2018 Nat Genet 50:1728) | Two-trait covariance dominated by a third confounder |
| LHC-MR | Bidirectional + heritable confounder | All SNPs | Joint likelihood over genome-wide markers; estimates both directions + confounder (Darrous 2021 Nat Commun 12:7274) | Computationally heavy; rare-variant trait |
| MRlap | Sample overlap + winner's curse + weak-IV jointly | Genome-wide sumstats | LDSC-scaffolded joint correction (Mounier & Kutalik 2023 Genet Epidemiol 47:314) | LDSC intercept poorly estimated (h^2 < 0.05); non-EUR without matched LD scores |
| Doubly-Ranked MR (DRMR) | Non-linear, non-parametric | 5+ strata | Non-parametric stratification (Tian 2023 PLoS Genet 19:e1010823); replaces residual stratification when linearity fails | Continuous exposures only; needs individual-level data; Hamilton 2023 medRxiv 23293658 shows stratum-specific bias from age/sex |

Methodology evolves; benchmark consensus shifts every 2-3 years. Verify against the current Slob & Burgess 2020 *Genet Epidemiol*, Burgess 2023 *Wellcome Open Res* "Guidelines for performing Mendelian randomization" (v3+), and STROBE-MR 2021 reporting standards before locking a method as primary.

## Decision Tree by Experimental Scenario

| Scenario | Primary method | Sensitivity battery | Why |
|----------|----------------|----------------------|-----|
| Standard two-sample, independent cohorts, polygenic exposure | IVW (random) | Egger + weighted median + MR-PRESSO + MR-RAPS + Steiger | Default; covers balanced, directional, outlier, weak-IV regimes |
| One-sample (e.g. UK Biobank both ends) | IVW with weak-IV-aware (MR-RAPS) | Egger + LCV + jackknife SE | One-sample bias formulas: Bowden 2019 IJE 48:728; F-stat floor shifts to F >= 20; jackknife SE preferred over analytic at one-sample scale; do NOT run exposure GWAS and outcome GWAS on the same individuals then claim two-sample (Hartwig 2021 collider bias); within-stratum MR (e.g. "MR among smokers") risks collider bias from the stratification variable |
| Partial sample overlap (UKB exposure + UKB outcome) | MR-RAPS with overlap correction | Sample-overlap-adjusted IVW (Burgess 2016 Genet Epidemiol 40:597) | Bias is intermediate, proportional to z-score correlation |
| Drug-target / cis-MR (cis-pQTL, cis-eQTL) | IVW restricted to cis window | Colocalization PP.H4 + LD-prune within window | Exclusion restriction relaxed because the protein/transcript directly mediates effect (Schmidt 2020 Nat Commun 11:3255) |
| MVMR for measured pleiotropy (e.g. LDL adjusted for HDL/TG) | `MVMR::ivw_mvmr` | Conditional F + Q_A heterogeneity | Required when exposures correlate via shared SNPs |
| Mediation MR (X -> M -> Y) | MVMR difference of total vs direct | Two-step MR + product-of-coefficients (Carter 2021 Eur J Epidemiol 36:465) | Network MR; quantifies indirect effect |
| Polygenic exposure with potential CHP (e.g. BMI -> CHD) | CAUSE (primary) + IVW (secondary) | Egger + MR-PRESSO + LCV | CAUSE explicitly models CHP via shared-factor; needs >=100 sig SNPs |
| Binary outcome (e.g. T2D) on linear scale | IVW on log-OR with log-additive coding | All sensitivity on log-OR; report exp(beta) | Linearity of MR estimating equation holds on log-OR not OR |
| Time-to-event (Cox) outcome | IVW on log-HR | Burgess & Labrecque 2018 Eur J Epidemiol 33:947 framework | Non-collapsibility caveats apply |
| Non-linear MR (e.g. alcohol J-curve) | DRMR (Tian 2023) + residual stratification side-by-side | Negative-control outcomes (genotype-vs-sex within strata); Hamilton 2023 limitation cited | Both methods produce stratum-specific bias from age/sex effects (Hamilton 2023 medRxiv 23293658); pre-specify the non-linear hypothesis, do not data-snoop the J-curve, report negative-control sanity checks |
| Single-patient rare disease | Not MR -- use FRASER/DROP outlier framework | See alternative-splicing/outlier-splicing-detection | MR requires summary stats; n=1 is wrong regime |

## One-Sample vs Two-Sample Bias Direction

| Design | Weak-IV bias direction | Reason |
|--------|------------------------|--------|
| One-sample, F<10 | Toward confounded observational estimate (overestimates causal effect if confounding is in same direction) | Sample correlation between IV-X and IV-Y residuals |
| Two-sample non-overlapping, F<10 | Toward null | Independent samples decouple residuals (Burgess 2011 IJE 40:755) |
| Two-sample with partial overlap | Intermediate; proportional to overlap fraction and z-score correlation | Burgess 2016 Genet Epidemiol 40:597; correction available |

**Operational rule:** Whenever both GWAS came from UK Biobank (or any single biobank), treat the analysis as one-sample-equivalent and prefer MR-RAPS as primary. Treating it as "two-sample because separate GWAS files" is a common error and produces overestimates.

### MRlap: unified correction for sample overlap + winner's curse + weak instruments

MRlap (Mounier & Kutalik 2023 Genet Epidemiol 47:314) jointly corrects three biases that previously required three separate tools: sample overlap, winner's curse, and weak-instrument bias. It builds on an LDSC scaffold (cross-trait LD-score regression intercept estimates the overlap-induced covariance) and reweights the IVW estimate against the analytical bias-correction formula.

```r
remotes::install_github('n-mounier/MRlap')   # never on CRAN; bioconductor unsuitable
library(MRlap)

fit <- MRlap(
    exposure = gwas_X_df, exposure_name = 'BMI',
    outcome = gwas_Y_df, outcome_name = 'T2D',
    ld = 'eur_w_ld_chr/', hm3 = 'w_hm3.snplist',   # LDSC reference files
    MR_threshold = 5e-8, MR_pruning_dist = 500, MR_pruning_LD = 0.05
)
fit$MRcorrection$corrected_effect       # overlap + winner's curse + weak-IV corrected
fit$MRcorrection$corrected_effect_se
fit$LDSC$h2_LDSC                         # heritability sanity check
fit$LDSC$lambda                          # sample-overlap-induced lambda; ~0 means no overlap
```

**Decision rule -- prefer MRlap when:** (a) any sample overlap is suspected, (b) only sumstats are available (no individual-level data for re-running GWAS on disjoint samples), (c) exposure discovery and outcome were both run inside the same biobank (UKB-on-UKB, FinnGen-on-FinnGen). MRlap returns NA / unstable estimates when h^2 < 0.05; in that regime, fall back to Burgess 2016 overlap-corrected IVW plus MR-RAPS for the weak-IV component.

## Drug-Target / cis-MR Framework

cis-MR restricts instruments to the cis-regulatory window of the gene encoding the protein/transcript exposure (Schmidt 2020 Nat Commun 11:3255), relaxing the exclusion-restriction assumption because the protein product directly mediates the SNP's effect on the outcome. Operational core: extract cis-pQTL/cis-eQTL within +/-500 kb of the gene; clump at r2 < 0.1 (looser than polygenic MR to retain power within a narrow window); require colocalization PP.H4 >= 0.7; flag protein-altering variants (PAV) which can break SomaScan/Olink aptamer/antibody binding rather than reflect biology.

Full drug-target cis-MR workflow including UKB-PPP / deCODE / Fenland pQTL panels, PAV flagging, Olink vs SomaScan replication (~15-30% cross-platform disagreement), and the operational claim ladder lives in causal-genomics/proteome-mr-drug-target. Use that skill for any drug-target nomination.

### Binary outcomes and non-collapsibility

MR with logistic-GWAS sumstats returns per-allele log-OR on the **population-averaged** (marginal) scale, NOT the conditional log-OR (Burgess & Labrecque 2018 Eur J Epidemiol 33:947; Burgess 2017 Stat Methods Med Res 26:2333). For rare disease (prevalence < 10%), OR ~= RR ~= HR and the distinction is harmless. For common disease, OR diverges from RR/HR and the MR estimate cannot be back-converted to a conditional effect without strong assumptions; report as "per 1-SD increase in genetically-predicted X, OR for Y = ..." rather than implying an individual-level intervention effect.

Collider bias in case-only or disease-progression cohorts (Hu 2022 IJE 51:1289): conditioning the outcome on disease status (e.g. studying progression among diabetics) opens a collider path between any cause of disease and any cause of progression. MR within affected subsets without explicit adjustment for selection probability is fragile; weight by inverse probability of selection or restrict claims to the unconditioned population.

## Per-Method Failure Modes

### IVW under directional pleiotropy

**Trigger:** Several SNPs affect the outcome through pathways not via the exposure, in a consistent direction.

**Mechanism:** IVW is a weighted regression through origin; non-zero mean pleiotropy shifts the slope.

**Symptom:** Egger intercept p < 0.05 with non-zero estimate; IVW differs from weighted median; MR-PRESSO global test p < 0.05.

**Fix:** Use Egger (if `I^2_GX >= 0.9` -- otherwise SIMEX-correct via the `simex` package applied to the Egger fit, treating `se.exposure` as measurement error in `beta.exposure`); cross-check with weighted median, MR-PRESSO, and CAUSE; report IVW only as one of a panel, never alone. The `MendelianRandomization::mr_egger()` function accepts `distribution='normal'` and reports NOME-corrected estimates internally but does NOT expose a SIMEX wrapper.

### Weak-instrument bias direction

**Trigger:** Mean per-instrument F-statistic < 10, or several individual F < 10.

**Mechanism:** Weak IVs amplify finite-sample correlation between IV-X and IV-Y errors; bias direction depends on overlap regime (see table above).

**Symptom:** Estimates shift markedly when removing the weakest instruments; one-sample MR estimates much larger than two-sample.

**Fix:** Compute F per instrument from the EXPOSURE GWAS, not the outcome; exclude F < 10; use MR-RAPS (handles weak IVs by design); for two-sample, also report unweighted IVW (less weak-IV-bias-inflated than weighted in some regimes).

### Winner's curse at P~5e-8

**Trigger:** Discovery GWAS is the source of both instrument selection and effect-size estimates.

**Mechanism:** SNPs that just cross 5e-8 in discovery have over-estimated effect sizes (regression toward the mean in independent replication); MR uses inflated `beta_X`, biasing causal estimate.

**Symptom:** MR effect shrinks substantially when using effect sizes from an independent replication GWAS.

**Fix:** (1) Three-sample design (discovery / replication-for-instrument-effect / outcome) where feasible. (2) When sumstats-only: MRlap (Mounier 2023), MR-SimSS (sample-splitting from sumstats), or RIVW (Ma 2023 JASA -- rerandomized IVW) jointly correct winner's curse + weak IVs + overlap. (3) Sadreev 2024 IJE 52:1209 empirical magnitude: 5-15% bias inflation at p < 5e-8; larger at relaxed thresholds.

### NOME violation invalidating Egger

**Trigger:** Running MR-Egger with `I^2_GX < 0.9`.

**Mechanism:** Egger assumes NO Measurement Error in exposure effect sizes (NOME); when violated, Egger slope is attenuated toward null with reciprocal bias on the intercept.

**Symptom:** `mr_pleiotropy_test()` Egger estimate disagrees with weighted median in magnitude but agrees in direction; `Isq()` function returns <0.9.

**Fix:** Compute `Isq(beta_X, se_X)` (Bowden 2016 IJE 45:1961); if <0.9, apply SIMEX correction via `simex` package or report Egger as exploratory only. The MendelianRandomization package's `mr_egger()` returns a NOME-aware corrected estimate.

### Steiger filter false flag under unmeasured confounding

**Trigger:** Applying `steiger_filtering()` on traits with unmeasured shared confounders (e.g. SES).

**Mechanism:** Steiger compares variance explained in exposure vs outcome per SNP; an unmeasured confounder upstream of both produces SNPs that explain more variance in the outcome than the exposure, falsely flagging "reverse causation" (Hemani Tilling 2022 Wellcome Open Res 7:14).

**Symptom:** Many SNPs flagged as wrong-direction yet biology and prior MR support forward causation.

**Fix:** Treat Steiger as a heuristic, not gospel; cross-validate direction with bidirectional MR (forward + reverse with independent instrument sets); for known-confounder-rich domains (psychiatric traits, SES proxies) use LCV or LHC-MR instead, which jointly model confounders.

### Palindromic SNP harmonization

**Trigger:** SNPs with alleles A/T or C/G near MAF 0.5.

**Mechanism:** Strand orientation is ambiguous for palindromic SNPs when allele frequencies are intermediate; flipping introduces sign errors that look like pleiotropy.

**Symptom:** `harmonise_data()` reports many palindromic SNPs dropped; remaining SNPs show heterogeneity from a handful.

**Fix:** Default `action = 2` (infer from allele frequencies) drops MAF~0.5 palindromes; `action = 3` drops ALL palindromes (most conservative); never use `action = 1` (assumes forward strand) unless both GWAS are guaranteed to use the same strand convention. Document choice in methods.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| F-statistic > 10 per instrument | Staiger & Stock 1997 (linear IV) | Heuristic; debated (Burgess 2011 IJE; Zhao 2020 argues 10 is too low for one-sample) |
| Conditional F > 10 per exposure (MVMR) | Sanderson 2019 IJE 48:713 | Total F can be high while conditional F low; per-exposure F is what matters |
| `I^2_GX >= 0.9` for Egger | Bowden 2016 IJE 45:1961 | NOME assumption; below this, SIMEX correction required |
| CAUSE >= 100 significant SNPs | Morrison 2020 Nat Genet 52:740 | Mixture model needs signal density for shared-factor estimation |
| Egger >= 10 instruments | Bowden 2015 IJE 44:512 | Power for slope test in weighted regression |
| Sample-overlap z-score correlation | Burgess 2016 Genet Epidemiol 40:597 | Use LDSC bivariate intercept as proxy; correct IVW SE accordingly |
| Clumping r2 < 0.001, 10 Mb window | TwoSampleMR default; matches GWAS LD norms | Polygenic MR; cis-MR uses r2 < 0.1 within window |
| Steiger p < 0.05 | Hemani 2017 PLoS Genet 13:e1007081 | Heuristic; subject to confounder caveat (Hemani Tilling 2022) |
| MR-PRESSO NbDistribution | 1000 exploratory; >= 5000 publication; >= 10000 stringent | Verbanck 2018 Nat Genet 50:693; precision of global p-value scales with NbDistribution |
| STROBE-MR all 20 items | Skrivankova 2021 JAMA 326:1614; BMJ 375:n2233 | Required since 2022 by most epidemiology journals |
| Bonferroni for pheWAS-MR | Standard | Many outcomes; FDR if exploratory |

## TwoSampleMR Standard Workflow

**Goal:** Produce a defensible primary IVW estimate plus a full sensitivity battery from two-sample summary statistics.

**Approach:** Extract genome-wide significant instruments -> clump (local plink preferred) -> extract outcome -> harmonise -> mr -> pleiotropy + heterogeneity + leave-one-out -> Steiger -> MR-PRESSO -> report.

```r
library(TwoSampleMR)
library(ieugwasr)

exposure_raw <- read_exposure_data(
    filename = 'exposure_gwas.tsv', sep = '\t',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2',
    eaf_col = 'EAF', pval_col = 'P'
)

exposure_sig <- subset(exposure_raw, pval.exposure < 5e-08)  # genome-wide significance

# F-statistic computed from EXPOSURE (Burgess 2011); ratio of squared effect to its variance
exposure_sig$f_stat <- (exposure_sig$beta.exposure / exposure_sig$se.exposure)^2
exposure_sig <- subset(exposure_sig, f_stat >= 10)  # Staiger-Stock 1997 weak-IV heuristic

clumped <- ld_clump(
    data.frame(rsid = exposure_sig$SNP, pval = exposure_sig$pval.exposure),
    clump_r2 = 0.001, clump_kb = 10000,  # polygenic MR convention
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = '1kg_EUR/EUR'
)
exposure_dat <- subset(exposure_sig, SNP %in% clumped$rsid)

outcome_dat <- read_outcome_data(
    filename = 'outcome_gwas.tsv', snps = exposure_dat$SNP, sep = '\t',
    snp_col = 'SNP', beta_col = 'BETA', se_col = 'SE',
    effect_allele_col = 'A1', other_allele_col = 'A2',
    eaf_col = 'EAF', pval_col = 'P'
)

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)  # infer from EAF; drops MAF~0.5 palindromes

primary <- mr(dat, method_list = c('mr_ivw', 'mr_egger_regression',
                                    'mr_weighted_median', 'mr_weighted_mode'))

heterogeneity <- mr_heterogeneity(dat)         # Cochran Q
pleiotropy <- mr_pleiotropy_test(dat)          # Egger intercept
loo <- mr_leaveoneout(dat)                     # influential-SNP check
steiger <- directionality_test(dat)            # variance-explained direction
```

## MR-PRESSO Outlier Detection

**Goal:** Detect horizontal-pleiotropy outliers, remove them, and test whether the corrected estimate differs from the uncorrected one (distortion test).

**Approach:** Three-step framework: global test (presence of pleiotropy), outlier test (per-SNP), distortion test (effect change after outlier removal).

```r
library(MRPRESSO)

presso <- mr_presso(
    BetaOutcome = 'beta.outcome', BetaExposure = 'beta.exposure',
    SdOutcome = 'se.outcome', SdExposure = 'se.exposure',
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
    data = dat, NbDistribution = 10000,  # >= 10000 for publication-grade p-value precision
    SignifThreshold = 0.05
)

print(presso$`MR-PRESSO results`$`Global Test`)         # any pleiotropy
print(presso$`MR-PRESSO results`$`Distortion Test`)     # change after outlier removal
outlier_snps <- which(presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05 / nrow(dat))
```

## CAUSE for Correlated Horizontal Pleiotropy

CAUSE (Morrison 2020 Nat Genet 52:740) fits a shared-factor mixture to genome-wide sumstats and compares causal vs sharing-only models by delta-ELPD. Workflow: `gwas_merge()` -> sample ~1M variants for `est_cause_params()` -> filter sig SNPs (P < 1e-3) and optionally LD-prune -> `cause(X, variants, param_ests)`. Needs >= 100 sig SNPs. Full annotated example and ELPD interpretation in causal-genomics/pleiotropy-detection.

## MVMR with Conditional F

**Goal:** Estimate the causal effect of exposure X1 on Y, adjusting for measured pleiotropy via X2.

**Approach:** Format exposures + outcome into MVMR object; compute conditional F per exposure (>10 required); run multivariable IVW; report Q_A heterogeneity.

```r
library(MVMR)

mvmr_dat <- format_mvmr(
    BXGs = cbind(dat$beta.x1, dat$beta.x2),
    BYG = dat$beta.y,
    seBXGs = cbind(dat$se.x1, dat$se.x2),
    seBYG = dat$se.y,
    RSID = dat$SNP
)

condF <- strength_mvmr(r_input = mvmr_dat, gencov = 0)  # per-exposure conditional F
# condF must be > 10 for EACH exposure (Sanderson 2019); total F is misleading

mv_ivw <- ivw_mvmr(r_input = mvmr_dat)
mv_qa <- qhet_mvmr(r_input = mvmr_dat, pcor = 0)  # Q_A heterogeneity
```

`gencov = 0` is valid ONLY if the exposure GWAS samples don't overlap; for overlapping exposures use the bivariate LDSC intercept matrix as `gencov`. If any conditional F < 10, the IVW point estimate is weak-IV-biased; switch to the Q-minimization estimator: `qhet_mvmr(r_input, pcor, CI = TRUE, iterations = 1000)` (Sanderson 2021 Stat Med 40:5434), which minimizes Q-statistic heterogeneity rather than weighting by inverse variance and is robust to weak conditional instruments.

## Bidirectional and Steiger

```r
exposure_rev <- format_data(outcome_raw, type = 'exposure')  # treat former outcome as exposure
outcome_rev <- format_data(exposure_raw, type = 'outcome')
dat_rev <- harmonise_data(exposure_rev, outcome_rev, action = 2)
results_rev <- mr(dat_rev, method_list = 'mr_ivw')

dat_filt <- steiger_filtering(dat)  # per-SNP; flags SNPs where variance(Y) > variance(X)
dir_test <- directionality_test(dat) # global; correct_causal_direction == TRUE if forward
```

**Report direction operationally:** forward p < 5e-8 with reverse p > 0.05; `directionality_test()` `correct_causal_direction == TRUE` with Steiger p < 0.05; point estimate in reverse direction has |effect| substantially smaller than forward (reflecting reverse-instrument-strength asymmetry). Example sentence: "Forward MR showed BMI -> T2D (IVW beta = 0.85, p = 2e-15); reverse MR was null (IVW beta = 0.02, p = 0.61); Steiger directionality test favored the forward direction (p = 3e-7)."

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| IVW sig, Egger null with non-zero intercept | Directional pleiotropy | Report Egger as primary if `I^2_GX >= 0.9`; otherwise SIMEX-correct |
| IVW sig, weighted median sig, mode null | Mode underpowered or bimodal pleiotropy | Trust the agreement of IVW + median |
| MR-PRESSO distortion test sig | Outliers materially shift estimate | Report PRESSO-adjusted estimate as primary |
| CAUSE sig, IVW sig, same direction | High-confidence causal claim | Report both; emphasize CAUSE rules out CHP |
| CAUSE null, IVW sig, same direction | CHP indistinguishable from causation | Downgrade to "consistent with causation but not separable from CHP" |
| Forward MR sig, reverse MR also sig | Bidirectional causation OR confounded | Use LHC-MR or LCV for genome-wide resolution |
| IVW sig but `mean F << 10` in one-sample design | Weak-IV bias toward observational | Re-run with MR-RAPS; report adjusted estimate |

**Operational rule for publication:** Primary IVW + concordant Egger (or weighted median if NOME violated) + non-significant MR-PRESSO global test + Steiger correct direction = publication-ready evidence. CAUSE concordance is required when the exposure is polygenic and the prior on CHP is high (e.g. BMI -> outcome, lipids -> outcome). Single-method "significant IVW" claims should be downgraded to exploratory.

## Cohort Gotchas

- UKB GWAS commonly include non-EUR participants; pan-UKB EUR/AFR/EAS/SAS subsets are separate releases. Mixing ancestries inflates instrument strength via population stratification rather than biology.
- FinnGen DF12 (2024) cohort: Finnish founder effects produce narrower LD blocks and higher winner's curse magnitude than UKB at matched sample size.
- MVP / GBMI / AoU: multi-ancestry meta-analyses; stratify by ancestry before MR or use ancestry-specific subsets.
- UKB-on-UKB MR creates one-sample-equivalent bias regardless of "different GWAS file" appearance; use MRlap or move the outcome to an external cohort (FinnGen, BBJ, MVP).

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Sample overlap between exposure and outcome GWAS?" | LDSC bivariate intercept reported; MRlap or sample-overlap-corrected IVW applied; document overlap fraction |
| "Weak instruments (F < 10)?" | F computed from EXPOSURE; per-instrument and mean F reported; MR-RAPS used as sensitivity if mean F borderline |
| "Horizontal pleiotropy?" | IVW + Egger + weighted median + weighted mode + MR-PRESSO; if rg > 0.3 also CAUSE (see pleiotropy-detection) |
| "Reverse causation?" | Steiger filter applied; bidirectional MR ran; LCV gcp reported if rg > 0.3 |
| "Pre-registered?" | OSF protocol filed; STROBE-MR all 20 items reported |
| "InSIDE assumption?" | INstrument Strength Independent of Direct Effect -- pleiotropic effects alpha uncorrelated with instrument-exposure effects gamma. Tested via Egger intercept + CHP-aware sensitivity (CAUSE) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| F-statistic computed from outcome | Reading off `beta.outcome / se.outcome` | Compute from EXPOSURE; outcome F is meaningless for IV strength |
| All instruments dropped at harmonise | All SNPs palindromic; or allele coding mismatch | Verify A1/A2 conventions match; try `action = 3` to drop, not assume |
| `Unauthorized` / `403` from OpenGWAS | OAuth deprecated May 2024; JWT token required | Generate at api.opengwas.io -> set `OPENGWAS_JWT=<token>` in `~/.Renviron` -> restart R -> verify with `ieugwasr::get_opengwas_jwt()`. For production, skip the API: use `ld_clump_local()` with local 1KG bfile (saves rate-limit + auth headaches). |
| Conditional F vs total F confusion in MVMR | Reporting `mean(F)` not `strength_mvmr()` per-exposure | Always use `MVMR::strength_mvmr()` and report each column |
| `install.packages("mr.raps")` fails | CRAN-archived 2025-03-01 | `remotes::install_github('qingyuanzhao/mr.raps')`; TwoSampleMR's `mr_raps()` wrapper internally calls this package |
| MR-PRESSO returns NA p-value | `NbDistribution` too small; signal too thin | Increase to >= 10000; check that >= 4 SNPs remain after harmonization |
| Egger intercept "highly significant" with 5 SNPs | Underpowered Egger over-fits the slope | Egger needs >= 10 SNPs; below that, intercept is unreliable |
| Sample-overlap correction ignored | Treating UKB-on-UKB as two-sample | Apply Burgess 2016 correction; or use MR-RAPS |
| `cause()` runs forever | Default model fit on too many SNPs | Filter to sig SNPs (P < 1e-3) before `cause()`; `est_cause_params` uses the random subset |
| MAF column missing -> harmonise action=2 silently downgrades | EAF unavailable | Provide EAF or use `action = 3` and document the loss |

## Tool Installation Notes

```r
# CRAN-stable
install.packages(c('remotes', 'MendelianRandomization', 'MVMR', 'coloc'))

# GitHub-only or recently archived
remotes::install_github('MRCIEU/TwoSampleMR')          # primary orchestrator
remotes::install_github('MRCIEU/ieugwasr')             # OpenGWAS client + local clumping
remotes::install_github('rondolab/MR-PRESSO')          # never on CRAN
remotes::install_github('qingyuanzhao/mr.raps')        # CRAN-archived 2025-03-01
remotes::install_github('jean997/cause')               # depends on mixsqp; suggests Rfast
remotes::install_github('cnfoley/mrclust')             # heterogeneity clusters
remotes::install_github('LizaDarrous/lhcMR')           # bidirectional + heritable confounder
remotes::install_github('HDTian/DRMR')                 # doubly-ranked stratification
remotes::install_github('n-mounier/MRlap')             # joint overlap + winner's-curse + weak-IV correction
```

`TwoSampleMR::mr_raps()` is a thin wrapper that calls `mr.raps::mr.raps()` under the hood; the GitHub `mr.raps` install above is therefore required. The `MendelianRandomization` package does NOT export `mr_raps()` (verify with `ls('package:MendelianRandomization')`); only TwoSampleMR offers a MR-RAPS entry point. For local clumping, install plink2 and download a 1KG EUR (or matched-ancestry) reference bfile.

## STROBE-MR Reporting

STROBE-MR (Skrivankova 2021 JAMA 326:1614; BMJ 375:n2233): 20-item checklist required by Eur J Epi / Nat Genet / JAMA / Diabetologia since 2022. See causal-genomics/pleiotropy-detection for the per-item table -- that skill is the consolidated reporting reference for the full MR + sensitivity battery.

## References

- Davey Smith G & Ebrahim S 2003 IJE 32:1 (foundational framework)
- Burgess S & Thompson SG 2017 SAGE textbook (MR canonical reference)
- Bowden J et al 2015 IJE 44:512 (MR-Egger)
- Bowden J et al 2016 Genet Epidemiol 40:304 (weighted median)
- Hartwig FP et al 2017 IJE 46:1985 (weighted mode)
- Zhao Q et al 2020 Ann Stat 48:1742 (MR-RAPS)
- Morrison J et al 2020 Nat Genet 52:740 (CAUSE)
- Zhu Z et al 2018 Nat Commun 9:224 (GSMR + HEIDI)
- Verbanck M et al 2018 Nat Genet 50:693 (MR-PRESSO)
- Sanderson E et al 2019 IJE 48:713 (MVMR conditional F)
- Foley CN et al 2020 Bioinformatics 37:531 (MR-Clust)
- Burgess S et al 2020 Nat Commun 11:376 (contamination mixture)
- O'Connor LJ & Price AL 2018 Nat Genet 50:1728 (LCV)
- Darrous L et al 2021 Nat Commun 12:7274 (LHC-MR)
- Tian H et al 2023 PLoS Genet 19:e1010823 (DRMR)
- Burgess S et al 2016 Genet Epidemiol 40:597 (sample-overlap correction)
- Schmidt AF et al 2020 Nat Commun 11:3255 (drug-target / cis-MR framework)
- Hemani G & Tilling K 2022 Wellcome Open Res 7:14 (Steiger filter caveat)
- Skrivankova VW et al 2021 JAMA 326:1614 (STROBE-MR statement)
- Mounier N & Kutalik Z 2023 Genet Epidemiol 47:314 (MRlap joint correction)
- Sanderson E et al 2021 Stat Med 40:5434 (qhet_mvmr Q-minimization estimator)
- Burgess S & Labrecque JA 2018 Eur J Epidemiol 33:947 (binary outcomes non-collapsibility)
- Hu Y et al 2022 IJE 51:1289 (collider bias in case-only / progression cohorts)
- Bowden J et al 2019 IJE 48:728 (one-sample MR bias formulas)
- Hamilton FW et al 2023 medRxiv 23293658 (NLMR stratum-specific bias critique)
- Sadreev II et al 2024 IJE 52:1209 (winner's curse empirical magnitude)
- Ma S et al 2023 JASA (RIVW rerandomized IVW)

## Related Skills

- causal-genomics/colocalization-analysis - Confirm shared causal variant for cis-MR drug-target work
- causal-genomics/pleiotropy-detection - Deep dive on MR-PRESSO, Egger, contamination mixture diagnostics
- causal-genomics/fine-mapping - Credible-set construction at instrument loci before cis-MR
- causal-genomics/mediation-analysis - Two-step MR and MVMR difference method for X -> M -> Y mediation
- causal-genomics/transcriptome-wide-association - TWAS as MR-adjacent framework for gene-level inference
- causal-genomics/proteome-mr-drug-target - Dedicated cis-pQTL drug-target MR workflow with UKB-PPP/deCODE/Fenland and coloc triangulation
- population-genetics/association-testing - GWAS source for instrument selection
- clinical-databases/clinvar-lookup - Annotate instrument SNPs for downstream interpretation
