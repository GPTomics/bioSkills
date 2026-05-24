---
name: bio-causal-genomics-pleiotropy-detection
description: Detect and adjust for horizontal pleiotropy in two-sample Mendelian randomization by distinguishing uncorrelated (UHP) from correlated (CHP) pleiotropy and choosing among Egger, MR-PRESSO, MR-RAPS, CAUSE, LHC-MR, LCV, MR-Clust, MR-Mix, and contamination-mixture methods. Use when validating an MR causal claim, running the STROBE-MR sensitivity battery, suspecting a shared heritable confounder, working under weak-instrument or polygenic-exposure regimes, or reconciling discordant estimates across robust methods.
tool_type: r
primary_tool: TwoSampleMR
---

## Version Compatibility

Reference examples tested with: TwoSampleMR 0.5.11+, MendelianRandomization 0.9.0+, MR-PRESSO 1.0+, CAUSE 1.2.0+, MR-Clust 0.1.0+, MRMix 0.1+, mr.raps 0.4.1+ (GitHub), LHC-MR 0.0.0.9000+ (GitHub), LCV (script-based, no version tag), simex 1.8+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- For GitHub-only packages, check the repo HEAD vs the local install date

If code throws errors, introspect the installed package and adapt the example rather than retrying.

# Pleiotropy Detection in Mendelian Randomization

**"Validate my MR result against pleiotropic bias"** -> Decompose violations of the exclusion-restriction assumption into uncorrelated horizontal pleiotropy (UHP, addressable by Egger / median / mode / MR-PRESSO) and correlated horizontal pleiotropy (CHP, addressable only by CAUSE / LHC-MR / LCV), then run a method battery whose assumptions span both regimes.

- R: `TwoSampleMR::mr()` (IVW + Egger + median + mode), `mr_pleiotropy_test()`, `mr_heterogeneity()`, `mr_leaveoneout()`, `directionality_test()`
- R: `MRPRESSO::mr_presso()` for UHP outlier removal + distortion test
- R: `cause::cause()` for CHP-aware estimation; `mrclust::mr_clust_em()` for mechanism-heterogeneous instruments
- R: `MendelianRandomization::mr_conmix()` for contamination mixture; `MRMix::MRMix()` for mixture-of-distributions

## UHP vs CHP: The Central Postdoc-Grade Distinction

Horizontal pleiotropy comes in two regimes, and most "standard" MR sensitivity methods address only one of them.

| Regime | Definition | InSIDE assumption | Methods that handle it |
|--------|------------|--------------------|------------------------|
| UHP (uncorrelated horizontal pleiotropy) | Pleiotropic effect alpha_j independent of instrument-exposure effect gamma_j | Holds | IVW (balanced UHP only), MR-Egger, weighted median, weighted mode, MR-PRESSO, MR-RAPS, MR-Mix, contamination mixture |
| CHP (correlated horizontal pleiotropy) | alpha_j correlates with gamma_j through a shared upstream factor (heritable confounder, network mediator) | Violated | CAUSE, LHC-MR, LCV, MR-Clust (partial), Steiger-filtered MR (partial) |

**InSIDE = INstrument Strength Independent of Direct Effect** (Bowden 2015 IJE 44:512). Plain English: across SNPs, the per-SNP pleiotropic effect alpha and per-SNP instrument-exposure effect gamma are treated as independent random variables. CHP is the case where they covary because both flow from a shared upstream genetic factor.

**The trap (Morrison 2020 Nat Genet 52:740):** IVW, MR-Egger, MR-PRESSO, and GSMR are all blind to CHP. Under a shared heritable confounder they each return a plausible-looking corrected causal estimate that is systematically biased in the direction of the confounder. The MR-PRESSO global test does not flag CHP because correlated pleiotropy is not an outlier pattern, it is a population mean shift in the alpha distribution conditional on gamma.

**Operational rule:** If genetic correlation rg(exposure, outcome) is high (LDSC `>= 0.3`) or biology strongly suggests a shared upstream factor, the IVW / Egger / PRESSO triple is insufficient. Add CAUSE (preferred when sig SNPs `>= 100`) or LHC-MR (preferred for polygenic genome-wide IVs).

## Operational Decision Flow (4 Steps)

1. **Compute genetic correlation (LDSC).** Run `ldsc.py --rg <exposure.sumstats.gz>,<outcome.sumstats.gz>` (see causal-genomics/genetic-correlation). If `|rg| > 0.3`, CHP is plausible -> flag for Step 3 escalation. If the LDSC rg standard error spans zero broadly, treat low-rg evidence as weak rather than confirming absence of CHP.
2. **Standard battery.** IVW (random-effects when Cochran Q p < 0.05) + MR-Egger (with NOME I^2_GX check) + weighted median + weighted mode + MR-PRESSO (NbDistribution `>= 10000` for stringent reporting). Report all five with point estimate, SE, p, 95% CI, and n_SNPs_used. Compute Egger I^2_GX; apply SIMEX if I^2_GX < 0.9 (see examples/simex_egger_correction.R).
3. **CHP escalation.** Trigger when rg > 0.3 OR PRESSO global p < 0.05 with > 50% nominal outliers OR Egger / median / mode disagree by > 2 SE. Run CAUSE (if `>= 100` significant SNPs after pruning) or LHC-MR (any N; uses genome-wide sumstats). Report ELPD delta + z + q (CHP fraction) + gamma (CHP-adjusted causal estimate).
4. **Triangulate.** Pre-MR Steiger filter; bidirectional MR (examples/bidirectional_mr.R); LCV gcp; LDSC rg report. Consensus across methods supports a publication-ready claim. Disagreement requires narrowing the scope (e.g., subgroup, cis-MR, time-varying analysis) rather than reporting a single point estimate.

## Algorithmic Taxonomy

| Method | Models | UHP-robust | CHP-robust | Min #SNPs | Fails when | Citation |
|--------|--------|-----------|-----------|-----------|------------|----------|
| Inverse-variance weighted (IVW) | Weighted regression through origin | Balanced UHP only | No | 3 | Directional UHP; CHP; weak IV bias; heterogeneity | Burgess 2013 Genet Epidemiol 37:658 |
| MR-Egger intercept + slope | IVW + free intercept | Directional UHP | No | >=10 for power | NOME violated (I^2_GX < 0.9); <10 SNPs; CHP | Bowden 2015 IJE 44:512 |
| Weighted median | Median of Wald ratios | Up to 50% invalid | No | >=10 | >50% invalid; CHP | Bowden 2016 Genet Epidemiol 40:304 |
| Weighted mode (MBE) | Mode of estimate density | Plurality valid | Partial | >=10 | Multimodal estimates from CHP clusters | Hartwig 2017 IJE 46:1985 |
| Cochran Q | Heterogeneity across Wald ratios | Total heterogeneity flag, not direction-specific | No | 3 | Cannot distinguish UHP from heterogeneity from CHP | Greco 2015 Stat Med 34:2926 |
| MR-PRESSO | Detect + remove UHP outliers via RSS-out | Yes (assumes majority valid) | No | >=4 | >50% pleiotropic; any CHP; small n | Verbanck 2018 Nat Genet 50:693 |
| GSMR + HEIDI-outlier | Outlier removal via single-instrument estimate heterogeneity | Yes | No | >=10 | CHP (HEIDI-outlier is heterogeneity-driven) | Zhu 2018 Nat Commun 9:224 |
| MR-RAPS | Profile likelihood with overdispersion + Huber/Tukey loss | Yes; weak-IV robust | Partial via overdispersion | >=10 | Strong CHP | Zhao 2020 Ann Stat 48:1742 |
| MR-Mix | Mixture-of-distributions over valid + invalid | Yes | Partial | >=20 | Few SNPs; very heterogeneous CHP | Qi & Chatterjee 2019 Nat Commun 10:1941 |
| Contamination mixture | Profile likelihood over contamination fraction | Yes | Partial | >=20 | Few SNPs | Burgess 2020 Nat Commun 11:376 |
| MR-Clust | k-means over Wald estimates with NULL cluster | Yes | Diagnostic for CHP via clusters | >=20 | Single-mechanism exposure (no clustering signal) | Foley 2020 Bioinformatics 37:531 |
| CAUSE | Bayesian mixture: shared causal + shared-factor (CHP) components | Yes | Yes (explicit) | >=100 sig SNPs at p<5e-8 | <100 sig SNPs; non-overlapping GWAS samples | Morrison 2020 Nat Genet 52:740 |
| LHC-MR | Latent heritable confounder + bidirectional + heritability | Yes | Yes | Genome-wide GWAS sumstats | Heritability mis-estimated; severe sample overlap | Darrous 2021 Nat Commun 12:7274 |
| LCV (latent causal variable) | gcp parameter on genome-wide rg | N/A (not an MR method) | Diagnostic | Genome-wide GWAS sumstats | Heritability low; non-Gaussian effect distribution | O'Connor & Price 2018 Nat Genet 50:1728 |

Methodology evolves; verify against the Burgess & Thompson textbook (2nd ed 2021), Hemani 2018 (basic four-method battery), and Sanderson 2022 Nat Rev Methods Primers 2:6 before locking a sensitivity battery.

## Decision Tree by Scenario

| Scenario | Primary estimator | Sensitivity / triangulation |
|----------|-------------------|----------------------------|
| Many strong IVs, no biological shared trait suspected | IVW + Egger + weighted median + weighted mode + PRESSO | Cochran Q; leave-one-out; F-stat; Steiger filtering |
| LDSC rg(exposure, outcome) `>= 0.3` or strong shared-factor biology | CAUSE (if sig SNPs `>= 100`) OR LHC-MR | LCV gcp; cross-check IVW after Steiger filter |
| Many weak IVs (mean F < 20) | MR-RAPS with overdispersion + robust Huber loss | MR-Mix or contamination mixture; report F-stat range |
| Suspected heterogeneous causal mechanisms (e.g. LDL on CHD via multiple lipoprotein pathways) | MR-Clust; report per-cluster IVW | Pathway annotation of cluster instruments; Bayesian mixture |
| Cis-MR drug target (single locus, few SNPs in LD) | Colocalization (causal-genomics/colocalization-analysis) + Steiger | PWCoCo; conditional analysis; not Egger (low SNP count) |
| Polygenic exposure (heritability spread genome-wide; few significant loci) | LHC-MR (uses all SNPs) | LDSC rg; genome-wide IVW with weak-IV-aware methods (RAPS) |
| Reverse causation suspected | Bidirectional MR with Steiger filter; LHC-MR (jointly estimates both directions) | directionality_test; effect-size r2 comparison |
| Population-level summary discordant with biology | Re-examine instrument selection; check Winner's curse; LD pruning settings | Triangulate with cis-MR; family-based MR if available |

## Per-Method Failure Modes

### MR-Egger NOME violation

**Trigger:** I^2_GX = (Q_GX - df) / Q_GX is below 0.9, indicating measurement-error attenuation of the Egger slope (NOME = "no measurement error" in the exposure GWAS effect sizes).

**Mechanism:** MR-Egger regresses outcome effects on exposure effects with a free intercept. Imprecise exposure effects (high beta.exposure SE relative to beta.exposure variability across instruments) introduce regression dilution that pulls the Egger slope toward the null and inflates the intercept.

**Symptom:** Egger slope much closer to zero than IVW, weighted median, and weighted mode estimates; large Egger SE.

**Fix:** Apply SIMEX correction (Bowden 2016 IJE 45:1961; Cook & Stefanski 1994 JASA 89:1314 SIMEX framework) using the `simex` package on the Egger regression, treating beta.exposure SE as measurement error. See examples/simex_egger_correction.R. Alternative: use MR-RAPS, which models the exposure-effect error explicitly via profile likelihood and does not suffer the NOME failure.

### MR-PRESSO majority-outlier breakdown

**Trigger:** More than 50% of instruments are pleiotropic (UHP), e.g. when instrument set was loosely selected (genome-wide significant but unfiltered).

**Mechanism:** MR-PRESSO's global RSS-out statistic and outlier detection both assume a majority-valid set; outliers are defined relative to that majority. With a pleiotropic majority, PRESSO removes the valid minority.

**Symptom:** PRESSO-corrected estimate is similar in magnitude (and sign) to the uncorrected estimate even after dropping nominally "outlier" SNPs; distortion-test p-value paradoxically non-significant; few or no outliers detected despite obvious global-test significance.

**Fix:** Do not trust PRESSO corrected estimate. Re-examine instrument selection (drop loose p-thresholds, prune LD harder); switch to CAUSE or LHC-MR; consider weighted-mode estimator which is plurality-valid rather than majority-valid.

### MR-PRESSO false negative under CHP

**Trigger:** Strong shared heritable confounder (high rg) producing CHP. Confirmed by significant LDSC rg or LCV gcp.

**Mechanism:** Correlated pleiotropy is a population-level mean shift in alpha conditional on gamma; it is not an outlier pattern. PRESSO's RSS-out distance is invariant under such a mean shift, so the global test is not powered against CHP.

**Symptom:** PRESSO global p > 0.05 (no detected pleiotropy) while a CHP-aware method (CAUSE, LHC-MR) returns a substantially different (often null) causal estimate.

**Fix:** When CHP is plausible, ALWAYS run CAUSE or LHC-MR in addition to PRESSO; do not rely on PRESSO global non-significance as evidence of no pleiotropy.

### MR-Egger underpowered with few SNPs

**Trigger:** Fewer than 10 instruments.

**Mechanism:** Egger's intercept variance is driven by the spread of beta.exposure across instruments; with few SNPs the intercept CI is so wide that even strongly pleiotropic data give non-significant intercepts.

**Symptom:** Non-significant Egger intercept p-value alongside obviously discordant IVW and weighted-median estimates.

**Fix:** Report intercept point estimate and CI rather than a binary "pleiotropy present / absent" verdict; do not use Egger as the only sensitivity method when SNP count is low; weight evidence toward weighted-median, weighted-mode, and CAUSE / LHC-MR.

### Steiger filter inverted by exposure measurement error (Hemani & Tilling 2022 IJE)

**Trigger:** Exposure is imprecisely measured (lower heritability ascertained in the exposure GWAS) and outcome is well-measured.

**Mechanism:** Steiger compares r^2_GX vs r^2_GY per SNP. Measurement error in the exposure underestimates r^2_GX; well-measured outcome captures r^2_GY accurately. Per-SNP, the inequality can flip even when the true causal direction is exposure -> outcome.

**Symptom:** A large fraction of instruments fail Steiger (`steiger_dir == FALSE`) in a direction that conflicts with biological plausibility.

**Fix:** Interpret Steiger as one signal among many, not a hard gate; cross-check with bidirectional MR; verify exposure GWAS heritability and sample size; switch to LHC-MR which models both directions jointly and accounts for heritability.

### CAUSE underpowered with few significant SNPs

**Trigger:** Fewer than 100 genome-wide-significant instruments (p < 5e-8) after harmonization and LD pruning.

**Mechanism:** CAUSE fits a Bayesian mixture model over a shared-factor (CHP) component, a shared-causal component, and a null component. Posterior identification of the mixture weights requires substantial signal across many SNPs.

**Symptom:** CAUSE delta_ELPD CI crosses zero; Pareto-k diagnostic flags unstable points; posterior intervals on q (CHP fraction) span [0, 1].

**Fix:** Use LCV gcp for genome-wide directional inference (does not require many significant SNPs); use LHC-MR if heritability and sumstats are available; or report CAUSE alongside an explicit caveat about its underpowered regime.

### LCV gcp under non-Gaussian effect distributions

**Trigger:** Highly polygenic trait with substantial sparsity in true effects (mixture of large-effect and zero-effect loci).

**Mechanism:** LCV assumes a bivariate normal model for effect sizes after LDSC adjustment. Sparse architectures (e.g. immune traits with HLA dominance) violate this and bias gcp estimates.

**Symptom:** LCV gcp point estimate appears extreme but heritability LDSC z-scores are modest; partitioned heritability shows extreme HLA enrichment.

**Fix:** Exclude HLA region from LDSC inputs; complement with CAUSE / LHC-MR; report gcp with awareness of the polygenicity caveat.

## Quantitative Thresholds

| Metric | Threshold | Source / rationale |
|--------|-----------|--------------------|
| F-statistic per SNP | `>=10` strong; `<10` weak | Burgess 2011 IJE 40:755 (rule of thumb); weak-IV bias toward observational confounded estimate |
| I^2_GX (NOME) | `>=0.9` Egger reliable | Bowden 2016 IJE 45:1961 |
| I^2_GX (NOME) intermediate | 0.6-0.9 SIMEX-corrected Egger | Bowden 2016 IJE |
| I^2_GX (NOME) severe | `<0.6` drop Egger; use MR-RAPS or CAUSE | Bowden 2016 IJE; SIMEX unreliable below 0.6 |
| Egger min SNP count for adequate power | `>=10` | Bowden 2015 IJE 44:512 |
| Cochran Q significance | p < 0.05 indicates heterogeneity | Greco 2015 Stat Med 34:2926 |
| MR-PRESSO NbDistribution | 1000 exploratory; `>=5000` publication; `>=10000` stringent | Verbanck 2018 Nat Genet 50:693 (Methods) |
| MR-PRESSO global test p | < 0.05 -> heterogeneity / outliers present | Verbanck 2018 Nat Genet |
| MR-PRESSO distortion test p | < 0.05 -> outliers materially shifted estimate; if `>= 0.05` report uncorrected IVW | Verbanck 2018 Nat Genet |
| MR-PRESSO min instruments | `>=4` to run; `>=10` for non-degenerate global test | Verbanck 2018 Nat Genet |
| MR-PRESSO SignifThreshold | 0.05 default | Verbanck 2018 Nat Genet |
| MR-PRESSO majority-valid breakdown | Fails when `>50%` instruments pleiotropic | Verbanck 2018 Nat Genet (theoretical limit) |
| Weighted median validity | Robust to `<=50%` invalid IVs | Bowden 2016 Genet Epidemiol 40:304 |
| Weighted mode validity | Plurality-valid (largest valid subset is most common estimate) | Hartwig 2017 IJE 46:1985 |
| CAUSE min #SNPs | `>=100` p < 5e-8 SNPs after pruning | Morrison 2020 Nat Genet 52:740 (Supplement) |
| CAUSE delta_ELPD criterion | one-sided p < 0.05; z = delta_elpd / se(delta_elpd); z > 1.96 standard; z > 3.0 stringent | Morrison 2020 Nat Genet |
| LDSC rg suggesting CHP | `>= 0.3` flags need for CAUSE / LHC-MR | Operational rule; see causal-genomics/genetic-correlation |
| Steiger r^2 difference | Reverse-causal flag at any per-SNP r2_GY > r2_GX | Hemani 2017 PLoS Genet 13:e1007081 |
| Standard sensitivity battery | IVW + Egger + median + mode + PRESSO + Steiger + LOO | Hemani 2018 eLife 7:e34408 / STROBE-MR 2021 |

LCV gcp interpretation thresholds (0, 0.5, 0.6, 1) are tabulated in usage-guide.md.

## Standard Sensitivity Battery (Working Reference)

**Goal:** Run the canonical UHP-focused MR sensitivity suite on harmonized two-sample data.

**Approach:** Compute IVW + Egger + median + mode side-by-side; test Egger intercept and heterogeneity; run MR-PRESSO with `>=5000` distributions for publication or `>=10000` for stringent reporting; apply Steiger filter; leave-one-out; report all estimates.

```r
library(TwoSampleMR)
library(MRPRESSO)

methods <- c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median', 'mr_weighted_mode')
res_mr <- mr(dat, method_list = methods)

het <- mr_heterogeneity(dat)
pleio <- mr_pleiotropy_test(dat)
loo <- mr_leaveoneout(dat)
steiger <- directionality_test(dat)

isq <- Isq(dat$beta.exposure, dat$se.exposure)
nome_pass <- isq >= 0.9

presso <- mr_presso(
    BetaOutcome='beta.outcome', BetaExposure='beta.exposure',
    SdOutcome='se.outcome', SdExposure='se.exposure',
    OUTLIERtest=TRUE, DISTORTIONtest=TRUE,
    data=dat, NbDistribution=10000, SignifThreshold=0.05)

global_p <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
outlier_p <- presso$`MR-PRESSO results`$`Outlier Test`$Pvalue
distortion_p <- presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
n_outliers <- sum(outlier_p < 0.05, na.rm=TRUE)
```

Full working pipeline incl SIMEX, MR-RAPS, contamination mixture, and STROBE-MR table: examples/sensitivity_battery.R.

## CAUSE for CHP-Aware Estimation

**Goal:** Distinguish causal from shared-factor (correlated horizontal pleiotropy) explanations of an exposure-outcome association.

**Approach:** Fit nuisance parameters (LD pruning + rho_GWAS sample-overlap correction) on a random SNP set; fit the sharing and causal posterior; compare ELPD (expected log predictive density) via Pareto-k smoothed importance sampling.

```r
library(cause)
params <- est_cause_params(dat_cause, variants = pruned_subset_snps)
res_cause <- cause(X = dat_cause, variants = pruned_snps, param_ests = params)
elpd <- summary(res_cause)$tab
```

Full posterior extraction + reporting: examples/cause_analysis.R.

**Interpreting CAUSE output:**

- `q`: posterior CHP fraction; 0 = no CHP, 1 = all instruments operate via the shared factor
- `eta`: shared-factor effect on Y (the "confounder pathway" magnitude)
- `gamma`: posterior causal effect of E on Y after partialling out CHP; report median + 95% credible interval
- `delta_ELPD` (sharing - causal): negative -> causal model preferred; z = delta_elpd / se(delta_elpd); z > 1.96 standard, z > 3.0 stringent; one-sided p reported alongside posterior gamma
- Pareto-k > 0.7 indicates unstable posterior on those points; if more than 10% of points are unstable, treat the posterior as unreliable; remediation: add more SNPs (loosen p-threshold one notch then re-prune in LD) or re-fit excluding flagged outliers

CAUSE requires sumstats from both exposure and outcome GWAS in matched effect-allele coding. The pruning step typically retains 100-5000 signature SNPs at LD r^2 < 0.01 in a 1 Mb window; nuisance estimation should use a larger random SNP subset (`>= 100,000` genome-wide SNPs) to fit rho (sample overlap) stably.

## MR-RAPS Loss Function and Overdispersion

**Trigger:** Weak instruments (mean F < 20) and/or suspected UHP requiring outlier-resistant estimation.

- `over.dispersion = TRUE` always for MR (horizontal-pleiotropy variance is real, not noise; turning this off underestimates SE)
- `loss.function = 'huber'` (default; outlier-resistant; suited to mild to moderate UHP)
- `loss.function = 'tukey'` (more aggressive; downweights extreme outliers more; choose when many obvious outliers suspected)
- `loss.function = 'l2'` (non-robust; equivalent to weighted least squares; do not use when UHP suspected)

Tukey is preferable when leave-one-out reveals 2+ SNPs single-handedly shifting the IVW estimate by > 1 SE.

## MR-Clust for Mechanism Heterogeneity

**Goal:** When a single causal estimate is misleading because instruments operate through multiple causal mechanisms (e.g. LDL on CHD via multiple lipoprotein subfractions), identify clusters of instruments with similar per-SNP Wald ratios.

```r
library(mrclust)
ratio_hat <- dat$beta.outcome / dat$beta.exposure
ratio_se <- abs(dat$se.outcome / dat$beta.exposure)
res_mc <- mr_clust_em(theta=ratio_hat, theta_se=ratio_se,
                     bx=dat$beta.exposure, by=dat$beta.outcome,
                     bxse=dat$se.exposure, byse=dat$se.outcome,
                     obs_names=dat$SNP)
per_cluster <- res_mc$results$best
```

Clusters with cluster_class = 'null' are pleiotropy-only instruments. Per-cluster IVW estimates may differ substantially; biological annotation of the SNPs in each cluster (pathway, target gene) is the interpretation step.

## LHC-MR Workflow

**Goal:** Jointly estimate forward causal effect, reverse causal effect, and the heritable-confounder contribution from genome-wide sumstats (not just significant SNPs).

```r
library(lhcMR)
ld <- list(ld_score_path='ldsc/eur_w_ld_chr/', hapmap_snp_path='ldsc/w_hm3.snplist')
sp_list <- sp_estimate(trait.df=list(X=df_x, Y=df_y), ld=ld)
res_lhc <- lhc_mr(sp_list, sample_overlap=overlap_matrix, n_cores=4)
```

LHC-MR is computationally heavy (hours on full sumstats) but among the most rigorous CHP-aware estimators when both GWAS are well-powered. Output includes axx, ayy, hxy (confounder effect on each trait), and bidirectional alpha_xy, alpha_yx.

**Choosing CAUSE vs LHC-MR (Darrous 2021):**

| Condition | Preferred method |
|-----------|------------------|
| `>= 100` genome-wide significant SNPs after pruning | CAUSE (Bayesian; CHP-explicit; mature posterior diagnostics) |
| Polygenic exposure with few significant loci | LHC-MR (uses genome-wide signal, not just significant SNPs) |
| Severe sample overlap between exposure and outcome GWAS | LHC-MR (jointly models overlap); CAUSE's rho correction is exposed to misspecification at high overlap |
| Bidirectionality of central interest | LHC-MR (jointly estimates alpha_xy and alpha_yx); CAUSE only models forward |
| Limited compute / quick turnaround | CAUSE (minutes to hours); LHC-MR may be > 24h on full sumstats |

When both apply, report both with the agreement / disagreement explicit in the discussion.

## Bidirectional MR Procedure

1. **Forward MR:** instrument exposure E, test effect on outcome Y (primary)
2. **Reverse MR:** instrument outcome Y, test effect on exposure E (using outcome-direction instruments)
3. **Steiger pre-filter both directions:** `steiger_filtering(dat)`; drop SNPs where outcome r^2 > exposure r^2 before primary IVW
4. **Compare estimates:** null reverse + significant forward strengthens the forward causal claim; bidirectional significance flags feedback / shared confounder / reciprocal causation
5. **LCV gcp orthogonal check:** genome-wide directional inference independent of the instrument set

Working code: examples/bidirectional_mr.R.

**Interpretation cheat-sheet:**

| Forward p | Reverse p | Reading |
|-----------|-----------|---------|
| significant | non-significant | Forward causal claim strengthened |
| non-significant | significant | Re-examine instrument-exposure assignment; the "outcome" may causally drive the "exposure" |
| significant | significant | Feedback loop, shared confounder, or reciprocal causation; resolve with LHC-MR |
| non-significant | non-significant | No evidence of causation in either direction |

When forward and reverse both clear Steiger and both IVW p < 0.05, run LHC-MR jointly rather than reporting two univariable estimates.

## LCV (Latent Causal Variable)

LCV uses LDSC-merged genome-wide sumstats and reports gcp (genetic causality proportion) on [-1, 1]. It is a complement to, not a replacement for, MR; gcp ~ 0 with high LDSC rg implies pure genetic correlation without partial causation.

```r
source('LCV/R/RunLCV.R')
res_lcv <- RunLCV(ldscores$L2, x$Z, y$Z)
# res_lcv$gcp; res_lcv$pval.gcpzero.2tailed
```

Full gcp interpretation table is in usage-guide.md.

## Required Supplementary Tables

**Instrument table** (one row per SNP retained for primary analysis):

| Column | Content |
|--------|---------|
| rsID, chr, pos | Variant identifier and genome position |
| EA, OA, EAF | Effect allele, other allele, effect allele frequency in exposure GWAS |
| beta_E, se_E, p_E, F | Exposure-side estimate, SE, p-value, per-SNP F-statistic |
| beta_Y, se_Y, p_Y | Outcome-side estimate, SE, p-value (harmonized to EA) |
| harmonization_action | 1=kept, 2=flipped, 3=dropped (palindromic ambiguity) |
| palindromic_flag | TRUE / FALSE; tracked per Hartwig 2016 IJE 45:1717 |
| Steiger_direction | forward / reverse / inconclusive |
| Steiger_p | per-SNP directionality p-value |

**Sensitivity-battery table** (one row per method):

| Column | Content |
|--------|---------|
| method | IVW / Egger / WM / WMode / PRESSO (raw + corrected) / RAPS / CAUSE |
| estimate, se, p, 95% CI | Point estimate and inference |
| n_SNPs_used | Post-harmonization, post-Steiger SNP count |
| heterogeneity_p | Q for IVW; Q' for Egger; global p for PRESSO |
| intercept_p | Egger only (directional UHP test) |
| ELPD_delta + z | CAUSE only (sharing - causal; negative + |z| > 1.96 -> causal) |

## Reconciliation Across Methods

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| IVW and Egger agree (small Egger intercept); median and mode agree | Likely true causal; minimal pleiotropy | Report all; STROBE-MR; emphasize agreement |
| IVW significant; Egger non-significant with similar slope; PRESSO global p > 0.05 | Egger underpowered (few SNPs) OR Egger NOME violated | Check I^2_GX; SIMEX-correct if 0.6 < I^2_GX < 0.9 |
| IVW shifted relative to Egger / median / mode; Egger intercept significant | Directional UHP | Trust Egger slope; PRESSO-corrected IVW; mode estimator |
| IVW + Egger + median + mode + PRESSO all agree but CAUSE delta_ELPD non-significant or in opposite direction | CHP via shared confounder | Trust CAUSE; report all five UHP methods alongside but flag the discordance; check LDSC rg |
| All UHP methods agree; LCV gcp ~ 0 | Genetic correlation only, not causation | Causation evidence is weak; report rg explicitly; consider colocalization (cis-MR) instead |
| MR-Clust shows >=2 distinct non-null clusters | Heterogeneous mechanisms | Report per-cluster estimates; do not summarize as a single effect |
| Steiger fails on a substantial fraction of instruments | Reverse causation OR exposure measurement error | Run bidirectional MR; check exposure GWAS heritability; LHC-MR |
| MR-PRESSO global p < 0.05 but corrected estimate similar to uncorrected | >50% pleiotropic OR CHP masquerading as UHP | Re-prune instruments; switch to weighted-mode / CAUSE / LHC-MR |

**Operational rule for publication:** Report IVW (primary), Egger slope + intercept, weighted median, weighted mode, MR-PRESSO global + distortion + corrected, Cochran Q, Steiger directionality, F-statistic distribution, I^2_GX (for Egger validity), and at least one CHP-aware method (CAUSE or LHC-MR) when rg `>= 0.3` or biology suggests shared upstream. Failure to report a CHP-aware result when CHP is plausible is a reviewer-flagged red flag since 2020.

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Did you check for CHP?" | LDSC rg reported (causal-genomics/genetic-correlation); if rg > 0.3, CAUSE or LHC-MR ran; q posterior reported |
| "Why CAUSE and not LHC-MR?" | CAUSE preferred when `>= 100` significant SNPs available (Morrison 2020). LHC-MR preferred when significant-SNP set is small or polygenic, using genome-wide sumstats (Darrous 2021) |
| "Egger NOME?" | I^2_GX computed; if 0.6 <= I^2_GX < 0.9, SIMEX correction applied; if < 0.6, Egger dropped in favor of MR-RAPS |
| "PRESSO doesn't catch CHP?" | Confirmed (Morrison 2020); CAUSE / LHC-MR reported alongside PRESSO for that reason |
| "Steiger filter applied pre-MR or post?" | Pre-MR: SNPs failing per-SNP Steiger directionality dropped before primary IVW |
| "Why no replication cohort?" | Two-sample design uses independent exposure and outcome cohorts; if same biobank, MRlap used or noted as a limitation |
| "Why not just trust the IVW?" | IVW assumes balanced UHP and no CHP; both violated routinely; sensitivity battery is the standard since STROBE-MR 2021 |
| "Effect size is implausibly large" | Re-examine F-statistic distribution for weak IV bias; check Winner's curse; consider Wald ratio at a single strong instrument as sanity check |

## STROBE-MR Reporting (Skrivankova 2021)

| Item | Required content |
|------|-----------------|
| 1-3 | Title / abstract / background indicates this is an MR study; pre-registered protocol |
| 4-7 | Study design, data sources, instrument selection criteria (p-threshold, LD clumping, MAF) |
| 8-11 | Harmonization, palindrome handling, allele alignment |
| 12-14 | F-statistic distribution; weak-instrument bias mitigation |
| 15-17 | Primary MR method + all sensitivity methods + CHP-aware method when relevant |
| 18-19 | Pleiotropy tests, Steiger, heterogeneity |
| 20 | Discussion of remaining assumption violations; limitations |

Sub-items (30 total) detail per-method reporting. The full statement (JAMA 326:1614) and explanation (BMJ 375:n2233) are now reviewer-required at most cardiovascular and psychiatric journals since 2022.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| MR-PRESSO crashes with `Not enough intrumental variables` | Fewer than 4 SNPs | Need >=4 for PRESSO; for cis-MR with few SNPs use colocalization |
| Egger intercept p < 0.05 but I^2_GX = 0.5 | NOME violated; intercept is artifactually inflated | SIMEX-correct or do not trust Egger; use MR-RAPS instead |
| `Isq()` not found | TwoSampleMR version where Isq is unexported | Compute manually: Q_GX = sum((beta_GX/se_GX)^2); I2 = (Q_GX - (n-1))/Q_GX, clipped to [0,1] |
| MR-RAPS `package not found` after CRAN install | CRAN-archived 2025-03-01 | `remotes::install_github('qingyuanzhao/mr.raps')`; call via `TwoSampleMR::mr_raps()` wrapper (MendelianRandomization does NOT export `mr_raps`) |
| CAUSE delta_ELPD CI spans zero; Pareto-k > 0.7 | <100 sig SNPs OR severe sample overlap | Use LHC-MR; or report CAUSE with the explicit caveat |
| Steiger labels most instruments reverse-causal | Exposure GWAS imprecise OR sample size mismatch | Treat as one signal; cross-check with bidirectional MR |
| LHC-MR runtime > 24h | Default n_cores=1 on full sumstats | Use n_cores >= 4; restrict to LDSC-overlapping SNPs first |
| MR-PRESSO outliers all on same chromosome | Genome-wide LD not properly pruned; clumping window too narrow | Re-clump at r^2 < 0.001 in 10 Mb window |
| MR-Mix returns NA | Few SNPs OR no variation in mixture support | Need >=20 SNPs; default mixture grid may need tuning |

## References

- Bowden J et al 2015 Int J Epidemiol 44:512 (MR-Egger; InSIDE)
- Bowden J et al 2016 Int J Epidemiol 45:1961 (NOME, I^2_GX, SIMEX)
- Cook JR & Stefanski LA 1994 JASA 89:1314 (SIMEX framework)
- Burgess S 2020 Nat Commun 11:376 (contamination mixture)
- Burgess S & Thompson SG 2021 (Mendelian Randomization 2nd ed)
- Darrous L et al 2021 Nat Commun 12:7274 (LHC-MR)
- Foley CN et al 2020 Bioinformatics 37:531 (MR-Clust)
- Hartwig FP et al 2017 Int J Epidemiol 46:1985 (weighted mode)
- Hemani G et al 2017 PLoS Genet 13:e1007081 (Steiger orientation)
- Hemani G et al 2018 eLife 7:e34408 (TwoSampleMR framework)
- Hemani G & Tilling K 2022 Int J Epidemiol (Steiger limitations under measurement error)
- Morrison J et al 2020 Nat Genet 52:740 (CAUSE; CHP)
- O'Connor LJ & Price AL 2018 Nat Genet 50:1728 (LCV)
- Sanderson E et al 2022 Nat Rev Methods Primers 2:6 (MR Primer)
- Skrivankova VW et al 2021 JAMA 326:1614 (STROBE-MR statement); BMJ 375:n2233 (explanation)
- Verbanck M et al 2018 Nat Genet 50:693 (MR-PRESSO)
- Zhao Q et al 2020 Ann Stat 48:1742 (MR-RAPS)

## Related Skills

- causal-genomics/mendelian-randomization - Primary causal estimation that this sensitivity battery validates
- causal-genomics/genetic-correlation - LDSC rg required for Step 1 of the decision flow; CHP escalation trigger
- causal-genomics/colocalization-analysis - Required for cis-MR drug-target signals where instruments are too few for Egger
- causal-genomics/fine-mapping - Identify causal variants underlying instrument loci
- causal-genomics/mediation-analysis - Multivariable MR for mediator-adjusted causal estimates
- population-genetics/association-testing - GWAS summary statistics underlying MR instruments
- clinical-biostatistics/effect-measures - Translate MR estimates to clinical effect measures
