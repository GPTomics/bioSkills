---
name: bio-clinical-biostatistics-trial-reporting
description: Prepares statistical reports for clinical trials following CONSORT 2025, SPIRIT 2025, ICH E9(R1) estimands, and FDA 2023 covariate adjustment guidance. Covers Table 1 generation, analysis populations (ITT/FAS/PP/Safety), the 5 ICH E9(R1) intercurrent-event strategies, MMRM under MAR (mmrm), reference-based MI (rbmi J2R/CR/CIR), Permutt tipping-point sensitivity, and Rubin's-rules vs frequentist variance debate. Use when preparing regulatory submissions, defining estimands, or implementing missing-data sensitivity analyses.
tool_type: python
primary_tool: tableone
---

## Version Compatibility

Reference examples tested with: tableone 0.9+, statsmodels 0.14+, scikit-learn 1.4+, pandas 2.1+, numpy 1.26+. R packages cited (essential for current regulatory work): mmrm 0.3+ (Roche/openpharma), rbmi 1.5+ (Roche/Bayer via insightsengineering), gMCP, RBesT.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Trial Reporting Under CONSORT 2025 + ICH E9(R1)

**"Prepare a clinical trial statistical report"** -> Define the estimand explicitly per ICH E9(R1); execute a covariate-adjusted primary analysis targeting the right summary measure; pre-specify the missing-data strategy and run regulatory-grade sensitivity analyses; structure the output per CONSORT 2025 and the new SPIRIT 2025 alignment.

## The Single Most Important Methodological Shift -- The Estimand Comes First

**Kahan, Cro, Li, Harhay 2023 *Am J Epidemiol* 192:987 ("Eliminating Ambiguous Treatment Effects Using Estimands"):** 98% of published trial reports do not describe what the reported treatment effect represents. 54% of trials: impossible to deduce the estimand from reported methods. In 74% of trials submitted for regulatory approval 1996-2017, "what-if" hypothetical effects were used but only 2 trials explained this.

**The framework:** ICH E9(R1) Addendum (November 2019, EMA effective Feb 2020, FDA May 2021) defines an estimand as the precise specification of what is being estimated, via five attributes:

1. **Treatment condition** -- what is being compared
2. **Population** -- which patients
3. **Endpoint** -- which variable
4. **Population-level summary measure** -- mean diff, OR, HR, RD
5. **Intercurrent-event (ICE) handling strategy** -- one of five

**The order is non-negotiable:** specify the estimand BEFORE choosing the statistical method. Choosing MMRM and retrofitting the estimand to match is the canonical error.

## The Five Intercurrent-Event Strategies

| Strategy | What it does | Typical implementation | Identification cost | Regulatory pattern |
|----------|--------------|------------------------|---------------------|---------------------|
| **Treatment policy** | Include all data regardless of ICE | ANCOVA on observed value (retrieved-dropout data); ITT-like | Trivially identified; needs full follow-up regardless of ICE | FDA-preferred default for cardio/HF/CV-safety; Fleming 2025 endorsement |
| **Hypothetical** | What would have been observed had ICE not occurred | MMRM under MAR; g-computation; IPCW; reference-based MI under MAR | Sequential ignorability (causal); MAR (missing-data shorthand) | Heavy use in CNS, diabetes, respiratory; EMA more accepting than FDA |
| **Composite** | Incorporate ICE into endpoint | Death = non-responder; PFS = composite of progression OR death; MACE | Identified from observed data; embeds a ranking choice | Standard in oncology PFS; acceptable when ICE has clinical signal |
| **While-on-treatment** | Use only pre-ICE values | Censor at discontinuation (for TTE); analyse last pre-ICE value (for repeated measures) | Estimates conditional quantity | Safety endpoints (AE rate per time on drug); FDA cautious for efficacy |
| **Principal stratum** | Confine to latent stratum (e.g. tolerators) | Bayesian estimation under monotonicity/principal ignorability | Latent membership; unverifiable assumptions | Rare as primary; some oncology/vaccine acceptance |

**Postdoc reading list:**

- Bornkamp, Rufibach, Lin, Liu, Mehrotra, Roychoudhury, Schmidli, Shentu, Wolbers 2021 *Pharm Stat* 20:737 — principal stratum critique
- Olarte Parra, Bartlett, Daniel 2022 *Stat Biopharm Res* — proves MMRM-MAR IS a hypothetical estimator under specific causal assumptions
- Permutt 2016 *Stat Med* 35:2865; 2020 *Stat Biopharm Res* 12:45 — FDA Missing Data Working Group taxonomy; "do covariates change the estimand?"
- Fleming, Carroll, Wittes et al 2025 *Stat Med* — argue treatment policy is the only strategy preserving randomisation; critique hypothetical
- Morris 2026 *Stat Med* — causal-perspective comment on Fleming
- Lipkovich, Ratitch, Mallinckrodt 2020 *Stat Biopharm Res* — Rubin causal model connection

## Decision Tree for Estimand Selection

| Trial scenario | Recommended estimand strategy | Why |
|----------------|-------------------------------|-----|
| Continuous endpoint, monotone missingness, MAR plausible | Hypothetical via MMRM (mmrm + KR) | Standard FDA-favoured MAR analysis; cite Mallinckrodt 2008/2014 |
| Continuous endpoint, ICE = treatment discontinuation, sponsor wants effectiveness | Treatment policy via retrieved-dropout MI | If post-ICE data available; ITT-respecting |
| Continuous endpoint, treatment-policy primary with ICE-related missingness | Hybrid: J2R for discontinuation ICEs, MMRM-MAR for other missingness | Aprocitentan precedent; de facto FDA standard 2024-2025 |
| Binary endpoint, RCT, FDA 2023-compliant | Marginal RD via g-computation; conditional OR supportive | See clinical-biostatistics/logistic-regression for g-computation |
| Oncology OS with crossover | Treatment policy as primary; hypothetical (RPSFT/IPCW) as sensitivity | Sotorasib CodeBreaK 200 precedent |
| Oncology PFS | Treatment policy with composite for subsequent therapy ICE | Lewis 2023 framework; Fleming 2025 |
| Weight management / chronic disease | Retrieved-dropout MI (Wegovy STEP precedent); J2R supportive | FDA 2025 obesity draft guidance explicitly endorses |
| Long-term safety endpoint | While-on-treatment for rate; treatment policy for cumulative incidence | Standard ICH E2A practice |
| AlloSCT in hematologic oncology | Composite "event-free survival" treating alloSCT as event | Rufibach 2020; "no alloSCT" hypothetical is clinically meaningless |
| Symptomatic palliative endpoint with high dropout | Composite with worst-rank for dropouts (Permutt trimmed means) | Permutt 2017 *Pharm Stat* 16:20 |

## MMRM -- The FDA-Favoured MAR Analysis

**Mallinckrodt 2008/2014, codified in DIA Scientific Working Group "three pillars" doctrine:** for continuous longitudinal endpoints under monotone (or near-monotone) MAR, an MMRM with treatment + visit + treatment-by-visit + baseline + baseline-by-visit, unstructured (UN) within-subject covariance, REML, contrast at the primary timepoint -- is the consistent and FDA-preferred analysis. LOCF is biased even under MCAR because it discards imputation uncertainty and assumes a flat post-withdrawal trajectory.

### The mmrm R package (Roche / openpharma)

```r
library(mmrm)
fit <- mmrm(
    formula = change_from_baseline ~ baseline + arm * visit + us(visit | subject),
    data = trial_data,
    method = "Kenward-Roger",  # or "Satterthwaite", "Kenward-Roger-Linear"
    reml = TRUE
)
summary(fit)  # treatment-by-visit contrast at primary timepoint
```

**The Kenward-Roger flavour question:** `method = "Kenward-Roger"` uses full second-order Kenward-Roger (Kenward-Roger 1997 *Biometrics* 53:983), which inflates SE for fixed-effect contrasts using an adjusted covariance estimator with second-order Taylor terms. **`method = "Kenward-Roger-Linear"` drops the second-order Cholesky-derivative term to match SAS PROC MIXED bit-for-bit.** Most submissions use Kenward-Roger-Linear to maintain SAS-R reproducibility.

### Convergence-vs-correctness trade-off

Unstructured (UN) covariance has p(p+1)/2 parameters for p visits. With ~30-50 patients per arm by week 12 and 6+ visits, UN can fail to converge. The industry-standard fallback hierarchy per pre-specified SAP:

1. UN with KR (preferred)
2. UN with Satterthwaite (if KR fails)
3. Heterogeneous Toeplitz (k+1 parameters)
4. AR(1) with heterogeneous variances
5. CS with heterogeneous variances (last resort)

Each step down imposes more structure and the structure can be wrong — biasing both SEs and point estimates. CS imposes equal correlation across time which is rarely true for treatment-ramp-up endpoints (HbA1c, BP). Pre-specify the fallback in the SAP, not at analysis time.

### MMRM = hypothetical estimator (Olarte Parra unification)

**Olarte Parra, Bartlett, Daniel 2022 *Stat Biopharm Res*:** under specific identifying assumptions, MMRM under MAR IS a causal hypothetical estimand via g-formula equivalence. The "issue" is articulation, not statistical machinery — MMRM-MAR implicitly answers a hypothetical estimand whose hypothetical scenario must be made explicit in the SAP (e.g., "what would the mean response at week 24 be had all patients continued randomised treatment and remained observable?").

## Reference-Based Multiple Imputation -- The rbmi Framework

**Carpenter, Roger, Kenward 2013 *J Biopharm Stat* 23:1352** — the canonical paper. Reference-based MI operationalises MNAR sensitivity not as a numeric delta but as a clinical narrative:

- **Jump-to-reference (J2R):** "after withdrawal the patient instantly resembles the placebo arm"
- **Copy-reference (CR):** "the entire post-baseline trajectory copied from placebo, with subject's baseline deviation preserved"
- **Copy-increments-in-reference (CIR):** "the patient retains the on-treatment increment but trends with the placebo arm thereafter, anchored at last observed value"
- **Last-mean-carried-forward (LMCF):** "patient stays at last on-treatment mean"

**rbmi R package (Wolbers et al 2022 *Pharm Stat* 21(6):1246-1257; CRAN; insightsengineering):**

```r
library(rbmi)
# Draws -> Impute -> Analyse -> Pool pipeline
draws <- draws(data = trial_data, vars = vars,
               method = method_bayes(n_samples = 100))
imputed <- impute(draws, references = c('Active' = 'Placebo', 'Placebo' = 'Placebo'))
analyses <- analyse(imputed, fun = ancova,
                    vars = list(outcome = 'change', visit = 'avisit',
                                group = 'arm', covariates = c('baseline')))
result <- pool(analyses)  # Rubin's rules pooling
```

Four inference engines:

1. **Bayesian MI + Rubin's rules** (historical default) — information-anchored variance
2. **Approximate Bayesian via REML + bootstrap** — frequentist variance from bootstrap
3. **Conditional mean imputation + jackknife** (Wolbers 2022 contribution) — single deterministic imputation per ANCOVA-linearity theorem, jackknife for SE; FDA-friendly because deterministic + frequentist
4. **BMLMI (bootstrapped MI of Lipkovich/Ratitch)** — within/between variance decomposition

## The Variance Debate -- Cro/Carpenter vs Bartlett/Wolbers

**The single most active methodological argument in current biostatistics.**

**Cro/Carpenter/Kenward 2019 *JRSS-A* 182:623 ("Information-Anchored Sensitivity Analysis"):** proved that Rubin's-rules variance applied to J2R/CR/CIR is *approximately information-anchored* — the relative loss of information from missingness in the sensitivity analysis matches the relative loss in the MAR primary analysis. True repeated-sampling variance is "information positive" because reference-based imputation borrows from the reference arm and reduces the marginal variance of the active arm BELOW what an MAR analysis with the same missingness would give.

**Philosophical position:** a sensitivity analysis should not import information the primary analysis did not have; if borrowing from placebo makes the active-arm CI narrower, the analysis is no longer "anchored" to the same information state.

**Bartlett 2021 *Stat Biopharm Res* 15(1):178 + Wolbers 2022 *Pharm Stat* counter:** if J2R is the actual sampling model under which inference is made, then the *correct* frequentist variance is the one that delivers nominal Type-I error and CI coverage under that model -- the jackknife/bootstrap variance, NOT Rubin's. Simulations in `rbmi` vignettes: Bayesian MI with Rubin's gives Type-I error 0.9-2.5% (over-conservative); CMI+jackknife gives 4.84-4.96% (nominal) under J2R; Bayesian MI loses real power.

**Regulatory practice 2024-2025 is bifurcating:** EMA tolerates either; FDA reviewers increasingly flag Rubin's-rules variance under reference-based MI as needing a frequentist sensitivity analysis in addition. **What postdocs argue about:** whether Type-I inflation under bootstrap is the price of correct inference, or evidence J2R was never coherent as a true sampling model.

## Permutt Tipping-Point Analysis -- The Analyst as Adversary

**Permutt 2016 *Stat Med* 35:2876** (Permutt was head of FDA Division of Biometrics IV): the regulator's question is not "what is a reasonable MNAR adjustment?" but "how bad would the missing data have to be in the active arm to overturn the significant primary result?"

**Delta-adjustment patterns:**

- One-arm shift (FDA preferred): add delta to imputed values in active arm only; vary delta from 0 to the value that nullifies the effect
- Symmetric shift: both arms worsened by delta (probes systematic optimism)
- Reverse shift: placebo improved by delta (more aggressive, rarely needed)

```r
# rbmi with delta adjustment
delta <- delta_template(imputed, delta = c(0, 5, 10, 15, 20), dlag = c(1, 1, 1, 1))
adjusted <- analyse(imputed, delta = delta, ...)
# Report: minimum delta that flips p-value below 0.05
```

The regulator then judges whether the tipping delta is clinically plausible — larger than the active-arm treatment effect itself? Larger than the MCID? FDA-preferred report: tipping delta in units of residual SD (for cross-trial comparison), not raw outcome units.

## Decisive Regulatory Cases -- The 2020-2025 Casebook

**Aducanumab (Biogen BLA 761178, 2021):** EMERGE and ENGAGE studies both stopped early for futility; EMERGE high-dose positive, ENGAGE negative. MMRM-MAR primary. FDA Office of Biostatistics (Tristan Massie review) argued futility-stop-induced missingness was not MAR (differential ARIA-driven unblinding); 6-1 AdCom against approval. Textbook case showing MAR-based primary in trial with high differential missingness is regulator-divisive.

**Aprocitentan (Idorsia PRECISION trial, FDA approval 2024):** documented in Mathur 2025 *Pharm Stat* (PMC12753554). FDA pushed back on sponsor's MMRM-MAR primary; MAR was not credible for treatment-discontinuers. Accepted compromise: stratified imputation — J2R for treatment-discontinuation ICEs, MAR-MMRM for other missingness. This hybrid is now de facto FDA standard for treatment-policy estimand.

**Wegovy/Ozempic STEP trials (Wilding 2021 *NEJM*; NDA 215256):** retrieved-dropout MI as primary for treatment-policy. Missing body weight at week 68 imputed by sampling from observed week-68 measurements among "retrieved dropouts" (patients who discontinued semaglutide but remained in follow-up). J2R-MI as supportive. RD-MI now standard for chronic weight management. FDA 2025 obesity guidance explicitly endorses MI as primary.

## Table 1 -- Baseline Characteristics

**CONSORT 2010 discouraged baseline significance tests** because randomisation is a known mechanism, not a hypothesis. Many journals still require them.

```python
from tableone import TableOne

columns = ['age', 'sex', 'race', 'bmi', 'baseline_score', 'disease_stage']
categorical = ['sex', 'race', 'disease_stage']

table1 = TableOne(df, columns=columns, categorical=categorical,
                  groupby='ARM', pval=True, smd=True,
                  missing=True, overall=True)
print(table1.tabulate(tablefmt='github'))
table1.to_excel('table1.xlsx')
```

**Use standardised mean differences (SMD) rather than p-values:** SMD > 0.1 suggests meaningful imbalance regardless of statistical significance.

**Senn's "balance testing is incoherent" (1994 *Stat Med* 13:1715; Altman 1985):** balancing via randomisation, testing balance, then adjusting only when the test fails is a selection rule that destroys nominal Type-I error. Pre-specify covariates in the SAP; do not condition adjustment on observed imbalance.

## Analysis Populations -- ITT vs FAS vs PP vs Safety

| Population | Definition | Bias direction | Primary use |
|-----------|------------|----------------|-------------|
| ITT | All randomised, as randomised | Conservative (toward null) | Primary efficacy per ICH E9 |
| FAS (Full Analysis Set) | ITT excluding eligibility failures + subjects with no post-baseline data | Middle ground; close to ITT | Common practical primary; ICH E9 "as complete as possible while remaining unbiased" |
| Per-Protocol | Completed treatment per protocol without major violations | Anti-conservative (inflates effect) | Sensitivity analysis only |
| Safety | All received at least one dose | n/a | AE analysis |
| mITT | Sponsor-defined modified ITT | Variable | Pre-specify and justify |

**FAS vs ITT distinction is critical for regulatory submissions** — FAS may exclude post-randomisation subjects (ineligibility, no post-baseline efficacy); ITT cannot. Sponsors often equate them on the SAP only to discover at submission that FDA expected stricter ITT. Pre-specification in protocol is essential.

```python
itt = dm.copy()
pp = dm[dm['USUBJID'].isin(completers) & ~dm['USUBJID'].isin(protocol_violators)]
dosed = ex[ex['EXDOSE'] > 0]['USUBJID'].unique()
safety = dm[dm['USUBJID'].isin(dosed)]
```

## Missing Data Mechanisms -- The Practical Framework

| Mechanism | Definition | Testable? | Valid method |
|-----------|------------|-----------|--------------|
| MCAR | Independent of all data | Partially (Little's test) | Complete-case unbiased but loses power |
| MAR | Depends on observed data only | NO (assumption) | MMRM under MAR; MI under MAR |
| MNAR | Depends on unobserved values | NO | Requires sensitivity analysis (J2R, CR, CIR, tipping point) |

**MAR vs MNAR cannot be distinguished from observed data alone — this is a fundamental limitation.** Pre-specify the assumed mechanism in the SAP; pre-specify the sensitivity analysis under MNAR (NRC 2010 Recommendation 15: "examining sensitivity to assumptions about the missing-data mechanism should be a mandatory component of reporting").

**Clinical reasoning beyond the abstraction:** examine the DS (Disposition) domain to tabulate reasons for discontinuation by treatment arm. If discontinuation rates or reasons differ between arms, missing data is likely informative and MNAR sensitivity analyses are mandatory.

## Standard MMRM-MAR Primary Analysis Code

**Goal:** Execute the FDA-preferred primary analysis for a continuous longitudinal endpoint under MAR with valid Type-I in small/moderate trials.

**Approach:** Fit MMRM with unstructured covariance + Kenward-Roger via the Roche/openpharma `mmrm` R package; in Python, use `statsmodels.mixedlm` as exploratory-only (lacks KR).

```python
# Python is weak for MMRM — current state of the art is R `mmrm`
# For Python users, statsmodels.mixedlm is the closest alternative but lacks
# Kenward-Roger; consider rpy2 to call R from Python for confirmatory work.

import statsmodels.formula.api as smf
import pandas as pd

# Random intercept LMM (suboptimal vs MMRM but Python-native)
model = smf.mixedlm(
    'change ~ baseline + C(ARM) * C(VISIT)',
    data=df_long,
    groups=df_long['USUBJID']
).fit(reml=True)
# WARNING: this is NOT FDA-equivalent to MMRM with UN+KR
# For confirmatory work, use R `mmrm` package via rpy2 or fit in R directly
```

## Multiple Imputation in Python (sklearn)

```python
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
import statsmodels.formula.api as smf
import numpy as np

n_imputations = 20  # rule: m >= 100 * FMI (fraction of missing info)
imputer = IterativeImputer(max_iter=10, random_state=0, sample_posterior=True)

results = []
for i in range(n_imputations):
    imputer.set_params(random_state=i)
    imputed = pd.DataFrame(imputer.fit_transform(df[numeric_cols]), columns=numeric_cols)
    for col in ['ARM', 'sex']:
        imputed[col] = df[col].values
    model = smf.logit(
        'outcome ~ C(ARM, Treatment(reference="Placebo")) + age', data=imputed
    ).fit(disp=0)
    results.append({'coef': model.params.iloc[1], 'se': model.bse.iloc[1]})

# Rubin's rules
pooled_coef = np.mean([r['coef'] for r in results])
within_var = np.mean([r['se']**2 for r in results])
between_var = np.var([r['coef'] for r in results], ddof=1)
total_var = within_var + (1 + 1/n_imputations) * between_var
pooled_se = np.sqrt(total_var)
pooled_or = np.exp(pooled_coef)
```

**Critical sklearn caveats:**

- `sample_posterior=True` is essential — without it all m imputations are nearly identical
- **`sample_posterior=True` only works with `BayesianRidge` (the default estimator)**. If estimator is changed (e.g., `RandomForestRegressor`), parameter is silently ignored and MI degenerates to single imputation
- IterativeImputer treats binary covariates as continuous; consider `miceforest` for mixed types or move to R `mice`/`rbmi`
- Only impute covariates and post-baseline outcomes, NOT treatment assignment (fully determined by randomisation)
- Include outcome in imputation model as predictor but exclude from imputed variables to avoid circular dependency
- IterativeImputer is **experimental** in sklearn — API may change without standard deprecation

**Congeniality (Meng 1994):** imputation model must be at least as flexible as analysis model. If analysis includes treatment-by-covariate interactions, imputation model should include them. Uncongenial imputation biases estimates and invalidates variance pooling.

## Co-Primary Endpoints and Multiplicity

| Method | Approach | Conservatism |
|--------|----------|--------------|
| Bonferroni | alpha / m | Most conservative |
| Hierarchical (gatekeeping) | Pre-specified order; proceed only if previous rejects | Moderate; full alpha for first |
| Graphical procedure (Bretz-Maurer) | Directed graph; alpha propagates on rejection | Flexible; standard in modern SAPs |
| Hochberg / Hommel step-up | Ordered p-values vs alpha/(m-k+1) | Less conservative than Bonferroni; requires PRDS |

See clinical-biostatistics/multiplicity-graphical for the full Bretz-Maurer-Hommel treatment with `gMCP`.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MMRM-MAR primary p < 0.05; reference-based MI sensitivity p > 0.05 | MNAR mechanism after treatment discontinuation; J2R imputes active arm toward placebo | Decide which estimand is regulatory primary; if treatment policy, switch to reference-based MI as primary (Aprocitentan precedent) |
| Retrieved-dropout MI vs J2R-MI differ | RD-MI uses observed post-ICE data; J2R uses reference-arm-based assumption | If post-ICE data available, RD-MI is empirically grounded (preferred); J2R as supportive |
| Cro information-anchored variance vs Wolbers frequentist variance differ | Information-anchored Rubin's pools within+between; frequentist jackknife/bootstrap captures borrowing-induced variance reduction | Report both; cite Cro 2019 and Wolbers 2022 (Bartlett 2021 frequentist critique); 2024-2025 FDA practice accepts both with one as supportive |
| ITT primary p < 0.05; PP secondary p > 0.05 | Per-protocol excludes non-completers (often differentially); PP-only effect inflated when significant; non-significance in PP is signal of fragility | ITT remains primary (ICH E9); PP as sensitivity flag; investigate completer pattern |
| Statistical significance achieved but effect estimate below MCID | Powered for δ << MCID; or δ = MCID with no precision margin | Pre-specify δ >= 1.5 × MCID in SAP (postdoc rule); report against pre-specified MCID, not 0 |
| Estimand strategies (hypothetical vs treatment-policy) give different effect sizes | Different ICE-handling strategies target different parameters | Report all pre-specified estimands; pick PRIMARY for the regulatory question (not the "winner"); cite Kahan 2023 |
| Subgroup effect estimate >2x main effect | Winner's curse (Sun 2010 documented median 2.4x inflation) | Apply Bayesian shrinkage (Dixon-Simon, RBesT) for corrected estimate; cite as discovery, not confirmatory |

## CONSORT 2025 (April 2025) -- What Changed

**Citation:** Hopewell, Chan, Collins et al 2025. CONSORT 2025 statement. *Lancet* 405:1633-1640. E&E in BMJ 2025;389:e081124.

**30-item checklist (was 25)**, 7 new items, 3 substantially revised, 1 deleted. Key new items relevant to statistical reporting:

- **Item 4** (data and code sharing — where/how de-identified IPD and analysis code can be accessed) — NEW
- **Item 15** (how harms were actually assessed — methods, attribution, intensity, monitoring rules) — NEW; absorbs CONSORT-Harms 2022
- **Item 21c** (missing-data handling methods) — NEW; reflects NRC 2010 / Carpenter-Kenward consensus
- **Item 24a/24b** (intervention/comparator actual delivery details) — absorbs TIDieR

**Estimands did NOT make consensus for mandatory inclusion** — they appear in Box 1 (terminology only). **For regulatory submissions, ICH E9(R1) is the operative estimand standard, NOT CONSORT 2025.** Sponsors should follow ICH E9(R1) directly for the 5-attribute estimand statement.

**No DOORS framework in CONSORT 2025** — diversity/equity additions are happening via a separate SAGER-SPIRIT-CONSORT alignment workstream (GENDRO/EASE).

## CONSORT Flow Diagram

```python
flow = {
    'screened': len(screening_log),
    'eligible': len(screening_log[screening_log['eligible']]),
    'randomized': len(dm),
    'allocated_drug': len(dm[dm['ARM'] == 'Drug']),
    'allocated_placebo': len(dm[dm['ARM'] == 'Placebo']),
    'completed_drug': len(dm[(dm['ARM'] == 'Drug') & dm['USUBJID'].isin(completers)]),
    'completed_placebo': len(dm[(dm['ARM'] == 'Placebo') & dm['USUBJID'].isin(completers)]),
    'analyzed_itt': len(itt),
    'analyzed_fas': len(fas),
    'analyzed_pp': len(pp),
    'analyzed_safety': len(safety),
}
```

CONSORT 2025 templates on consort-spirit.org now explicitly accommodate non-1:1 allocation, cluster, multi-arm, and crossover variants.

## Per-Method Failure Modes

### MMRM-MAR with high differential missingness

- **Trigger:** Differential dropout between arms (e.g., toxic active vs tolerated placebo).
- **Mechanism:** MAR assumption requires that conditional on observed covariates, missingness is unrelated to outcome — implausible when patients drop out *because* the treatment isn't working.
- **Symptom:** Discontinuation reasons differ qualitatively between arms; primary p-value sensitive to model specification.
- **Fix:** Treatment-policy estimand with retrieved-dropout MI as primary; J2R as sensitivity. Tipping-point delta in active arm (Permutt 2016).

### Reference-based MI with Rubin's variance

- **Trigger:** J2R/CR/CIR with Rubin's-rules variance reported as primary.
- **Mechanism:** Rubin's variance is information-anchored but over-conservative (Cro 2019); under-rejects relative to nominal.
- **Symptom:** Power loss vs frequentist variance; Type-I 0.9-2.5% vs nominal 5% in J2R simulations.
- **Fix:** Report frequentist variance via CMI+jackknife as supportive (Wolbers 2022); cite both Cro and Bartlett.

### sample_posterior=False in IterativeImputer

- **Trigger:** Default behaviour or estimator changed away from BayesianRidge.
- **Mechanism:** Imputations are point predictions, near-identical across draws; between-imputation variance approximately zero.
- **Symptom:** Artificially narrow Rubin's-pooled CIs.
- **Fix:** Always set `sample_posterior=True`; verify estimator is BayesianRidge; consider `miceforest` or R `mice` for mixed types.

### Per-protocol as primary

- **Trigger:** SAP specifies PP analysis as primary.
- **Mechanism:** Excludes non-completers who may have dropped due to treatment failure; inflates effect.
- **Symptom:** PP effect much larger than ITT; reviewers raise post-randomisation bias concern.
- **Fix:** ITT or FAS as primary per ICH E9; PP as sensitivity only.

### Method-before-estimand

- **Trigger:** SAP describes MMRM, MI, or LOCF without specifying the estimand the method targets.
- **Mechanism:** Retrofits the estimand to match the chosen method.
- **Symptom:** Reviewer asks "what is the estimand?" and sponsor cannot articulate.
- **Fix:** Per ICH E9(R1), define the 5 attributes BEFORE choosing method; cite Kahan 2023.

### Missing the FAS vs ITT distinction

- **Trigger:** Treating FAS and ITT as synonymous.
- **Mechanism:** FAS may exclude eligibility failures or no-post-baseline; ITT cannot.
- **Symptom:** Submission counts differ between FAS and ITT; reviewer asks for ITT.
- **Fix:** Pre-specify both populations explicitly in protocol; report both with reconciliation.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| m >= 100 * FMI imputations | von Hippel 2020 *J Roy Stat Soc A* | Adequate for stable pooled SE; with 40% missingness and FMI ~0.3, m=30 needed |
| SMD > 0.1 = meaningful imbalance | Austin 2009 *Stat Med* 28:3083 | Beyond what randomisation would normally produce |
| Missing data > 40% on key variable | NRC 2010 | Above this, MI under MAR is unreliable; treat as hypothesis-generating |
| Kenward-Roger DF correction for MMRM with UN | Kenward-Roger 1997 | Without it, MMRM-REML under-covers in small/moderate trials; Type-I inflates 1-2 pp |
| Information-anchored vs frequentist variance for reference-based MI | Cro 2019 vs Wolbers 2022 | Active regulatory debate; report both for safety |
| Tipping delta in residual SD units, not raw | FDA Division of Biometrics preference | Cross-trial comparison; report adjacent to raw |
| Treatment policy default for cardio/HF/CV | Fleming 2025 *Stat Med* | Only strategy preserving randomisation; FDA-favoured |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Reviewer: "what is the estimand?" with no answer | Method-before-estimand | Pre-specify 5 attributes in protocol; cite Kahan 2023 finding 98% don't articulate |
| LOCF used as primary | Inertia or "conservative" misconception | LOCF is biased even under MCAR (Mallinckrodt 2008); switch to MMRM or MI |
| MMRM with CS covariance treated as equivalent to UN | Convergence forced fallback without pre-specification | Pre-specify fallback hierarchy; document deviation if invoked |
| Rubin's variance for J2R with no sensitivity | Cro 2019 information-anchored argument applied without acknowledgement | Cite Bartlett 2021 + Wolbers 2022; report frequentist variance via CMI+jackknife |
| Tipping delta in raw outcome units only | Hard to compare cross-trial | Also report in residual SD units (FDA preference) |
| ITT and FAS conflated | SAP vague on distinction | Pre-specify both with explicit criteria for FAS exclusion |
| Per-protocol significant, ITT not — sponsor highlights PP | Post-randomisation bias inflation | ITT as primary; PP as sensitivity with explicit caveat |
| Imputation only of outcome | Throws away covariates | Joint imputation of outcome + covariates; cite Carpenter-Roger 2013 |
| `sample_posterior` silently ignored | Non-default estimator | Verify with `imputer.estimator` is BayesianRidge or use rbmi/mice |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "What is the estimand?" | Articulate 5 ICH E9(R1) attributes; cite Kahan 2023; state ICE strategy with mechanism (not just label). |
| "Is MAR plausible?" | Examine DS for differential discontinuation patterns; if differential, switch to treatment-policy with retrieved-dropout MI or J2R sensitivity. |
| "Why MMRM not LOCF?" | Mallinckrodt 2008/2014 case for MMRM; LOCF is biased even under MCAR. Cite DIASWG three-pillars doctrine. |
| "Why Rubin's not frequentist variance for J2R?" | Cite Cro 2019 information-anchored argument as rationale; report frequentist (CMI+jackknife) as supportive per Wolbers 2022. |
| "Tipping point analysis result?" | Minimum delta in active arm that flips p > 0.05; report in residual SD units; judge plausibility against MCID and treatment effect. |
| "FAS vs ITT?" | Pre-specified in protocol; FAS excludes [explicit criteria]; both populations analysed; reconciliation table provided. |
| "Has the SAP been registered?" | Yes — clinicaltrials.gov NCTxxxxx with full SAP appended; EU CTR EUCTxxxx; SPIRIT 2025 compliant. |
| "How is subsequent therapy handled?" | Pre-specified per ICH E9(R1); composite strategy (e.g., subsequent therapy = treatment failure) or treatment-policy (include all data). |
| "Where is the multiplicity adjustment for co-primary / key secondary?" | Bretz-Maurer graphical procedure via gMCP pre-specified in SAP; alpha allocation diagram in CSR appendix; cite CONSORT 2025 item 20 + FDA Multiple Endpoints Final Oct 2022; see clinical-biostatistics/multiplicity-graphical. |
| "How are missing baseline covariates handled in ANCOVA?" | Complete-case ANCOVA is unbiased under MCAR baseline missingness (Kahan-Morris 2012); proportion with complete baseline reported in flow diagram; if >5% missing, multiple imputation of baseline as sensitivity analysis. |

## References

- Bornkamp B, Rufibach K, Lin J, Liu Y, Mehrotra DV, Roychoudhury S, Schmidli H, Shentu Y, Wolbers M. 2021. Principal stratum strategy: potential role in drug development. *Pharm Stat* 20:737-751.
- Carpenter JR, Roger JH, Kenward MG. 2013. Analysis of longitudinal trials with protocol deviation: a framework for relevant, accessible assumptions, and inference via multiple imputation. *J Biopharm Stat* 23:1352-1371.
- Cro S, Carpenter JR, Kenward MG. 2019. Information-anchored sensitivity analysis: theory and application. *JRSS-A* 182:623-645.
- EMA. 2010. Guideline on Missing Data in Confirmatory Clinical Trials. EMA/CPMP/EWP/1776/99 Rev.1.
- Fleming TR, Carroll KJ, Wittes JT, Emerson SS, Rothmann M, Collins R, Levin G. 2025. Perspectives on appropriate implementation of the ICH E9(R1) estimand framework. *Stat Med* e70104.
- Hopewell S, Chan AW, Collins GS, et al. 2025. CONSORT 2025 statement. *Lancet* 405:1633-1640.
- ICH. 2019. E9(R1) Addendum on Estimands and Sensitivity Analysis in Clinical Trials. Step 4.
- Kahan BC, Cro S, Li F, Harhay MO. 2023. Eliminating ambiguous treatment effects using estimands. *Am J Epidemiol* 192:987-994.
- Kenward MG, Roger JH. 1997. Small sample inference for fixed effects from restricted maximum likelihood. *Biometrics* 53:983-997.
- Lipkovich I, Ratitch B, Mallinckrodt CH. 2020. Causal inference and estimands in clinical trials. *Stat Biopharm Res* 12:54-67.
- Mallinckrodt CH, Lane PW, Schnell D, Peng Y, Mancuso JP. 2008. Recommendations for the primary analysis of continuous endpoints in longitudinal clinical trials. *Drug Information Journal* 42:303-319.
- Mathur AR et al. 2025. Regulatory experiences with the use of multiple imputation for missing data in a phase 3 confirmatory trial. *Pharm Stat*.
- Morris TP. 2026. Comment on Fleming et al. *Stat Med* e70455.
- NRC. 2010. The Prevention and Treatment of Missing Data in Clinical Trials. National Academies Press.
- Olarte Parra C, Bartlett JW, Daniel RM. 2022. Hypothetical estimands: a unification of causal inference and missing data methods. *Stat Biopharm Res*.
- Permutt T. 2016. Sensitivity analysis for missing data in regulatory submissions. *Stat Med* 35:2876-2879.
- Permutt T. 2020. Do covariates change the estimand? *Stat Biopharm Res* 12:45-53.
- von Hippel PT. 2020. How many imputations are needed (a two-stage calculation using a quadratic rule). *Sociol Methods Res* 49:699-718.
- Wolbers M, Noci A, Delmar P, Gower-Page C, Yiu S, Bartlett JW. 2022. Standard and reference-based conditional mean imputation. *Pharm Stat* 21(6):1246-1257.

## Related Skills

- clinical-biostatistics/cdisc-data-handling - SDTM/ADaM data preparation feeding the analysis dataset
- clinical-biostatistics/logistic-regression - Primary analysis models (binary) with marginal-vs-conditional
- clinical-biostatistics/effect-measures - Effect-measure CIs under MI pooling
- clinical-biostatistics/subgroup-analysis - Pre-specified subgroup analyses for the CSR
- clinical-biostatistics/missing-data-sensitivity - MMRM, reference-based MI, tipping point (in depth)
- clinical-biostatistics/multiplicity-graphical - Bretz-Maurer graphs for co-primary and key secondary
- clinical-biostatistics/survival-analysis - Time-to-event estimand framing per Lewis 2023
- clinical-biostatistics/power-and-sample-size - Sample size justification per CONSORT 2025 item 7
- reporting/rmarkdown-reports - Formatted statistical report generation
- experimental-design/multiple-testing - General multiplicity correction
