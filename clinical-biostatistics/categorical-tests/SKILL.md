---
name: bio-clinical-biostatistics-categorical-tests
description: Tests associations between categorical variables in clinical data using chi-square, Fisher's exact, Boschloo, Cochran-Mantel-Haenszel, and modern McNemar variants with calibrated confidence intervals (Wilson, Newcombe, Miettinen-Nurminen). Use when analyzing categorical outcomes, paired binary endpoints, or testing treatment-outcome independence in confirmatory or exploratory clinical trials.
tool_type: python
primary_tool: scipy
---

## Version Compatibility

Reference examples tested with: scipy 1.12+ (Boschloo and Barnard added in 1.7), statsmodels 0.14+, pingouin 0.5+, exact2x2 (R) 1.6+, pandas 2.1+, numpy 1.26+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R packages cited for reference (exact2x2, Exact, ratesci): use `packageVersion()` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Categorical Association Tests for Clinical Data

**"Test association between categorical variables"** -> Determine whether treatment and a categorical clinical outcome are statistically independent (or that marginal proportions agree, for paired data) using a test calibrated to the design, the sample size, and the regulatory question.

## Algorithmic Taxonomy

| Test | Design | Asymptotic / exact | Conditioning | Strength | Fails when |
|------|--------|--------------------|--------------|----------|------------|
| Pearson chi-square (no continuity correction) | Independent groups, any RxC | Asymptotic | None | Standard for n>=40 with all expected counts >=5; matches Miettinen-Nurminen score CI | Any expected cell <1; >20% of cells with expected <5 (Cochran 1954) |
| Fisher's exact (conditional) | Independent 2x2 | Exact | Conditions on BOTH margins | Exact small-sample guarantee on level | Conservative (true alpha << nominal); discards information by double conditioning (Mehta-Senchaudhuri 2003) |
| Boschloo's exact | Independent 2x2 | Exact unconditional | Conditions on ONE margin only | Uniformly more powerful than Fisher (Boschloo 1970; Mehta-Senchaudhuri 2003); preserves nominal alpha exactly | Computationally heavier; RxC extensions limited |
| Barnard's exact | Independent 2x2 | Exact unconditional | Conditions on ONE margin only | Maximises nuisance parameter; well-calibrated | Slightly less powerful than Boschloo on average; compute scales O(n^2) |
| CMH (Mantel-Haenszel) | Stratified independent groups | Asymptotic | Conditions within strata | Tests common-OR null across strata; pooled OR estimator | Assumes no qualitative interaction; misleading when ORs reverse direction across strata |
| Breslow-Day | Stratified independent groups | Asymptotic | Within strata | Tests homogeneity of stratum ORs | Underpowered with few strata or sparse strata; non-significance does NOT prove homogeneity |
| McNemar (asymptotic, no continuity correction) | Paired binary | Asymptotic | Conditions on discordant pairs | Fagerland 2013 default; outperforms exact conditional | Discordant pair count b+c < 25 (chi-square approximation breaks) |
| Mid-p McNemar | Paired binary | Quasi-exact | Discordant pairs | Fagerland-Lydersen-Laake 2013 recommended default; less conservative than exact conditional | Slight under-coverage tolerable at small b+c |
| Exact conditional McNemar (Liddell 1983) | Paired binary | Exact | Discordant pairs only | Guaranteed coverage | Over-conservative; loses power vs mid-p or unconditional |
| Suissa-Shuster exact unconditional | Paired binary | Exact unconditional | All N pairs | Uniformly more powerful than exact conditional McNemar; 20-40% smaller n for same power | Implementation only in R `exact2x2::mcnemarExactDP` and SAS macros |

**Postdoc reading:** Lydersen, Fagerland & Laake 2009 *Stat Med* 28:1159 ("Recommended tests for association in 2x2 tables") argues **Fisher's exact should be retired from routine use** in favour of Boschloo or asymptotic Pearson; Fagerland-Lydersen-Laake 2013 *BMC Med Res Methodol* 13:91 makes the parallel case for mid-p or asymptotic McNemar over exact conditional. Regulatory practice (FDA reviewers) is moving in this direction but Fisher and exact conditional McNemar remain entrenched in many SAPs by inertia.

## Decision Tree by Experimental Scenario

| Scenario | Recommended test | Why |
|----------|------------------|-----|
| Independent 2x2, all expected >=5, n>=40 | Pearson chi-square, `correction=False` | Standard asymptotic; Yates' continuity correction is overly conservative and now discouraged |
| Independent 2x2, expected <5 in any cell OR n<40 | Boschloo's exact (`scipy.stats.boschloo_exact`) | Uniformly more powerful than Fisher; preserves exact Type-I control |
| Independent RxC, expected <5 in >20% of cells | Permutation chi-square or Fisher-Freeman-Halton (R `coin::chisq_test(distribution = approximate())`) | Exact RxC asymptotic invalid; permutation preserves level |
| Stratified design (multi-site, multi-stratum randomisation) | CMH for the pooled test + Breslow-Day for homogeneity + per-stratum ORs | Stratification factor in randomisation MUST appear in analysis (Kahan-Morris 2012 *Stat Med* 31:328: ignoring inflates Type-I error by up to 30%) |
| Stratified with sign-reversing effect (Simpson's paradox suspected) | Always report stratum-specific ORs + visual diagnostic; consider logistic regression with interaction | CMH pooled OR can mask sign reversal; homogeneity test underpowered |
| Paired binary (pre/post on same subjects; matched case-control) | Asymptotic McNemar (`mcnemar(table, exact=False, correction=False)`) when b+c >= 25; mid-p when b+c < 25 | Fagerland 2013 simulations show mid-p and asymptotic outperform exact conditional |
| Matched-pair non-inferiority (especially diagnostics) | Suissa-Shuster exact unconditional via R `exact2x2` | 20-40% smaller n than exact conditional for same power |
| Composite endpoint (any of several events) | Logistic regression with covariate adjustment, not chi-square | Composite changes the estimand under ICH E9(R1); see clinical-biostatistics/effect-measures |

## Chi-Square Test (Pearson, no continuity correction)

**Goal:** Test whether treatment group and outcome category are independent under asymptotic Type-I control.

**Approach:** Construct a contingency table, verify expected cell counts, compute the Pearson chi-square statistic without Yates' continuity correction.

```python
from scipy.stats import chi2_contingency
import pandas as pd

table = pd.crosstab(df['treatment'], df['outcome'])
chi2, p, dof, expected = chi2_contingency(table, correction=False)
if (expected < 5).any():
    print('WARNING: switch to Boschloo (2x2) or permutation chi-square (RxC)')
```

**Cochran 1954 rule (the precise version, not the textbook caricature):** "no expected cell should be <1 AND no more than 20% of cells should have expected <5." The textbook "all >=5" rule is the conservative simplification. With well-balanced 2x2 trials this rarely matters; with sparse RxC tables it materially expands the asymptotic range. The R `chisq.test` issues a warning under the strict Cochran rule; Python users must check manually.

**Why Yates' correction is now discouraged:** continuity correction was introduced to approximate the exact distribution under H0 but inflates Type-II error by ~10% (D'Agostino, Casagrande, Pike 1988 *Stat Med* 7:347). Modern computing makes Boschloo's exact test cheap; the correct fix for sparse 2x2 is Boschloo, not chi-square + continuity.

## Fisher's Exact -- and why Boschloo is usually better

**Goal:** Test 2x2 association with exact Type-I control.

**Approach:** Use Fisher's exact only when historical SAP requires it; otherwise prefer Boschloo's test.

```python
from scipy.stats import fisher_exact, boschloo_exact

odds_ratio, p_fisher = fisher_exact(table.values, alternative='two-sided')

# Boschloo (uniformly more powerful than Fisher):
# n= controls Sobol sampling resolution for the null distribution (scipy 1.12+; default 32);
# higher = more precise p-value at higher CPU cost. NOT the sample size per arm.
result = boschloo_exact(table.values, alternative='two-sided', n=64)
p_boschloo = result.pvalue
```

**The conditioning critique (Mehta-Senchaudhuri 2003):** Fisher's exact conditions on both margins of the 2x2 table, discarding information about the marginal totals. Boschloo conditions on one margin only and treats the second as a nuisance to be maximised over -- recovering the discarded information. Power gain at n=10/arm is 16-20 percentage points for moderate effects. Boschloo *uses Fisher's p-value as its test statistic*, then computes the exact unconditional null distribution of that p-value -- so it is automatically at least as powerful as Fisher.

**Since scipy 1.10**, `fisher_exact` returns the sample (unconditional) odds ratio, not the conditional MLE. For the conditional MLE matching R's `fisher.test`, use `scipy.stats.contingency.odds_ratio(table, kind='conditional')`.

## Cochran-Mantel-Haenszel (Stratified)

**Goal:** Test treatment-outcome association while controlling for a stratification variable; quantify the common odds ratio across strata.

**Approach:** Construct per-stratum 2x2 tables, compute MH pooled OR and CMH test of H0: common-OR = 1; test homogeneity via Breslow-Day.

```python
from statsmodels.stats.contingency_tables import StratifiedTable
import pandas as pd

tables = []
for stratum in df['site'].unique():
    stratum_data = df[df['site'] == stratum]
    t = pd.crosstab(stratum_data['treatment'], stratum_data['outcome']).values
    if t.shape == (2, 2) and t.min() > 0:
        tables.append(t)

st = StratifiedTable(tables)
print(st.test_null_odds())          # CMH H0: common OR = 1
print(st.oddsratio_pooled)          # MH pooled OR
print(st.oddsratio_pooled_confint(method='normal'))
print(st.test_equal_odds())         # Breslow-Day H0: equal stratum ORs
```

### Per-method failure modes

**CMH -- Simpson's paradox masking**

- **Trigger:** Stratum-specific ORs reverse direction across strata while the pooled MH OR appears null or modestly different from 1.
- **Mechanism:** CMH pools weighted log-ORs; equal-magnitude opposite-sign ORs cancel.
- **Symptom:** Breslow-Day p < 0.05 with stratum ORs visually reversing.
- **Fix:** Report stratum-specific ORs as primary; the MH pooled estimate is not a valid summary. Move to logistic regression with treatment-by-stratum interaction.

**Breslow-Day -- low-power false reassurance**

- **Trigger:** Few strata (k<5) or sparse strata (mean cell count <10).
- **Mechanism:** Breslow-Day chi-square has k-1 df; with k=3 and modest heterogeneity, power can be <40%.
- **Symptom:** Breslow-Day p > 0.5 with stratum ORs visually heterogeneous on a forest plot.
- **Fix:** Always supplement with a forest plot of stratum-specific ORs. Use likelihood-ratio interaction test from logistic regression as a second check.

**CMH -- ignoring randomisation stratification factors**

- **Trigger:** Randomisation was stratified (sex, region, baseline severity) but the primary analysis pools across strata.
- **Mechanism:** Stratified randomisation removes between-stratum variability that the unstratified SE still counts.
- **Symptom:** Type-I error inflation up to 30% (Kahan-Morris 2012).
- **Fix:** Strata variables from randomisation MUST appear in analysis -- either CMH, logistic regression with strata, or stratified log-rank.

## McNemar's Test for Paired Binary Data

**Goal:** Test the null of marginal homogeneity (P(positive at time 1) = P(positive at time 2)) for paired binary observations.

**Approach:** Default to asymptotic McNemar without continuity correction when discordant pairs >=25; switch to mid-p when discordant pairs <25; reserve exact conditional only when regulator-mandated.

```python
from statsmodels.stats.contingency_tables import mcnemar
import numpy as np

# table[i,j] = count with outcome i at time 1 and j at time 2
table = np.array([[45, 15], [5, 35]])  # b=15, c=5 discordant

# Asymptotic, no continuity correction -- the Fagerland 2013 recommended default
result = mcnemar(table, exact=False, correction=False)
print(result.statistic, result.pvalue)

# Exact conditional (Liddell 1983) -- only when b+c is very small or required by SAP
result_exact = mcnemar(table, exact=True)
```

**Fagerland-Lydersen-Laake 2013 *BMC Med Res Methodol* 13:91 simulation findings:** mid-p McNemar and asymptotic McNemar (no continuity correction) outperform exact conditional McNemar across small-to-moderate samples. The exact conditional is *too conservative* because it conditions on a discrete margin (the discordant pair count). Their title is the methodological provocation -- "The McNemar test: asymptotic and mid-p are better than exact conditional."

**Suissa-Shuster 1991 *Biometrics* 47:361 exact unconditional** uses *all* N pairs (not just discordant) -- uniformly more powerful than exact conditional McNemar; sample sizes 20-40% smaller for the same power. Available in R `exact2x2::mcnemarExactDP`. Practically essential for matched-pair non-inferiority in diagnostic device trials.

## Effect Sizes for Categorical Data

**Goal:** Quantify association strength beyond p-values.

**Approach:** Phi for 2x2, Cramer's V for RxC; bias-corrected variants in pingouin.

```python
import numpy as np
import pingouin as pg

n = table.values.sum()
phi = np.sqrt(chi2 / n)
k = min(table.shape) - 1
cramers_v = np.sqrt(chi2 / (n * k))

# Pingouin with multiple test variants and bias correction:
expected, observed, stats = pg.chi2_independence(df, x='treatment', y='outcome')
# stats columns: test, lambda, chi2, dof, pval, cramer, power
```

**Cohen 1988 effect-size benchmarks:**

| df | Small | Medium | Large |
|----|-------|--------|-------|
| 1 | 0.10 | 0.30 | 0.50 |
| 2 | 0.07 | 0.21 | 0.35 |
| 3 | 0.06 | 0.17 | 0.29 |

Phi equals Cramer's V for 2x2 (k=1). For RxC, only Cramer's V is valid because phi can exceed 1.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Fisher exact p > 0.05 but Boschloo p < 0.05 | Fisher over-conservative via double-margin conditioning (Mehta-Senchaudhuri 2003); Boschloo recovers power | Cite Boschloo as primary; provide Fisher in appendix for transparency |
| Pearson chi-square p < 0.05 but Fisher exact p > 0.05 | Asymptotic approximation breaking down at small expected counts (Cochran rule violation) | Use Boschloo (exact unconditional, more powerful than Fisher); document expected-count diagnostic in SAP |
| CMH pooled OR not significant, stratum-specific ORs strongly differ | Simpson's paradox -- opposite-sign cancellation OR effect modification | Report stratum-specific ORs as primary; switch to logistic regression with treatment-by-stratum interaction; forest plot stratum ORs |
| Breslow-Day non-significant but stratum-OR forest plot visually heterogeneous | Low power of Breslow-Day with few/sparse strata | Cite low-power caveat; report LR interaction test from logistic as secondary; do NOT claim homogeneity |
| McNemar exact conditional p > 0.05 but mid-p McNemar p < 0.05 | Exact conditional over-conservative due to discrete-margin conditioning | Cite Fagerland-Lydersen-Laake 2013; mid-p or asymptotic recommended; exact conditional only when SAP-mandated |
| Suissa-Shuster unconditional p < exact conditional McNemar p | Unconditional uses all N pairs; conditional discards concordant pairs | Suissa-Shuster preferred for matched-pair NI (esp. diagnostic devices) due to 20-40% smaller n |
| Wald CI excludes null but Wilson/Newcombe CI overlaps null | Wald has poor coverage near 0 and 1 (Brown-Cai-DasGupta 2001) | Wilson/Newcombe/MN preferred; cite as regulatory standard |
| Multiple categorical secondary endpoints, correlation structure unclear | Bonferroni overly conservative; Hochberg requires PRDS (Sarkar 1998) | See clinical-biostatistics/multiplicity-graphical for Bretz-Maurer graphical procedures and PRDS check |

## Confidence Intervals for Proportions and Differences

For a single proportion, **Wald is bad** for small samples and extreme p (Brown-Cai-DasGupta 2001 *Stat Sci* 16:101 documents "chaotic" coverage with coverage dropping to 0.0 in extreme cells). Use Wilson score or Jeffreys.

For a 2x2 risk difference or risk ratio, the regulatory standard for CI is **Miettinen-Nurminen score-based** (1985 *Stat Med* 4:213) -- consistent with the Pearson chi-square test and accepted by FDA/EMA for noninferiority margins.

```python
from statsmodels.stats.proportion import proportion_confint, proportions_ztest

# Single proportion: Wilson is the modern default
ci = proportion_confint(45, 60, alpha=0.05, method='wilson')
# Also available: 'jeffreys', 'agresti_coull', 'beta' (Clopper-Pearson exact)

# Difference of proportions: Newcombe-Wilson hybrid / MOVER
from statsmodels.stats.proportion import confint_proportions_2indep
ci_diff = confint_proportions_2indep(45, 60, 30, 60, method='newcomb', alpha=0.05)
# 'wald' is discouraged; 'newcomb' (Newcombe-Wilson hybrid) and 'agresti-caffo' are calibrated
```

For Miettinen-Nurminen CIs (the regulatory standard for stratified RD or RR), use R `ratesci::scoreci(contrast='RD'|'RR', distrib='bin', stratified=TRUE)` -- there is no production-grade Python implementation as of 2026.

## Post-Hoc Pairwise Comparisons

```python
from statsmodels.stats.multitest import multipletests
from itertools import combinations

categories = df['outcome'].unique()
pvalues, comparisons = [], []
for cat1, cat2 in combinations(categories, 2):
    subset = df[df['outcome'].isin([cat1, cat2])]
    sub_table = pd.crosstab(subset['treatment'], subset['outcome'])
    _, p_val, _, _ = chi2_contingency(sub_table, correction=False)
    pvalues.append(p_val)
    comparisons.append(f'{cat1} vs {cat2}')

reject, adjusted_p, _, _ = multipletests(pvalues, method='holm')
```

`method='holm'` (FWER) for confirmatory; `method='fdr_bh'` for exploratory. **Critical bug:** `multipletests` default is `method='hs'` (Holm-Sidak), NOT Holm or Bonferroni -- always specify explicitly. The FDA Multiple Endpoints Final Guidance (October 2022) requires FWER control for key secondary endpoints in regulatory submissions; FDR is acceptable for exploratory subgroup screens only.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| n >= 40 for 2x2 chi-square | Cochran 1954 *Biometrics* 10:417 | Below this, asymptotic chi-square distribution approximation degrades regardless of expected counts |
| All expected >=5 OR <=20% with expected <5 AND none <1 | Cochran 1954 (strict) | Textbook "all >=5" is overconservative; the strict rule expands chi-square's valid range |
| Yates' correction discouraged | D'Agostino, Casagrande, Pike 1988 *Stat Med* 7:347 | Overly conservative; correct fix for sparse 2x2 is Boschloo's exact, not continuity-corrected chi-square |
| Discordant pairs >=25 for asymptotic McNemar | Fagerland-Lydersen-Laake 2013 *BMC MRM* 13:91 | Below this, chi-square approximation breaks; switch to mid-p, not exact conditional |
| Newcombe-Wilson / Miettinen-Nurminen for RD CI | Newcombe 1998a *Stat Med* 17:873 | Wald CI for RD has poor coverage and can produce limits outside [-1, 1] |
| Boschloo > Fisher for small 2x2 | Mehta-Senchaudhuri 2003; Lydersen-Fagerland-Laake 2009 *Stat Med* 28:1159 | Boschloo uniformly more powerful at same Type-I |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `multipletests(p)` returns Holm-Sidak adjusted p | Default method is 'hs' not Holm/Bonferroni | Always specify `method='holm'` or `'bonferroni'` explicitly |
| `Table2x2(crosstab.values)` gives reciprocal OR | `pd.crosstab` orders columns alphabetically; statsmodels expects event-positive column first | Reorder: `cross[[1, 0]]` or `cross[['Yes', 'No']]` |
| Fisher's exact in published paper, Boschloo missing | SAP inertia; reviewers unfamiliar with Boschloo | Cite Mehta-Senchaudhuri 2003 in the SAP; use Boschloo as primary with Fisher in appendix |
| `fisher_exact` returns "wrong" OR vs R | Since scipy 1.10, scipy returns sample (unconditional) OR; R returns conditional MLE | Use `scipy.stats.contingency.odds_ratio(table, kind='conditional')` to match R |
| CMH significant but stratum ORs reverse direction | Simpson's paradox; Breslow-Day underpowered | Forest plot stratum ORs; report stratum-specific as primary; switch to logistic with interaction |
| Yates' correction enabled by default in `chi2_contingency` | scipy default is `correction=True` for 2x2 | Always pass `correction=False` for Pearson chi-square |
| McNemar p-value much larger than expected | Default may be exact conditional in some packages; over-conservative | Use asymptotic without continuity correction (Fagerland 2013) |
| Stratified randomisation ignored in primary analysis | Common SAP error | Include strata in CMH, logistic, or stratified log-rank; ignoring inflates Type-I (Kahan-Morris 2012) |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Why not Fisher's exact?" | Cite Lydersen-Fagerland-Laake 2009; Boschloo is uniformly more powerful at same alpha. Provide Fisher p in appendix for direct comparison. |
| "Why no continuity correction?" | D'Agostino-Casagrande-Pike 1988 -- Yates' inflates Type-II by ~10%. The correct fix for sparse 2x2 is Boschloo's exact, not Yates'. |
| "Are these ORs collapsible?" | OR is non-collapsible (see clinical-biostatistics/effect-measures); marginal and conditional ORs differ even without confounding. Cite Permutt 2020. |
| "Why mid-p McNemar over exact conditional?" | Fagerland-Lydersen-Laake 2013 simulations show exact conditional is over-conservative; mid-p and asymptotic maintain nominal Type-I with better power. |
| "Adjustment for stratification factors?" | Per ICH E9 and FDA 2023 covariate adjustment guidance, strata from randomisation must appear in analysis. CMH or logistic with strata as covariates. |
| "What is the estimand?" | Per ICH E9(R1), categorical-test analyses target a specific estimand (treatment policy is implicit if all randomised analysed). Articulate explicitly. |

## References

- Boschloo RD. 1970. Raised conditional level of significance for the 2x2 table when testing the equality of two probabilities. *Stat Neerl* 24:1.
- Brown LD, Cai TT, DasGupta A. 2001. Interval estimation for a binomial proportion. *Stat Sci* 16:101-117.
- Cochran WG. 1954. Some methods for strengthening the common chi-squared tests. *Biometrics* 10:417-451.
- D'Agostino RB, Chase W, Belanger A. 1988. The appropriateness of some common procedures for testing the equality of two independent binomial populations. *Am Stat* 42:198-202.
- Fagerland MW, Lydersen S, Laake P. 2013. The McNemar test: asymptotic and mid-p are better than exact conditional. *BMC Med Res Methodol* 13:91.
- FDA. 2022. Multiple Endpoints in Clinical Trials -- Guidance for Industry. Federal Register Oct 2022.
- Kahan BC, Morris TP. 2012. Improper analysis of trials randomised using stratified blocks or minimisation. *Stat Med* 31:328-340.
- Liddell FDK. 1983. Simplified exact analysis of case-referent studies: matched pairs. *Appl Stat* 32:127-132.
- Lydersen S, Fagerland MW, Laake P. 2009. Recommended tests for association in 2x2 tables. *Stat Med* 28:1159-1175.
- Mehta CR, Senchaudhuri P. 2003. Conditional vs unconditional exact tests for comparing two binomials. (Cytel technical report; widely cited in subsequent literature.)
- Miettinen O, Nurminen M. 1985. Comparative analysis of two rates. *Stat Med* 4:213-226.
- Newcombe RG. 1998a. Interval estimation for the difference between independent proportions: comparison of eleven methods. *Stat Med* 17:873-890.
- Newcombe RG. 1998b. Improved confidence intervals for the difference between binomial proportions based on paired data. *Stat Med* 17:2635-2650.
- Permutt T. 2020. Do covariates change the estimand? *Stat Biopharm Res* 12:45-53.
- Suissa S, Shuster JJ. 1991. The 2x2 matched-pairs trial: exact unconditional design and analysis. *Biometrics* 47:361-372.

## Related Skills

- clinical-biostatistics/effect-measures - Detailed OR/RR/RD with modern CI methods (Wilson, Newcombe, Miettinen-Nurminen)
- clinical-biostatistics/logistic-regression - Regression alternative with covariate adjustment; modified Poisson for RR
- clinical-biostatistics/subgroup-analysis - Stratified analysis with interaction terms and HTE methods
- clinical-biostatistics/multiplicity-graphical - Bretz-Maurer graphical procedures for confirmatory multiplicity
- clinical-biostatistics/trial-reporting - CONSORT 2025 and ICH E9(R1) reporting of categorical analyses
- experimental-design/multiple-testing - General multiple testing correction methods
