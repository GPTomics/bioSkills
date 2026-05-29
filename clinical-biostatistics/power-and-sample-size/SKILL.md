---
name: bio-clinical-biostatistics-power-sample-size
description: Computes sample size and power for clinical trials including continuous, binary, and time-to-event endpoints; superiority, non-inferiority, and equivalence designs; FDA 2016 non-inferiority margin selection with M1/M2 framework; Schoenfeld 1981 and Lakatos 1988 for survival; Schuirmann TOST and 80-125% bioequivalence; minimum clinically important difference (MCID) vs δ distinction. Use when justifying trial size in protocol or SAP per CONSORT 2025 item 7.
tool_type: mixed
primary_tool: statsmodels
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: statsmodels 0.14+, scipy 1.12+, numpy 1.26+, pandas 2.1+. R packages cited: pwr, gsDesign (Anderson/Merck), gsDesign2, rpact (Wassmer/Brannath), presize, npsurvSS, nph, simtrial.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Power and Sample Size for Clinical Trials

**"Justify the trial's sample size"** -> Compute the n needed to detect a pre-specified alternative δ with power 1-β at significance α, accounting for endpoint distribution, design (superiority/NI/equivalence), expected dropout, multiplicity, and stratification — and distinguish δ (the effect the trial is powered to detect) from MCID (the clinically meaningful difference).

## The Foundational Distinction -- δ vs MCID

**δ (the alternative effect):** what the trial is *powered to detect*. Usually set above the MCID because sponsors want a strong signal that exceeds noise + design uncertainty.

**MCID (Minimum Clinically Important Difference):** the smallest effect size considered clinically meaningful. Jaeschke-Singer-Guyatt 1989 *Control Clin Trials* 10:407 (anchor-based) and Norman-Sloan-Wyrwich 2003 *Med Care* 41:582 ("the remarkable universality of half a standard deviation") established the modern conventions.

**Confusing the two has produced both:**

- Underpowered trials where sponsor sets δ = MCID and gets a CI straddling zero
- Overgenerous NI margins where sponsor sets M2 = full MCID (NI margin should be a fraction of MCID)

**Postdoc rule of thumb:** for superiority, δ >= 1.5 × MCID; for NI, M2 <= 0.5 × MCID.

## Algorithmic Taxonomy

| Design | Formula / approach | Software | Strength | Fails when |
|--------|--------------------|----------|----------|------------|
| Two-sample t-test, continuous | Cohen's d; n = 2 × (z_α/2 + z_β)² / d² | `pwr::pwr.t.test` (R); `statsmodels.power.tt_ind_solve_power` (Py) | Standard | Heteroscedasticity; non-normal outcomes |
| Two-sample proportions (Fleiss) | Asymptotic normal approximation with/without continuity correction | `power.prop.test` (R) -- uncorrected; `pwr::pwr.2p.test`; statsmodels | Standard | n < 100/arm: continuity correction debate (D'Agostino 1988) |
| Survival (Schoenfeld 1981) | events ≈ 4(z_α/2 + z_β)² / (log HR)² for 1:1 | `gsDesign::nSurv`; `npsurvSS::size_two_arm` | Standard PH-conformant | PH violated (immuno-oncology) under-estimates by 20-50% |
| Survival under non-PH (Lakatos 1988) | Markov chain accommodating time-varying HR, accrual, dropout | `gsDesign::nSurv`; `npsurvSS`; `simtrial` | Handles immuno-oncology delayed effects | Requires explicit specification of HR(t) and accrual |
| MaxCombo SS under NPH | Simulation-based; pre-specify weight family | `nphRCT`; `simtrial` | Robust to NPH pattern | Computationally heavier |
| Non-inferiority fixed-margin | n = (z_α + z_β)² × variance / M² | `pwr::pwr.t2n.test` adapted; `rpact::getSampleSizeMeans` | Pre-discounted M | Constancy assumption violation invisible |
| Non-inferiority synthesis | Pool historical control-vs-placebo + current test-vs-control | `gsDesign::ssTwoArmTest` | More efficient than fixed-margin | Constancy assumption MUST hold exactly |
| Equivalence TOST | Two one-sided tests at α each | `pwr::pwr.t.test` adapted; `presize` | No multiplicity adjustment needed | Wrong question when superiority/NI is intended |
| Group-sequential | Lan-DeMets spending function | `rpact`; `gsDesign` | Interim analyses; early stopping | More complex SAP |
| Sample-size re-estimation (Mehta-Pocock) | Promising-zone conditional power | `rpact::getSampleSizeMeans` with reestimation | Recovers power if interim shows promise | Unblinded SSR scares FDA |
| Cluster-randomised | Adjust for design effect = 1 + (m-1)ICC | `clusterPower`; `pwr` adapted | Standard | ICC misspecification |

**Postdoc reading list:**

- Fleiss JL 1981 *Statistical Methods for Rates and Proportions* (with/without continuity correction tables)
- Schoenfeld DA 1981 *Biometrika* 68:316 (canonical TTE formula)
- Lakatos E 1988 *Biometrics* 44:229 (Markov chain SS for complex survival)
- Schuirmann DJ 1987 *J Pharmacokinet Biopharm* 15:657 (TOST)
- Jaeschke R, Singer J, Guyatt GH 1989 *Control Clin Trials* 10:407 (MCID anchor-based)
- Norman GR, Sloan JA, Wyrwich KW 2003 *Med Care* 41:582 (0.5 SD heuristic)
- Snapinn SM 2000 *Curr Control Trials Cardiovasc Med* 1:19 (NI biocreep)
- Temple R, Ellenberg SS 2000 *Ann Intern Med* 133:455 + 133:464 (NI assay sensitivity; **Ann Intern Med NOT NEJM** — common citation error)
- Hung HMJ, Wang SJ, O'Neill RT 2005 *Biom J* (NI within-trial Type-I not guaranteed)
- Mehta CR, Pocock SJ 2011 *Stat Med* 30:3267 (promising zone)
- Lin RS, Lin J, Roychoudhury S et al 2020 *Stat Biopharm Res* (NPH Working Group)

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Continuous outcome, two-arm parallel, superiority | Fleiss / Cohen's d formula via `power.t.test` or statsmodels | Standard; cite Fleiss 1981 |
| Binary outcome, n > 100/arm | Uncorrected normal approximation via `power.prop.test`; statsmodels `power_proportions_2indep` | Adequate for n > 100 |
| Binary outcome, n < 100/arm | Fisher exact-based simulation OR Fleiss with continuity correction | Continuity correction debate; cite D'Agostino 1988 if uncorrected |
| Time-to-event with PH | Schoenfeld 1981 via `gsDesign::nSurv` | Standard; pre-specify hazards |
| Time-to-event under expected NPH (immuno-oncology) | Lakatos 1988 OR simulation via `simtrial::simtrial`/`nphRCT` | Schoenfeld under-estimates by 20-50%; cite Lin 2020 |
| Non-inferiority continuous | Fixed-margin with M2 = 0.5 × historical M1 (discounted) | FDA 2016 NI guidance; cite Temple-Ellenberg 2000 |
| Non-inferiority binary | Fixed-margin with Miettinen-Nurminen CI for RD | Cite EMA NI guidance; MN-CI standard |
| Bioequivalence (Cmax/AUC) | TOST with 80-125% margins on geometric mean ratio | FDA 1992 BE guidance; Schuirmann 1987 |
| Cluster-randomised | n × design effect = 1 + (m-1)ICC | Cite Murray 1998 cluster RCT methodology |
| Group-sequential | Lan-DeMets spending function | `rpact` or `gsDesign` |
| Pilot for sample size re-estimation | Blinded SSR (Friede-Kieser 2006) | No Type-I inflation; safer than unblinded |
| Promising-zone reestimation | Mehta-Pocock 2011 with pre-specified increase rule | Recovers power; cite caveat re Jennison-Turnbull 2015 critique |

## Continuous Outcomes -- Two-Sample t-Test

```python
from statsmodels.stats.power import tt_ind_solve_power

# Solve for n per group
n = tt_ind_solve_power(
    effect_size=0.5,  # Cohen's d = (mu1 - mu2) / sigma
    alpha=0.05,
    power=0.80,
    alternative='two-sided'
)
print(f'n per arm = {np.ceil(n):.0f}')

# Solve for power given n
power = tt_ind_solve_power(effect_size=0.5, alpha=0.05, nobs1=100, alternative='two-sided')
```

**Cohen's d benchmarks:** small = 0.2, medium = 0.5, large = 0.8. Choose d to detect based on prior literature or clinically meaningful effect, NOT post-hoc to fit the affordable n.

**Inflate for dropout:** if dropout rate is q, multiply final n by 1/(1-q). For q=0.20, n_total = n/0.80.

**Stratified randomisation efficiency:** if randomisation stratifies on prognostic factors with combined R² = r against outcome, the effective n is n/(1-r²). Senn 2013 *Stat Med* 32:1439 makes the precision argument explicit.

## Binary Outcomes -- The Continuity Correction Debate

```python
from statsmodels.stats.power import NormalIndPower
import math

# Without continuity correction
p1, p2 = 0.30, 0.20
pbar = (p1 + p2) / 2
z_alpha = 1.96  # two-sided alpha=0.05
z_beta = 0.84   # power=0.80
n = ((z_alpha * math.sqrt(2 * pbar * (1 - pbar)) +
      z_beta * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2))) / (p1 - p2)) ** 2
print(f'n per arm (uncorrected) = {math.ceil(n)}')

# Continuity-corrected (Fleiss formula)
correction = (1 + math.sqrt(1 + 4 / (n * abs(p1 - p2)))) ** 2 / 4
n_corrected = math.ceil(n * correction)
print(f'n per arm (continuity corrected) = {n_corrected}')
```

**The continuity correction debate** (D'Agostino, Casagrande, Pike 1988 *Stat Med* 7:347):

- Pro (Fleiss, Yates): better matches the exact distribution under H0
- Con: overly conservative, wastes ~10% sample size; modern computing makes exact (Fisher-Irwin, Boschloo) tests feasible

**R's `power.prop.test` uses uncorrected normal approximation** — can over-state power 5-10% in small samples (n<100). For confirmatory work with n<100, use simulation-based SS via `simr` or `presize`.

## Survival (Time-to-Event)

### Schoenfeld 1981 -- the canonical formula

```python
import math

# Events needed for two-sample log-rank under PH (Schoenfeld 1981)
from scipy.stats import norm

def schoenfeld_events(hr, alpha=0.05, power=0.80, allocation_ratio=1.0, two_sided=True):
    """
    events = (z_crit + z_beta)^2 / (p*(1-p) * log(HR)^2)

    where z_crit = z_{1-alpha/2} for two-sided test OR z_{1-alpha} for one-sided.
    """
    z_crit = norm.ppf(1 - alpha/2) if two_sided else norm.ppf(1 - alpha)
    z_beta = norm.ppf(power)
    p = 1 / (1 + allocation_ratio)
    return math.ceil((z_crit + z_beta) ** 2 / (p * (1 - p) * math.log(hr) ** 2))

# For HR=0.70, alpha=0.05 two-sided, 80% power, 1:1 -> approx 247 events
events_needed = schoenfeld_events(hr=0.70)
print(f'Events needed for HR=0.70: {events_needed}')
```

**Events drive power, not subjects.** Convert to n via expected event rate and follow-up:

```python
# n needed = events / overall event probability over follow-up
expected_overall_event_prob = 0.40  # from pilot or historical
n_per_arm = math.ceil(events_needed / (2 * expected_overall_event_prob))
# Inflate for dropout
n_per_arm_with_dropout = math.ceil(n_per_arm / (1 - 0.15))  # 15% dropout
```

### Lakatos 1988 -- complex survival

Use R `gsDesign::nSurv()` or `npsurvSS::size_two_arm()` for complex scenarios with:

- Time-varying hazard ratios (delayed effect)
- Non-uniform accrual
- Different dropout rates per arm
- Cure fractions

```r
library(gsDesign)
n_lakatos <- nSurv(
    lambdaC = 0.04,     # control hazard per month
    hr = 0.70,           # treatment HR
    eta = 0.005,         # dropout hazard per month
    T = 24,              # total study duration in months
    minfup = 12,         # minimum follow-up
    accrualTime = 12,    # accrual duration
    alpha = 0.025,       # one-sided
    beta = 0.10          # power = 90%
)
print(n_lakatos)
```

### Under non-PH (immuno-oncology with delayed effect)

**Schoenfeld 1981 ASSUMES PH.** Under immuno-oncology with delayed separation, the formula under-estimates required events by 20-50%. The fix:

```r
# Simulate under expected hazard pattern (Lin 2020 NPH Working Group recommendation)
library(simtrial)
sim_result <- simtrial(
    n_per_arm = 250,
    enroll_rate = piecewise_enroll(),
    fail_rate = piecewise_fail(  # delayed-effect specification
        duration = c(6, 18),
        fail_rate = c(0.04, 0.04),  # control
        hr = c(1.0, 0.5),            # treatment HR is 1 for first 6 months, then 0.5
        dropout_rate = c(0.005, 0.005)
    ),
    total_duration = 36
)
# Compute empirical power via MaxCombo or weighted log-rank
```

## Non-Inferiority Designs -- The FDA 2016 Framework

**FDA NI Trials Guidance (November 2016)** establishes the canonical M1/M2 framework:

- **M1** = entire effect of active control vs placebo (lower bound of CI from historical meta-analysis)
- **M2** = clinically acceptable loss, typically M1/2 (50% effect retention) but as conservative as M1/4 for mortality endpoints

**The "double discount" trap:** M2 = 0.5 × (lower CI bound of M1) effectively requires 50% retention of a discounted historical estimate, yielding a much tighter margin than naive point-estimate retention. This has effectively halted new antibiotic NDA development under stringent NI margins.

```python
# Fixed-margin NI sample size (binary endpoint)
import math

p_control = 0.85         # historical/observed control success rate
p_test = 0.85            # null assumption: test = control (NI null is test - control < -M2)
margin_M2 = 0.05         # pre-specified clinically acceptable loss
z_alpha = 1.645          # one-sided alpha=0.025
z_beta = 0.84            # power=0.80

# Approximate: n per arm
variance_term = p_control * (1 - p_control) + p_test * (1 - p_test)
n = math.ceil((z_alpha + z_beta) ** 2 * variance_term / margin_M2 ** 2)
print(f'NI sample size per arm: {n}')
```

**Use Miettinen-Nurminen score CI for RD** (the regulatory standard) when computing the NI margin from observed proportions. Cite `ratesci::scoreci` (R) for production.

**Constancy assumption:** the historical placebo-vs-active effect must persist in the NI trial. If standard of care drifted, M1 over-estimates and margin is too liberal. **Hung-Wang-O'Neill 2005** *Biom J* critique: NI trials have NO within-trial Type-I error guarantee — alpha is conditional on constancy.

**Assay sensitivity (Temple-Ellenberg 2000 *Ann Intern Med* 133:455 + 133:464):** NI trials have only *external* assay sensitivity (inferred from historical data + constancy); placebo-controlled trials have *internal* assay sensitivity. This is why NI trials are second-best when placebo is ethical.

## Equivalence Designs and Bioequivalence

**TOST (Two One-Sided Tests; Schuirmann 1987):** reject H0 of inequivalence iff both one-sided tests reject at α each — equivalent to the 1-2α CI lying within (-δ, +δ). Closed under intersection-union -> no multiplicity adjustment needed despite two tests.

```python
# TOST sample size for continuous outcome
import math

def tost_sample_size(mu_diff, sigma, margin, alpha=0.05, power=0.80):
    """Sample size per arm for TOST equivalence test."""
    z_alpha = 1.645  # one-sided
    z_beta = 0.84
    effect = abs(margin - abs(mu_diff))
    if effect <= 0:
        raise ValueError('|mu_diff| must be < margin')
    n = math.ceil(2 * sigma ** 2 * (z_alpha + z_beta) ** 2 / effect ** 2)
    return n

n_eq = tost_sample_size(mu_diff=0.0, sigma=10, margin=5)
print(f'Equivalence sample size per arm: {n_eq}')
```

**FDA 1992 bioequivalence (80-125% rule):** geometric mean ratio of Cmax and AUC must have 90% CI within (0.80, 1.25). Implemented in `PowerTOST` R package for crossover BE designs.

## Crossover Designs -- Bioequivalence and Repeated-Measures

**Crossover trials** randomise each subject to receive both treatments in different periods. Power calculation is fundamentally different from parallel: the within-subject SD (sigma_w) drives power, not between-subject SD (sigma_b). With reasonable carryover-free designs, crossover requires ~25-50% the n of parallel for the same precision.

```python
# Two-period crossover sample size (continuous endpoint)
import math

def crossover_n_continuous(mean_diff, sd_within, alpha=0.05, power=0.80):
    """n per sequence (so total n = 2*n) for two-period crossover."""
    z_alpha = 1.96  # two-sided
    z_beta = 0.84   # power=0.80
    n_per_seq = math.ceil(2 * (z_alpha + z_beta)**2 * sd_within**2 / mean_diff**2)
    return n_per_seq, 2 * n_per_seq  # n_per_seq, n_total
```

### Bioequivalence (FDA 1992 / EMA 2010 framework)

**Average bioequivalence:** geometric mean ratio (GMR) of test/reference Cmax and AUC must have 90% CI within (0.80, 1.25). Log-transform pharmacokinetic parameters; analyse with mixed-effects ANOVA (period, sequence, treatment, subject random).

```r
library(PowerTOST)
# Sample size for 2x2 crossover bioequivalence with CV (within-subject)
n_be <- sampleN.TOST(
    alpha = 0.05,
    targetpower = 0.80,
    theta0 = 0.95,        # expected GMR
    theta1 = 0.80,        # lower BE bound
    theta2 = 1.25,        # upper BE bound
    CV = 0.25,             # within-subject CV
    design = '2x2x2'      # standard 2-period 2-sequence
)
print(n_be)
```

**Highly variable drugs (CV > 30%):** standard 80-125% fails with feasible n. Use **scaled average bioequivalence (SABE)** per FDA 2010:
- Reference-scaled limits: theta = exp(0.760 * sigma_WR) when sigma_WR > 0.294
- Widens BE bounds proportional to within-subject variability of reference
- `PowerTOST::sampleN.scABEL()` for sample size

### Carryover assessment (Grizzle 1965)

**Grizzle test** for carryover: compares baseline-adjusted period 1 vs period 2 effects. **Now controversial:** Senn 2002 *Statistics in Medicine* showed Grizzle test has poor power; significant carryover should drive design choice (washout extension), not analysis switch.

**Modern operational rule:** pre-specify adequate washout (>=5 half-lives); do NOT routinely test for carryover in primary analysis. If carryover is biologically plausible, use parallel design or extend washout.

### Period effect

Period effects (calendar/learning) are estimable in 2x2 crossover. Standard analysis includes period as fixed effect:

```r
library(nlme)
fit <- lme(response ~ treatment + period, random = ~1 | subject, data = df)
```

If period significant, treatment effect is still unbiased (orthogonal in balanced 2x2); report period effect for transparency.

### Decision tree -- when crossover is appropriate

| Scenario | Use crossover? | Why |
|----------|----------------|-----|
| Bioequivalence of pharmacokinetic parameters | YES | Standard regulatory; within-subject precision much higher |
| Chronic stable disease (HTN, GERD, asthma) | YES | Reversible response; ~50% sample size savings |
| Acute disease (sepsis, MI) | NO | Cannot retreat same patient |
| Curative intent (oncology, surgery) | NO | Treatment alters disease state irreversibly |
| Long-half-life drug (>1 week) | NO | Washout impractical |
| Patient-reported outcomes with strong period bias | CAUTION | Period bias may dominate |

## MCID -- Anchor-Based vs Distribution-Based

| Method | Source | Description |
|--------|--------|-------------|
| Anchor-based | Jaeschke 1989 | Patient-reported anchor question; MCID = mean change in subjects who report "small important" change |
| 1 SEM (Wyrwich 1999) | *J Clin Epidemiol* | MCID ≈ standard error of measurement |
| 0.5 SD (Norman 2003) | *Med Care* 41:582 | "Remarkable universality of half a standard deviation" across PROMs |
| 0.5 SD baseline | Common rule of thumb | Quick approximation for unknown PROMs |

**Postdoc rule:** δ in power calculation should EXCEED MCID (often 1.5-2× MCID) so that trials are powered to detect effects clearly larger than the noise threshold.

## Multiple Endpoints -- See multiplicity-graphical

For co-primary endpoints, power-adjust per FDA Multiple Endpoints Guidance (October 2022):

- **Co-primary (all-must-win):** each endpoint at full alpha but joint power = product -> inflate n
- **Multiple primary (any-wins):** alpha split (e.g., Bonferroni or graphical)
- **Hierarchical:** test in order; if any fail, downstream cannot be claimed

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Schoenfeld n much smaller than Lakatos simulation n | PH violated under expected HR(t) pattern (immuno-oncology delayed effect) | Use Lakatos OR simulation under expected HR(t); cite Lin 2020 NPH Working Group |
| Continuity-corrected vs uncorrected proportion SS differ by ~10% | Asymptotic approximation debate (D'Agostino 1988) | For confirmatory n < 100/arm, use simulation-based SS; document correction approach in SAP |
| NI margin per FDA double discount (50% of M1 lower CI) vs sponsor proposed margin (50% of M1 point estimate) | Sponsor uses point estimate; FDA expects lower-CI discount | FDA 2016 NI guidance: M2 from lower CI bound is the regulatory expectation |
| Bayesian predictive probability of success vs frequentist power give different n | Bayesian integrates prior uncertainty; frequentist conditions on hypothesised δ | Bayesian PPoS more honest about pre-trial uncertainty; frequentist preferred for confirmatory regulatory; report both for adaptive contexts |
| TOST and CI-inclusion approach give different equivalence conclusions | TOST is two-one-sided at α; CI is 1-2α — mathematically equivalent so should NOT disagree | If disagree, check α conventions (Schuirmann 1987: TOST at α, CI at 1-2α) |
| Sponsor's MCID vs published anchor-based MCID differ | Different anchor questions or population | Cite primary MCID source; cross-check with Norman 2003 0.5 SD heuristic; declare in SAP |
| Stratified-randomisation SS calculation vs unstratified differ substantially | Stratification efficiency gain via Senn 2013 formula | Use stratified formula; cite efficiency = 1/(1-r²) where r² is R² of strata against outcome |
| Cluster-RCT SS underpowered after enrollment | ICC misspecified at design (most common SS failure in cluster RCTs) | Use historical ICC + sensitivity analysis at ±50% of estimate; cite Murray 1998 |

## Per-Method Failure Modes

### Schoenfeld under non-PH

- **Trigger:** Immunotherapy or other delayed-effect mechanism
- **Mechanism:** Schoenfeld assumes constant log-HR
- **Symptom:** Trial under-powered; observed events insufficient
- **Fix:** Lakatos 1988 or simulation under expected HR(t); cite Lin 2020 NPH Working Group

### Continuity correction debate

- **Trigger:** Small-sample binary design
- **Mechanism:** Uncorrected normal approx over-states power; corrected wastes ~10%
- **Symptom:** Trial reach significance at lower n than expected
- **Fix:** Use simulation-based SS for n<100; pre-specify correction approach in SAP

### NI margin too liberal

- **Trigger:** M2 set to MCID rather than 50% of historical M1
- **Mechanism:** Allows clinically unacceptable inferiority within "non-inferior" CI
- **Symptom:** Regulators flag margin as inappropriate at scientific advice
- **Fix:** M2 <= 0.5 × MCID; M2 <= 0.5 × historical M1 lower bound (FDA double discount)

### Constancy assumption violated

- **Trigger:** NI trial with new standard of care differing from historical
- **Mechanism:** M1 derived from historical data no longer reflects current control efficacy
- **Symptom:** Active control performs worse than expected; "successful" NI may mask true inferiority
- **Fix:** Synthesis approach with formal sensitivity to constancy; consider three-arm trial with small placebo

### MCID confused with δ

- **Trigger:** Sponsor sets δ = MCID in power calculation
- **Mechanism:** Trial powered to detect smallest meaningful effect, with no margin for noise/dropout
- **Symptom:** Effect estimate near MCID with wide CI straddling zero
- **Fix:** δ >= 1.5 × MCID; document rationale in SAP

### Promising-zone "stealth alpha inflation"

- **Trigger:** Mehta-Pocock SSR applied without sufficient simulation
- **Mechanism:** Jennison-Turnbull 2015 critique: increase rule + statistical test combination may inflate alpha beyond ~0.001
- **Symptom:** Independent reanalysis finds Type-I 5.3% vs nominal 5%
- **Fix:** Pre-specify increase rule transparently; report simulation operating characteristics

### Unblinded SSR -- DMC firewall failure

- **Trigger:** Interim effect estimate leaked beyond IDMC
- **Mechanism:** Sponsor inference from sample-size increase decision reveals direction of interim effect
- **Symptom:** Regulator audit reveals unblinding
- **Fix:** Strict firewall; only "increase / no increase" communicated to sponsor; document SOP

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| δ >= 1.5 × MCID for superiority | Postdoc rule; Norman 2003 implication | δ = MCID gives no margin for sampling variation |
| M2 <= 0.5 × historical M1 (FDA double discount) | FDA NI 2016 guidance | Conservative retention against biocreep |
| M2 <= 0.5 × MCID for NI | Standard regulatory expectation | Acceptable loss < clinically meaningful difference |
| Schoenfeld under non-PH under-estimates by 20-50% | Lin 2020 NPH WG | Switch to Lakatos or simulation |
| Continuity correction wastes ~10% sample size | D'Agostino 1988 *Stat Med* 7:347 | Modern computing makes exact tests cheap |
| 90% CI within (0.80, 1.25) for BE | FDA 1992 BE guidance | Geometric mean ratio of Cmax/AUC |
| Unblinded SSR Type-I requires CHW weights | Cui-Hung-Wang 1999 | Naive increase inflates Type-I |
| Promising-zone CP range ~30-80% | Mehta-Pocock 2011 *Stat Med* 30:3267 | Mathematical calibration for Type-I preservation |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Trial under-powered for immunotherapy | Schoenfeld formula used under non-PH | Lakatos or simulation; cite Lin 2020 |
| δ = MCID, CI straddles zero | δ-MCID confusion | Pre-specify δ >= 1.5 × MCID |
| NI "success" but active much worse than historical | Constancy violation | Cite Hung 2005; consider synthesis or 3-arm |
| Unblinded SSR with naive sample increase | Type-I inflation | CHW-weighted test; cite Cui-Hung-Wang 1999 |
| `power.prop.test` over-states power n < 100 | Uncorrected normal approximation | Simulation-based SS or Fleiss continuity correction |
| BE failure after large n with point estimate near 1.0 | Variability higher than assumed | Re-estimate from pilot; consider replicate design for highly variable drugs |
| Promising-zone tipping to favourable not pre-specified | Operational confusion | Document increase rule in SAP; pre-specify cap |
| Cluster RCT n based on individual-level | Design effect ignored | n × (1 + (m-1)ICC); cite Murray 1998 |
| TOST sample size with non-zero true difference | Treating null as abs(mu_diff)=0 | Include realistic mu_diff in calculation |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Is δ realistic?" | Justified from prior phase 2; δ exceeds MCID by factor X; pre-specified |
| "Why this NI margin?" | M2 = 50% of historical lower-CI M1 per FDA 2016 double discount; M2 < MCID |
| "Constancy assumption?" | Reviewed historical placebo-vs-active effect over time; stable; sensitivity to lower-end constancy provided |
| "Schoenfeld vs Lakatos?" | PH plausible from phase 2 KM; pre-specified Schoenfeld; Lakatos simulation under expected delayed-effect as sensitivity |
| "Continuity correction?" | n > 200/arm; uncorrected normal approximation valid (Fleiss 1981) |
| "Bioequivalence variability assumption?" | Conservative (upper bound of historical CV); pilot will trigger re-estimation if CV >X% |
| "Cluster ICC source?" | Historical ICC from prior cluster trial in similar setting; sensitivity at ICC ±50% provided |
| "MCID source?" | Anchor-based MCID from validation study (Jaeschke 1989); supported by 0.5 SD heuristic (Norman 2003) |
| "Promising zone or just blinded SSR?" | Promising zone with pre-specified CP boundaries; CHW weights for Type-I control |

## References

- Cui L, Hung HMJ, Wang SJ. 1999. Modification of sample size in group sequential clinical trials. *Biometrics* 55:853-857.
- D'Agostino RB, Chase W, Belanger A. 1988. The appropriateness of some common procedures for testing the equality of two independent binomial populations. *Am Stat* 42:198-202.
- FDA. 2016. Non-Inferiority Clinical Trials to Establish Effectiveness. Final Guidance.
- FDA. 2022. Multiple Endpoints in Clinical Trials. Final Guidance.
- Fleiss JL. 1981. *Statistical Methods for Rates and Proportions* (2nd ed). Wiley.
- Friede T, Kieser M. 2006. Sample size recalculation in internal pilot study designs. *Biom J* 48:537-555.
- Hung HMJ, Wang SJ, O'Neill RT. 2005. A regulatory perspective on choice of margin and statistical inference issue in non-inferiority trials. *Biom J* 47:28-36.
- Jaeschke R, Singer J, Guyatt GH. 1989. Measurement of health status: ascertaining the minimal clinically important difference. *Control Clin Trials* 10:407-415.
- Lakatos E. 1988. Sample sizes based on the log-rank statistic in complex clinical trials. *Biometrics* 44:229-241.
- Lin RS, Lin J, Roychoudhury S, Anderson KM, Hu T, Huang B, Leon LF, Liao JJZ, Liu R, Luo X, Mukhopadhyay P, Qin R, Tatsuoka K, Wang X, Wang Y, Zhu J, Chen TT, Iacona R. 2020. Alternative analysis methods for time to event endpoints under nonproportional hazards: a comparative analysis. *Stat Biopharm Res*.
- Mehta CR, Pocock SJ. 2011. Adaptive increase in sample size when interim results are promising. *Stat Med* 30:3267-3284.
- Norman GR, Sloan JA, Wyrwich KW. 2003. Interpretation of changes in health-related quality of life: the remarkable universality of half a standard deviation. *Med Care* 41:582-592.
- Schoenfeld DA. 1981. The asymptotic properties of nonparametric tests for comparing survival distributions. *Biometrika* 68:316-319.
- Schuirmann DJ. 1987. A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. *J Pharmacokinet Biopharm* 15:657-680.
- Senn S. 2013. Seven myths of randomisation in clinical trials. *Stat Med* 32:1439-1450.
- Snapinn SM. 2000. Noninferiority trials. *Curr Control Trials Cardiovasc Med* 1:19-21.
- Temple R, Ellenberg SS. 2000. Placebo-controlled trials and active-control trials in the evaluation of new treatments. *Ann Intern Med* 133:455-470.
- Wyrwich KW, Tierney WM, Wolinsky FD. 1999. Further evidence supporting an SEM-based criterion for the identification of meaningful intra-individual changes in health-related quality of life. *J Clin Epidemiol* 52:861-873.

## Related Skills

- clinical-biostatistics/survival-analysis - TTE-specific sample size (Schoenfeld, Lakatos, simulation)
- clinical-biostatistics/effect-measures - δ on OR/RR/RD scales
- clinical-biostatistics/categorical-tests - Binary endpoint test selection
- clinical-biostatistics/multiplicity-graphical - Power adjustment for co-primary endpoints
- clinical-biostatistics/adaptive-designs - SSR, promising zone, group-sequential
- clinical-biostatistics/bayesian-trials - Bayesian SS via predictive probability of success
- clinical-biostatistics/trial-reporting - CONSORT 2025 item 7 (sample size justification)
- experimental-design/sample-size - General sample-size methods
- experimental-design/power-analysis - General power methods
