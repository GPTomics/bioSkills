---
name: bio-clinical-biostatistics-adaptive-designs
description: Designs adaptive clinical trials including group-sequential (O'Brien-Fleming, Pocock, Lan-DeMets spending), sample-size re-estimation (blinded Friede-Kieser, unblinded Cui-Hung-Wang, Mehta-Pocock promising zone), seamless Phase 2/3 with treatment-arm selection, population enrichment, and response-adaptive randomisation. Covers FDA 2019 Final Adaptive Designs Guidance, FDA 2022 Master Protocols, and ICH E20 Step 2b/3 draft (June 2025, NOT final). Use when planning interim analyses, sample-size re-estimation, or master/platform-trial designs.
tool_type: r
primary_tool: rpact
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: R `rpact` 4.2+ (Wassmer/Brannath), `gsDesign` 3.6+ and `gsDesign2` 1.1+ (Anderson/Merck), `adaptr`, `simtrial`. Commercial: East/EastHorizon (Cytel), ADDPLAN (ICON), FACTS (Berry Consultants).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name`
- Python adaptive packages are limited; R is the regulatory de facto standard

If code throws an error, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Adaptive Clinical Trial Designs

**"Design an adaptive trial"** -> Pre-specify a design with one or more interim adaptations (early stopping, sample-size re-estimation, treatment selection, population enrichment, randomisation ratio changes) that strongly controls Type-I error at the trial-wide level via combination tests or the Conditional Rejection Probability principle.

## Regulatory Status -- The 2024-2026 Landscape

**FDA 2019 Final Adaptive Designs Guidance** (Federal Register 2019-25986, Dec 2 2019) finalised the 2010 and 2018 drafts. Recognises 5 design types: group-sequential, blinded SSR, unblinded SSR, adaptive enrichment, adaptive randomisation.

**FDA 2022 Final Master Protocols Guidance** (March 2022, NOT 2018 — common citation error): basket (one drug, many diseases), umbrella (multiple drugs, one disease), platform (perpetual, drugs enter/exit).

**ICH E20 Adaptive Clinical Trials**: **Step 2b draft June 25 2025; Step 3 public consultation (EU deadline Nov 30 2025; FDA Federal Register Sept 30 2025); Step 4 final targeted late 2026.** As of May 2026, ICH E20 is NOT final. The EFPIA/PhRMA position paper preceded the formal ICH work; Berry Consultants public comment letter is one of the more important submissions.

**FDA CDER Bayesian Methodology Draft (Jan 2026)** (FDA-2025-D-3217): first-ever drug-side Bayesian guidance; permits Bayesian primary inference in pivotals with simulation-based Type-I error calibration.

**Project Optimus (FDA OCE, 2021-2024)**: rewrites Phase I/II oncology by requiring randomised dose comparison before registration, replacing MTD-and-go. Made BOIN, mTPI-2, and multi-arm dose-finding the default.

## Algorithmic Taxonomy

| Design type | Adaptation | Type-I preservation | Software | Strength | Fails when |
|-------------|-----------|---------------------|----------|----------|------------|
| Group-sequential (O'Brien-Fleming) | Early stopping for efficacy/futility | Boundary calculation; very conservative early, near-nominal at end | rpact, gsDesign | FDA's preferred adaptive design | More complex SAP; IDMC firewall essential |
| Group-sequential (Pocock) | Early stopping | Constant nominal alpha at each look | rpact, gsDesign | Easy early stopping | Large penalty at final analysis |
| Wang-Tsiatis power family | Early stopping | Parameterised by Delta | rpact | Tunable conservatism | Δ choice matters |
| Lan-DeMets spending function | Early stopping (flexible timing) | Alpha-spending function | rpact, gsDesign | Operational flexibility; analyses don't need pre-specified number | FDA's de facto preferred framework |
| Blinded SSR (Friede-Kieser 2006) | Re-estimate variance/event-rate; recompute n | No Type-I inflation; agency-uncontroversial | rpact | EMA/FDA endorsed | Variance estimate must be blinded |
| Unblinded SSR (Cui-Hung-Wang 1999) | Increase n based on interim effect estimate | Requires CHW weights for control; or Mehta-Pocock promising zone | rpact | Recovers power if interim promising | IDMC firewall must be perfect; Jennison-Turnbull 2015 critique |
| Mehta-Pocock promising zone (2011) | Increase n if conditional power in (0.3, 0.8) | Calibrated so Type-I inflation negligible (~0.001) | rpact | Operational simplicity | "Stealth alpha inflation" critique (Jennison 2015) |
| Bauer-Köhne 1994 combination | Combine stagewise p-values via Fisher product | Any pre-specified design modification | rpact | Most flexible; theoretical foundation | Power loss vs designed group-sequential |
| Müller-Schäfer 2001 CRP principle | Preserve null conditional rejection probability | Any adaptation at any time | rpact | Modern theoretical bedrock | Implementation complexity |
| Adaptive enrichment | Drop sub-populations failing futility | Closed-test stage-wise | rpact, adaptr | Recovers power on responders | Selection bias on enriched population |
| Response-adaptive randomisation | Update allocation probabilities | Stratification + time-trend covariates required | adaptr, FACTS | Patient-welfare; learn-and-confirm | Drift bias, estimator bias; controversial (Hey-Kimmelman 2015 ethics) |
| Bayesian platform (I-SPY 2 style) | RAR + biomarker stratification + graduation criterion | Frequentist OCs via simulation | FACTS, custom Stan/JAGS | Modern oncology adaptive | Operational complexity; requires IDMC sophistication |

**Postdoc reading list:**

- Bauer P, Köhne K 1994 *Biometrics* 50:1029 (combination test; original adaptive)
- Cui L, Hung HMJ, Wang SJ 1999 *Biometrics* 55:853 (CHW weighted test for unblinded SSR)
- Müller HH, Schäfer H 2001 *Biometrics* 57:886 (CRP principle — theoretical bedrock)
- Mehta CR, Pocock SJ 2011 *Stat Med* 30:3267 (promising zone)
- Jennison C, Turnbull BW 2015 *Stat Med* 34(29):3793-3810 (Mehta-Pocock critique)
- Friede T, Kieser M 2006 *Biom J* 48:537 (blinded SSR)
- Lan KKG, DeMets DL 1983 *Biometrika* 70:659 (alpha spending function)
- O'Brien PC, Fleming TR 1979 *Biometrics* 35:549 (OBF boundary)
- Pocock SJ 1977 *Biometrika* 64:191 (Pocock boundary)
- Hey SP, Kimmelman J 2015 *Clin Trials* 12:102 (RAR ethics critique)
- Berry DA 2015 commentary *Clin Trials* 12:107 (counter)
- Robertson DS, Lee KM, López-Kolkovska BC, Villar SS 2023 *Stat Sci* 38:185 (canonical modern RAR review)
- Wassmer G, Brannath W 2016 *Group Sequential and Confirmatory Adaptive Designs in Clinical Trials* (Springer)
- Jennison C, Turnbull BW 2000 *Group Sequential Methods with Applications to Clinical Trials* (CRC)

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Confirmatory trial wanting interim early stopping | Group-sequential with O'Brien-Fleming boundaries via gsDesign | FDA-preferred; near-nominal final alpha |
| Group-sequential with flexible look timing | Lan-DeMets spending function | Operational flexibility; FDA de facto preferred |
| Phase 3 with uncertain nuisance parameter (variance, event rate) | Blinded SSR (Friede-Kieser) | No Type-I inflation; agency-uncontroversial |
| Phase 3 wanting to increase n if interim shows promise | Mehta-Pocock promising zone with CHW weights | Recovers power; calibrated Type-I |
| Seamless Phase 2/3 with arm selection | Bauer-Köhne combination test + closed testing | Most flexible; cite Müller-Schäfer CRP |
| Adaptive enrichment (drop subpopulation) | Adaptive enrichment with closed-test stage-wise | Recovers power on responders |
| Multi-arm oncology platform | Bayesian platform with RAR (I-SPY 2 model) | Patient-welfare argument strong for multi-arm |
| 2-arm phase 3 oncology with potential RAR | Avoid RAR; group-sequential preferred | Hey-Kimmelman 2015 ethics + drift bias |
| Continuous endpoint, treatment discontinuation, follow-up data available | Hybrid: J2R imputation for treatment-discontinuation ICEs, MMRM-MAR for other missingness | Aprocitentan PRECISION precedent (2024); FDA de facto standard 2024-2025 for treatment-policy estimands |
| Phase 1 dose-finding | BOIN (FDA Fit-for-Purpose qualified 2021) | Transparent, tabulated decisions; no bedside Bayesian software |
| Phase 1b/2 dose-optimisation (Project Optimus) | Multi-arm BOIN-12 or multi-dose randomised | FDA Aug 2024 final dose-optimisation guidance |
| Basket trial (one drug, multiple diseases) | EXNEX or robust MAP via RBesT | Borrows across baskets while permitting one to detach |
| Umbrella trial (one disease, multiple drugs) | Bayesian platform with shared control | FDA Master Protocols 2022 |
| Pediatric extrapolation borrowing from adults | Power prior with discount γ in 0.3-0.6 | FDA Bayesian Jan 2026 draft endorses |

## Group-Sequential Designs

### O'Brien-Fleming -- the regulatory default

```r
library(gsDesign)

# OBF boundaries; 3 interim looks at 33%, 67%, 100% information
design <- gsDesign(
    k = 4,                  # total analyses including final
    test.type = 1,          # 1-sided efficacy
    alpha = 0.025,
    beta = 0.10,            # power = 0.90
    sfu = sfLDOF,           # Lan-DeMets approximation of OBF
    timing = c(0.25, 0.50, 0.75, 1.0)
)
print(design)
plot(design)
```

OBF is **very conservative at early looks** (nominal alpha approximately 0.0001 at 25% info) and **near-nominal at final analysis** (~0.024 of 0.025). Preferred by FDA because the final-analysis penalty is small.

### Pocock -- constant nominal

Constant nominal alpha at each look. Easy early stopping but large final-analysis penalty (~0.018 of 0.025 with k=4). Rarely used in confirmatory.

### Lan-DeMets spending function -- the modern flexibility

```r
# Lan-DeMets OBF-like spending function (sfLDOF)
# Allows analysis timing to differ from pre-specified
design_flex <- gsDesign(
    k = 3,
    sfu = sfLDOF,      # OBF-like spending
    alpha = 0.025,
    beta = 0.10
)

# Actual analyses can occur at different information fractions
# Spending function returns alpha to spend at each look based on actual timing
```

**The flexibility:** sponsor can perform analyses at different information fractions than originally planned. FDA's de facto preferred framework.

### Sample-size for group-sequential

```r
# Time-to-event group-sequential
library(gsDesign)
n_gs <- gsSurv(
    k = 3,
    test.type = 2,    # 2-sided
    alpha = 0.025,
    beta = 0.10,
    sfu = sfLDOF,
    lambdaC = 0.04,   # control hazard per month
    hr = 0.70,        # treatment HR
    eta = 0.005,      # dropout hazard
    T = 24,           # total study duration
    minfup = 12       # minimum follow-up
)
print(n_gs)
```

## Sample-Size Re-Estimation

### Blinded SSR (Friede-Kieser 2006)

Re-estimate nuisance parameter (variance σ² for continuous, control event rate p_0 for binary, overall event rate for survival) from blinded interim data. **No Type-I error inflation when test statistic ignores the SSR.**

```r
library(rpact)
# Blinded SSR for continuous outcome
design_blinded_ssr <- getDesignGroupSequential(
    kMax = 2,
    alpha = 0.025,
    beta = 0.20,
    sided = 1,
    informationRates = c(0.5, 1)
)

# At interim, re-estimate variance and recompute n
# (manual implementation; rpact has built-in support via getDesignInverseNormal for unblinded)
```

EMA Reflection Paper 2007 and FDA 2019 explicitly endorse blinded SSR. **Uncontroversial.**

### Unblinded SSR (Cui-Hung-Wang 1999)

Interim effect estimate triggers sample-size change. **Type-I inflation if naive:** Cui-Hung-Wang showed 8% Type-I vs 2.5% target.

The Cui-Hung-Wang weighted test uses pre-specified weights from the original design:

```
Z_weighted = w_1 * Z_1 + w_2 * Z_2_residual
```

where w_1, w_2 are the pre-specified weights (based on original n_1, n_2) and Z_2_residual is the test statistic on the data after the interim. **Pre-specified weights preserve alpha** even if the actual n at stage 2 differs.

```r
library(rpact)
design_unblinded_ssr <- getDesignInverseNormal(
    kMax = 2,
    alpha = 0.025,
    beta = 0.20,
    sided = 1,
    informationRates = c(0.5, 1),
    typeOfDesign = 'WT',  # Wang-Tsiatis power family
    deltaWT = 0.25
)

# Use inverse normal combination for adaptive SSR
analysis_result <- getAnalysisResults(
    design_unblinded_ssr,
    dataInput = getDataMeans(...)
)
```

### Mehta-Pocock Promising Zone (2011)

At interim, compute conditional power (CP) given observed effect:

- **Unfavourable zone** (CP < ~30%): stop or continue without modification
- **Promising zone** (CP in 30-80%): increase n to recover power; NO Type-I penalty if increase rule pre-specified and uses original test statistic with original weights
- **Favourable zone** (CP > 80%): continue without change

```r
# rpact implementation
# Sample size recalculation in promising zone
n_increased <- getSampleSizeMeans(
    design_unblinded_ssr,
    alternative = 5,        # detect mean diff of 5
    stDev = 12,
    groups = 2
)
```

**The mathematical sleight:** promising zone is constructed so unconditional Type-I error inflation is negligible (~0.001) even WITHOUT CHW weighting. **Jennison-Turnbull 2015 critique:** stealth alpha inflation in unpublished simulation assumptions; inefficient relative to CHW-weighted GSD. Mehta defends on operational grounds.

**Hsiao et al 2020 *Trials* 21:1003** is the systematic review.

## Combination Tests and CRP Principle

**Bauer-Köhne 1994** *Biometrics* 50:1029: combine stagewise p-values via Fisher's product test. Permits design modifications post-interim while controlling Type-I error.

**Müller-Schäfer 2001** *Biometrics* 57:886: **Conditional Rejection Probability (CRP) principle** — preserve the null conditional rejection probability at every adaptation, and unconditional Type-I is preserved. The theoretical bedrock of all post-2001 confirmatory adaptive designs.

Müller-Schäfer 2004 *Stat Med* 23:2497 extended to ANY design change at ANY time.

```r
# rpact natively supports combination tests
design_comb <- getDesignFisher(
    kMax = 3,
    alpha = 0.025,
    sided = 1
)

# Or inverse normal combination
design_inv_norm <- getDesignInverseNormal(
    kMax = 3,
    alpha = 0.025,
    informationRates = c(0.33, 0.67, 1.0)
)
```

## Adaptive Enrichment

Drop sub-populations failing futility; re-power on responders. **Closed-test stage-wise** to control familywise error across full and enriched populations.

```r
# rpact: enrichment design via getDesignEnrichmentSubgroup
# Standard implementation requires explicit definition of full population (F)
# and enriched population (S)
```

**Postdoc concern:** selection bias on the enriched population — the observed treatment effect on the enriched subgroup is biased upward by selection. Bias-correction via simulation or hierarchical Bayesian.

## Response-Adaptive Randomisation -- The Ethics Fight

**Hey & Kimmelman 2015 *Clin Trials* 12:102 "Are outcome-adaptive allocation trials ethical?"** Argued RAR's purported ethical advantage (equipoise, sub-optimal exposure minimisation) **fails in two-arm and early-phase settings** because:

1. Drift bias inflates Type-I error / biases estimates (time trends confounded with allocation)
2. Consent dynamics confused — patients believe allocation is "personalised" when stochastic
3. Marginal patient-welfare benefit is statistical and small while operational risks real

**Counter-arguments:**

- **Berry DA 2015 commentary** *Clin Trials* 12:107: RAR enables learn-and-confirm, multi-arm platforms (I-SPY 2 model) where the patient-welfare argument IS the point and equal allocation would be unethical given accumulating evidence.
- **Saville & Berry 2016 *Clin Trials* 13:358:** RAR's operating characteristics are competitive in multi-arm; bias and inflation problems disappear with stratification, time-trend covariates, and proper analysis weights.
- **Buyse 2015 *Clin Trials* 12:108:** Hey-Kimmelman correct for 2-arm but wrong for multi-arm.

**Consensus position (2020s; ICH E20):** RAR appropriate when (a) multi-arm (>=3 arms), (b) rare disease / limited pool, (c) strong PoC of differential biomarker response, (d) robust drift-bias adjustment and pre-specified analysis weights. **Inappropriate for confirmatory 2-arm trials.**

**Robertson, Lee, López-Kolkovska, Villar 2023 *Stat Sci* 38:185 ("Response-adaptive randomization: from myths to practical considerations")** is the canonical modern review settling the debate.

## Bayesian Platform Trials

**I-SPY 2 (Barker-Sigman 2009 *Clin Pharmacol Ther*; Park-Liu 2016 *NEJM* 375:11):** neoadjuvant breast cancer; 10 biomarker-defined subtypes × multiple arms; Bayesian RAR; graduation criterion = posterior predictive probability of success in 300-patient Phase 3 ≥ 85%. Berry Consultants designed the engine. Multiple drugs graduated (neratinib, veliparib, pembrolizumab).

**GBM AGILE (Alexander 2018; published readouts beginning 2024):** glioblastoma; response-adaptive Bayesian; first global registrational platform in neuro-oncology. Regorafenib readout 2025 *JCO* JCO-25-01137.

**REMAP-CAP (Angus 2020 *JAMA*):** severe pneumonia, repurposed for COVID-19 in 2020; **Bayesian factorial multi-domain design** — multiple intervention domains tested simultaneously and combinatorially. Generated corticosteroid signal in COVID independently of RECOVERY.

### Drop-the-loser vs promising-the-winner

- **Adaptive arm-dropping (futility):** Bayesian posterior probability of beating control drops below threshold → arm closes. Mathematically straightforward; FDA-acceptable.
- **"Promising-the-winner" (graduate to Phase 3):** introduces selection bias. Bias-adjusted estimators (Robertson 2023; conditional MLE) now standard in I-SPY 2 reports.

## Phase I Dose-Finding -- BOIN, mTPI, CRM

| Design | Citation | Idea | Where it wins |
|--------|----------|------|---------------|
| CRM | O'Quigley-Pepe-Fisher 1990 | Single-parameter logistic/power model; updates posterior MTD probability after each cohort | Statistically efficient; skeleton calibration needed |
| EWOC | Babb-Rogatko-Zacks 1998 | CRM-like with explicit overdose-control constraint (P(dose > MTD) <= 0.25) | Safer than CRM in small trials |
| mTPI | Ji et al 2010 *Clin Trials* 7:653 | Beta-binomial; UPM decision rule on under/proper/over-dosing intervals | Pre-tabulated decisions; documented over-shoot bias |
| mTPI-2 / Keyboard | Guo-Wang-Chen-Ji 2017 | Fixes mTPI Ockham bias by equal-width intervals | Default mTPI replacement |
| BOIN | Liu-Yuan 2015 *J R Stat Soc C* 64:507 | Pre-tabulated escalation interval bounds optimised to minimise incorrect-decision probability | **FDA Fit-for-Purpose qualified Dec 2021**; near-CRM with no bedside software |

**Why FDA prefers BOIN operationally:** qualified as Fit-for-Purpose under FDA Drug Development Tools program (review document FDA-2020-X-XXXX, posted 2021). Investigator uses pre-printed escalation table — no real-time Bayesian software at the bedside.

R packages: `BOIN`, `dfcrm` (Cheung — author of CRM textbook), `trialr` (Brock — includes EffTox), `escalation` (Brock — unified framework).

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Blinded SSR n vs unblinded SSR n differ substantially | Unblinded SSR uses interim effect estimate; blinded uses nuisance parameter only | Blinded is Type-I-clean; unblinded requires CHW weighting; pre-specify the approach in SAP |
| Group-sequential rejects at interim; Cui-Hung-Wang weighted final test does not | Naive interim rejection used original test statistic; CHW weights downweight late data | Pre-specify boundary and weights; do NOT switch tests mid-stream |
| Mehta-Pocock promising-zone vs CHW-weighted GSD give different n increases | Promising zone calibrated for Type-I (~0.001 inflation); CHW more efficient under known effect | Jennison-Turnbull 2015 critique: promising zone "stealth alpha"; pre-specify with simulation OCs |
| Adaptive enrichment selects subgroup at interim; replication shows smaller effect | Selection bias on enriched population (Sun 2010 winner's curse) | Bias-correction via conditional MLE or hierarchical Bayesian; cite Robertson 2023 |
| RAR posterior allocation favours active in 2-arm trial; randomisation drift bias suspected | Time trends confounded with allocation changes | Pre-specify time-trend covariates in analysis; use proper analysis weights; cite Robertson 2023 RAR consensus (RAR INAPPROPRIATE for 2-arm confirmatory) |
| BOIN vs CRM choose different MTD on same data | CRM uses model; BOIN uses tabulated boundaries; differ when skeleton mis-calibrated | BOIN Fit-for-Purpose qualified (Dec 2021); CRM more efficient under correct skeleton; report OCs over both |
| I-SPY 2 graduation criterion met but Phase 3 replication fails | Selection bias on graduated arm; PP threshold not bias-corrected | Apply conditional MLE; cite Robertson 2023; report both raw and bias-corrected estimates |
| Müller-Schäfer CRP preserved but ad hoc rule appears Type-I-inflated in simulation | Implementation deviation from formal CRP | Verify CRP equation precisely; report OCs via simulation; cite Müller-Schäfer 2001 |

## Per-Method Failure Modes

### Unblinded SSR with naive sample increase

- **Trigger:** Sponsor increases n based on interim effect without CHW weighting.
- **Mechanism:** Type-I inflation up to 8% (Cui-Hung-Wang 1999 simulation).
- **Symptom:** Independent reanalysis finds Type-I > 5%.
- **Fix:** Pre-specify CHW weights from original design; use combination test in rpact.

### Mehta-Pocock promising-zone "stealth alpha"

- **Trigger:** Promising-zone applied without sufficient simulation.
- **Mechanism:** Jennison-Turnbull 2015 critique — Type-I inflation hidden in unpublished simulation assumptions.
- **Symptom:** Independent reanalysis finds Type-I 5.3% vs nominal 5%.
- **Fix:** Pre-specify increase rule transparently; report simulation OCs.

### RAR drift bias

- **Trigger:** RAR in trial with time trends (calendar effects, learning curves).
- **Mechanism:** Time trends confounded with allocation changes; biased effect estimate.
- **Symptom:** Effect estimate sensitive to time-trend adjustment.
- **Fix:** Pre-specify time-trend covariates in analysis; use proper analysis weights; cite Robertson 2023.

### Schoenfeld formula under immunotherapy delayed effect

- **Trigger:** Sample size calculated via Schoenfeld 1981 assuming PH.
- **Mechanism:** Delayed effect violates PH; events under-estimated by 20-50%.
- **Symptom:** Trial under-powered; observed events insufficient.
- **Fix:** Lakatos 1988 or simulation under expected HR(t); cite Lin 2020 NPH Working Group.

### IDMC firewall failure in unblinded SSR (IDMC = Independent Data Monitoring Committee; the regulatory-standard term)

- **Trigger:** Interim effect estimate leaks beyond IDMC.
- **Mechanism:** Sponsor inference from increase decision reveals direction of interim effect.
- **Symptom:** Regulator audit reveals unblinding.
- **Fix:** Strict firewall SOP; only "increase / no increase" communicated to sponsor; cite ICH E20.

### Adaptive enrichment selection bias

- **Trigger:** Enriched population effect reported without bias correction.
- **Mechanism:** Selection on subgroup with promising interim effect inflates estimate.
- **Symptom:** Independent replication on enriched subgroup gives smaller effect.
- **Fix:** Bias-correction via simulation or hierarchical Bayesian; cite Sun 2010 winner's curse.

### RAR in 2-arm confirmatory

- **Trigger:** RAR applied to confirmatory 2-arm trial.
- **Mechanism:** Hey-Kimmelman 2015 critique — drift bias, consent confusion, marginal benefit.
- **Symptom:** Reviewer rejects RAR as inappropriate for setting.
- **Fix:** Group-sequential with futility/efficacy boundaries instead; cite Robertson 2023 consensus.

### CRM skeleton mis-calibration

- **Trigger:** CRM applied with default skeleton without simulation.
- **Mechanism:** Skeleton dictates target dose; mis-calibration biases MTD.
- **Symptom:** MTD selection differs systematically from clinical expectation.
- **Fix:** Calibrate skeleton via Lee-Cheung 2009 indifference-interval method; or switch to BOIN.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| FDA Fit-for-Purpose BOIN qualification (Dec 2021) | FDA Drug Development Tools program | First dose-finding design with formal FDA endorsement |
| Mehta-Pocock promising zone CP 30-80% | Mehta-Pocock 2011 | Mathematical calibration for Type-I preservation |
| RAR appropriate >= 3 arms | Robertson 2023 consensus | Multi-arm patient-welfare argument |
| OBF nominal alpha ~0.024 at final / 0.025 | gsDesign | Small final penalty preferred by FDA |
| Schoenfeld under non-PH under-estimates 20-50% | Lin 2020 NPH WG | Use Lakatos or simulation |
| I-SPY 2 graduation: PP success in Phase 3 >= 85% | Barker 2009 | Bayesian platform standard |
| Power prior discount γ 0.3-0.6 for pediatric extrapolation | FDA Bayesian Jan 2026 draft | Partial borrowing default |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Unblinded SSR with naive sample increase | No CHW weighting | Pre-specify CHW weights; cite Cui-Hung-Wang 1999 |
| Mehta-Pocock without simulation OCs | Stealth alpha inflation | Report simulation OCs; cite Jennison 2015 |
| RAR in 2-arm confirmatory | Misapplication | Group-sequential instead; cite Robertson 2023 |
| Schoenfeld for immuno-oncology | PH assumption violated | Lakatos or simulation; cite Lin 2020 |
| Adaptive enrichment effect reported uncorrected | Selection bias | Bias-correction; cite Sun 2010 |
| CRM with default skeleton | Mis-calibration | Calibrate via Lee-Cheung 2009 or switch to BOIN |
| ICH E20 cited as "finalised April 2024" | Confusion with EFPIA position paper | ICH E20 is Step 2b/3 draft (June 2025); not final |
| FDA Master Protocols "2018" | 2018 was draft | March 2022 was the final |
| Bauer-Köhne combination test as "old-fashioned" | Misunderstanding | Foundational; cited in modern combination-test implementations |
| Stop-for-efficacy at first interim with OBF | OBF nominal alpha ~0.0001 at 25% info | Trial must show very strong evidence to stop early; expected |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "How is Type-I error controlled?" | Closed testing via Müller-Schäfer CRP principle; specific implementation is inverse normal combination test in rpact |
| "Why these boundaries?" | OBF via Lan-DeMets sfLDOF spending function; preserves final-analysis power; pre-specified in SAP |
| "Pre-specification of SSR rule?" | Promising zone CP in (0.3, 0.8) triggers increase to n_max via CHW-weighted statistic; pre-specified n_max in SAP |
| "IDMC firewall?" | IDMC receives interim effect estimate; sponsor receives only "increase / no increase" decision; SOP documented; pre-specified |
| "RAR ethics?" | Multi-arm (4 arms) setting; Berry 2015 consensus that patient-welfare argument valid; drift-bias adjustment in primary analysis |
| "Promising zone vs CHW-weighted GSD?" | Operational simplicity preferred; OCs from simulation confirm Type-I ~5%; supportive Cui-Hung-Wang analysis |
| "Adaptive enrichment bias?" | Bias-correction via simulation; conditional MLE for enriched-population effect; cite Sun 2010 |
| "Phase 1 BOIN vs CRM?" | BOIN Fit-for-Purpose qualified by FDA Dec 2021; tabulated decisions; no bedside Bayesian software |

## References

- Babb J, Rogatko A, Zacks S. 1998. Cancer Phase I clinical trials: efficient dose escalation with overdose control. *Stat Med* 17:1103-1120.
- Bauer P, Köhne K. 1994. Evaluation of experiments with adaptive interim analyses. *Biometrics* 50:1029-1041.
- Berry DA. 2015. Commentary on Hey & Kimmelman. *Clin Trials* 12:107-110.
- Cui L, Hung HMJ, Wang SJ. 1999. Modification of sample size in group sequential clinical trials. *Biometrics* 55:853-857.
- FDA. 2019. Adaptive Designs for Clinical Trials of Drugs and Biologics. Final Guidance.
- FDA. 2022. Master Protocols: Efficient Clinical Trial Design Strategies to Expedite Development of Oncology Drugs and Biologics. Final Guidance, March 2022.
- FDA. 2026. Use of Bayesian Methodology in Clinical Trials. Draft Guidance, January 2026.
- Friede T, Kieser M. 2006. Sample size recalculation in internal pilot study designs. *Biom J* 48:537-555.
- Hey SP, Kimmelman J. 2015. Are outcome-adaptive allocation trials ethical? *Clin Trials* 12:102-106.
- Jennison C, Turnbull BW. 2015. Adaptive sample size modification in clinical trials: start small then ask for more? *Stat Med* 34(29):3793-3810.
- Lan KKG, DeMets DL. 1983. Discrete sequential boundaries for clinical trials. *Biometrika* 70:659-663.
- Liu S, Yuan Y. 2015. Bayesian optimal interval designs for phase I clinical trials. *JRSS-C* 64:507-523.
- Mehta CR, Pocock SJ. 2011. Adaptive increase in sample size when interim results are promising. *Stat Med* 30:3267-3284.
- Müller HH, Schäfer H. 2001. Adaptive group sequential designs for clinical trials: combining the advantages of adaptive and of classical group sequential approaches. *Biometrics* 57:886-891.
- O'Brien PC, Fleming TR. 1979. A multiple testing procedure for clinical trials. *Biometrics* 35:549-556.
- O'Quigley J, Pepe M, Fisher L. 1990. Continual reassessment method: a practical design for phase 1 clinical trials in cancer. *Biometrics* 46:33-48.
- Pocock SJ. 1977. Group sequential methods in the design and analysis of clinical trials. *Biometrika* 64:191-199.
- Robertson DS, Lee KM, López-Kolkovska BC, Villar SS. 2023. Response-adaptive randomization in clinical trials: from myths to practical considerations. *Stat Sci* 38:185-208.
- Wassmer G, Brannath W. 2016. *Group Sequential and Confirmatory Adaptive Designs in Clinical Trials*. Springer.

## Related Skills

- clinical-biostatistics/power-and-sample-size - Sample size for adaptive designs
- clinical-biostatistics/multiplicity-graphical - Closed testing in adaptive contexts
- clinical-biostatistics/bayesian-trials - Bayesian platform trials, BOIN/CRM/EWOC
- clinical-biostatistics/trial-reporting - Reporting adaptive trial results per CONSORT 2025
- clinical-biostatistics/survival-analysis - Adaptive designs for TTE endpoints
- experimental-design/sample-size - General sample-size methods
