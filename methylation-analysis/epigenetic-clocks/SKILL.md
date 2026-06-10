---
name: bio-methylation-epigenetic-clocks
description: Computes DNA methylation age (DNAm age) and pace of aging by applying frozen elastic-net epigenetic clocks to a clean beta matrix with methylclock, dnaMethyAge, or methylCIPHER. Covers the clock menu by question (chronological Horvath/Hannum/skin&blood; health-mortality PhenoAge/GrimAge; DunedinPACE pace; pediatric/gestational; mitotic epiTOC), age acceleration (EAA/IEAA/EEAA) as the real endpoint, the principal-component (PC) clock fix for the per-CpG reliability crisis, and EPICv2 clock-CpG dropout with missing-CpG imputation bias. Use when estimating epigenetic age, computing age acceleration, choosing a clock for an outcome, assessing clock reliability, or porting a clock to EPICv2. A clock is a frozen predictor: do not GO-enrich its CpGs and do not train it here. For cell-count adjustment (IEAA) see cell-type-deconvolution; for predictor training/validation/leakage see machine-learning/model-validation; for survival modeling of age acceleration see clinical-biostatistics/survival-analysis.
tool_type: r
primary_tool: methylclock
---

## Version Compatibility

Reference examples tested with: methylclock 1.8+, dnaMethyAge (GitHub yiluyucheng), methylCIPHER (GitHub MorganLevineLab), DunedinPACE (GitHub danbelsky).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The clock coefficient sets are FIXED and version-pinned (a clock is a frozen list of CpGs and weights), so the package version mostly controls which clocks ship and what the clock-name strings are. Verify accepted names live: `methylclock` via `checkClocks(beta)`; `dnaMethyAge` via `availableClock()`. The ARRAY PLATFORM is the version that matters most: EPICv2 drops a clock-specific fraction of CpGs, so always report how many of each clock's CpGs were actually present.

# Epigenetic Clocks

**"How old is this sample epigenetically?"** -> Apply a frozen elastic-net clock to the beta matrix, then report the age ACCELERATION (residual vs chronological age), not the raw age - because a clock is a predictor, and the residual is the signal.
- R: `DNAmAge(beta, clocks = c('Horvath', 'Hannum', 'Levine'), age = pheno$age)`

Scope: applying pre-trained clocks (DNAm age, pace, mitotic) and computing age acceleration from a clean beta/M-value matrix. Cell-count adjustment for IEAA -> cell-type-deconvolution. Clean beta matrix / EPICv2 replicate-probe collapse -> array preprocessing. Training a predictor / cross-validation / leakage -> machine-learning/model-validation. Survival/mortality modeling of EAA -> clinical-biostatistics/survival-analysis. Per-CpG and region testing -> differential-cpg-testing, dmr-detection.

## The Single Most Important Modern Insight -- A Clock Is a Predictor, Not a Mechanism, and Not Even a Reliable One Per-CpG

DNAm age is a frozen elastic-net weighted sum over CpGs that a penalty chose for out-of-sample prediction. The CpGs are prediction features, never an aging pathway. Four corollaries every common misuse violates:

1. **The CpG set is not biology.** Do not GO-enrich clock CpGs. Two clocks for the same outcome can share almost zero CpGs (the elastic net arbitrarily keeps one of many correlated predictors), so non-overlap is expected, not a contradiction.
2. **The endpoint is age ACCELERATION, not the raw age.** Raw DNAm age just recapitulates chronological age (r often > 0.9). The signal is the residual of DNAm age on chronological age (EAA); IEAA additionally residualizes on cell counts - which is exactly where deconvolution meets this skill.
3. **First-gen per-CpG reliability can be smaller than the effect being chased.** Many first-gen clock CpGs have low test-retest ICC (Sugden 2020 *Patterns* 1:100014), so the same sample can age several years between technical replicates. PC clocks (Higgins-Chen 2022 *Nat Aging* 2:644) exist specifically to fix this for longitudinal and trial use.
4. **Association is not causation is not transfer.** EAA associating with an exposure does not make the clock causal; a blood clock does not automatically work in another tissue or ancestry without a recalibration check.

Organize the analysis around defending these four, not around listing clock names.

## The Clock Menu by Question

The question dictates the clock; there is no single best clock. Pick by what is being predicted, not by popularity.

| Generation | Clock | Citation | Predicts | Tissue | Note |
|------------|-------|----------|----------|--------|------|
| 1st (chronological) | Horvath multi-tissue | Horvath 2013 *Genome Biol* 14:R115 | chronological age | 51 tissues | 353 CpGs; works cross-tissue; log-linear age transform for <20y |
| 1st | Hannum | Hannum 2013 *Mol Cell* 49:359 | chronological age | whole blood | 71 CpGs; tight in blood, poor cross-tissue |
| 1st | skin & blood | Horvath 2018 *Aging* 10:1758 | chronological age | skin, blood, fibroblasts | 391 CpGs; for in-vitro/fibroblast/skin work |
| 2nd (health-mortality) | PhenoAge | Levine 2018 *Aging* 10:573 | morbidity/mortality composite | blood | 513 CpGs; trained on a 9-biomarker phenotypic age |
| 2nd | GrimAge | Lu 2019 *Aging* 11:303 | lifespan/healthspan | blood | composite of DNAm protein surrogates; strongest mortality predictor |
| pace | DunedinPACE | Belsky 2022 *eLife* 11:e73420 | RATE of aging | blood | 173 CpGs; ~1.0 = one biological year per calendar year; NOT an age |
| pediatric | PedBE | McEwen 2020 *PNAS* 117:23329 | age 0-20 | buccal | buccal-specific |
| gestational | Knight / Bohlin | Knight 2016 *Genome Biol* 17:206 / Bohlin 2016 *Genome Biol* 17:207 | gestational age | cord blood | newborn GA estimation |
| mitotic | epiTOC | Yang 2016 *Genome Biol* 17:205 | cumulative stem-cell divisions | normal tissue | tracks mitotic, not chronological, age; cancer-risk relevant |

DunedinPACE is reported as the raw PACE value (already a rate); never residualize it like an age clock and never compare its number to Horvath years.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Chronological-age accuracy in blood | Hannum or Horvath | trained on chronological age; Hannum tighter in blood |
| Cross-tissue or non-blood sample | Horvath multi-tissue or skin&blood | only the multi-tissue clocks transfer; report a known-age check |
| Morbidity / mortality / healthspan endpoint | GrimAge (or PhenoAge) | second-gen; trained on health outcomes, not just age |
| Pace of aging / intervention sensitivity | DunedinPACE | a rate; sensitive to caloric-restriction-style trials |
| Longitudinal or clinical-trial endpoint | PC clocks (methylCIPHER) | first-gen test-retest noise can exceed the intervention effect |
| Pediatric buccal / newborn cord blood | PedBE / Knight or Bohlin | age-and-tissue-matched clocks |
| Cancer / mitotic-age question | epiTOC / epiTOC2 | estimates cell divisions, a different aging axis |
| EPICv2 data | clock that retains its CpGs (Horvath/PhenoAge) | GrimAge/Hannum/DunedinPACE lose >10% of CpGs on EPICv2 |
| Adjust EAA for cell composition (IEAA) | -> cell-type-deconvolution | residualize the clock on estimated cell counts |
| Train or validate a new predictor | -> machine-learning/model-validation | clocks here APPLY frozen models; they are not trained here |
| Survival/mortality model of EAA | -> clinical-biostatistics/survival-analysis | EAA-to-outcome modeling lives there |

## Apply a Clock and Extract Age Acceleration

**Goal:** Compute DNAm age for several clocks and turn the raw ages into age acceleration, the actual endpoint.

**Approach:** Run `DNAmAge` with chronological age supplied so it returns acceleration columns directly; `ageAcc` is the raw DNAm-minus-chronological difference and `ageAcc2` is the residual of DNAm age on chronological age (the EAA to test). Always check clock-CpG coverage first.

```r
library(methylclock)

cpg_report <- checkClocks(beta)        # which clock CpGs are missing BEFORE estimating
ages <- DNAmAge(beta, clocks = c('Horvath', 'Hannum', 'Levine', 'skinHorvath'),
                age = pheno$age,       # supplying age yields ageAcc and ageAcc2 columns
                cell.count = FALSE,    # set TRUE only when adjusting toward IEAA-style estimates
                min.perc = 0.8)        # refuse a clock missing >20% of its CpGs (default 0.8)
# ageAcc  = DNAm age - chronological age (raw difference)
# ageAcc2 = residual of DNAm age on chronological age = the EAA endpoint with cell.count=FALSE
#   (methylclock labels ageAcc2 "similar to IEAA"; confirm the exact column semantics in the
#    installed vignette, and the cell-count-adjusted residual when cell.count=TRUE)
```

The `dnaMethyAge` package returns acceleration in one call and exposes author-year clock IDs:

```r
library(dnaMethyAge)
availableClock()                                  # confirm the installed clock-name strings
phenoage <- methyAge(beta, clock = 'LevineM2018',
                     age_info = pheno,            # data.frame with Sample, Age (Sex for GrimAge variants)
                     fit_method = 'Linear')       # adds an Age_Acceleration column
```

## Pace of Aging Is Not an Age

**Goal:** Compute DunedinPACE as a rate and keep it on its own scale.

**Approach:** Use the dedicated package; the output is a per-sample pace (~1.0 = normal). Do not regress it on chronological age and do not merge it with age-clock acceleration.

```r
library(DunedinPACE)
pace <- PACEProjector(beta)   # returns the DunedinPACE pace values (~1.0 = normal, >1 = faster aging); report as-is
# Never residualize PACE on chronological age and never compare its value to Horvath years.
```

## Reliability and PC Clocks

For longitudinal, interventional, or clinical-trial endpoints, first-gen per-CpG noise (Sugden 2020) can swamp a small intervention effect. PC clocks (Higgins-Chen 2022) train the elastic net on principal components across thousands of CpGs, averaging out per-CpG noise and lifting test-retest ICC toward ~0.9. They live in methylCIPHER (MorganLevineLab), not in methylclock. Use a PC clock, or at minimum document an ICC/reliability assessment, for any repeated-measures design.

## Per-Method Failure Modes

### GO-enriching clock CpGs
**Trigger:** running pathway enrichment on a clock's CpG list to "explain aging." **Mechanism:** the CpGs are penalty-selected prediction features, one arbitrary representative per correlated cluster. **Symptom:** a plausible-looking enrichment that is an artifact of feature selection. **Fix:** do not enrich clock CpGs; treat them as predictors only.

### Reporting raw DNAm age instead of acceleration
**Trigger:** correlating raw DNAm age with an exposure. **Mechanism:** raw age is dominated by chronological age (r > 0.9). **Symptom:** every clock "associates" with age-correlated variables. **Fix:** test the residual (ageAcc2 / Age_Acceleration / IEAA), not the raw age.

### Missing clock CpGs mean-imputed
**Trigger:** a platform or failed probes drop clock CpGs; the tool imputes to the training mean. **Mechanism:** mean-imputation pulls the prediction toward the training population age and shrinks variance. **Symptom:** age acceleration biased toward zero; attenuated associations; on EPICv2 Hannum can return NEGATIVE ages. **Fix:** report the fraction of clock CpGs present (`checkClocks`); flag/refuse samples with high missingness; prefer a clock that retains its CpGs on the platform.

### EPICv2 replicate probes not collapsed
**Trigger:** feeding a raw EPICv2 matrix with suffixed replicate probe IDs. **Mechanism:** EPICv2 carries multiple beads per CpG with suffixed names, so the clock cannot find its CpGs. **Symptom:** huge apparent CpG dropout, nonsensical ages. **Fix:** collapse replicate probes to one value per CpG upstream before any clock call.

### First-gen clock used for a longitudinal endpoint
**Trigger:** detecting a small intervention effect with Horvath/Hannum across timepoints. **Mechanism:** per-CpG test-retest noise (Sugden 2020) rivals the effect. **Symptom:** unstable EAA between replicates; the effect is inside the noise band. **Fix:** PC clocks (methylCIPHER) or a documented reliability assessment.

### Cross-tissue or cross-ancestry application
**Trigger:** a blood-trained clock (Hannum) on saliva, or a European-cohort clock applied elsewhere. **Mechanism:** clocks do not automatically transfer. **Symptom:** a systematic age offset vs known age. **Fix:** use a tissue-appropriate clock (skin&blood, PedBE, gestational) and report a known-age calibration check.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Report fraction of clock CpGs present | Higgins-Chen 2022 *Nat Aging* 2:644 | high imputed fraction invalidates the estimate |
| `min.perc` >= 0.8 of clock CpGs | methylclock docs | below ~80% coverage mean-imputation dominates the prediction |
| EPICv2 dropout 3.5-32.6% per clock; GrimAge/Hannum/DunedinPACE > 10% | EPICv2 clock benchmarks | platform-specific; pick a CpG-retaining clock on EPICv2 |
| First-gen clock CpG ICC often < 0.5 | Sugden 2020 *Patterns* 1:100014 | technical noise rivals signal; use PC clocks for repeated measures |
| PC clock ICC ~0.9+ | Higgins-Chen 2022 *Nat Aging* 2:644 | the reliability bar for longitudinal/trial designs |
| EAA = residual of DNAm age on chronological age | Horvath 2013 *Genome Biol* 14:R115 | the endpoint; raw age is uninformative (r > 0.9 with chronological age) |
| Horvath age transform applied for age < 20 | Horvath 2013 *Genome Biol* 14:R115 | the clock is log-linear below 20y; do not compare pre/post-transform values |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Every clock correlates with an age-linked variable | testing raw DNAm age | test age acceleration (residual), not raw age |
| Hannum returns negative ages | EPICv2 clock-CpG dropout | use a CpG-retaining clock; report coverage; collapse replicate probes |
| Age acceleration shrunk toward zero | clock CpGs mean-imputed | report `checkClocks` coverage; refuse high-missingness samples |
| Many clock CpGs "missing" on EPICv2 | replicate probes not collapsed | collapse suffixed probes to one value per CpG upstream |
| DunedinPACE compared to Horvath years | mixing a rate with an age | keep PACE on its own scale; never residualize it |
| Unstable EAA across replicates | first-gen per-CpG noise | PC clocks (methylCIPHER) or a reliability assessment |
| `methyAge` clock name not found | wrong clock-ID string | `availableClock()` for installed author-year IDs |

## References

- Horvath S. 2013. DNA methylation age of human tissues and cell types. *Genome Biol* 14:R115.
- Hannum G, Guinney J, Zhao L, et al. 2013. Genome-wide methylation profiles reveal quantitative views of human aging rates. *Mol Cell* 49:359-367.
- Horvath S, Oshima J, Martin GM, et al. 2018. Epigenetic clock for skin and blood cells applied to Hutchinson Gilford Progeria Syndrome and ex vivo studies. *Aging (Albany NY)* 10:1758-1775.
- Levine ME, Lu AT, Quach A, et al. 2018. An epigenetic biomarker of aging for lifespan and healthspan. *Aging (Albany NY)* 10:573-591.
- Lu AT, Quach A, Wilson JG, et al. 2019. DNA methylation GrimAge strongly predicts lifespan and healthspan. *Aging (Albany NY)* 11:303-327.
- Belsky DW, Caspi A, Corcoran DL, et al. 2022. DunedinPACE, a DNA methylation biomarker of the pace of aging. *eLife* 11:e73420.
- McEwen LM, O'Donnell KJ, McGill MG, et al. 2020. The PedBE clock accurately estimates DNA methylation age in pediatric buccal cells. *PNAS* 117:23329-23335.
- Knight AK, Craig JM, Theda C, et al. 2016. An epigenetic clock for gestational age at birth based on blood methylation data. *Genome Biol* 17:206.
- Bohlin J, Haberg SE, Magnus P, et al. 2016. Prediction of gestational age based on genome-wide differentially methylated regions. *Genome Biol* 17:207.
- Yang Z, Wong A, Kuh D, et al. 2016. Correlation of an epigenetic mitotic clock with cancer risk. *Genome Biol* 17:205.
- Sugden K, Hannon EJ, Arseneault L, et al. 2020. Patterns of reliability: assessing the reproducibility and integrity of DNA methylation measurement. *Patterns (N Y)* 1:100014.
- Higgins-Chen AT, Thrush KL, Wang Y, et al. 2022. A computational solution for bolstering reliability of epigenetic clocks. *Nat Aging* 2:644-661.
- Pelegi-Siso D, de Prado P, Ronkainen J, et al. 2021. methylclock: a Bioconductor package to estimate DNA methylation age. *Bioinformatics* 37:1759-1760.

## Related Skills

- array-preprocessing - Provides the clean beta matrix clocks consume
- cell-type-deconvolution - IEAA: adjust age acceleration for cell composition
- machine-learning/model-validation - Predictor training, cross-validation, leakage (clocks apply frozen models)
- clinical-biostatistics/survival-analysis - Survival/mortality modeling of age acceleration
- ewas-design - Age-acceleration association study design
- workflows/methylation-pipeline - End-to-end methylation pipeline
