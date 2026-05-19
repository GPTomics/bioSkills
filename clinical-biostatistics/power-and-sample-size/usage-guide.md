# Power and Sample Size - Usage Guide

## Overview

Computes sample size and statistical power for clinical trials covering superiority, non-inferiority, equivalence, and bioequivalence designs across continuous, binary, and time-to-event endpoints. Implements FDA 2016 NI margin selection (M1/M2 with double discount), Schoenfeld 1981 events formula with Lakatos 1988 extension for non-proportional hazards, Schuirmann TOST for equivalence, and the critical δ-vs-MCID distinction.

## Prerequisites

```bash
pip install statsmodels scipy numpy pandas
```

R is strongly recommended for production sample-size work, especially complex survival:

```r
install.packages(c('pwr', 'gsDesign', 'rpact', 'presize', 'npsurvSS', 'simtrial', 'nphRCT'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Compute the sample size for a two-arm trial to detect a 5-point difference in HbA1c with 80% power at α=0.025 one-sided"
- "Calculate events needed for OS in oncology with HR=0.70 and 90% power"
- "Justify the NI margin for my antibiotic trial using FDA 2016 double discount with historical M1"
- "Run simulation-based sample size for immunotherapy with delayed-effect HR(t)"
- "TOST equivalence sample size with margin ±5 and assumed mean difference 0"

## Example Prompts

### Continuous endpoints

> "Two-arm parallel trial; primary endpoint is change in SBP at 24 weeks; expected effect 5 mmHg; SD 12; 80% power; α=0.05 two-sided; expected dropout 15%. Compute n per arm."

> "Re-power the same trial under stratified randomisation by region with R² of 0.10 against outcome. Compute the efficiency gain."

### Time-to-event

> "Schoenfeld events formula for OS with HR=0.65, α=0.025 one-sided, 90% power, 2:1 allocation. Convert to n per arm with expected event rate 35% over 24-month follow-up."

> "Immunotherapy with expected delayed effect: HR=1.0 for first 6 months then HR=0.50 thereafter. Schoenfeld formula will under-estimate. Use simulation via R simtrial."

> "Lakatos 1988 sample size in R gsDesign for survival with non-uniform accrual (linear ramp), 15% annual dropout, control hazard 0.04/month, treatment HR 0.70."

### Non-inferiority

> "Antibiotic NI trial. Historical placebo-vs-active effect 25% (lower CI bound 20%). Apply FDA 2016 double discount: M2 = 50% of 20% = 10% margin. Compute n per arm for 80% power."

> "My NI margin is 5% RD. Use Miettinen-Nurminen score CI for the RD per regulatory expectation. Compute n with control rate 85%."

### Equivalence / Bioequivalence

> "Bioequivalence study of generic vs reference. Geometric mean ratio range (0.80, 1.25). CV from pilot 20%. Compute n for two-period crossover."

> "TOST for clinical equivalence in PRO outcome. Margin ±5 on a 100-point scale; assume true difference 0; SD 10; 80% power."

### Multiplicity / Co-primary

> "Two co-primary endpoints; per FDA 2022 Multiple Endpoints Guidance, both must reach significance. Compute n at each endpoint's per-marginal power = 90% so joint power is 0.81. Compare to inflating n for joint 80%."

### Adaptive / SSR

> "Pre-specify Mehta-Pocock promising-zone SSR. Initial n=300/arm; if interim CP ∈ (0.3, 0.8), increase to maximum n=450/arm. Use CHW weights for Type-I control."

## What the Agent Will Do

1. Verify the endpoint type (continuous, binary, TTE) and the design (superiority, NI, equivalence)
2. Distinguish δ (alternative effect to detect) from MCID (clinical meaningful difference)
3. Compute sample size with appropriate formula (Fleiss, Schoenfeld, Lakatos, TOST)
4. For NI: validate margin per FDA 2016 (M2 <= 0.5 × historical M1 lower CI; M2 <= 0.5 × MCID)
5. For non-PH survival: switch from Schoenfeld to Lakatos or simulation
6. Inflate for expected dropout
7. Adjust for stratified randomisation efficiency (Senn 2013 precision gain)
8. Report n per arm with explicit assumptions and sensitivity analyses

## Tips

- **δ is what you want to detect; MCID is what is clinically meaningful.** Postdoc rule: δ >= 1.5 × MCID for superiority. Setting δ = MCID produces underpowered trials.
- **Schoenfeld 1981 assumes proportional hazards.** Under non-PH (immuno-oncology), it under-estimates events by 20-50%. Switch to Lakatos 1988 or simulation under expected HR(t) pattern.
- **For NI margins, apply the FDA 2016 double discount:** M2 = 50% of the *lower CI bound* of historical M1, not the point estimate. M2 should also be <= 0.5 × MCID.
- **Constancy assumption for NI is unverifiable** but critical. Hung-Wang-O'Neill 2005: NI trials have NO within-trial Type-I error guarantee -- alpha is conditional on constancy.
- **Temple-Ellenberg 2000 is *Ann Intern Med* 133:455-463, NOT NEJM** (common citation error). Established the assay sensitivity framework.
- **`power.prop.test` in R uses uncorrected normal approximation** -- can over-state power 5-10% for n<100/arm. For confirmatory work with small n, use simulation-based SS.
- **TOST is closed under intersection-union** -- no multiplicity adjustment despite two tests.
- **MCID has two derivation traditions:** anchor-based (Jaeschke 1989) and distribution-based (0.5 SD heuristic, Norman 2003). Report both when available.
- **Unblinded SSR scares FDA** because of DMC firewall risk. Blinded SSR (Friede-Kieser 2006) is uncontroversial; unblinded with CHW weights (Cui-Hung-Wang 1999) is conditionally acceptable.
- **Mehta-Pocock promising zone has been called "stealth alpha inflation"** by Jennison-Turnbull 2015. Pre-specify carefully and report operating characteristics.
- **Cluster-randomised trials need design-effect adjustment**: multiply individual-level n by 1 + (m-1)ICC. ICC misspecification is the single largest source of underpowering.

## Related Skills

- clinical-biostatistics/survival-analysis - TTE sample size (Schoenfeld, Lakatos)
- clinical-biostatistics/effect-measures - δ scales (OR, RR, RD)
- clinical-biostatistics/categorical-tests - Binary endpoint test selection
- clinical-biostatistics/multiplicity-graphical - Power for co-primary
- clinical-biostatistics/adaptive-designs - SSR, promising zone, group-sequential
- clinical-biostatistics/bayesian-trials - Predictive probability of success
- clinical-biostatistics/trial-reporting - CONSORT 2025 sample-size justification
- experimental-design/sample-size - General methods
- experimental-design/power-analysis - General power methods
