---
name: bio-metabolomics-isotope-tracing
description: Designs and analyzes stable-isotope-resolved metabolomics (SIRM / isotope tracing / fluxomics) experiments that measure metabolic ACTIVITY via 13C/15N/2H tracers, distinct from steady-state pool profiling. Covers tracer choice, isotopologue vs isotopomer, mass-isotopomer distributions (MID), fractional enrichment, the mandatory natural-abundance + tracer-purity correction (IsoCor, AccuCor), and the metabolic/isotopic steady-state vs non-stationary (INST-MFA) distinction. Use when feeding a labeled tracer and interpreting labeling patterns, correcting raw isotopologue intensities, computing or plotting an MID, or deciding tracing vs abundance profiling. For absolute pool concentration and MRM mechanics see metabolomics/targeted-analysis; for constraint-based genome-scale flux (FBA, not empirical tracing) see systems-biology/flux-balance-analysis; for feature detection see metabolomics/xcms-preprocessing; for pathway enrichment that ignores the pool-vs-flux caveat see metabolomics/pathway-mapping.
tool_type: python
primary_tool: isocor
---

## Version Compatibility

Reference examples tested with: isocor 2.2+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Isotope Tracing / Stable-Isotope-Resolved Metabolomics

**"What is this pathway actually doing, not just how much metabolite is there?"** -> Feed a labeled tracer, then measure how label propagates into downstream metabolites as a mass-isotopomer distribution (MID) over time.
- Python: `isocor.mscorrectors.MetaboliteCorrectorFactory().correct()` (IsoCor) for natural-abundance correction
- R: `accucor::natural_abundance_correction()` (AccuCor) for high-resolution correction
- Modeling layer (separate discipline): `INCA`, `13CFLUX2`, `OpenFLUX` for 13C-MFA / INST-MFA flux fitting

## The Single Most Important Insight -- Labeling Reports Flux; Pool Size Does Not

A metabolite's concentration is *how much* is there; its labeling pattern (MID) is *where the carbon came from and how fast it got there*. These are independent measurements and frequently move in OPPOSITE directions: block a downstream-consuming enzyme and the intermediate pool rises (it backs up) while the labeling of downstream products falls (flux through them dropped). Reading the pool alone reports the opposite of the biology. An isotope-tracing experiment therefore answers a fundamentally different question than untargeted or targeted abundance profiling, and its analysis is dominated by two mandatory corrections that fabricate flux if skipped: natural-abundance / tracer-purity correction (raw isotopologue areas are NOT the labeling), and the steady-state assumption (a single MID is a snapshot whose meaning depends on whether labeling has plateaued). Fractional enrichment is concentration-independent (it is a ratio within one pool), which is why it survives the recovery/matrix problems that plague absolute quant -- but it says nothing about amount.

## Core Concepts

| Term | Meaning | Why it matters |
|---|---|---|
| Tracer / tracee | The labeled substrate fed (tracer, e.g. U-13C6-glucose) vs the unlabeled endogenous pool (tracee) | The experiment measures how tracer atoms replace tracee atoms over time |
| Isotopologue | A molecule differing only in number of heavy atoms (M+0, M+1, M+2 ...) | Resolved by MASS; this is what MS measures and what an MID counts |
| Isotopomer | Same number of heavy atoms but at different POSITIONS (e.g. 1-13C vs 6-13C lactate) | Resolved by POSITION; needs NMR or positional tracers, NOT mass spectra alone |
| MID (mass-isotopomer distribution) | The fractional vector of M+0, M+1, ... for one metabolite | The primary readout; its shape encodes which route carbon took |
| Fractional / mean enrichment | Weighted-mean labeled-atom fraction = sum(i * MID_i) / n_atoms | One-number summary of how labeled a pool is; concentration-independent |
| Atom transitions | The map of which substrate carbons land on which product carbons per reaction | Defines the expected MID for each pathway; the basis of flux models |
| Metabolic steady state | Pool sizes constant in time | Required for classical MFA; if pools drift, plateau MIDs do not give fluxes |
| Isotopic steady state | Labeling has equilibrated to a stable plateau | Classical MFA reads fluxes from the plateau; sampling before it is invalid |

## Decision Tree by Scenario

| Goal / situation | Do | Why |
|---|---|---|
| Want amount/concentration, units, biomarker level | Use abundance profiling -> metabolomics/targeted-analysis | Pool size is not flux; tracing cannot give a concentration |
| Want pathway ACTIVITY/route, central carbon metabolism | 13C tracing (U-13C6-glucose, 13C5-glutamine); measure MIDs | Labeling reports flux through the route the carbon took |
| Trace nitrogen handling (transamination, urea, nucleotides) | 15N tracer (e.g. 15N2-glutamine, 15N-ammonia) | N-flux is invisible to a 13C tracer |
| Distinguish two carbon entry points into one pool | Positional / partially-labeled tracer (e.g. 1,2-13C2-glucose) | The M+1 vs M+2 split of products separates PPP from glycolysis |
| Fast-labeling small pools, cultured cells, clear metabolic steady state | Steady-state 13C-MFA from plateau MIDs (INCA, 13CFLUX2) | Plateau labeling + atom transitions -> flux estimates |
| Slow labeling, large pools, autotrophs, primary/quiescent cells | INST-MFA from the labeling TIME COURSE (INCA) | Drops the isotopic-steady-state assumption; fits transient + pool sizes |
| Have raw isotopologue areas (low-res QqQ / high-res Orbitrap) | Natural-abundance + purity correction FIRST (IsoCor / AccuCor) | Uncorrected MID is wrong by construction; see below |
| Want genome-scale predicted flux without a tracer | systems-biology/flux-balance-analysis | FBA is constraint-based prediction, NOT empirical label measurement |

## Natural-Abundance + Tracer-Purity Correction (mandatory)

**Goal:** Turn raw measured isotopologue areas into a true MID that reflects only tracer-derived label.

**Approach:** Even a fully unlabeled molecule shows an M+1, M+2 ladder because ~1.07% of carbon is naturally 13C (plus 15N, 2H, 18O, 34S, and derivatization Si). Build the natural-abundance ladder from the molecular (and derivative) formula, deconvolve it out, then correct for the tracer not being 100% isotopically pure. Feeding uncorrected areas to a flux model is the equivalent of reporting an uncalibrated peak area as a concentration.

```python
import isocor

# corrector knows the formula's natural-abundance ladder and the tracer
corrector = isocor.mscorrectors.MetaboliteCorrectorFactory(
    'C6H12O6', tracer='13C',
    correct_NA_tracer=True,           # also strip the labeled element's own natural abundance
    tracer_purity=[0.01, 0.99])       # [unlabeled, labeled] per-position purity of the tracer

# raw measured areas M+0..M+6 for a partially labeled glucose pool
corrected_area, iso_fraction, residuum, mean_enrichment = corrector.correct(
    [50000., 8000., 12000., 3000., 1500., 6000., 25000.])
# iso_fraction is the corrected MID; mean_enrichment is fractional enrichment
```

High-resolution Orbitrap data resolves 13C from 15N/2H by exact mass, enabling a different (often simpler) correction; AccuCor (R) is tuned for that case:

```r
library(accucor)
# El-MAVEN / MAVEN isotopologue table; Resolution is the instrument resolving power
corrected <- natural_abundance_correction(path = 'elmaven_export.xlsx',
                                          resolution = 100000, purity = 0.99)
```

Pick the corrector by tracer count and resolution: IsoCor handles any tracer at any resolution; AccuCor (single tracer) and AccuCor2 (dual 13C-15N / 13C-2H) target high-res. Verify the chosen tool's current argument names before running -- both APIs drift across versions.

## Computing and Plotting an MID / Fractional Enrichment

**Goal:** Summarize a corrected isotopologue vector as an MID and one fractional-enrichment number, comparably across conditions.

**Approach:** Normalize corrected areas to sum 1 (the MID), then take the atom-weighted mean over isotopologue index divided by the number of tracer atoms.

```python
import numpy as np

corrected = np.array([26000., 2200., 5600., 1200., 500., 2300., 12500.])
mid = corrected / corrected.sum()                              # M+0..M+n fractions
fractional_enrichment = np.sum(np.arange(len(mid)) * mid) / (len(mid) - 1)
# stacked-bar MID per condition is the standard visualization; never plot raw (uncorrected) areas
```

## Steady-State Check (does the labeling number mean anything yet?)

**Goal:** Decide whether a measured MID may be read as flux-informative or is still a kinetic transient.

**Approach:** Sample labeling at several timepoints; isotopic steady state is reached when the MID stops changing (plateau). Only plateau MIDs license classical-MFA flux inference; a rising MID is kinetic data requiring INST-MFA.

```python
import numpy as np

# fractional enrichment per timepoint (minutes) for one metabolite
t = np.array([0, 5, 15, 30, 60, 120])
fe = np.array([0.00, 0.18, 0.31, 0.39, 0.42, 0.43])
reached_plateau = abs(fe[-1] - fe[-2]) < 0.02     # <2% change between last points = plateau
# if not reached_plateau: the pool is still labeling -> use the full time course (INST-MFA), not one point
```

## Per-Method Failure Modes

### Skipping natural-abundance correction
- **Trigger:** Reporting or modeling raw isotopologue areas straight from El-MAVEN / Skyline.
- **Mechanism:** ~1.07% natural 13C (plus 15N, 2H, derivatization Si) creates an M+1/M+2 ladder on every molecule independent of the tracer; raw M+1 is mostly natural abundance for short-chain metabolites.
- **Symptom:** Apparent labeling in unlabeled controls; inflated M+1; flux fits with tight CIs that are simply wrong.
- **Fix:** Always run IsoCor/AccuCor with the correct formula (and derivative formula for GC-MS), tracer element, and tracer purity before any interpretation.

### Assuming steady state when it is not reached
- **Trigger:** Inferring flux from a single early-timepoint MID.
- **Mechanism:** Classical MFA assumes both metabolic AND isotopic steady state; a transient MID encodes kinetics, not the flux plateau.
- **Symptom:** Fluxes that change with sampling time; large residuals; biologically implausible splits.
- **Fix:** Verify plateau across a time course, or switch to INST-MFA (INCA) which fits the transient and estimates pool sizes too.

### Tracer impurity ignored
- **Trigger:** Treating a "U-13C6" tracer as 100% labeled.
- **Mechanism:** Per-position purity is ~99%, so a fraction of tracer molecules carry a 12C, distorting the fully-labeled isotopologue; the error compounds with atom count.
- **Symptom:** Fully-labeled isotopologue (M+n) systematically under-counted; enrichment biased low.
- **Fix:** Supply the measured tracer purity to the corrector (`tracer_purity` / `purity`).

### Pool-size-vs-labeling confound
- **Trigger:** Concluding "flux changed" from a changed pool concentration (or vice versa).
- **Mechanism:** Pool and labeling are independent and can anticorrelate; a rising intermediate can mean LESS downstream flux.
- **Symptom:** Pool-based and label-based conclusions disagree; "activation" that is actually a backup.
- **Fix:** Interpret MID/enrichment for flux and concentration for amount separately; report both, never substitute one for the other.

### Quench/extraction continuing turnover
- **Trigger:** Slow quench between harvest and metabolism arrest.
- **Mechanism:** High-turnover metabolites keep reacting post-harvest, scrambling labeling before extraction.
- **Symptom:** Variable, sample-dependent MIDs; collapsed nucleotide/energy-charge metabolites.
- **Fix:** Fast cold quench (-40 to -80 C aqueous methanol/acetonitrile); standardize and minimize harvest-to-quench time.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| 13C natural abundance ~1.07% | IUPAC isotopic composition | Sets the natural-abundance ladder corrected out of every MID |
| Tracer purity ~99% per position | Vendor U-13C specs | Must be supplied to correction; compounds with atom count |
| Isotopic-steady-state = <~2% MID change between timepoints | Convention | Below this, plateau reached; classical MFA licensed |
| Quench at -40 to -80 C aqueous organic | Quenching literature (convention) | Arrests metabolism fast enough for high-turnover pools |
| INST-MFA when labeling is slow / pools large / autotrophic | Cheah & Young 2018 | Isotopic steady state is unreachable in time, so fit the transient |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `correct()` length mismatch in IsoCor | Measurement vector is not n_tracer_atoms + 1 long | Pass M+0..M+n with n = count of tracer-element atoms in the formula |
| Labeling appears in unlabeled control | No natural-abundance correction | Run IsoCor/AccuCor before interpreting |
| M+n isotopologue under-reported | Tracer purity left at 1.0 | Set `tracer_purity` / `purity` to the measured value |
| GC-MS MID still wrong after correction | Derivatization atoms (TMS/TBDMS Si, extra C) omitted | Provide the derivative formula to the corrector |
| Flux estimates shift with sampling time | Isotopic steady state not reached | Use a time course + INST-MFA, not a single MID |
| `ValueError` half-defined resolution in IsoCor | Gave `mz_of_resolution`/`charge` without `resolution` | Provide all high-res parameters together or none |

## References

- Millard P, Delepine B, Guionnet M, Heuillet M, Bellvert F, Letisse F. 2019. IsoCor: isotope correction for high-resolution MS labeling experiments. *Bioinformatics* 35:4484-4487.
- Su X, Lu W, Rabinowitz JD. 2017. Metabolite Spectral Accuracy on Orbitraps. *Analytical Chemistry* 89:5940-5948.
- Clasquin MF, Melamud E, Rabinowitz JD. 2012. LC-MS Data Processing with MAVEN: A Metabolomic Analysis and Visualization Engine. *Current Protocols in Bioinformatics* 37:14.11.1-14.11.23.
- Cheah YE, Young JD. 2018. Isotopically nonstationary metabolic flux analysis (INST-MFA): putting theory into practice. *Current Opinion in Biotechnology* 54:80-87.
- Young JD. 2014. INCA: a computational platform for isotopically non-stationary metabolic flux analysis. *Bioinformatics* 30:1333-1335.
- Antoniewicz MR. 2018. A guide to 13C metabolic flux analysis for the cancer biologist. *Experimental & Molecular Medicine* 50:1-13.

## Related Skills

- metabolomics/targeted-analysis - Absolute pool quantification and MRM/SRM mechanics
- metabolomics/xcms-preprocessing - Upstream LC-MS feature detection
- metabolomics/pathway-mapping - Pathway enrichment that interprets pools, not flux
- systems-biology/flux-balance-analysis - Constraint-based predicted flux, distinct from empirical tracing
