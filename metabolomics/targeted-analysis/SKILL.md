---
name: bio-metabolomics-targeted-analysis
description: Designs and validates quantitative targeted metabolomics assays (MRM/SRM on triple-quadrupole, PRM on high-resolution instruments) to report absolute concentrations. Covers the internal-standard strategy (external cal -> global IS -> standard addition -> stable-isotope-labeled IS), weighted calibration judged by back-calculated %RE not R-squared, ion-ratio quantifier/qualifier confirmation, matrix-effect/recovery characterization, and ICH M10 method validation. Use when quantifying a closed panel of known metabolites with units, building or validating an LC-MS/MS assay, choosing an IS or calibration weighting, or judging whether a reported concentration is trustworthy. For untargeted feature detection see metabolomics/xcms-preprocessing; for group statistics see metabolomics/statistical-analysis; for flux/MID/tracing see metabolomics/isotope-tracing.
tool_type: mixed
primary_tool: skyline
---

## Version Compatibility

Reference examples tested with: R 4.3+, ggplot2 3.5+, Skyline 23.1+, pandas 2.2+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

An absolute concentration requires three inputs the code cannot supply: an authentic reference standard (its certificate-of-analysis purity scales every reported number), a stable-isotope-labeled internal standard that co-elutes with the analyte, and a per-analyte validation record. Without these, the workflow below produces relative peak-area ratios dressed as concentrations.

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Targeted Metabolomics Analysis

**"Quantify these specific metabolites and give me concentrations with units"** -> Fit weighted calibration curves on authentic standards, normalize each analyte to a co-eluting stable-isotope-labeled internal standard, confirm identity by ion ratio, and report concentrations only within the validated range.
- CLI: Skyline builds the small-molecule transition list, integrates peaks, fits weighted calibration, and exports via the Document Grid.
- R / Python: post-export curve fitting, back-calculated %RE checks, IS normalization, ion-ratio confirmation, validation metrics.

## The Single Most Important Modern Insight -- A Concentration Is a Chain of Cancellations, and Matrix Effects Are the Term That Fails to Cancel

Ion suppression is competition for charge and droplet surface in the electrospray source: a co-eluting matrix component (phospholipids late in a reversed-phase gradient, salts at the void) steals ionization from the analyte. Suppression is a property of the co-elution, not of the analyte, so it is retention-time-dependent and lot-dependent. The only mechanism that truly removes it is a stable-isotope-labeled internal standard (SIL-IS) that sits in the identical droplet at the identical instant: the suppression cancels in the analyte/IS area ratio. An IS that elutes even half a minute away samples a different point on the suppression landscape and injects new error rather than removing it. Every other safeguard in this skill -- weighting, ion ratios, validation -- assumes this cancellation is working; the gap between a solvent calibration curve and a matrix-matched curve is a direct readout of how badly the IS is failing.

## Targeted vs Untargeted -- Different Experiments, Not Two Settings

| Axis | Untargeted (discovery) | Targeted (quantification) |
|---|---|---|
| Analyte set | Open -- everything ionizable | Closed -- a panel defined before acquisition |
| Output | Relative fold-change; often putative IDs | Absolute concentration for confirmed analytes |
| Instrument | High-res full-scan / DDA (Orbitrap, QTOF) | Triple-quad SRM/MRM, or high-res PRM |
| Validation | QA/QC framework (Broadhurst, mQACC) | Full bioanalytical validation possible (ICH M10) |
| Question | "What changed?" | "How much is there?" |

Targeted buys sensitivity and absolute quant by spending scope (only what is on the list is seen) and up-front method development. Common pattern: untargeted discovery -> targeted validation of the hits. Feature detection upstream is metabolomics/xcms-preprocessing; this skill begins once the panel and transitions are defined.

## Acquisition Mechanics -- MRM/SRM and PRM

A transition is a precursor-m/z -> product-m/z pair plus a tuned collision energy. On a triple quadrupole, Q1 isolates the precursor, the collision cell fragments it, Q3 isolates one product: the double mass filter is the source of MRM sensitivity. SRM monitors one transition; MRM multiplexes many. Each analyte should carry at least two transitions -- a quantifier (most intense/cleanest, used for concentration) and one or more qualifiers (orthogonal confirmation). PRM replaces Q3 with a high-resolution analyzer that records the full product spectrum in parallel, so transitions are chosen post hoc and isobaric interferences are resolved by exact mass; MRM still wins on absolute sensitivity and very large panels. Dwell time is the signal-accumulation time per transition; cycle time must stay short enough for at least 10-15 points across each chromatographic peak (convention). Scheduled MRM monitors each transition only within a retention-time window so a large panel keeps adequate dwell -- but a peak that drifts out of its window vanishes with no error message, the classic scheduled-MRM failure.

## Decision Tree -- Quant Goal -> IS + Calibration + Validation Depth

| Goal / situation | Internal standard | Calibration | Validation depth | Why |
|---|---|---|---|---|
| Clinical / regulated / PK number | One SIL-IS per analyte (13C/15N) | Multi-level weighted curve, judged by %RE | Full ICH M10 (accuracy, precision, MF, recovery, carryover, stability, ISR) | A number driving a decision must carry its evidence |
| Cross-study quantitative claim | SIL-IS per analyte or per RT/chemical cluster | Multi-level weighted | Accuracy/precision + matrix-factor on QCs | Comparability across runs demands characterized bias |
| Exploratory research, ranking | Few global IS, or per-class | 1/x^2 weighted, low-end %RE checked | Broadhurst/mQACC QC discipline (pooled QC, blanks, RSD filtering) | Relative comparison tolerates residual matrix bias |
| Dirty matrix, isobaric interferences | SIL-IS + high-res | PRM, post-hoc transitions | Selectivity dominated | Exact-mass product resolves co-eluters a unit-resolution Q3 cannot |
| Large standardized panel (600+) | Kit-supplied class IS | Single/limited-point (vendor) | Vendor + bridging study before pooling sites | Kit buys comparability and throughput, not per-analyte full-validation accuracy |
| Carbon source / pathway rate | (tracer, not IS) | -- | -- | Flux question: hand off to metabolomics/isotope-tracing; MID measures rate, not pool size |

The IS rule of thumb: ask how far (in retention time and chemistry) each analyte is from its assigned IS -- that distance is the size of the uncorrected matrix error. 13C/15N at non-exchangeable positions are preferred over deuterium: deuterium causes a small reversed-phase retention shift (the deuterium isotope effect) that can chromatographically separate the IS from its analyte so it stops correcting suppression, and labile deuteriums back-exchange to H. If forced to a deuterated IS, verify co-elution by overlaying analyte and IS chromatograms.

## Calibration and Weighting

| Weighting | When | Effect |
|---|---|---|
| Unweighted (OLS) | Narrow range, near-constant variance | High points dominate; low-end bias on heteroscedastic MS data -- usually wrong |
| 1/x | Moderate range (1-2 orders) | Down-weights high concentrations; restores low-end fit |
| 1/x^2 | Wide range (3+ orders), the common LC-MS default | Aggressively down-weights the top; can over-weight the low end -- still compare, do not reflex |
| Quadratic | Genuine, mechanism-explained curvature (detector saturation) | Never to paper over a bad linear fit |

MS detector response is heteroscedastic -- absolute variance grows with concentration -- so weighting models the variance structure (1/x and 1/x^2 are parametric stand-ins for 1/variance). Select empirically: fit candidate weightings, then pick the one minimizing the sum of absolute back-calculated relative error (%RE) across levels, especially the bottom two or three. R-squared is the wrong instrument: it is dominated by high-leverage top points, so a curve with R-squared 0.999 can be +40% biased at the LLOQ. Use a fitted (non-zero) intercept; forcing the line through the origin re-introduces low-end bias. The blank (matrix only) and zero (matrix + IS) are diagnostic, not calibration points.

### Build a Weighted Calibration Curve With a Back-Calculated %RE Check

**Goal:** Fit a calibration curve on the analyte/IS response ratio and accept it by per-level back-calculation accuracy, not by R-squared.

**Approach:** Fit 1/x^2-weighted linear regression of response ratio on nominal concentration, back-calculate every standard, flag any non-LLOQ level outside +/-15% and the LLOQ outside +/-20%, and set the LLOQ to the lowest passing level.

```r
standards <- data.frame(
  conc = c(1, 5, 10, 25, 50, 100, 250, 500, 1000),
  analyte_area = c(480, 2500, 4900, 12100, 24500, 49000, 121000, 245000, 488000),
  istd_area = c(100000, 98000, 99000, 101000, 100000, 98000, 101000, 99000, 100000)
)
standards$ratio <- standards$analyte_area / standards$istd_area

fit <- lm(ratio ~ conc, data = standards, weights = 1 / standards$conc^2)
standards$back_calc <- (standards$ratio - coef(fit)[1]) / coef(fit)[2]
standards$re_pct <- (standards$back_calc - standards$conc) / standards$conc * 100

# ICH M10: each calibrator within +-15%, +-20% at the LLOQ (lowest level)
tol <- ifelse(standards$conc == min(standards$conc), 20, 15)
standards$pass <- abs(standards$re_pct) <= tol
lloq <- min(standards$conc[standards$pass])
```

### Internal-Standard Normalization

**Goal:** Convert raw analyte area to a matrix-corrected response that the calibration curve maps to concentration.

**Approach:** Divide analyte area by co-eluting SIL-IS area per sample, then invert the same response-ratio calibration; matrix effect and extraction recovery cancel in the ratio.

```r
samples$ratio <- samples$analyte_area / samples$istd_area
samples$conc <- (samples$ratio - coef(fit)[1]) / coef(fit)[2]
samples$conc[samples$conc < lloq] <- NA   # below validated range -> not reportable
```

### Ion-Ratio Confirmation

**Goal:** Guard against quantifying an isobaric co-eluter as the analyte.

**Approach:** Compute the qualifier/quantifier area ratio per sample, compare to the mean calibrator ratio, and flag samples outside the tolerance window -- a drifted ratio means the quantifier peak is partly something else.

```r
cal_ratio <- mean(standards$qualifier_area / standards$quantifier_area)
samples$ion_ratio <- samples$qualifier_area / samples$quantifier_area
# SANTE/2020/12830 uses +-30% relative for LC-MS/MS qualifier/quantifier ratios
samples$id_confirmed <- abs(samples$ion_ratio - cal_ratio) / cal_ratio <= 0.30
```

MRM gives mass selectivity, not identity: two compounds can share a precursor->product transition (many acylcarnitines share m/z 85; lipids share head-group fragments). Identity needs retention time plus the ion ratio plus an authentic standard. Near the LLOQ the qualifier may fall below its own detection limit, so ion-ratio confirmation is usually only enforceable above a few times the LLOQ -- state that limit rather than hiding it. A single-transition method has no defense against isobaric interference and is a documented compromise, not a default.

### LOD and LLOQ

**Goal:** Set the lowest reliably quantifiable concentration from noise and accuracy, not from an extrapolated curve.

**Approach:** Estimate LOD from blank-signal scatter (S/N ~3) and confirm the LLOQ as the lowest calibrator meeting the +/-20% back-calculation and precision criteria; never anchor the curve below the real noise floor to claim sensitivity.

```r
blank_areas <- c(100, 120, 95, 110, 105)
slope <- coef(fit)[2]
lod <- (mean(blank_areas) + 3 * sd(blank_areas)) / slope   # S/N~3 convention
# LLOQ is the lowest calibrator passing +-20% %RE AND precision -- not 10*SD/slope alone
```

## Per-Method Failure Modes

### Matrix suppression unaccounted
- **Trigger:** Neat-solvent calibration, or an IS that does not co-elute with the analyte.
- **Mechanism:** Co-eluting phospholipids/salts suppress analyte ionization; with no co-eluting IS the suppression does not cancel.
- **Symptom:** Solvent and matrix-matched curves disagree; IS-normalized matrix factor far from 1; lot-to-lot drift.
- **Fix:** Co-eluting SIL-IS per analyte; map suppression zones by post-column infusion and move peaks off them; report IS-normalized matrix factor across at least six matrix lots.

### One IS shared across chemically diverse analytes
- **Trigger:** A single global IS used to correct a heterogeneous panel.
- **Mechanism:** The IS corrects suppression and recovery only for analytes co-eluting and chemically near it; distant analytes carry the difference of two suppressions.
- **Symptom:** Excellent CVs but biased group means -- precision (set by the IS correcting injection/drift) and accuracy (set by the per-analyte residual) are decoupled.
- **Fix:** SIL-IS per analyte, or per RT/chemical cluster; treat low CV as no evidence of correctness.

### Unweighted calibration over a wide range
- **Trigger:** OLS fit on heteroscedastic data; acceptance judged by R-squared.
- **Mechanism:** High-concentration points dominate least squares; the low end is fit poorly.
- **Symptom:** R-squared 0.999 yet +30-40% bias at the LLOQ.
- **Fix:** Compare 1/x and 1/x^2, pick by minimizing low-end |%RE|, judge by per-level back-calculation.

### Isotopic crosstalk between analyte and IS
- **Trigger:** IS-analyte mass gap below ~3-4 Da, or high IS:analyte ratio at the LLOQ.
- **Mechanism:** The analyte natural-isotope envelope bleeds into the IS channel (bends the high end, mis-read as saturation); IS isotopic impurity bleeds into the analyte quantifier channel (inflates the low end, biases the LLOQ badly).
- **Symptom:** Non-linear high end; a matrix-blank-plus-IS sample shows signal in the analyte channel.
- **Fix:** Choose an IS mass gap of at least 3-4 Da or a less-abundant SIL isotopologue transition; always run a zero (matrix + IS only) and require its analyte-channel signal below 20% of the LLOQ.

### Pre-analytical degradation
- **Trigger:** Delayed quench, freeze-thaw, slow time-to-freezer.
- **Mechanism:** Metabolism continues post-collection (glycolysis, esterases, redox auto-oxidation); labile metabolites collapse in seconds to minutes.
- **Symptom:** No chromatographic error -- the true value is biased before injection; low/variable adenylate energy charge across samples is the tell that quenching, not biology, drove the numbers.
- **Fix:** Cold (-40 to -80 C) aqueous-organic quench matched to the metabolite, measure freeze-thaw and long-term stability per labile analyte, control collection-to-freeze time as a study variable.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| Calibrator back-calc within +/-15% (+/-20% at LLOQ), >=75% of >=6 levels pass | ICH M10 (Step 4, 2022) | Per-level accuracy, not correlation, defines a usable curve |
| QC accuracy +/-15% (+/-20% at LLOQ); precision CV <=15% (<=20% at LLOQ) | ICH M10 | Intra- and inter-day acceptance at >=4 levels |
| IS-normalized matrix factor CV <=15% across >=6 lots | ICH M10 / Matuszewski 2003 | Proof the IS cancels matrix effect; raw MF may be poor while IS-normalized MF ~1 |
| Carryover <=20% of LLOQ (analyte), <=5% (IS) | ICH M10 | Measured in a blank after the ULOQ; concentration-dependent, must be quantified not eyeballed |
| Selectivity: interference at LLOQ <=20% of analyte, <=5% of IS response | ICH M10 | Across >=6 individual matrix lots |
| ISR: >=2/3 of reanalyzed study samples within +/-20% | ICH M10 | Only test that catches incurred-sample-specific problems spiked QCs cannot |
| Ion-ratio tolerance +/-30% relative (LC-MS/MS) | SANTE/2020/12830 | Illustrative codified window; enforce only above a few times the LLOQ |
| >=10-15 points across a chromatographic peak | Convention | Reliable integration; sets the cycle-time ceiling |
| S/N ~3 = LOD, ~5-10 = LLOQ | Convention | Detection vs reliable quantification; LLOQ also bounded by accuracy/precision |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| Curve accepted on R-squared, biased at LLOQ | Unweighted heteroscedastic fit | Weight (1/x, 1/x^2); accept by per-level back-calculated %RE |
| Deuterated IS gives lot-dependent ratios | Deuterium isotope effect separates IS from analyte; lost matrix correction | Use 13C/15N at non-exchangeable positions, or verify co-elution explicitly |
| High-end curvature mis-read as detector saturation | Analyte natural isotopes bleed into a too-close IS channel | Widen IS-analyte mass gap to >=3-4 Da; use nonlinear isotopic-crosstalk correction |
| Beautiful CVs, wrong group means | One global IS across diverse analytes -- precision/accuracy decoupled | SIL-IS per analyte or per RT/chemical cluster |
| Low samples after a high sample read high | Concentration-dependent carryover | Inject a blank after the ULOQ, randomize run order, report measured carryover |
| Skyline never ratios analyte to IS | IS not tagged Label Type = heavy and paired to its light analyte | Set Label Type heavy in the transition list; pair by molecule name |
| Validated assay, study numbers still wrong | Pre-analytical degradation (no error message) | Quench fast, measure stability, monitor adenylate energy charge |

## References

- MacLean B, Tomazela DM, Shulman N, Chambers M, Finney GL, Frewen B, Kern R, Tabb DL, Liebler DC, MacCoss MJ. 2010. Skyline: an open source document editor for creating and analyzing targeted proteomics experiments. *Bioinformatics* 26(7):966-968.
- Peterson AC, Russell JD, Bailey DJ, Westphall MS, Coon JJ. 2012. Parallel reaction monitoring for high resolution and high mass accuracy quantitative, targeted proteomics. *Molecular & Cellular Proteomics* 11(11):1475-1488.
- Matuszewski BK, Constanzer ML, Chavez-Eng CM. 2003. Strategies for the assessment of matrix effect in quantitative bioanalytical methods based on HPLC-MS/MS. *Analytical Chemistry* 75(13):3019-3030.
- ICH M10 Bioanalytical Method Validation and Study Sample Analysis. ICH Harmonised Guideline, Step 4, adopted 24 May 2022 (FDA implemented November 2022; EMA effective January 2023).
- Broadhurst D, Goodacre R, Reinke SN, Kuligowski J, Wilson ID, Lewis MR, Dunn WB. 2018. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. *Metabolomics* 14(6):72.
- Wang S, Cyronak M, Yang E. 2007. Does a stable isotopically labeled internal standard always correct analyte response? A matrix effect study on a LC/MS/MS method for the determination of carvedilol enantiomers in human plasma. *Journal of Pharmaceutical and Biomedical Analysis* 43(2):701-707.
- Teo G, Chew WS, Burla BJ, Herr DR, Tai ES, Wenk MR, Torta F, Choi H. 2020. MRMkit: automated data processing for large-scale targeted metabolomics analysis. *Analytical Chemistry* 92(20):13677-13682.

## Related Skills

- metabolomics/xcms-preprocessing - Upstream feature detection for untargeted discovery before targeted validation
- metabolomics/statistical-analysis - Group comparison and multivariate analysis of quantified concentrations
- metabolomics/isotope-tracing - Stable-isotope tracing and flux (MID), the adjacent discipline this skill hands off to
- metabolomics/normalization-qc - QC-sample-driven drift correction and RSD filtering
- clinical-biostatistics/cdisc-data-handling - Regulated-trial bioanalysis data handling when targeted numbers feed a clinical study
