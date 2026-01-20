# Targeted Metabolomics Usage Guide

## Overview

Targeted metabolomics focuses on quantifying a predefined set of metabolites using selected reaction monitoring (SRM/MRM). This approach provides absolute quantification with high sensitivity and reproducibility.

## When to Use This Skill

- Absolute quantification required
- Known metabolites of interest
- Biomarker validation studies
- Pharmacokinetic studies
- Clinical assay development

## Comparison: Targeted vs Untargeted

| Aspect | Targeted | Untargeted |
|--------|----------|------------|
| Metabolites | Pre-defined panel | Discovery-based |
| Quantification | Absolute (with standards) | Relative |
| Sensitivity | Higher (optimized) | Variable |
| Throughput | Lower setup, faster run | Higher discovery |
| Method development | Required per analyte | Generic |

## Standard Curve Requirements

### Calibration Levels
- Minimum 6 concentration points
- Span expected sample range
- Include blank and LLOQ

### Regression Models
- Linear: `y = mx + b`
- Weighted linear: `1/x` or `1/x²` weighting
- Quadratic: For wide dynamic ranges

### Acceptance Criteria
- R² > 0.99 (preferably > 0.995)
- Back-calculated accuracy 85-115%

## Quality Control Samples

| QC Level | Concentration | Purpose |
|----------|---------------|---------|
| LLOQ | At quantification limit | Sensitivity check |
| Low | 3× LLOQ | Low range accuracy |
| Medium | Mid-range | Mid accuracy |
| High | 75% of ULOQ | High range accuracy |

## Validation Parameters

### Accuracy
- Measure of trueness
- (Measured / Nominal) × 100%
- Accept: 85-115% (80-120% at LLOQ)

### Precision
- CV% of replicate measurements
- Accept: CV < 15% (< 20% at LLOQ)

### Linearity
- R² of calibration curve
- Accept: R² > 0.99

### LOD/LOQ
- LOD: 3.3 × σ / slope
- LOQ: 10 × σ / slope

## Internal Standards

### Selection Criteria
- Stable isotope labeled (preferred)
- Similar ionization to analyte
- Similar retention time
- Not present in samples

### Common IS Types
- ¹³C-labeled
- D-labeled (deuterated)
- Structural analogs

## Common Issues

### Matrix effects
- Use matrix-matched standards
- Evaluate with post-extraction spike

### Carryover
- Include blank after high samples
- Consider wash steps

### Stability
- Assess autosampler stability
- Freeze-thaw cycles

## Software Options

- **Skyline** - Free, comprehensive
- **TraceFinder** - Thermo
- **MassHunter** - Agilent
- **MultiQuant** - SCIEX

## References

- FDA Bioanalytical Method Validation Guidance (2018)
- EMA Guideline on bioanalytical method validation
- Skyline: doi:10.1093/bioinformatics/btq054
