# Targeted Metabolomics - Usage Guide

## Overview

Targeted metabolomics quantifies a predefined set of metabolites using selected reaction monitoring (SRM/MRM). This approach provides absolute quantification with high sensitivity and reproducibility.

## Prerequisites

```bash
# Skyline (free, comprehensive)
# Download from: https://skyline.ms/

# R packages
install.packages(c("calibrate", "ggplot2"))

# Python
pip install pandas numpy scipy scikit-learn
```

## Quick Start

Tell your AI agent what you want to do:
- "Build calibration curves and quantify metabolites from my MRM data"
- "Validate my targeted assay with QC sample statistics"

## Example Prompts

### Calibration
> "Build calibration curves for my standard dilution series with 1/x weighting"
> "Calculate R-squared and back-calculated accuracy for each analyte"
> "Determine LOD and LOQ from calibration curve residuals"

### Quantification
> "Calculate absolute concentrations using my calibration curves"
> "Normalize to internal standards before quantification"
> "Apply dilution factors and report final concentrations"

### Validation
> "Calculate accuracy and precision from QC samples at low, medium, and high levels"
> "Check if QC accuracy is within 85-115% acceptance criteria"
> "Report CV% for replicate measurements"

### Quality Assessment
> "Flag samples with concentrations outside the calibration range"
> "Check for carryover using blank samples after high concentration samples"
> "Assess matrix effects using post-extraction spike"

## What the Agent Will Do

1. Import peak areas/heights and standard concentrations
2. Build calibration curves with appropriate weighting
3. Calculate regression statistics and LOD/LOQ
4. Quantify unknowns using calibration
5. Validate with QC sample statistics
6. Export concentrations with quality flags

## Tips

- Use weighted regression (1/x or 1/x^2) for wide concentration ranges
- Include at least 6 calibration points spanning expected sample range
- QC samples at 3 levels (low, medium, high) track assay performance
- Accept calibrators with 85-115% back-calculated accuracy (80-120% at LLOQ)
- Use stable isotope-labeled internal standards when available

## Acceptance Criteria

| Parameter | Threshold |
|-----------|-----------|
| Calibration R^2 | > 0.99 |
| Accuracy (calibrators) | 85-115% |
| Accuracy at LLOQ | 80-120% |
| Precision (CV%) | < 15% |
| Precision at LLOQ | < 20% |

## LOD/LOQ Calculation

- LOD: 3.3 x (residual SD) / slope
- LOQ: 10 x (residual SD) / slope

## Standard Curve Levels

| Level | Purpose |
|-------|---------|
| Blank | No analyte, check contamination |
| LLOQ | Lower limit of quantification |
| Low | 3x LLOQ |
| Medium | Mid-range |
| High | 75% of ULOQ |
| ULOQ | Upper limit |

## Software Options

- **Skyline** - Free, comprehensive
- **TraceFinder** - Thermo
- **MassHunter** - Agilent
- **MultiQuant** - SCIEX

## References

- FDA Bioanalytical Method Validation Guidance (2018)
- EMA Guideline on bioanalytical method validation
- Skyline: doi:10.1093/bioinformatics/btq054
