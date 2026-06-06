# machine-learning

## Overview

Decision-grade machine learning for omics and biomedical data: feature selection, classification, survival prediction, model interpretation, validation, and single-cell reference mapping, with the p>>n, leakage, calibration, and correlated-feature traps made explicit.

**Tool type:** python | **Primary tools:** scikit-learn, scikit-survival, shap, scvi-tools, boruta

## Skills

| Skill | Description |
|-------|-------------|
| biomarker-discovery | All-relevant vs minimal-optimal selection (Boruta, mRMR, elastic-net, stability), and why most signatures do not replicate |
| omics-classifiers | Regularized logistic, random forest, and gradient-boosted classifiers for p>>n, with batch-shortcut, imbalance, and calibration handling |
| model-validation | Nested CV, the full leakage taxonomy, calibration vs discrimination, decision-curve net benefit, and TRIPOD+AI reporting |
| prediction-explanation | SHAP, LIME, and permutation importance with the conditional-vs-interventional and attribution-is-not-causation boundaries |
| survival-analysis | Predictive time-to-event models (penalized Cox, RSF, deep survival) with Uno's C, time-dependent AUC, integrated Brier, and competing risks |
| atlas-mapping | Single-cell reference mapping and label transfer (scArches, Symphony, CellTypist) with out-of-distribution uncertainty gating |

## Example Prompts

- "Select an all-relevant biomarker set with Boruta and report its stability index"
- "Build an elastic-net classifier and check whether my high AUC is a batch artifact"
- "Run nested cross-validation and report calibration, not just AUC"
- "Explain my classifier with interventional SHAP and detect shortcut learning"
- "Fit a penalized Cox prognostic signature and evaluate it with Uno's C and integrated Brier"
- "Map my single-cell query to a reference atlas and flag cells that do not belong"

## Requirements

```bash
pip install scikit-learn scikit-survival xgboost shap lime lifelines Boruta mrmr-selection scvi-tools celltypist imbalanced-learn
```

## Related Skills

- **single-cell/preprocessing** - QC and normalization before atlas mapping
- **differential-expression/de-results** - Pre-filter genes and pseudobulk validation
- **experimental-design/multiple-testing** - FDR control and sample-size planning for discovery
- **clinical-biostatistics/survival-analysis** - Confirmatory Kaplan-Meier, log-rank, and classical Cox inference
- **pathway-analysis/go-enrichment** - Functional enrichment of selected features
- **workflows/biomarker-pipeline** - End-to-end orchestration of selection, validation, and interpretation
