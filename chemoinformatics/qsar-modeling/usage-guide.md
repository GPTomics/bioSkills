# QSAR Modeling Usage Guide

## Overview

Build target-specific QSAR / QSPR models from in-house assay data using chemprop D-MPNN (modern default), RandomForest + ECFP4 (small data baseline), or transformer-based MolFormer / Uni-Mol / ChemBERTa for large datasets. Apply OECD 5 principles with scaffold-balanced splits, ensemble uncertainty, applicability domain, and conformal prediction.

## Prerequisites

```bash
pip install chemprop rdkit scikit-learn mapie shap
```

chemprop 2.0+ (note major API change from 1.x).

## Quick Start

Tell the AI agent what to do:
- "Train chemprop D-MPNN classifier on this hERG dataset with scaffold split"
- "Build RandomForest QSAR on 150 compounds (small data); use ECFP4 + nested CV"
- "Add conformal prediction intervals to my chemprop classifier"
- "Apply applicability domain filter (kNN distance) to new predictions"
- "Compute SHAP atomic attribution to interpret hERG predictions"

## Example Prompts

### Small-data baseline
> "Build RF + ECFP4 binary classifier on 150-compound dataset (hERG_blocker column). 5-fold scaffold-split cross-validation. Report AUC, sensitivity, specificity, applicability domain by kNN Tanimoto."

### chemprop classifier
> "Train chemprop 2.0 D-MPNN ensemble on hERG_data.csv. Use scaffold_balanced split, 5 folds x 5 ensemble. Include rdkit_2d_normalized features. Report test AUC + calibration."

### Conformal prediction
> "Wrap my trained chemprop model in MAPIE conformal predictor. Use 90% coverage. For each new SMILES, report prediction + 90% interval."

### Multi-task CYP inhibition
> "Train chemprop multitask classifier on CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4. Scaffold split. Compare to per-target single-task models."

## What the Agent Will Do

1. Standardize input data (recommend `chemoinformatics/molecular-standardization`).
2. Apply scaffold split (Bemis-Murcko) for train / val / test.
3. Featurize with ECFP4 + RDKit 2D descriptors (for chemprop hybrid).
4. Train ensemble (5 folds × 5 seeds for chemprop; bootstrap for RF).
5. Compute test metrics + calibration plots.
6. Build applicability domain assessment (kNN Tanimoto distance or conformal prediction).
7. Optionally compute SHAP feature importance.

## Tips

- Always scaffold-split; never random split.
- For < 200 compounds, RF + ECFP4 is competitive and more interpretable.
- For > 10k compounds, transformer-based methods may outperform chemprop.
- Calibration is critical: deep learning probabilities rarely calibrated; Platt-correct.
- Conformal prediction gives distribution-free 90% coverage interval.
- Use ensemble variance as applicability domain measure.

## Related Skills

- chemoinformatics/molecular-descriptors - Featurization choices
- chemoinformatics/molecular-standardization - Mandatory upstream
- chemoinformatics/scaffold-analysis - Scaffold split implementation
- chemoinformatics/admet-prediction - ADMET-specific QSAR
- chemoinformatics/generative-design - QSAR as scoring component
- machine-learning/model-validation - General ML validation
- machine-learning/biomarker-discovery - Adjacent ML approaches
