# Model Interpretation Usage Guide

## Overview

Explain machine learning predictions on omics data using SHAP values for global and local feature attribution, and LIME for model-agnostic explanations.

## Prerequisites

```bash
pip install shap lime scikit-learn xgboost matplotlib pandas
```

## Quick Start

Tell your AI agent what you want to do:
- "Explain my classifier predictions with SHAP"
- "Which genes are most important for my model's predictions?"
- "Show me a SHAP summary plot for my classifier"
- "Use LIME to explain individual predictions"

## Example Prompts

### Global Feature Importance

> "I trained a random forest on gene expression data. Use SHAP to show me which genes are most important across all predictions."

> "Create a SHAP beeswarm plot for my XGBoost classifier. Show the top 20 features."

### Individual Predictions

> "Use SHAP to explain why my classifier predicted 'disease' for sample X. Show a waterfall plot."

> "Generate LIME explanations for the 5 samples with the highest prediction confidence."

### Feature Ranking

> "Extract the mean absolute SHAP values for all features and save to a CSV file."

### Interactions

> "Create a SHAP dependence plot for BRCA1 to see how its expression level affects the prediction."

## What the Agent Will Do

1. Load trained model and test data
2. Create appropriate SHAP explainer (TreeExplainer for tree models)
3. Calculate SHAP values for test samples
4. Generate summary visualizations (beeswarm, bar)
5. Save top features to file

## Tips

- Use `explainer(X)` not `.shap_values(X)` for SHAP v0.47+
- TreeExplainer is fastest for tree-based models (RF, XGBoost, LightGBM)
- LIME is model-agnostic but slower; use for non-tree models
- Beeswarm plots show both importance and direction of effect
- For multi-class, SHAP values have an extra dimension for classes
- SHAP values sum to the difference between prediction and expected value
