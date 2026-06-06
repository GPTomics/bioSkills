# Model Interpretation Usage Guide

## Overview

Explain ML predictions on omics data with SHAP, LIME, and permutation importance, while respecting the boundaries that make attributions trustworthy: an attribution describes the model, not biology; the conditional-vs-interventional Shapley choice changes which correlated genes get credit; and a SHAP ranking is a debugging and hypothesis tool, never a validated biomarker panel.

## Prerequisites

```bash
pip install shap lime scikit-learn xgboost pandas
```

Conceptual prerequisites: a model that already generalizes (interpreting an over-fit or batch-confounded model explains an artifact), a representative background dataset for interventional SHAP, and a co-expression module map for aggregating attributions over correlated genes.

## Quick Start

Tell your AI agent what you want to do:
- "Explain my classifier with interventional SHAP and a representative background"
- "Use SHAP to check whether my model is keying on batch instead of biology"
- "Aggregate SHAP over co-expression modules before ranking genes"
- "Explain one prediction locally, clearly labeled as model-internal"

## Example Prompts

### Choosing the estimand

> "Which genes is my random forest actually keying on? Use interventional SHAP with a background, not the path-dependent default, and tell me the estimand."

> "I want a descriptive view of which genes are informative here. Use path-dependent SHAP and explain why it is not the same as model reliance."

### Debugging shortcuts

> "Use SHAP to check whether my classifier's top features are batch indicators or housekeeping genes tracking sequencing depth."

### Correlated genes

> "My top SHAP genes are co-expressed. Aggregate the attributions over co-expression modules so I do not over-interpret the within-module order."

### Local explanation

> "Explain why the model predicted disease for sample 12 with a SHAP waterfall, and note it is background-relative and model-internal."

## What the Agent Will Do

1. Confirm the model generalizes before interpreting it
2. Choose the Shapley conditioning for the question (interventional for reliance, path-dependent for description) and state it
3. Supply a representative background and report the attribution scale (log-odds vs probability)
4. Aggregate attributions over correlated modules before ranking
5. Use attribution to audit for batch/shortcut learning
6. Route any "select a biomarker panel" request to biomarker-discovery with independent validation

## Tips

- A high-SHAP gene can be a correlate of a batch shortcut the model exploited -- attribution explains the model, not biology
- `feature_perturbation='auto'` (shap 0.47+) silently flips the estimand on whether you pass `data=`; set it explicitly
- Path-dependent SHAP can give a gene nonzero credit even if the model never uses it (correlation leak); interventional gives unused genes zero
- The split of credit among correlated genes is mode-dependent and not identifiable from biology -- aggregate over modules before ranking
- LIME is non-reproducible across seeds and kernel widths; use it only to eyeball one local prediction, never for global ranking
- Permutation importance breaks under correlation the same way SHAP does; cluster features or use conditional permutation
- The best, most trustworthy use of attribution is catching shortcut/batch learning, not selecting biomarkers

## Related Skills

- machine-learning/omics-classifiers - Train the model being explained; debug batch shortcuts
- machine-learning/biomarker-discovery - Validated feature selection (SHAP ranking is not selection)
- machine-learning/model-validation - Confirm the model generalizes before interpreting it
- data-visualization/heatmaps-clustering - Visualize module-aggregated attributions
