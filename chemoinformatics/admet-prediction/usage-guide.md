# ADMET Prediction - Usage Guide

## Overview
Predict absorption, distribution, metabolism, excretion, and toxicity properties for drug discovery. Includes drug-likeness filtering and structural alerts.

## Prerequisites
```bash
pip install rdkit requests
pip install deepchem  # For ML-based predictions
```

## Quick Start
Tell your AI agent what you want to do:
- "Predict ADMET properties for my lead compounds"
- "Filter my library for drug-like compounds"
- "Check for PAINS alerts in my hit list"
- "Predict hERG liability for my compounds"

## Example Prompts

### ADMET Prediction
> "Use ADMETlab 3.0 to predict ADMET properties for these SMILES."

> "Predict CYP inhibition profiles for my lead series."

### Drug-Likeness
> "Calculate Lipinski violations and QED scores for my compounds."

> "Annotate Lipinski and Veber results, then rank using thresholds justified for this project."

### Safety Filtering
> "Check my hits for PAINS and other structural alerts."

> "Identify compounds with potential hERG liability."

## What the Agent Will Do
1. Calculate drug-likeness properties (Lipinski, QED)
2. Use the current official ADMETlab 3.0 API workflow or web service for predictions
3. Flag PAINS and structural alerts for review
4. Rank compounds only after calibrating model outputs and project decision rules
5. Generate summary report

## Tips
- ADMETlab 3.0 reports 119 platform features: 77 prediction models, 34 computed properties, and 8 rules
- SwissADME has no public API and its terms restrict automated data retrieval; use its documented manual batches
- Chemprop ensembles require explicit prediction-time uncertainty estimation and separate calibration when calibrated outputs are needed
- DeepChem supports both PyTorch and TensorFlow (TF not deprecated)
- QED is a ranking descriptor, not a universal pass/fail definition of drug-likeness
- Interpret hERG potency relative to exposure, assay conditions, and safety margin; no single IC50 cutoff establishes safety
- PAINS matches flag possible assay interference for orthogonal testing; they are not categorical exclusions

## Related Skills
- molecular-descriptors - Calculate descriptors for ML
- substructure-search - Filter reactive groups
- virtual-screening - Screen after ADMET filtering
