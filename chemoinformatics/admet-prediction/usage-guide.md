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

> "Filter for compounds passing both Lipinski and Veber rules."

### Safety Filtering
> "Check my hits for PAINS and other structural alerts."

> "Identify compounds with potential hERG liability."

## What the Agent Will Do
1. Calculate drug-likeness properties (Lipinski, QED)
2. Call ADMETlab 3.0 API for predictions
3. Filter for PAINS and structural alerts
4. Rank compounds by safety profile
5. Generate summary report

## Tips
- ADMETlab 3.0 provides 119 endpoints (use this, not ADMETlab 2.0)
- SwissADME has NO API - it is web-only, do not try programmatic access
- DeepChem supports both PyTorch and TensorFlow (TF not deprecated)
- QED > 0.5 is generally drug-like
- hERG IC50 > 10 μM is typically considered safe
- PAINS filter removes promiscuous compounds that cause assay interference

## Related Skills
- molecular-descriptors - Calculate descriptors for ML
- substructure-search - Filter reactive groups
- virtual-screening - Screen after ADMET filtering
