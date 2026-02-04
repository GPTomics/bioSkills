# Molecular Descriptors - Usage Guide

## Overview
Calculate molecular fingerprints and physicochemical properties for compound characterization, drug-likeness assessment, and machine learning featurization.

## Prerequisites
```bash
pip install rdkit
pip install numpy pandas
```

## Quick Start
Tell your AI agent what you want to do:
- "Calculate ECFP4 fingerprints for my compound library"
- "Check Lipinski rule of 5 compliance for these molecules"
- "Calculate QED drug-likeness scores"
- "Generate descriptors for machine learning"

## Example Prompts

### Fingerprints
> "Generate Morgan fingerprints with radius 2 and 2048 bits for my molecules."

> "Calculate MACCS keys for similarity searching."

### Drug-Likeness
> "Check which compounds pass Lipinski's rule of 5."

> "Calculate QED scores and filter for drug-like compounds (QED > 0.5)."

### Full Descriptor Set
> "Calculate all available RDKit descriptors for my molecules."

> "Generate 3D conformers and calculate shape descriptors."

## What the Agent Will Do
1. Load molecules from provided structures
2. Calculate requested descriptors/fingerprints
3. Compile results into a DataFrame
4. Apply filters based on thresholds if requested
5. Export results in requested format

## Tips
- ECFP4 = radius 2, ECFP6 = radius 3 (diameter = 2 * radius + 2)
- Include `useChirality=True` for stereo-sensitive fingerprints
- QED > 0.5 is generally considered drug-like
- ETKDGv3 is now the default conformer generator in RDKit
- Lipinski thresholds: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10
- 3D descriptors require conformer generation first

## Related Skills
- molecular-io - Load molecules for descriptor calculation
- similarity-searching - Use fingerprints for similarity
- admet-prediction - Predict ADMET from descriptors
- machine-learning/biomarker-discovery - ML on molecular features
