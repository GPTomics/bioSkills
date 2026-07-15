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

> "Calculate QED scores and apply QED > 0.5 as a documented project triage heuristic; show the unfiltered distribution too."

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
- ECFP4 = radius 2 and ECFP6 = radius 3 (the ECFP number denotes diameter, `2 * radius`)
- Include `useChirality=True` for stereo-sensitive fingerprints
- QED > 0.5 can be used as a repository/project triage heuristic, but calibrate it for the library and do not treat it as a universal drug-likeness boundary
- RDKit defaults to ETKDG; select `AllChem.ETKDGv3()` explicitly when the v3 small-ring and macrocycle improvements are intended
- Lipinski thresholds: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10
- 3D descriptors require conformer generation first

## Related Skills
- molecular-io - Load molecules for descriptor calculation
- similarity-searching - Use fingerprints for similarity
- admet-prediction - Predict ADMET from descriptors
- machine-learning/biomarker-discovery - ML on molecular features
