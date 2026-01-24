# Modern Structure Prediction

## Overview

Predict protein structures using state-of-the-art ML models including ESMFold, AlphaFold3, Chai-1, and Boltz-1. These tools enable rapid structure prediction for single proteins and complexes, including protein-ligand interactions.

## Prerequisites

```bash
# ESMFold (local)
pip install fair-esm

# Chai-1
pip install chai-lab

# Boltz-1
pip install boltz

# ColabFold
pip install colabfold

# For API-only usage, only requests is needed
pip install requests
```

## Quick Start

Tell your AI agent what you want to do:
- "Predict the structure of this protein sequence using ESMFold"
- "Run AlphaFold3 on my protein complex"
- "Compare predictions from multiple structure prediction methods"
- "Predict a protein-ligand complex with Chai-1"

## Example Prompts

### Single Protein Prediction
> "Predict the structure of this protein sequence using the fastest method"

> "Run ESMFold on my sequence and analyze the confidence scores"

> "Get an AlphaFold3 prediction for this protein"

### Complex Prediction
> "Predict the structure of this protein-protein complex"

> "Run Boltz-1 on my heterodimer sequences"

> "Model my protein binding to this small molecule ligand"

### Comparison and Validation
> "Compare ESMFold, AlphaFold3, and Chai-1 predictions for this sequence"

> "Calculate RMSD between different structure predictions"

> "Which regions have high confidence across all prediction methods?"

### Batch Processing
> "Predict structures for all sequences in this FASTA file"

> "Run ESMFold on my list of protein sequences"

## What the Agent Will Do

1. Determine the best prediction method based on your needs (speed, accuracy, complex support)
2. Prepare the input sequence(s) in the appropriate format
3. Run the prediction via API or local installation
4. Extract and interpret confidence metrics (pLDDT, pTM, PAE)
5. Save the predicted structure in PDB or mmCIF format
6. Provide guidance on interpreting low-confidence regions

## Tips

- ESMFold is fastest for single chains and doesn't require MSA
- AlphaFold3 server is best for complexes but has usage limits
- Chai-1 and Boltz-1 are open-source alternatives for complex prediction
- pLDDT > 70 indicates confident predictions; < 50 suggests disorder
- For complexes, check ipTM (interface pTM) to assess binding prediction quality
- Compare multiple methods for critical applications
- Low-confidence regions may indicate intrinsically disordered regions
