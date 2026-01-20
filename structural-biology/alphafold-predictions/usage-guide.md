# AlphaFold Predictions - Usage Guide

## Overview

Download and analyze AI-predicted protein structures from the AlphaFold database, including confidence scores (pLDDT) and predicted aligned error (PAE) for assessing prediction quality.

## Prerequisites

```bash
pip install biopython requests numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Download the AlphaFold structure for UniProt P04637"
- "Check the confidence scores for this AlphaFold prediction"
- "Identify low-confidence regions in this predicted structure"

## Example Prompts

### Downloading Structures
> "Download the AlphaFold model for UniProt ID P53_HUMAN"

> "Get the AlphaFold structure for P04637 in mmCIF format"

> "Fetch AlphaFold predictions for these 5 UniProt IDs"

### Analyzing Confidence
> "Show me the pLDDT scores for this AlphaFold structure"

> "Which regions have low confidence (pLDDT < 70)?"

> "Plot the per-residue confidence scores"

### PAE Analysis
> "Download and visualize the PAE matrix for this protein"

> "Identify domain boundaries from the PAE data"

> "Which domain-domain interactions are reliable?"

### Quality Assessment
> "Compare this AlphaFold prediction with the experimental structure"

> "Highlight the confident vs disordered regions"

> "Is this prediction reliable enough for docking?"

## What the Agent Will Do

1. Construct the AlphaFold database URL from UniProt ID
2. Download the structure file (PDB or mmCIF format)
3. Parse confidence scores from B-factor column (pLDDT)
4. Optionally fetch PAE matrix for inter-residue error estimates
5. Analyze and report on prediction quality

## Database Access

**Direct download URLs:**
- Structure: `https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-model_v4.pdb`
- mmCIF: `https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-model_v4.cif`
- PAE: `https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-predicted_aligned_error_v4.json`

**API endpoint:** `https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT_ID}`

## Confidence Score Interpretation

| pLDDT Range | Interpretation | Color (AlphaFold DB) |
|-------------|----------------|----------------------|
| 90-100 | Very high confidence | Blue |
| 70-90 | Confident | Cyan |
| 50-70 | Low confidence | Yellow |
| <50 | Very low (likely disordered) | Orange |

## Tips

- **Check pLDDT before using** - Regions with pLDDT < 70 should be treated cautiously
- **Use PAE for domains** - Low PAE between residues indicates reliable relative positioning
- **B-factor column holds pLDDT** - Parse with Bio.PDB and read `atom.bfactor`
- **Single conformation only** - AlphaFold predicts one state, not conformational ensembles
- **Compare when possible** - Validate against experimental structures if available
- **200M+ structures available** - Coverage includes most UniProt sequences
