# AlphaFold Predictions Usage Guide

The AlphaFold Protein Structure Database provides AI-predicted structures for most known proteins.

## Database Coverage

- **200+ million** structures predicted
- Covers proteomes of major model organisms
- Predictions based on UniProt sequences

## Access Methods

### Direct Download URLs

```
https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-model_v4.pdb
https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-model_v4.cif
https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-predicted_aligned_error_v4.json
```

### API Endpoint

```
https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT_ID}
```

## Quick Start

```python
import requests

# Download structure
uniprot_id = 'P04637'
url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
response = requests.get(url)
with open(f'AF-{uniprot_id}.pdb', 'w') as f:
    f.write(response.text)
```

## Understanding Confidence Scores

### pLDDT (per-residue)
- Stored in B-factor column
- 0-100 scale
- Higher = more confident

### PAE (inter-residue)
- Matrix showing expected position error
- Lower = better
- Identifies domain boundaries

## Best Practices

1. **Check pLDDT scores** before using predictions
2. **Use PAE** to identify reliable domain-domain contacts
3. **Compare with experimental** structures when available
4. **Be cautious** with low-confidence regions (<70 pLDDT)
5. **Consider ensemble** - AlphaFold predicts one conformation

## Limitations

- Single static conformation predicted
- May miss alternative conformations
- Disordered regions have low confidence
- Oligomeric states not always correct
- Ligand binding may be absent

## When to Use AlphaFold

**Good for:**
- Proteins without experimental structures
- Homology modeling starting points
- Domain architecture analysis
- Identifying folded vs disordered regions

**Caution needed for:**
- Active site details
- Protein-protein interfaces
- Conformational changes
- Membrane protein orientations
