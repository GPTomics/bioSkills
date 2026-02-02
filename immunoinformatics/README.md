# immunoinformatics

## Overview

Computational immunology including MHC binding prediction, neoantigen identification, epitope prediction, and TCR-antigen matching for vaccine design and cancer immunotherapy.

**Tool type:** mixed | **Primary tools:** mhcflurry, pVACtools, BepiPred, IEDB

## Skills

| Skill | Description |
|-------|-------------|
| mhc-binding-prediction | Peptide-MHC binding with MHCflurry, NetMHCpan |
| neoantigen-prediction | Tumor neoantigens with pVACtools |
| epitope-prediction | B/T-cell epitopes with BepiPred, IEDB |
| immunogenicity-scoring | Neoantigen prioritization and ranking |
| tcr-epitope-binding | TCR-epitope specificity with ERGO-II |

## Example Prompts

- "Predict MHC binding for these peptides with HLA-A*02:01"
- "Find neoantigens from my somatic VCF"
- "Identify B-cell epitopes in this spike protein"
- "Rank these neoantigens by immunogenicity"
- "Predict what antigens this TCR sequence recognizes"
- "Design a personalized cancer vaccine from this mutation list"

## Requirements

```bash
pip install mhcflurry
mhcflurry-downloads fetch  # Download models

# pVACtools (complex dependencies - use conda)
conda create -n pvactools python=3.8
conda activate pvactools
pip install pvactools
```

## Related Skills

- **clinical-databases** - HLA typing, variant annotation
- **variant-calling** - Somatic mutation calling
- **tcr-bcr-analysis** - TCR repertoire analysis
