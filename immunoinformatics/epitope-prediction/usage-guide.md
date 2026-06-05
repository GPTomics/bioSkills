# Epitope Prediction - Usage Guide

## Overview

Predict B-cell and T-cell epitopes for vaccine antigen design and epitope mapping, with the discipline to apply the right method and the right confidence to each: T-cell epitope prediction is mature (it reduces to MHC presentation), while B-cell epitope prediction is unreliable for native antigens unless a 3D structure is used.

## Prerequisites

```bash
pip install bepipred3   # linear B-cell; auto-downloads ESM-2 weights on first run
# DiscoTope-3.0 (conformational B-cell), NetMHCpan/NetMHCIIpan (T-cell), ElliPro, SEPPA,
# NetChop/NetCTLpan are standalone/web (IEDB or DTU Health Tech). AlphaFold (or a PDB)
# is needed to run DiscoTope-3.0 on antigens without an experimental structure.
```

## Quick Start

Tell your AI agent what you want to do:
- "Predict conformational B-cell epitopes on this antigen from its AlphaFold model"
- "Map T-cell (CD8) epitopes in this protein for common HLA alleles"
- "Find linear B-cell epitopes for a peptide-ELISA reagent"
- "Which epitopes are conserved across these viral strains?"

## Example Prompts

### B-Cell Epitopes

> "Fold this antigen and run DiscoTope-3.0, keeping calls only in high-pLDDT regions"

> "Predict linear B-cell epitopes with BepiPred-3.0 for a denatured-target ELISA, and caveat the conformational limitation"

### T-Cell Epitopes

> "Tile this antigen and score CD8 epitopes with EL-mode presentation, without stacking NetChop"

> "Predict CD4 epitopes against the patient's DR alleles"

### Conservation and Coverage

> "Check epitope conservancy across these strains and report HLA population coverage"

> "Which high-scoring epitopes fall in hypervariable loops and should be dropped?"

## What the Agent Will Do

1. Determine whether the request is a T-cell (mature) or B-cell (unreliable) problem and say so
2. For T-cell: tile the antigen and score with EL-mode MHC presentation; skip explicit NetChop by default
3. For B-cell with a structure: run DiscoTope-3.0 and gate calls by pLDDT
4. For B-cell sequence-only: use BepiPred-3.0 at its real default (0.1512) and flag that it misses most native epitopes
5. Treat classical propensity scales as obsolete decoration
6. Add conservation and HLA population-coverage analysis for vaccine selection, then recommend functional validation

## Tips

- **Two maturity levels** - trust T-cell predictions to prioritize peptides; treat B-cell predictions as hypotheses
- **Conformational dominance** - ~90% of B-cell epitopes are conformational; fold a structure and use DiscoTope-3.0
- **BepiPred threshold** - the default is 0.1512, not 0.5; or use top-X% mode
- **pLDDT gating** - DiscoTope is least reliable in low-confidence flexible loops (where antibodies often bind)
- **Skip NetChop** - EL-trained MHC models already encode proteasomal cleavage; explicit cleavage is usually redundant
- **Propensity scales** - Kolaskar/Parker/Emini are obsolete; do not report them as data
- **Immunodominance is unsolved** - presentation is necessary, not sufficient; validate by ELISpot/tetramer

## Related Skills

- immunoinformatics/mhc-binding-prediction - T-cell (CD8) epitope prediction reduces to class I presentation
- immunoinformatics/mhc-class-ii-prediction - T-cell (CD4) epitopes
- immunoinformatics/immunogenicity-scoring - ranking epitope candidates by likely T-cell response
- structural-biology/alphafold-predictions - fold an antigen to enable DiscoTope-3.0
