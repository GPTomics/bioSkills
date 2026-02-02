# Epitope Prediction - Usage Guide

## Overview

Predict B-cell and T-cell epitopes from protein sequences using IEDB tools for vaccine design and antibody development.

## Prerequisites

```bash
pip install requests pandas
# IEDB tools accessed via API (no local install needed)
```

## Quick Start

Tell your AI agent what you want to do:
- "Identify B-cell epitopes in this spike protein"
- "Find immunogenic regions in my antigen"
- "Predict both B-cell and T-cell epitopes for vaccine design"

## Example Prompts

### B-Cell Epitopes

> "Predict B-cell epitopes using BepiPred-2.0"

> "Find antibody binding sites in this sequence"

### T-Cell Epitopes

> "Identify MHC-I epitopes for HLA-A*02:01"

> "Find helper T-cell epitopes (MHC-II)"

### Vaccine Design

> "Map all immunogenic regions in this antigen"

> "Combine predictions from multiple methods"

## What the Agent Will Do

1. Submit sequence to IEDB prediction tools
2. Parse prediction scores for each residue
3. Identify continuous epitope regions
4. Filter by score threshold
5. Return epitope locations and sequences

## Tips

- **Threshold** - BepiPred >0.5 default; increase for stringency
- **Epitope length** - B-cell: 5-15aa; T-cell: 8-11aa (MHC-I)
- **Consensus** - Multiple methods increase confidence
- **Conformational** - Need structure for ~90% of B-cell epitopes
- **Validation** - Experimental validation always recommended
