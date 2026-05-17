# Scaffold Analysis Usage Guide

## Overview

Analyze chemical libraries by underlying scaffolds (Bemis-Murcko, generic frameworks), cluster compounds by chemotype, perform R-group decomposition, run matched molecular pair analysis (mmpdb), and produce scaffold-balanced QSAR train/test splits to prevent data leakage.

## Prerequisites

```bash
pip install rdkit mmpdb scikit-learn datamol
```

## Quick Start

Tell the AI agent what to do:
- "Compute Bemis-Murcko scaffolds for my compound library and report cluster sizes"
- "Generate a scaffold-balanced 80/20 train/test split for QSAR"
- "Run MMPA on the bioassay data to find activity-improving transformations"
- "Decompose this 50-compound series into R-groups around the central scaffold"
- "Identify analog series of size >= 5 in this library"

## Example Prompts

### Library chemotype analysis
> "Compute Bemis-Murcko scaffolds for all 10k compounds in library.csv. Output a table of scaffold SMILES with member count. Report the top 20 most-populated scaffolds."

### Scaffold-balanced ML split
> "Apply Bemis-Murcko scaffold split (80% train, 20% test) to qsar_data.csv. Ensure no scaffold appears in both train and test. Output train.csv and test.csv."

### MMPA mining
> "Run mmpdb on hERG_data.csv (smiles + pIC50 columns). Report top 10 transformations with N >= 10 pairs and |delta(pIC50)| >= 0.5. Format as table."

### R-group decomposition
> "Given scaffold 'c1ccc(-[*:1])cc1-[*:2]' and 50 analogs, decompose into R-group table. Output CSV with columns R1, R2, activity for downstream Free-Wilson analysis."

## What the Agent Will Do

1. Parse SMILES library; standardize where needed.
2. Compute Bemis-Murcko scaffolds (and generic framework if requested).
3. Group compounds by scaffold; report chemotype distribution.
4. For ML splits: assign whole scaffolds to train or test; output partitioned data.
5. For MMPA: fragment compounds, index by mmpdb, run transform queries.
6. For R-decomposition: align compounds to scaffold, extract R-group SMILES.

## Tips

- Bemis-Murcko gives **empty scaffold** for linear molecules; supplement with linear features.
- Generic framework loses heteroatom info; use only for topology comparison.
- mmpdb needs >= 1000 compounds for reliable transformation statistics.
- Always use scaffold split for QSAR; random splits inflate metrics.
- R-decomposition fails if compounds match multiple scaffold templates; use FMCS first.

## Related Skills

- chemoinformatics/molecular-io - Parse compounds
- chemoinformatics/molecular-standardization - Standardize before scaffold extraction
- chemoinformatics/reaction-enumeration - Free-Wilson on R-decomp output
- chemoinformatics/similarity-searching - 2D scaffold-hopping
- chemoinformatics/qsar-modeling - Mandatory scaffold split for QSAR
- chemoinformatics/generative-design - Scaffold-decoration generation
