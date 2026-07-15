# Molecular Standardization Usage Guide

## Overview

Standardize molecular structures consistently for ML training, deduplication, and database joining. The skill distinguishes ChEMBL representation/parent standardization from optional RDKit tautomer canonicalization and the canSARchem workflow, which canonicalizes tautomers before parent extraction.

## Prerequisites

```bash
pip install rdkit chembl_structure_pipeline pandas
```

Examples target RDKit 2024.09+ and the maintained `rdkit.Chem.MolStandardize.rdMolStandardize` API. The standalone MolVS package is a legacy dependency, not a prerequisite.

## Quick Start

Tell the AI agent what to do:
- "Standardize my SMILES library using the ChEMBL pipeline"
- "Strip salts and canonicalize tautomers for compound deduplication"
- "Build a standardized + deduplicated training set for QSAR from this assay data"
- "Apply canSARchem-style standardization (canonical tautomer before parent extraction)"
- "Standardize but preserve charged quaternary ammoniums"

## Example Prompts

### ChEMBL-style standardization
> "Apply the ChEMBL structure pipeline to compound_library.csv. Run standardize_mol then get_parent_mol on each SMILES. Output deduplicated, canonical SMILES + InChIKey + activity mean per duplicate group."

### Custom rdMolStandardize pipeline
> "Run a multi-step RDKit standardization: sanitize -> largest fragment -> normalize functional groups -> neutralize (preserve quaternary N) -> canonicalize tautomer -> remove isotopes. Compare outputs to ChEMBL pipeline."

### ML data preparation
> "Standardize hERG.csv (column smiles, label hERG_blocker). Deduplicate on InChIKey, average activity, remove inorganics, output train.csv with 1 row per unique compound."

### Tautomer-aware canonicalization
> "For these 100 SMILES, use `TautomerEnumerator` with its documented scoring rules and report the selected canonical tautomer per input; do not describe it as a physical lowest-energy calculation."

### Database join cleanup
> "Two datasets, chembl_data.csv and zinc_data.csv, should be joined on standardized SMILES. Standardize both with the ChEMBL pipeline; compare InChIKey overlap; report duplicates."

## What the Agent Will Do

1. Read the input file (CSV / SDF / SMI) and parse SMILES with `Chem.MolFromSmiles`.
2. Apply the chosen standardization pipeline (ChEMBL default unless specified).
3. Apply the stages defined by that pipeline; ChEMBL does not canonicalize tautomers, while the custom RDKit and canSARchem-style workflows shown here can.
4. Generate canonical SMILES and InChIKey for each standardized compound.
5. Deduplicate by InChIKey; aggregate activity if multiple records.
6. Output standardized CSV with `smiles`, `inchikey`, `activity`, `n_replicates` columns.
7. Report parse failures, ChEMBL exclusion flags, fragments stripped, and any tautomer changes when a tautomer step was requested. Do not silently admit `exclude_flag=True` structures into registration or QSAR output.

## Tips

- Always standardize before fingerprinting, similarity searching, or ML training.
- Use InChIKey for cross-database identity (more robust than canonical SMILES).
- Preserve permanent charges with the selected uncharging policy (normally `force=False`) and inspect charge-sensitive structures; `canonicalOrder=True` controls site ordering, not whether a charge is chemically permanent.
- Tautomer canonicalization is the most contentious step; document choice and apply consistently across train + test.
- For natural products / peptides, default tautomer rules may not apply; manual review.
- Standardize entire library before deduplication; tautomer differences cause duplicate misses.

## Related Skills

- chemoinformatics/molecular-io - Parse and write molecular files
- chemoinformatics/molecular-descriptors - Featurize standardized molecules
- chemoinformatics/similarity-searching - Compare standardized molecules
- chemoinformatics/qsar-modeling - QSAR training requires standardization upstream
- chemoinformatics/scaffold-analysis - Scaffold extraction after standardization
