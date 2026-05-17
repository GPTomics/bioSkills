# Molecular Standardization Usage Guide

## Overview

Standardize molecular structures into a single canonical form for ML training, deduplication, and database joining. Skipping standardization causes silent ML data leakage and database join misses; this skill applies industry-validated pipelines (ChEMBL `chembl_structure_pipeline`, RDKit `rdMolStandardize`, canSARchem) covering sanitization, salt stripping, neutralization, tautomer canonicalization, and stereo standardization.

## Prerequisites

```bash
pip install rdkit chembl_structure_pipeline pandas
```

RDKit version 2024.09+ is required (uses `rdMolStandardize` C++ implementation; older Python `MolStandardize` deprecated Q1 2024).

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
> "For these 100 SMILES, enumerate tautomers and pick the lowest-energy canonical form for each. Use TautomerEnumerator with default scoring rules and report the canonical tautomer per input."

### Database join cleanup
> "Two datasets, chembl_data.csv and zinc_data.csv, should be joined on standardized SMILES. Standardize both with the ChEMBL pipeline; compare InChIKey overlap; report duplicates."

## What the Agent Will Do

1. Read the input file (CSV / SDF / SMI) and parse SMILES with `Chem.MolFromSmiles`.
2. Apply the chosen standardization pipeline (ChEMBL default unless specified).
3. For each compound: sanitize, strip salts, neutralize charges, canonicalize tautomers.
4. Generate canonical SMILES and InChIKey for each standardized compound.
5. Deduplicate by InChIKey; aggregate activity if multiple records.
6. Output standardized CSV with `smiles`, `inchikey`, `activity`, `n_replicates` columns.
7. Report parse failures, fragments stripped, and tautomer changes.

## Tips

- Always standardize before fingerprinting, similarity searching, or ML training.
- Use InChIKey for cross-database identity (more robust than canonical SMILES).
- Quaternary ammoniums require `Uncharger(canonicalOrdering=True)` to preserve charge.
- Tautomer canonicalization is the most contentious step; document choice and apply consistently across train + test.
- For natural products / peptides, default tautomer rules may not apply; manual review.
- Standardize entire library before deduplication; tautomer differences cause duplicate misses.

## Related Skills

- chemoinformatics/molecular-io - Parse and write molecular files
- chemoinformatics/molecular-descriptors - Featurize standardized molecules
- chemoinformatics/similarity-searching - Compare standardized molecules
- chemoinformatics/qsar-modeling - QSAR training requires standardization upstream
- chemoinformatics/scaffold-analysis - Scaffold extraction after standardization
