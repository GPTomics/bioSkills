# Molecular I/O - Usage Guide

## Overview
Read, write, and convert molecular file formats including SMILES, SDF, MOL2, and PDB. Inspect representation-sensitive details such as salts, charges, and stereochemistry during conversion, then use the molecular-standardization skill when a registration policy is required.

## Prerequisites
```bash
pip install rdkit
pip install openbabel-wheel  # For additional format support
```

## Quick Start
Tell your AI agent what you want to do:
- "Load my compound library from this SDF file"
- "Convert these SMILES to an SDF file with 3D coordinates"
- "Inspect my molecule database for salts, mixtures, charge states, and missing stereochemistry before conversion"
- "Read MOL2 files and convert to SMILES"

## Example Prompts

### Reading Molecules
> "Load all molecules from compounds.sdf and show me how many were successfully parsed."

> "Read SMILES from this CSV file where the SMILES column is named 'structure'."

### Writing Molecules
> "Save my filtered compounds to an SDF file with properties included."

> "Export canonical SMILES for my molecule list."

### Checks Before Standardization
> "Parse these molecules and report salts, mixtures, charge states, and stereochemistry issues that need a documented standardization policy."

> "Convert this compound library without silently changing charge, stereochemistry, or tautomer state, then identify records that need the molecular-standardization workflow."

## What the Agent Will Do
1. Parse molecular files using RDKit or Open Babel
2. Handle parsing errors gracefully
3. Preserve representation-sensitive details during conversion and route registration-grade standardization to the dedicated skill
4. Convert between formats as needed
5. Preserve molecular properties during conversion

## Tips
- Use the maintained `rdkit.Chem.MolStandardize.rdMolStandardize` module for custom standardization; standalone MolVS is a legacy package
- For Open Babel 3.x, use `from openbabel import pybel` not `import pybel`
- Standardization order is policy-specific; use the molecular-standardization skill and document whether tautomer canonicalization occurs before or after parent selection
- Use rdMolDraw2D when direct control over molecular drawing is needed; `Draw.MolToImage` remains a supported convenience API
- Always check for None when loading molecules (invalid structures return None)

## Related Skills
- molecular-descriptors - Calculate properties after loading
- similarity-searching - Compare loaded molecules
- virtual-screening - Prepare ligands for docking
