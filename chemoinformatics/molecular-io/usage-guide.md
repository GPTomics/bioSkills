# Molecular I/O - Usage Guide

## Overview
Read, write, and convert molecular file formats including SMILES, SDF, MOL2, and PDB. Includes structure standardization pipeline for consistent molecular representations.

## Prerequisites
```bash
pip install rdkit
pip install openbabel-wheel  # For additional format support
```

## Quick Start
Tell your AI agent what you want to do:
- "Load my compound library from this SDF file"
- "Convert these SMILES to an SDF file with 3D coordinates"
- "Standardize the structures in my molecule database"
- "Read MOL2 files and convert to SMILES"

## Example Prompts

### Reading Molecules
> "Load all molecules from compounds.sdf and show me how many were successfully parsed."

> "Read SMILES from this CSV file where the SMILES column is named 'structure'."

### Writing Molecules
> "Save my filtered compounds to an SDF file with properties included."

> "Export canonical SMILES for my molecule list."

### Standardization
> "Standardize these molecules by removing salts, neutralizing charges, and canonicalizing tautomers."

> "Clean up my compound library for consistent representations."

## What the Agent Will Do
1. Parse molecular files using RDKit or Open Babel
2. Handle parsing errors gracefully
3. Apply standardization pipeline if requested
4. Convert between formats as needed
5. Preserve molecular properties during conversion

## Tips
- Use rdMolStandardize module (Python MolStandardize was removed in Q1 2024)
- For Open Babel 3.x, use `from openbabel import pybel` not `import pybel`
- Standardization order: Sanitize, Normalize, Neutralize, Canonicalize tautomer, Strip salts
- Use rdMolDraw2D for molecular drawing (legacy Draw.MolToImage is deprecated)
- Always check for None when loading molecules (invalid structures return None)

## Related Skills
- molecular-descriptors - Calculate properties after loading
- similarity-searching - Compare loaded molecules
- virtual-screening - Prepare ligands for docking
