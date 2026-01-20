# Structure Navigation - Usage Guide

## Overview

This skill covers navigating the SMCRA (Structure-Model-Chain-Residue-Atom) hierarchy in Biopython Bio.PDB. Use it to access specific parts of structures, iterate over components, extract sequences, and handle disordered atoms.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "List all chains in this PDB file"
- "Extract the sequence from chain A"
- "Find all ligands in this structure"

## Example Prompts

### Accessing Structure Parts
> "Get chain A from model 0"

> "Show me all atoms in residue 100"

> "What is the B-factor of the CA atom in residue 50?"

### Iterating
> "Print all residue numbers in chain B"

> "List all C-alpha coordinates"

> "Count atoms in each chain"

### Extracting Information
> "Extract the protein sequence from this structure"

> "Find all ARG residues"

> "List all hetero groups (ligands)"

### Handling Disorder
> "Show alternative conformations for residue 42"

> "List all disordered atoms in chain A"

## What the Agent Will Do

1. Parse the structure file
2. Navigate to requested level using SMCRA hierarchy
3. Extract requested information (coordinates, properties, sequences)
4. Handle disordered atoms/residues appropriately
5. Return organized results

## SMCRA Hierarchy

```
Structure (whole PDB entry)
    Model (NMR conformer or asymmetric unit)
        Chain (polypeptide or nucleic acid)
            Residue (amino acid, nucleotide, ligand, water)
                Atom (individual atom with coordinates)
```

## Tips

- **Use `get_` methods** for iteration: `get_chains()`, `get_residues()`, `get_atoms()`
- **Residue ID is a tuple** - `(hetfield, resseq, icode)`, not just a number
- **hetfield values**: `' '` = amino acid, `'W'` = water, `'H_xxx'` = hetero
- **PPBuilder** extracts polypeptide sequences from coordinates
- **CaPPBuilder** works better for structures with missing atoms
- **Selection.unfold_entities** extracts all entities at a level (R=residue, A=atom)
