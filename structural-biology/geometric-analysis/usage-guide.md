# Geometric Analysis - Usage Guide

## Overview

This skill covers geometric calculations on protein structures: measuring distances, angles, and dihedrals; finding neighbor atoms; superimposing structures; calculating RMSD; and computing solvent accessible surface area (SASA).

## Prerequisites

```bash
pip install biopython numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Measure the distance between residue 50 CA and residue 100 CA"
- "Calculate the RMSD between these two structures"
- "Find all residues within 5 Angstroms of the ligand"

## Example Prompts

### Distance Measurements
> "What is the distance between these two atoms?"

> "Create a CA-CA distance matrix for chain A"

> "Find the closest residue to the ligand"

### Angles
> "Calculate the phi/psi angles for all residues"

> "What is the N-CA-C bond angle in residue 42?"

> "Calculate chi1 angles for all aromatic residues"

### Superimposition
> "Superimpose these two structures and calculate RMSD"

> "Align the mobile structure onto the reference"

> "What is the RMSD between these NMR models?"

### Neighbor Searches
> "Find all atoms within 4A of the active site"

> "List residue contacts within the protein"

> "Which residues are in contact with chain B?"

### Surface Analysis
> "Calculate the solvent accessible surface area"

> "Which residues are buried vs exposed?"

> "What is the SASA of chain A?"

## What the Agent Will Do

1. Parse structure files
2. Extract requested atoms/coordinates
3. Perform geometric calculations
4. For superimposition: find optimal rotation/translation
5. Return measurements with appropriate units

## Key Functions

| Function | Purpose |
|----------|---------|
| `atom1 - atom2` | Distance between atoms |
| `calc_angle()` | Angle between 3 atoms |
| `calc_dihedral()` | Dihedral angle (4 atoms) |
| `NeighborSearch` | Find atoms within radius |
| `Superimposer` | Structural alignment + RMSD |
| `CEAligner` | Align dissimilar structures |
| `ShrakeRupley` | Solvent accessible surface area |

## Tips

- **Atom subtraction** returns distance directly: `distance = atom1 - atom2`
- **Use `get_vector()`** for angle/dihedral calculations, not `coord`
- **NeighborSearch levels**: A=atom, R=residue, C=chain for contact finding
- **Superimposer** requires equal-length atom lists
- **CEAligner** works for structures with different sequences
- **RMSD units** are always Angstroms
