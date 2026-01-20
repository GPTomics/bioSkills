# Structure Modification - Usage Guide

## Overview

This skill covers modifying protein structures: transforming coordinates, removing/adding atoms and residues, modifying B-factors and occupancies, and building structures programmatically.

## Prerequisites

```bash
pip install biopython numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Remove all water molecules from this structure"
- "Translate the structure by 10 Angstroms in X"
- "Set B-factors based on this conservation score data"

## Example Prompts

### Transformations
> "Center this structure at the origin"

> "Rotate the structure 90 degrees around the Z axis"

> "Translate chain A by [10, 5, 0]"

### Removing Entities
> "Remove hydrogens from this structure"

> "Delete chain B"

> "Remove all hetero atoms (keep only protein)"

### Modifying Properties
> "Set all B-factors to 20"

> "Color by conservation score in the B-factor column"

> "Set occupancy to 0.5 for chain A"

### Structure Building
> "Renumber residues starting from 1"

> "Rename chain A to X"

> "Merge these two PDB files"

## What the Agent Will Do

1. Parse input structure(s)
2. Navigate to target entities
3. Apply requested modifications
4. Save modified structure to new file

## Key Operations

| Operation | Method |
|-----------|--------|
| Remove atom | `residue.detach_child(atom_id)` |
| Remove residue | `chain.detach_child(residue_id)` |
| Remove chain | `model.detach_child(chain_id)` |
| Transform | Modify `atom.coord` directly |
| Modify B-factor | Set `atom.bfactor` |
| Modify occupancy | Set `atom.occupancy` |
| Rename chain | Set `chain.id` |
| Renumber residue | Set `residue.id` |

## Tips

- **Use `detach_child()`** to remove entities, not `del`
- **Modify coord directly** - `atom.coord` is a numpy array
- **B-factors for visualization** - Use B-factor column to encode custom data
- **Chain IDs are single characters** - typically A-Z
- **Residue ID is a tuple** - `(hetfield, resnum, icode)`
- **Copy before modifying** if you need to preserve original
- **StructureBuilder** for building from scratch
