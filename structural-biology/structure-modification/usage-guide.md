# Structure Modification - Usage Guide

## Overview

This skill modifies protein structures in place with Biopython Bio.PDB: transforming coordinates, stripping waters and heteroatoms, overloading the B-factor column for coloring, renumbering, and building or copying entities. Every edit mutates the parsed object directly, so the two decisions that govern correctness are which rotation convention a transform matrix uses (row-convention `Entity.transform`/`Superimposer` versus column-convention biological-assembly operators) and what to preserve before a destructive edit. Getting the convention or the HETFLAG filter wrong corrupts geometry or deletes functional atoms silently, with no error. For mmCIF-fidelity edits (entities, anisotropic B-factors, assembly generation) Bio.PDB is the wrong tool - reach for gemmi.

## Prerequisites

```bash
pip install biopython numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Remove the waters but keep the zinc and heme"
- "Apply this REMARK 350 assembly operator to build the dimer"
- "Color by conservation in the B-factor column but keep the originals"
- "Renumber chain A starting from 1"
- "Merge these two PDB files without chain-id clashes"

## Example Prompts

### Transformations
> "Center this structure at the origin"

> "Rotate the structure 90 degrees around the Z axis about its center of mass"

> "Apply this rotation-plus-translation matrix from the assembly record"

### Removing Entities
> "Strip crystallographic waters but keep all metals and cofactors"

> "Remove hydrogens from this structure"

> "Extract chain A residues 50-100 into a new file"

### Modifying Properties
> "Write these pLDDT values into the B-factor column for PyMOL coloring"

> "Set occupancy to 1.0 for chain A"

### Structure Building
> "Renumber residues starting from 1, preserving insertion codes"

> "Duplicate chain A as chain B"

> "Merge these two PDB files"

## What the Agent Will Do

1. Parse the input structure(s), taking a deep copy first if the original must be preserved
2. Identify the transform convention (row for Bio.PDB/Superimposer, column for assembly operators) and transpose if needed
3. Apply the requested edit - transform, HETFLAG-filtered strip, B-factor overload with an originals snapshot, renumber, or build
4. Wire SMCRA parent-child links and check for id collisions when building, copying, or merging
5. Write the modified structure to a new file with PDBIO

## Key Operations

| Operation | Method |
|-----------|--------|
| Transform (row convention) | `Entity.transform(rot, tran)` |
| Assembly operator (column convention) | `Entity.transform(R.T, t)` |
| Remove atom / residue / chain | `parent.detach_child(id)` |
| Filtered write | `PDBIO.save(path, Select())` |
| Overload B-factor | Set `atom.bfactor` after snapshotting originals |
| Renumber residue | Set `residue.id = (hetflag, seq, icode)` |
| Build tree | `StructureBuilder` |
| Copy subtree | `copy.deepcopy(entity)` |

## Tips

- Every edit mutates in place; `copy.deepcopy` the structure before modifying if the original is still needed
- Bio.PDB `Superimposer`/`Entity.transform` use the row convention (`coords @ rot`); REMARK 350 and `_pdbx_struct_oper_list` operators are column-convention (`R @ x + t`) and must be transposed before `Entity.transform`
- Prefer `Entity.transform`/`atom.transform` over hand-rolled `np.dot` so the convention is handled automatically
- Overloading the B-factor column destroys the real temperature factors, and for AlphaFold models it overwrites the pLDDT stored there - snapshot originals and never send the file to refinement
- Strip solvent on the HETFLAG (`r.id[0] == 'W'`), not the residue name; a blanket `r.id[0] != ' '` deletes catalytic metals, cofactors, and modified residues like MSE
- Collect child ids into a list before `detach_child` to avoid mutating the dict mid-iteration
- Renumbering must carry the full `(hetflag, resseq, icode)` tuple; renumber through a temporary offset to avoid mid-loop id collisions
- For mmCIF-fidelity edits (anisotropic B-factors, entities, assembly generation) use gemmi - Bio.PDB does not round-trip them

## Related Skills

- structure-io - Parse and write structure files; mmCIF vs PDB format ceilings
- structure-navigation - Walk chains/residues/atoms and the SMCRA id tuple
- geometric-analysis - Superimpose structures and read back the row-convention rotation
- interface-analysis - Analyze interfaces after generating the biological assembly
- sequence-manipulation/seq-objects - Generate sequences from modified structures
