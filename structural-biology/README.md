# Structural Biology

Protein structure analysis using Biopython's Bio.PDB module for parsing, navigating, analyzing, and modifying 3D macromolecular structures.

## Overview

This category covers protein structure operations: reading/writing PDB and mmCIF files, navigating the SMCRA hierarchy (Structure-Model-Chain-Residue-Atom), performing geometric calculations (distances, angles, RMSD, superimposition), and modifying structures programmatically.

**Tool type:** `python`
**Primary tools:** Bio.PDB, Bio.PDB.NeighborSearch, Bio.PDB.Superimposer

## Skills

| Skill | Description |
|-------|-------------|
| [structure-io](structure-io/) | Parse PDB/mmCIF/MMTF files, download from RCSB, write structures |
| [structure-navigation](structure-navigation/) | Navigate SMCRA hierarchy, extract sequences, handle disorder |
| [geometric-analysis](geometric-analysis/) | Distances, angles, dihedrals, neighbor search, superimposition, RMSD, SASA |
| [structure-modification](structure-modification/) | Transform coordinates, remove/add entities, modify properties |

## Workflow

```
RCSB PDB / Local File
    |
    v
[structure-io] -----> Parse PDB/mmCIF, download structures
    |
    v
[structure-navigation] --> Access chains, residues, atoms
    |                      Extract sequences, find ligands
    |
    +---> [geometric-analysis]
    |         |
    |         +--> Measure distances/angles
    |         +--> Find contacts (NeighborSearch)
    |         +--> Superimpose structures, calc RMSD
    |
    +---> [structure-modification]
              |
              +--> Remove water/hydrogens
              +--> Transform coordinates
              +--> Modify B-factors
              +--> Build new structures
              |
              v
          [structure-io] --> Write modified structure
```

## SMCRA Hierarchy

Bio.PDB organizes structures in a five-level hierarchy:

```
Structure (PDB entry)
    |
    +-- Model (0, 1, ...)      # NMR ensemble or crystal asymmetric unit
          |
          +-- Chain (A, B, ...) # Polypeptide chains
                |
                +-- Residue     # Amino acids, ligands, water
                      |
                      +-- Atom  # Individual atoms with coordinates
```

## Example Prompts

### Structure I/O
- "Download PDB structure 4HHB"
- "Parse this mmCIF file and show the chains"
- "Convert this PDB to mmCIF format"

### Structure Navigation
- "List all chains and their lengths"
- "Extract the protein sequence from chain A"
- "Find all ligands in this structure"

### Geometric Analysis
- "Measure the distance between CA atoms of residues 50 and 100"
- "Calculate the RMSD between these two structures"
- "Find all residues within 5 Angstroms of the ligand"

### Structure Modification
- "Remove all water molecules"
- "Center the structure at the origin"
- "Set B-factors based on conservation scores"

## Requirements

```bash
pip install biopython numpy
```

## Notes

- **mmCIF preferred** - PDB format is frozen; mmCIF is the modern standard for full metadata
- **MMCIF2Dict for metadata** - Use for accessing any mmCIF field not in structure object
- **NeighborSearch for contacts** - KD-tree based fast neighbor finding
- **Superimposer for alignment** - Returns RMSD and applies transformation
- **CEAligner for dissimilar structures** - When sequences differ significantly
- **SASA built-in** - Use `Bio.PDB.SASA.ShrakeRupley` for solvent accessible surface area (no external dependency)
- **DSSP requires external program** - Secondary structure assignment via `Bio.PDB.DSSP` requires DSSP installation; not covered in these skills

## Deprecations Avoided

- **Bio.PDB.Vector module** - Import `calc_dihedral`, `calc_angle` from `Bio.PDB` directly
- **Bio.PDB.Polypeptide.three_to_one** - Use `Bio.Data.PDBData.protein_letters_3to1` instead
- **Bio.PDB.Residue.get_atom()** - Removed; use `residue[atom_name]` or `residue.has_id()`
- **Bio.PDB.QCPSuperimposer** - Use standard `Superimposer` or `Bio.PDB.qcprot`

## Related Skills

- **sequence-io** - Read sequences to compare with structure-derived sequences
- **sequence-manipulation** - Analyze extracted protein sequences
- **database-access** - Fetch structure metadata from NCBI/UniProt
- **alignment** - Sequence alignment for structure-based analysis

## References

- [Bio.PDB documentation](https://biopython.org/docs/latest/api/Bio.PDB.html)
- [Bio.PDB Tutorial](https://biopython.org/docs/dev/Tutorial/chapter_pdb.html)
- [Structural Bioinformatics FAQ](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ)
