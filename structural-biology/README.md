# structural-biology

## Overview

Protein structure analysis using Biopython's Bio.PDB module. Covers reading/writing PDB and mmCIF files, navigating the SMCRA hierarchy (Structure-Model-Chain-Residue-Atom), geometric calculations, superimposition, and working with AlphaFold predictions.

**Tool type:** python | **Primary tools:** Bio.PDB, Bio.PDB.Superimposer

## Skills

| Skill | Description |
|-------|-------------|
| structure-io | Parse PDB/mmCIF/MMTF files, download from RCSB, write structures |
| structure-navigation | Navigate SMCRA hierarchy, extract sequences, handle disorder |
| geometric-analysis | Distances, angles, dihedrals, neighbor search, superimposition, RMSD, SASA |
| structure-modification | Transform coordinates, remove/add entities, modify properties |
| alphafold-predictions | Download and analyze AlphaFold predicted structures, pLDDT confidence |

## Example Prompts

- "Download PDB structure 4HHB"
- "Parse this mmCIF file and show the chains"
- "Convert this PDB to mmCIF format"
- "List all chains and their lengths"
- "Extract the protein sequence from chain A"
- "Find all ligands in this structure"
- "Measure the distance between CA atoms of residues 50 and 100"
- "Calculate the RMSD between these two structures"
- "Find all residues within 5 Angstroms of the ligand"
- "Superimpose these two structures"
- "Remove all water molecules"
- "Center the structure at the origin"
- "Set B-factors based on conservation scores"
- "Download the AlphaFold structure for this UniProt ID"
- "Analyze the pLDDT confidence scores"
- "Plot the predicted aligned error (PAE)"

## Requirements

```bash
pip install biopython numpy
```

## Related Skills

- **sequence-io** - Read sequences to compare with structure-derived sequences
- **sequence-manipulation** - Analyze extracted protein sequences
- **database-access** - Fetch structure metadata from NCBI/UniProt
- **alignment** - Sequence alignment for structure-based analysis
