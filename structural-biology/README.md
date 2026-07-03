# structural-biology

## Overview

Protein structure analysis using Biopython's Bio.PDB module and modern prediction tools. A coordinate file is an interpreted model, not ground truth: reliability varies per-atom, mmCIF (PDBx) is the canonical format while legacy PDB is a frozen lossy container, and a prediction's pLDDT/PAE report per-residue confidence, not correctness. These skills cover reading/writing structures, navigating the SMCRA hierarchy, geometry and comparison, interface analysis, validation, and predicted-structure interpretation with those principles built in.

**Tool type:** python | **Primary tools:** Bio.PDB, ESMFold, AlphaFold DB

## Skills

| Skill | Description |
|-------|-------------|
| structure-io | Parse and write PDB/mmCIF/BinaryCIF, download from RCSB, auth-vs-label numbering, asymmetric-unit-vs-biological-assembly (MMTF retired) |
| structure-navigation | Navigate the SMCRA hierarchy, handle altlocs/insertion codes/disorder, extract observed-vs-SEQRES sequences, select NMR models |
| geometric-analysis | Distances, angles, dihedrals, superimposition, RMSD-vs-TM-score-vs-lDDT metric choice, SASA with the probe-radius caveat, radius of gyration |
| structure-modification | Transform coordinates (row-vs-column rotation convention), strip waters/hetero by HETFLAG, edit B-factors/occupancies, build structures |
| interface-analysis | Protein-protein and protein-ligand contacts, buried surface area, biological interface vs crystal-packing artifact, epitope/binding-site residues |
| structure-validation | Judge whether a structure/region is trustworthy: resolution/R-free/B-factor, MolProbity clash/Ramachandran/rotamer outliers, validate predicted models before docking |
| structure-preparation | Make a structure simulation-ready: add hydrogens, assign protonation/tautomer/flip states at a stated pH, fill missing atoms and short loops (PDBFixer, reduce, PROPKA/PDB2PQR) |
| binding-site-detection | Detect and rank cavities/pockets on apo structures, geometric vs ML ligandability, cryptic-pocket and predicted-model limits (fpocket, P2Rank, CASTp) |
| alphafold-predictions | Retrieve AlphaFold DB models and read pLDDT as per-residue confidence (low = disorder), PAE for inter-domain placement, downstream-use limits |
| modern-structure-prediction | Predict structures/complexes with ESMFold, AlphaFold3, Chai-1, Boltz; choose a predictor by input and reconcile via pLDDT/PAE/pTM/ipTM |

## Example Prompts

- "Download PDB structure 4HHB"
- "Parse this mmCIF file and show the chains"
- "Convert this PDB to mmCIF format"
- "List all chains and their lengths"
- "Extract the protein sequence from chain A"
- "Find all ligands in this structure"
- "Measure the distance between CA atoms of residues 50 and 100"
- "Calculate the RMSD between these two structures, and tell me if they share a fold"
- "Superimpose these two structures"
- "Map the interface residues between chains A and B and the buried surface area"
- "Is this dimer interface biological or a crystal-packing artifact?"
- "Find all residues within 5 Angstroms of the ligand"
- "Add hydrogens and assign protonation states at pH 7.4 for docking"
- "Fill in the missing side-chain atoms and short loops before simulation"
- "Find the druggable pockets on this apo structure and rank them"
- "Remove all water molecules but keep the metals and cofactors"
- "Set B-factors based on conservation scores"
- "Is this structure reliable enough in the loop I care about? Check the resolution and B-factors"
- "Validate my structure: flag clashes and Ramachandran outliers"
- "Download the AlphaFold structure for this UniProt ID"
- "Analyze the pLDDT confidence scores and tell me which regions are disordered"
- "Plot the predicted aligned error (PAE)"
- "Is my AlphaFold model good enough to dock into?"
- "Predict the structure of this sequence with ESMFold"
- "Run AlphaFold3 on my protein complex and judge the interface by ipTM"
- "Compare predictions from ESMFold, Chai-1, and Boltz-1"

## Requirements

```bash
pip install biopython numpy requests freesasa

# Structure preparation (add hydrogens, protonation, fill atoms)
pip install pdbfixer openmm

# For modern structure prediction (choose per tool)
pip install fair-esm chai-lab boltz

# External CLI tools installed separately, per skill:
# validation: MolProbity / Phenix, DSSP (mkdssp)
# preparation: reduce, PROPKA (propka3), PDB2PQR
# binding-site detection: fpocket, P2Rank
```

## Related Skills

- **alignment** - structural-alignment for cross-protein fold comparison and Foldseek/TM-align search
- **sequence-io** - Read sequences to compare with structure-derived sequences
- **sequence-manipulation** - Analyze extracted protein sequences
- **database-access** - Fetch structure metadata and UniProt cross-references
- **chemoinformatics** - Dock and screen against predicted or experimental structures
- **immunoinformatics** - Structural epitope mapping from interface analysis
