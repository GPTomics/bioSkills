# chemoinformatics

## Overview
Computational chemistry for drug discovery covering molecular representations, standardization, conformer generation, descriptor calculation, similarity / shape / pharmacophore search, virtual screening (classical + ML), pose validation, free-energy calculations, QSAR, generative design, retrosynthesis, covalent inhibitor design, and PROTAC degrader design.

**Tool type:** mixed | **Primary tools:** RDKit, AutoDock Vina, GNINA, ADMETlab 3.0, DiffDock-L, Boltz-2, REINVENT 4, AiZynthFinder, chemprop, OpenFE

## Skills
| Skill | Description |
|-------|-------------|
| molecular-io | Read, write, convert molecular formats (SMILES, InChI, SDF, MOL2, PDB) with aromaticity / stereo / charge handling |
| molecular-standardization | ChEMBL pipeline standardization and parent handling; optional separately configured tautomer canonicalization |
| conformer-generation | ETKDGv3 + MMFF94 / CREST + GFN2-xTB 3D conformer ensembles |
| molecular-descriptors | ECFP4/6, MACCS, MAP4, MHFP6, RDKit FP, AtomPair, TopTorsion, physchem descriptors |
| similarity-searching | Tanimoto / Tversky / Dice / cosine; Butina clustering; activity cliffs; MCS |
| substructure-search | SMARTS pattern matching; PAINS / BRENK / REOS / ZINC alert catalogs |
| scaffold-analysis | Bemis-Murcko scaffolds; MMPA via mmpdb; R-group decomposition; scaffold-balanced splits |
| shape-similarity | 3D shape similarity (USRCAT, Open3DAlign, ROCS, ShaEP) |
| pharmacophore-modeling | RDKit Pharm3D ligand-based + apo2ph4 receptor-based pharmacophore |
| reaction-enumeration | Reaction SMARTS enumeration; RECAP / BRICS; R-group decomposition |
| retrosynthesis | AiZynthFinder 4.4+ planning and route scoring; maintained or pinned template-free comparisons (Chemformer is archived) |
| qsar-modeling | chemprop 2.2.x D-MPNN + RDKit 2D; OECD 5 principles; applicability domain |
| generative-design | REINVENT 4 generator types (Reinvent, Libinvent, Linkinvent, Mol2Mol); MPO scoring |
| admet-prediction | ADMETlab 3.0 (119 platform features, including 77 prediction models), ensemble hERG prediction, PAINS filters |
| virtual-screening | AutoDock Vina, SMINA, GNINA with CNN scoring; ZINC22 / Enamine REAL ultralarge |
| pose-validation | PoseBusters physical-validity, strain energy, geometric checks |
| ml-docking-rescoring | DiffDock-L + GNINA hybrid; Boltz-1/Boltz-2; Chai-1; AlphaFold3 ligand |
| free-energy-calculations | OpenFE RBFE / ABFE; FEP+; MBAR/BAR; MM/GBSA endpoint |
| covalent-design | Cys-selective warheads; GSH stability; DOCKovalent / HCovDock |
| protac-degraders | E3 ligase choice; linker design; ternary complex prediction (PRosettaC, DeepTernary) |

## Example Prompts
- "Standardize my SMILES library using the ChEMBL pipeline and deduplicate"
- "Build a 50-conformer ETKDGv3 + MMFF94 ensemble for descriptor calculation"
- "Find scaffold hops with high 3D shape similarity but low ECFP4 Tanimoto"
- "Run AiZynthFinder retrosynthesis on top REINVENT 4 generated molecules"
- "Train chemprop hERG classifier with scaffold-balanced split + conformal prediction"
- "Hybrid VS: DiffDock-L pose + GNINA rescore + PoseBusters validation"
- "Predict ternary complex for kinase + PROTAC + CRBN using PRosettaC"
- "Run OpenFE RBFE on 6-compound lead series; report cycle closure"
- "Design covalent BTK inhibitor with acrylamide warhead and predict reactivity"
- "Compute Boltz-2 affinity prediction for 1000 docking candidates"

## Requirements
```bash
pip install rdkit chembl_structure_pipeline deepchem chemprop aizynthfinder posebusters
# OpenFE is distributed through conda-forge; the PyPI package named openfe is unrelated:
conda install -c conda-forge openfe=1.7
# REINVENT 4: follow the official environment and installation instructions:
# https://github.com/MolecularAI/REINVENT4
# For semi-empirical (CREST + GFN2-xTB):
conda install -c conda-forge xtb crest
# For GNINA:
# https://github.com/gnina/gnina
# For OpenMM / FEP:
conda install -c conda-forge openmm openff-toolkit gromacs
```

## Related Skills
- **structural-biology** - Protein structure prep for docking; AlphaFold3 / Boltz-1 for receptor prediction
- **machine-learning** - Adjacent ML approaches (biomarker discovery, model validation)
- **clinical-databases** - Drug bioactivity / pharmacogenomics overlay
