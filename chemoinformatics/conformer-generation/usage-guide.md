# Conformer Generation Usage Guide

## Overview

Generate 3D conformer ensembles for molecules from 2D structures. Default is ETKDGv3 (Hawkins 2010, knowledge-enhanced distance geometry) + MMFF94 optimization for drug-like molecules; CREST + GFN2-xTB for macrocycles / peptides; ML methods (GeoMol) for ultralarge libraries.

## Prerequisites

```bash
pip install rdkit
# For CREST + GFN2-xTB (semi-empirical):
conda install -c conda-forge xtb crest
```

## Quick Start

Tell the AI agent what to do:
- "Generate 50 ETKDGv3 conformers for this SMILES with MMFF94 optimization"
- "Build a 3D conformer ensemble suitable for docking input"
- "Sample macrocycle conformers using CREST + GFN2-xTB"
- "Generate Boltzmann-weighted descriptor averages over a 100-conformer ensemble"
- "Prune conformers below 0.5 RMSD and within a 5 kcal/mol energy window"

## Example Prompts

### Drug-like conformer ensemble
> "Generate 20 ETKDGv3 conformers for each compound in library.csv. MMFF94 optimize. Prune at RMSD < 0.5 A and energy window 10 kcal/mol. Output SDF with conformer index in property."

### Macrocycle high-quality sampling
> "Sample 200 conformers of cyclosporine (CAS 59865-13-3) with CREST + GFN2-xTB. Filter to 5 kcal/mol energy window. Output cluster centroids."

### 3D QSAR descriptor input
> "For each compound in active_set.sdf, generate 50 ETKDGv3 conformers and compute Boltzmann-averaged 3D descriptors (asphericity, eccentricity, PMI). Output descriptor table."

### Docking input
> "Generate a single low-energy ETKDGv3 + MMFF94 conformer for each SMILES; write to multi-SDF for downstream Vina docking."

## What the Agent Will Do

1. Parse SMILES, add hydrogens, embed initial 3D coordinates with ETKDGv3.
2. Generate the requested number of conformers (using `EmbedMultipleConfs`).
3. Optimize each conformer with MMFF94 (fallback to UFF if MMFF parameters unavailable).
4. Compute energy of each conformer.
5. Prune by RMSD and energy window.
6. For macrocycles or complex molecules: optionally pipe to CREST + GFN2-xTB.
7. Output: SDF / XYZ with conformer ensemble + energies + properties.

## Tips

- Use n_conf = max(10, 5 * NumRotatableBonds + 10) heuristic (Hawkins 2017).
- For docking: 1 conformer usually OK if MM-optimized; up to 50 for cross-docking.
- For QSAR descriptors: 50-200 conformers + Boltzmann averaging.
- For macrocycles: use `useMacrocycleTorsions=True`; CREST gold standard.
- Set `randomSeed` for reproducibility.
- Always `AddHs()` before embedding.

## Related Skills

- chemoinformatics/molecular-io - Parse molecules
- chemoinformatics/molecular-descriptors - 3D descriptors from ensembles
- chemoinformatics/shape-similarity - Multi-conformer 3D shape searching
- chemoinformatics/virtual-screening - Generate 3D for docking
- chemoinformatics/free-energy-calculations - Conformers for MD/FEP setup
- chemoinformatics/pharmacophore-modeling - 3D pharmacophore from ensembles
