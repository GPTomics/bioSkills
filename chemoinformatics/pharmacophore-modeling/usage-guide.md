# Pharmacophore Modeling Usage Guide

## Overview

Build 3D pharmacophore models from ligands or protein pockets. apo2ph4 derives LigandScout PML models from an apo pocket; PLIP can identify interactions in a co-crystal. Apply a pharmacophore only after representing it in the format required by Pharmer or Pharmit. PharmacoForge can generate candidate pharmacophores conditioned on a protein pocket, followed by compound-library retrieval.

## Prerequisites

```bash
pip install rdkit
# For PLIP interaction analysis (receptor-based):
pip install plip
```

## Quick Start

Tell the AI agent what to do:
- "Derive a pharmacophore from co-crystal of receptor + ligand"
- "Build a ligand-based pharmacophore from 5 active compounds"
- "Search library for compounds matching pharmacophore query"
- "Identify pharmacophore-equivalent bioisosteres"

## Example Prompts

### Receptor-based pharmacophore
> "From co-crystal complex.pdb, identify the bound ligand and compute pharmacophore features (donor, acceptor, hydrophobe, aromatic) at contact points. Use PLIP for interaction analysis."

### Ligand-based common pharmacophore
> "Given 5 active compounds (SMILES list), align and extract common 3D pharmacophore features. Use RDKit Pharm3D framework."

### Pharmacophore screening
> "Configure this four-feature pharmacophore in Pharmit and screen the selected library. Report compounds matching at least three features within 1.5 A tolerance and record the exact query representation."

### Bioisostere expansion
> "Replace each pharmacophore feature with bioisosteric alternatives (carboxylate -> tetrazole, amine -> guanidine). Search for matches."

## What the Agent Will Do

1. Identify pharmacophore features (donor, acceptor, hydrophobe, aromatic) from input.
2. For a co-crystal: use PLIP or PoseView to map ligand-residue contacts; for an apo pocket, use the documented apo2ph4 workflow and retain its PML output.
3. For ligand-based: align defensible active conformers, derive conserved feature correspondences and distance bounds with a documented workflow, then use RDKit or another search engine to apply the resulting model. RDKit `EmbedPharmacophore` embeds against an existing model; it does not derive the consensus.
4. Set geometric tolerances from aligned-feature variability, coordinate uncertainty, and retrospective validation.
5. Convert the model to the search engine's documented query format, then search by feature-distance constraints.
6. Output: ranked hits + retrospective enrichment.

## Tips

- Compare ligand- and receptor-based models on target-relevant validation; neither is universally more reliable.
- Use bioactive conformer when possible; not first-generated conformer.
- Calibrate geometric tolerances per feature and model; any tabulated ranges in the skill are repository starting heuristics. Flexible ligands or uncertain coordinates may justify wider starting values, but there is no universal drug-like tolerance.
- Measure specificity and recall on the project dataset; they depend on feature count, tolerances, conformers, and the fingerprint baseline.
- Combine pharmacophore + 2D fingerprint for hybrid search.
- Validate on a retrospective active/decoy set; enrichment >5x is only a repository starting heuristic and must be calibrated for the dataset and decoy construction.

## Related Skills

- chemoinformatics/molecular-io - Parse PDB / SDF
- chemoinformatics/conformer-generation - Generate 3D conformer ensembles
- chemoinformatics/shape-similarity - 3D shape adjacent to pharmacophore
- chemoinformatics/virtual-screening - Pharmacophore as docking pre-filter
- chemoinformatics/scaffold-analysis - 2D scaffold context
- chemoinformatics/generative-design - Generate molecules after pharmacophore-based retrieval
- structural-biology/structure-io - PDB handling
