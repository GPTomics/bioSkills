# Pharmacophore Modeling Usage Guide

## Overview

Build 3D pharmacophore models from ligands (ligand-based) or pockets (receptor-based, apo2ph4). Apply pharmacophore queries to library search via Pharmer / Pharmit. Supports scaffold hopping, bioisosteric replacement, and cross-target SAR transfer.

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
> "Apply the pharmacophore from query.ph4 to screen library.sdf. Report compounds matching at least 3 of 4 features within 1.5 A tolerance."

### Bioisostere expansion
> "Replace each pharmacophore feature with bioisosteric alternatives (carboxylate -> tetrazole, amine -> guanidine). Search for matches."

## What the Agent Will Do

1. Identify pharmacophore features (donor, acceptor, hydrophobe, aromatic) from input.
2. For receptor-based: use PLIP or PoseView to map ligand-residue contacts.
3. For ligand-based: align actives, find conserved features across set.
4. Apply geometric tolerances (1-2 A typical).
5. Search library by matching feature distance constraints.
6. Output: ranked hits + retrospective enrichment.

## Tips

- Receptor-based (with co-crystal) more reliable than ligand-based.
- Use bioactive conformer when possible; not first-generated conformer.
- Geometric tolerance 1.0-1.5 A for drug-like; up to 2 A for flexible.
- Pharmacophore matches are MORE specific than fingerprint matches but lower recall.
- Combine pharmacophore + 2D fingerprint for hybrid search.
- Validate pharmacophore on retrospective active set; enrichment > 5x.

## Related Skills

- chemoinformatics/molecular-io - Parse PDB / SDF
- chemoinformatics/conformer-generation - Generate 3D conformer ensembles
- chemoinformatics/shape-similarity - 3D shape adjacent to pharmacophore
- chemoinformatics/virtual-screening - Pharmacophore as docking pre-filter
- chemoinformatics/scaffold-analysis - 2D scaffold context
- chemoinformatics/generative-design - PharmacoForge for de novo
- structural-biology/structure-io - PDB handling
