# Shape Similarity Usage Guide

## Overview

3D shape-based similarity searching using USRCAT (ultrafast), Open3DAlign (RDKit), ROCS (commercial), or ShaEP. Find scaffold-hopped compounds that 2D fingerprints miss; identify bioisosteric replacements via shape + color (Tanimoto-Combo).

## Prerequisites

```bash
pip install rdkit usrcat
```

## Quick Start

Tell the AI agent what to do:
- "Find compounds with similar 3D shape to my query molecule"
- "Score library with USRCAT shape descriptors; top 100 hits"
- "Open3DAlign rescoring on top USRCAT hits"
- "Find scaffold-hopped compounds (shape > 0.7, ECFP4 < 0.5)"

## Example Prompts

### USRCAT pre-filter
> "Compute USRCAT descriptors for query.sdf and library.sdf. Rank library by similarity; return top 500."

### Open3DAlign rescore
> "Take top 500 USRCAT hits and rescore with Open3DAlign for accurate shape alignment. Return top 50."

### Scaffold-hopping
> "Find scaffold hops: shape similarity > 0.7 AND ECFP4 < 0.5 to query. Output 20 candidate scaffold-hopped molecules."

### Conformer-aware shape search
> "For each library compound, generate 20 conformers; find best-shape conformer match to query. Use Open3DAlign."

## What the Agent Will Do

1. Generate 3D conformers for query and library (ETKDGv3 + MMFF94).
2. Compute shape descriptors (USRCAT) or align (Open3DAlign).
3. Rank hits by shape Tanimoto / Open3DAlign score.
4. Optionally compute Tanimoto-Combo (shape + color/pharmacophore).
5. Output ranked list with shape score + ECFP4 similarity for diversity check.

## Tips

- USRCAT 100k mols/sec; Open3DAlign 100 mols/sec; ROCS commercial GPU 1k mols/sec.
- Always use conformer ensembles (≥20); single conformer misses ~30% of true hits.
- Shape > 0.7 AND ECFP4 < 0.5 = scaffold-hopping gold quadrant.
- ROCS Tanimoto-Combo > 0.7 = strong hit; > 1.0 rare.
- Validate with docking on top shape hits; not all shape matches dock well.

## Related Skills

- chemoinformatics/molecular-io - Parse molecules
- chemoinformatics/conformer-generation - Generate conformer ensembles
- chemoinformatics/similarity-searching - 2D similarity comparison
- chemoinformatics/pharmacophore-modeling - Pharmacophore alternative
- chemoinformatics/virtual-screening - Shape as pre-filter
