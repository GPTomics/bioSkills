# Shape Similarity Usage Guide

## Overview

3D shape-based similarity searching using USRCAT (ultrafast), Open3DAlign (RDKit), ROCS (commercial), or ShaEP. Find scaffold-hopped compounds that 2D fingerprints miss; identify bioisosteric replacements via shape + color (Tanimoto-Combo).

## Prerequisites

```bash
conda install -c conda-forge rdkit
```

RDKit provides USRCAT through `rdMolDescriptors`; no separate `usrcat` package is required for the documented workflow.
ShaEP is a separate binary; verify the official current release and use its documented `-q`, positional input/output, `-s`, and `--output-file` syntax.

## Quick Start

Tell the AI agent what to do:
- "Find compounds with similar 3D shape to my query molecule"
- "Score library with USRCAT shape descriptors; top 100 hits"
- "Open3DAlign rescoring on top USRCAT hits"
- "Find scaffold-hopped compounds using shape-high and ECFP4-low cutoffs calibrated on my reference set"

## Example Prompts

### USRCAT pre-filter
> "Compute USRCAT descriptors for query.sdf and library.sdf. Rank library by similarity; return top 500."

### Open3DAlign rescore
> "Take top 500 USRCAT hits and rescore with Open3DAlign for accurate shape alignment. Return top 50."

### Scaffold-hopping
> "Calibrate shape-similarity and ECFP4-dissimilarity cutoffs on my reference set, then output 20 candidate scaffold hops."

### Conformer-aware shape search
> "For each library compound, generate 20 conformers; find best-shape conformer match to query. Use Open3DAlign."

## What the Agent Will Do

1. Generate 3D conformers, verify embedding and MMFF coverage/convergence, and record failures.
2. Compute shape descriptors (USRCAT) or align (Open3DAlign).
3. Rank hits by shape Tanimoto / Open3DAlign score.
4. Optionally compute Tanimoto-Combo (shape + color/pharmacophore).
5. Output ranked list with shape score + ECFP4 similarity for diversity check.

## Tips

- Benchmark throughput on the actual conformer count, hardware, and library; published or vendor rates are not interchangeable.
- ROCS TanimotoCombo is shape Tanimoto plus color Tanimoto (range 0-2), whereas RDKit O3A's raw score is unnormalized.
- Use conformer ensembles sized from a convergence check; 20 conformers is a repository starting budget, not a universal minimum.
- Calibrate the shape-high/ECFP4-low quadrant on a task-relevant reference set rather than applying universal 0.7/0.5 cutoffs.
- ROCS TanimotoCombo ranges from 0 to 2; choose follow-up cutoffs from a query- and library-matched benchmark.
- Validate with docking on top shape hits; not all shape matches dock well.
- Open3DAlign mutates probe coordinates; preserve or copy the best-aligned conformer when coordinates are part of the output.

## Related Skills

- chemoinformatics/molecular-io - Parse molecules
- chemoinformatics/conformer-generation - Generate conformer ensembles
- chemoinformatics/similarity-searching - 2D similarity comparison
- chemoinformatics/pharmacophore-modeling - Pharmacophore alternative
- chemoinformatics/virtual-screening - Shape as pre-filter
