# Lineage Tracing - Usage Guide

## Overview

Lineage tracing reads heritable marks (CRISPR/Cas9 scars, static expressed barcodes, or somatic mtDNA mutations) across single cells and asks which cells share which marks to recover an ontogenetic phylogeny or clonal grouping. The central idea: transcriptomic state does not fully predict fate (Weinreb 2020), so lineage information is orthogonal to expression and reconstruction is phylogenetics on error-prone scars where homoplasy and dropout dominate the result.

## Prerequisites

```bash
pip install cassiopeia-lineage cospar scanpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a lineage tree from my CRISPR scar character matrix"
- "Compare solvers and tell me how robust the topology is"
- "Group clonally related cells from mtDNA mutations"
- "Recover early fate bias from sparse clonal barcodes"

## Example Prompts

### Tree Reconstruction
> "Reconstruct a lineage tree with Cassiopeia and collapse mutationless edges"
> "My data has heavy homoplasy and dropout - which solver should I use?"
> "Run a panel of solvers and compare them with Robinson-Foulds and triplets-correct"

### Clonal Dynamics and State
> "Integrate clones with transcriptomic state using CoSpar"
> "Compute early fate bias for monocyte vs neutrophil from clonal data"
> "Which clones expanded over time and what is their signature?"

### Mitochondrial and Native Tracing
> "Use mtDNA mutations to group clonally related cells in this human sample"
> "Why is my mtDNA giving clonal blobs instead of a deep tree?"

### Lineage vs State
> "Does transcriptomic state predict fate in my system, or is lineage adding new information?"
> "Validate my state-based trajectory against the lineage tree"

## What the Agent Will Do

1. Choose the recording assay by the question (deep topology, clonal state->fate, or retrospective native tissue)
2. Build a character matrix keeping missing (-1) distinct from unedited (0)
3. Assess missingness, informativeness, and barcode-collision risk before reconstruction
4. Reconstruct with a solver matched to scale and to homoplasy/dropout severity
5. Weight indels by formation probability and run a panel of solvers for robustness
6. Compare trees with Robinson-Foulds and depth-stratified triplets-correct
7. Integrate clone with state (CoSpar) to recover hidden fate bias and test state->fate

## Tips

- **State does not predict fate** - Weinreb 2020 showed sisters in the same state diverge; lineage data is orthogonal to expression, not redundant.
- **Missing is not unedited** - keep the -1 missing state distinct from 0; collapsing them is the most consequential preprocessing error because heritable dropout deletes whole clades and biases topology.
- **Homoplasy breaks parsimony** - non-uniform Cas9 indels mean unrelated cells share frequent scars; weight indels by probability and use Startle's star-homoplasy model under heavy convergence.
- **Run a solver panel** - VanillaGreedy/Hybrid/ILP/NJ rarely agree everywhere; agreement is the practical certainty signal. HybridSolver is the scalable default.
- **Deep splits are least certain** - early near-root splits rest on the fewest, most-overwritten characters yet matter most; report branch support and trust leaf structure more.
- **mtDNA gives clonal grouping, not deep trees** - low mutation rate, hotspot homoplasy, heteroplasmy drift, and selection limit it; blacklist NUMTs, RNA-edit, and hotspot positions.
- **Prospective vs retrospective matters** - only prospective barcodes installed before the process can test state->fate; retrospective mtDNA recovers ancestry but cannot establish the preceding state.
- **Static barcodes are flat clones** - LARRY/CellTag give clonal membership, not division-order topology; use library complexity far above founder number to avoid collisions.
- **Integrate clone with state** - CoSpar recovers fate bias from sparse clones rather than assuming the manifold encodes fate, but does not build a phylogeny.

## Related Skills

single-cell/trajectory-inference - state-based pseudotime/velocity that lineage data tests and corrects
single-cell/preprocessing - QC, doublet handling, and normalization upstream of barcode and clone calls
single-cell/clustering - cell-type labels annotated onto tree leaves and clones
phylogenetics/modern-tree-inference - general phylogenetic inference, parsimony vs ML, and branch support
