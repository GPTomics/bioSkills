# Trajectory Inference - Usage Guide

## Overview

Trajectory inference orders cells along a continuous differentiation manifold and assigns a pseudotime, fate probabilities, or RNA velocity direction. The central caution: a snapshot is not a movie. Pseudotime is geometry, not a clock; the existence of a continuum is a biological judgment made before ordering; the root choice flips every gene trend; and near bifurcations the honest output is a distribution over fates, not a hard branch label.

## Prerequisites

```r
# R packages
install.packages(c('BiocManager', 'Seurat', 'remotes'))
BiocManager::install(c('monocle3', 'slingshot', 'tradeSeq'))
remotes::install_github('satijalab/seurat-wrappers')   # SeuratWrappers: as.cell_data_set() for Seurat -> monocle3
```

```bash
# Python packages
pip install scanpy scvelo cellrank palantir
```

## Quick Start

Tell your AI agent what you want to do:
- "Check whether my cells form a real trajectory or are discrete types"
- "Order cells by pseudotime from the stem-cell population"
- "Compute fate probabilities at the branch point"
- "Run RNA velocity and tell me whether it is trustworthy here"

## Example Prompts

### Topology and Continuum Test
> "Run PAGA and tell me whether these clusters are actually connected"
> "Is this a real continuum or a mixture of discrete cell types?"
> "Build a PAGA-initialized UMAP so the global topology is faithful"

### Pseudotime and Rooting
> "Order cells by diffusion pseudotime rooted at the HSC marker, not by eye"
> "Show how the gene trends change if the root cluster changes"
> "Which method fits a bifurcating tree best here?"

### Fate Probabilities
> "Give me fate probabilities at the branch point, not hard branch labels"
> "Run CellRank 2 with a pseudotime kernel and find terminal states"
> "Compute differentiation potential (entropy) across the manifold"

### RNA Velocity
> "Run scVelo dynamical mode and check velocity confidence first"
> "The velocity arrows point backward in my mature cells - is velocity valid here?"
> "Re-run velocity with a second quantifier and check the direction is stable"

## What the Agent Will Do

1. Decide whether a continuum exists at all (biology + PAGA connectivity) before ordering anything
2. Fix topology first with PAGA, pruning weak edges by threshold
3. Choose a method by expected topology (Slingshot for trees, Palantir/CellRank for branching fate)
4. Anchor the root with orthogonal evidence (marker, real time, velocity, stemness), never UMAP aesthetics
5. Compute pseudotime, fate probabilities, or velocity with the matched method
6. For velocity, inspect confidence and phase portraits and check direction is stable across quantifiers
7. Report fate as probabilities near bifurcations and cross-validate across methods

## Tips

- **Topology before pseudotime** - run PAGA first; isolated clusters with no surviving connectivity edges are discrete types, not branches, and must not be ordered.
- **Root choice is not cosmetic** - it flips the sign of every gene trend; anchor it with a marker, real time, velocity, or a stemness score.
- **Pseudotime is not real time** - intervals are not durations and rates differ across lineages; only the RealTimeKernel/WOT exploit actual time.
- **Use fate probabilities near branches** - a progenitor's fate is genuinely undetermined, so hard branch labels are biologically false (Palantir branch_probs, CellRank fate probabilities).
- **State does not fully predict fate** - Weinreb 2020 showed sister cells in the same state diverge; treat state-based fate calls as predictions and validate with lineage data.
- **Velocity failure is the default expectation** - Bergen 2021; velocity is unreliable in mature/terminal/non-dividing systems and a clean stream plot is not evidence. Inspect phase portraits and velocity_confidence.
- **Quantifier choice is first-order** - velocyto vs kb-nac vs alevin-fry vs STARsolo can flip velocity sign; a direction not stable across two quantifiers is a pipeline artifact (Soneson 2021).
- **Do not consume raw arrows** - feed velocity into CellRank 2 as one kernel, combine with connectivity, and check robustness to dropping it.
- **CellRank macrostates assume metastability** - justify n_states with a Schur/eigenvalue gap or they may be coarse-graining artifacts of a smooth flow.

## Related Skills

single-cell/clustering - Leiden clusters and the kNN graph PAGA, DPT, and velocity moments depend on
single-cell/preprocessing - normalization, HVG, and PCA choices the inferred axis inherits
single-cell/lineage-tracing - orthogonal lineage ground truth that tests state-based fate calls
single-cell/cell-communication - downstream signaling analysis along the trajectory
differential-expression/deseq2-basics - pseudobulk DE between trajectory endpoints or branches
