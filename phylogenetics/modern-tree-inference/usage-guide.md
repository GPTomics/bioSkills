# Modern ML Tree Inference - Usage Guide

## Overview

Maximum-likelihood phylogenetics with IQ-TREE2 and RAxML-NG. An ML tree maximizes the likelihood under an assumed substitution model on a fixed alignment, so it inherits every flaw of both, and the support it reports measures repeatability under resampling, not whether a clade is true. This skill covers ModelFinder model selection, UFBoot2 and SH-aLRT support with their correct cutoffs, gene and site concordance factors that are the honest measure at genome scale, partitioning, AU topology tests, and long-branch-attraction control. The recurring lesson: high support is consistent with being wrong, and the cure for confident-wrong trees is better models, alignments, and diagnostics, not more bootstrap replicates.

## Prerequisites

```bash
# IQ-TREE2 (model selection, UFBoot2, SH-aLRT, concordance factors, AU test, PMSF all built in)
conda install -c bioconda iqtree

# RAxML-NG (very large trees, transfer bootstrap, precise branch lengths, checkpointing)
conda install -c bioconda raxml-ng
```

- Input is a multiple-sequence alignment (FASTA, PHYLIP, or NEXUS). Its column homology is assumed correct and is the single largest source of confident-wrong trees, so trim and inspect it first (alignment/alignment-trimming, alignment/alignment-io).
- For concordance factors you also need per-locus alignments or gene trees.
- Conceptual prerequisite: support is not accuracy: a node can carry UFBoot 100 or PP 1.00 and still be wrong, especially on a contested deep branch. The branch-support and concordance-factor sections below cover how to tell.

## Quick Start

Tell your AI agent what you want to do:
- "Build a maximum-likelihood tree from alignment.fasta with model selection and dual support"
- "Run IQ-TREE2 with ModelFinder, 1000 ultrafast bootstrap replicates, and SH-aLRT, and interpret the support"
- "Find the best substitution model and partition scheme for my concatenated multi-gene alignment"
- "Compute gene and site concordance factors for my phylogenomic tree"
- "Test whether my constrained topology is rejected against the ML tree with the AU test"
- "My two longest branches group together at 100% bootstrap: check for long-branch attraction"

## Example Prompts

### Tree Inference and Model Selection
> "Infer an ML tree from alignment.fasta with IQ-TREE2, using ModelFinder Plus and reporting both UFBoot2 and SH-aLRT support."

> "Should this concatenated alignment use +G or FreeRate, and should I partition it? Run ModelFinder with partition merging."

> "Build a deep protein tree under a site-heterogeneous C60/PMSF model to guard against long-branch attraction."

### Support Interpretation
> "My tree has UFBoot 100 almost everywhere. Is it well resolved, and how do I tell which nodes are genuinely supported?"

> "Explain the difference between the UFBoot >=95 cutoff and the bootstrap >=70 rule, and which applies to my output."

> "Compute gene and site concordance factors and flag any node with high bootstrap but low concordance."

### Partitioning and Large Trees
> "Run a partitioned analysis with separate models per gene and proportional branch lengths."

> "I have 4000 taxa: infer the tree in RAxML-NG with transfer bootstrap and explain why TBE not the standard bootstrap."

### Topology Testing and Troubleshooting
> "Test with the AU test whether a tree constraining these two genera to be sisters is rejected."

> "A fast-evolving taxon is attracted to the outgroup. Detect and treat the long-branch attraction."

## What the Agent Will Do

1. Inspect the alignment format and sequence type, and flag that column homology and trimming gate everything downstream.
2. Run ModelFinder (`-m MFP`, BIC) for a single locus, or `-m MFP+MERGE -rcluster 10` for a partitioned dataset; choose a C60/PMSF mixture for deep / LBA-prone protein data.
3. Search the ML tree with dual support (`-B 1000 -bnni -alrt 1000`) and a fixed seed for reproducibility.
4. Interpret support with the correct cutoffs: a branch is strongly supported only if SH-aLRT >=80 AND UFBoot >=95, never the bootstrap-70 rule on UFBoot.
5. For phylogenomic data, compute gene and site concordance factors and flag UFBoot-100/gCF-~33 nodes as effectively unresolved (ILS or introgression).
6. Run AU topology tests when an a-priori hypothesis is being tested against the ML tree.
7. Diagnose and treat long-branch attraction with site-heterogeneous models, fast-site/taxon removal, and recoding cross-checks.
8. Route out-of-scope work: model-free distance trees, posteriors, species trees under ILS, or dating to the sibling skills.

## Tips

- Use `-m MFP` (not the legacy `-m TEST`) so FreeRate models are tested; expect `+R` to win on concatenated data and `+G` on single short genes.
- Always add `-bnni` with `-B`; it reins in the UFBoot inflation that model violation causes.
- UFBoot >=95 is strong support; do NOT apply the standard-bootstrap >=70 threshold to it; they are different scales.
- Report at least two support measures (UFBoot + SH-aLRT) and, for any phylogenomic tree, concordance factors as well.
- A node with UFBoot 100 and gCF ~35 is essentially unresolved; report it as contested, not as a clade, and consider a coalescent species tree.
- Prefer `-p` (edge-linked proportional) for partitions; reserve `-Q` (edge-unlinked) for genuine heterotachy, and merge with `-m MFP+MERGE` to avoid over-partitioning.
- For LBA, change the model first (C60/PMSF), then remove fast sites/taxa; believe a deep node only when it survives a site-heterogeneous model, fast-site removal, and recoding.
- Reach for RAxML-NG on very large trees, for transfer bootstrap, for the most precise branch lengths, or for robust checkpointing.

## Related Skills

- distance-calculations - model-corrected distances and fast NJ trees as a model-free alternative
- bayesian-inference - posteriors, MCMC convergence, and CAT-GTR site-heterogeneous models
- species-trees - coalescent species-tree estimation when concordance factors reveal ILS
- divergence-dating - time-scaled trees from the ML topology
- tree-manipulation - rooting, pruning, and collapsing low-support nodes
- tree-visualization - mapping support and concordance factors onto branches
- alignment/alignment-io - reading and writing the alignment the ML tree trusts as fixed
