# Modern Tree Inference - Usage Guide

## Overview

Build publication-quality maximum likelihood phylogenetic trees using IQ-TREE2 and RAxML-NG. Covers automatic model selection with ModelFinder, ultrafast bootstrap with proper interpretation, concordance factors for phylogenomic datasets, partitioned analyses, topology testing, and long branch attraction awareness.

## Prerequisites

```bash
# Install IQ-TREE2
conda install -c bioconda iqtree

# Install RAxML-NG
conda install -c bioconda raxml-ng

# Or download binaries:
# IQ-TREE2: http://www.iqtree.org/
# RAxML-NG: https://github.com/amkozlov/raxml-ng
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a maximum likelihood tree from my alignment with bootstrap support"
- "Find the best substitution model for my sequences and infer a tree"
- "Run IQ-TREE2 with ultrafast bootstrap and SH-aLRT on my FASTA file"
- "Compare multiple gene partitions with different evolutionary rates"
- "Compute concordance factors for my phylogenomic dataset"

## Example Prompts

### Basic Tree Inference
> "Build a phylogenetic tree from alignment.fasta using maximum likelihood with proper model selection"

> "Run IQ-TREE2 with ModelFinder, 1000 ultrafast bootstrap replicates, and SH-aLRT"

> "Infer a tree with RAxML-NG using GTR+G and 100 bootstrap replicates"

### Model Selection
> "What's the best substitution model for my DNA alignment? Use ModelFinder with FreeRate models"

> "Run ModelFinder on my protein alignment and explain why it chose this model"

> "Should I use AIC or BIC for model selection?"

### Branch Support
> "Add both UFBoot and SH-aLRT support to my tree and explain the thresholds"

> "Compute gene and site concordance factors for my phylogenomic dataset"

> "My tree has low support throughout. What does that mean and what should I try?"

### Partitioned Analysis
> "Analyze my concatenated multi-gene alignment with separate models per gene"

> "Should I use edge-linked or edge-unlinked partition models?"

> "Run partitioned analysis with automatic partition merging"

### Topology Testing
> "Compare these two candidate tree topologies using the AU test"

> "Is my constrained topology significantly worse than the unconstrained ML tree?"

### Troubleshooting
> "I suspect long branch attraction in my tree. How do I detect and fix it?"

> "My two longest branches are grouping together. Is this an artifact?"

## What the Agent Will Do

1. Check the alignment file format and sequence type (DNA/protein)
2. Select appropriate model selection strategy (MFP for most cases)
3. Run ML tree search with proper flags (-bnni for UFBoot2)
4. Generate both UFBoot2 and SH-aLRT support if appropriate
5. Interpret support values with correct thresholds (UFBoot >= 95, not >= 70)
6. Compute concordance factors for phylogenomic datasets
7. Flag potential issues (LBA, low support, model adequacy)
8. Set random seed for reproducibility

## Tips

- Always use `-m MFP` (not `-m TEST`) to include FreeRate models in model testing
- The `-bnni` flag reduces UFBoot overestimation, so always include it
- UFBoot >= 95 means strong support; do NOT apply the standard bootstrap >= 70 threshold
- For publication, report at minimum two support measures (UFBoot + SH-aLRT)
- Add concordance factors for phylogenomic datasets to distinguish genuine discordance from noise
- Set `--seed` for reproducibility; use `-T AUTO` for automatic thread detection
- For large datasets (>500 taxa), consider `-fast` mode
- If support is low throughout, consider coalescent methods (ASTRAL) or Bayesian approaches
- When LBA is suspected, try site-heterogeneous models (C60, CAT-GTR) or add taxa to break long branches

## Related Skills

- bayesian-inference - Bayesian alternative when ML is insufficient
- species-trees - Coalescent methods for gene tree discordance
- divergence-dating - Molecular clock and dating analyses
- tree-io - Read and convert tree output files
- tree-visualization - Visualize trees with support annotations
- distance-calculations - Distance-based tree methods
