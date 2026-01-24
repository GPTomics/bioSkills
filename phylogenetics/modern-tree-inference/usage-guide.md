# Modern Tree Inference

## Overview

Build publication-quality maximum likelihood phylogenetic trees using IQ-TREE2 and RAxML-ng. These tools provide automatic model selection, ultrafast bootstrap, and partitioned analyses for multi-gene datasets.

## Prerequisites

```bash
# Install IQ-TREE2
conda install -c bioconda iqtree

# Install RAxML-ng
conda install -c bioconda raxml-ng

# Or download binaries from:
# IQ-TREE2: http://www.iqtree.org/
# RAxML-ng: https://github.com/amkozlov/raxml-ng
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a maximum likelihood tree from my alignment with bootstrap support"
- "Find the best substitution model for my sequences and infer a tree"
- "Run IQ-TREE2 with ultrafast bootstrap on my FASTA file"
- "Compare multiple gene partitions with different evolutionary rates"

## Example Prompts

### Basic Tree Inference
> "Build a phylogenetic tree from alignment.fasta using maximum likelihood"

> "Run IQ-TREE2 with automatic model selection and 1000 bootstrap replicates"

> "Infer a tree with RAxML-ng using the GTR+G model"

### Model Selection
> "What's the best substitution model for my DNA alignment?"

> "Compare GTR, HKY, and K2P models for my sequences"

> "Run ModelFinder on my protein alignment"

### Bootstrap Analysis
> "Add ultrafast bootstrap support to my tree"

> "Run both UFBoot and SH-aLRT tests for branch support"

> "Generate a consensus tree from bootstrap replicates"

### Partitioned Analysis
> "Analyze my concatenated alignment with separate models per gene"

> "Create a partition file for my multi-gene dataset"

> "Run edge-linked partition analysis with IQ-TREE2"

### Advanced
> "Compare multiple candidate tree topologies with AU test"

> "Infer a tree with a monophyly constraint"

> "Resume my interrupted IQ-TREE2 run"

## What the Agent Will Do

1. Check that your alignment file exists and is properly formatted
2. Select an appropriate substitution model (or use automatic model selection)
3. Run the ML tree search with specified parameters
4. Generate bootstrap support values if requested
5. Report the best tree location and key statistics
6. Provide guidance on interpreting branch support values

## Tips

- Use ultrafast bootstrap (-B 1000) for quick analyses; standard bootstrap (-b 100) for publication
- For multi-gene datasets, partitioned models often provide better results
- IQ-TREE2's ModelFinder (-m MFP) automatically selects the best-fit model
- UFBoot values >= 95 and SH-aLRT >= 80 indicate strong support
- Set a random seed (--seed) for reproducible results
- Use -T AUTO to automatically detect available CPU threads
- For large datasets (>500 taxa), consider -fast mode to reduce runtime
