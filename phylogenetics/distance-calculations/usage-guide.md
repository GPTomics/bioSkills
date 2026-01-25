# Distance Calculations - Usage Guide

## Overview

This skill computes evolutionary distances from sequence alignments and builds phylogenetic trees using distance-based methods (Neighbor Joining, UPGMA), parsimony, and bootstrap consensus approaches.

## Prerequisites

```bash
pip install biopython numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a neighbor joining tree from this alignment"
- "Calculate the distance matrix for these sequences"
- "Create a bootstrap consensus tree with 100 replicates"

## Example Prompts

### Distance Matrices
> "Calculate pairwise distances from this protein alignment"

> "Create a distance matrix using BLOSUM62 scoring"

> "Show me the distance matrix for my sequences"

### Tree Building
> "Build a UPGMA tree from my alignment"

> "Create a neighbor joining tree from these sequences"

> "Build a parsimony tree from this alignment"

### Bootstrap Analysis
> "Generate 1000 bootstrap trees and find the consensus"

> "What is the bootstrap support for each clade?"

> "Create a majority rule consensus tree"

### Tree Metrics
> "What is the total branch length of this tree?"

> "Calculate the distance between Human and Mouse in this tree"

> "Find the tree height (maximum root-to-tip distance)"

## What the Agent Will Do

1. Read alignment file
2. Select appropriate distance model
3. Compute distance matrix
4. Build tree using specified method
5. Optionally perform bootstrap analysis
6. Return tree and/or distance matrix

## Tree Building Methods

| Method | Best For |
|--------|----------|
| Neighbor Joining | General use, fast |
| UPGMA | Molecular clock assumption |
| Parsimony | Small datasets, morphology |
| Bootstrap Consensus | Confidence assessment |

## Tips

- Use `identity` model for quick exploratory analysis
- Use `blosum62` for protein alignments
- NJ is generally preferred over UPGMA (no clock assumption)
- 100-1000 bootstrap replicates is typical
- Check alignment quality before tree building
