# phylogenetics

## Overview

Phylogenetic tree analysis using Biopython's Bio.Phylo module. Covers reading/writing tree files, manipulating tree structure (rooting, pruning, ladderizing), visualizing trees with matplotlib, and building trees from sequence alignments.

**Tool type:** python | **Primary tools:** Bio.Phylo, Bio.Phylo.TreeConstruction

## Skills

| Skill | Description |
|-------|-------------|
| tree-io | Read, write, convert tree files (Newick, Nexus, PhyloXML, NeXML) |
| tree-visualization | Draw trees with matplotlib, customize labels and colors, export figures |
| tree-manipulation | Root, prune, ladderize, collapse, and modify tree structure |
| distance-calculations | Compute distance matrices, build NJ/UPGMA/parsimony trees, bootstrap consensus |

## Example Prompts

- "Read this Newick tree file and show the taxa"
- "Convert my Nexus tree to Newick format"
- "Parse all trees from this MrBayes output"
- "Draw this tree and save as PDF"
- "Show bootstrap values on the tree"
- "Create a tree figure with branch lengths labeled"
- "Root this tree using Mouse as outgroup"
- "Remove all bacterial sequences from the tree"
- "Ladderize the tree for cleaner visualization"
- "Collapse clades with bootstrap < 70"
- "Build a neighbor joining tree from this alignment"
- "Create a bootstrap consensus with 1000 replicates"
- "Calculate pairwise distances between all taxa"
- "Build UPGMA tree from distance matrix"

## Requirements

```bash
pip install biopython matplotlib numpy
```

## Related Skills

- **alignment** - Prepare MSAs for tree building
- **sequence-io** - Read sequences for alignment and tree building
- **database-access** - Fetch sequences from NCBI for phylogenetic analysis
