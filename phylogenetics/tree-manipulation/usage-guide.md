# Tree Manipulation - Usage Guide

## Overview

This skill modifies phylogenetic tree structure, including rooting, pruning, ladderizing, and extracting subtrees. Use it to prepare trees for analysis or visualization.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Root this tree using Mouse as outgroup"
- "Remove all bacteria from this tree"
- "Ladderize the tree for better visualization"

## Example Prompts

### Rooting
> "Root this tree with TaxonA as the outgroup"

> "Reroot the tree at the midpoint"

> "Check if this tree is rooted"

### Pruning
> "Remove all sequences except humans and primates"

> "Drop taxa with incomplete data from the tree"

> "Keep only the top 10 most divergent taxa"

### Ladderizing
> "Sort the branches so larger clades are at the bottom"

> "Ladderize the tree before visualization"

### Subtrees
> "Extract the primate subtree"

> "Find the common ancestor of Human and Mouse"

> "Get all taxa in the mammal clade"

### Other Modifications
> "Rename taxa using this mapping table"

> "Collapse branches with support below 70%"

> "Scale branch lengths to percentages"

## What the Agent Will Do

1. Read the tree from file or string
2. Apply requested modifications
3. Verify tree integrity after changes
4. Save or return the modified tree

## Common Operations

| Operation | Use Case |
|-----------|----------|
| Rooting | Required before some analyses |
| Ladderizing | Cleaner visualization |
| Pruning | Focus on taxa of interest |
| Collapsing | Simplify poorly-supported nodes |
| Renaming | Use readable names |
| Random trees | Testing, simulations |
| Splitting | Expand terminal nodes |

## Tips

- Root trees before comparing topologies
- Ladderize before drawing for consistent appearance
- Check `is_monophyletic()` before rooting with multiple outgroups
- Prune iteratively when removing many taxa
- Save backup before destructive operations
