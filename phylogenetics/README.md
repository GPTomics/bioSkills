# Phylogenetics

Phylogenetic tree analysis using Biopython's Bio.Phylo module for parsing, manipulating, visualizing, and constructing evolutionary trees.

## Overview

This category covers phylogenetic tree operations: reading/writing tree files in various formats, manipulating tree structure (rooting, pruning, ladderizing), visualizing trees with matplotlib, and building trees from sequence alignments using distance-based and parsimony methods.

**Tool type:** `python`
**Primary tools:** Bio.Phylo, Bio.Phylo.TreeConstruction, Bio.Phylo.Consensus

## Skills

| Skill | Description |
|-------|-------------|
| [tree-io](tree-io/) | Read, write, convert tree files (Newick, Nexus, PhyloXML, NeXML) |
| [tree-visualization](tree-visualization/) | Draw trees with matplotlib, customize labels and colors, export figures |
| [tree-manipulation](tree-manipulation/) | Root, prune, ladderize, collapse, and modify tree structure |
| [distance-calculations](distance-calculations/) | Compute distance matrices, build NJ/UPGMA/parsimony trees, bootstrap consensus |

## Workflow

```
Multiple Sequence Alignment
    |
    v
[distance-calculations] - Compute distances, build trees
    |
    +---> NJ/UPGMA tree
    |         |
    |         v
    |     [tree-manipulation] - Root, ladderize, prune
    |         |
    |         v
    |     [tree-visualization] - Draw, export figures
    |
    +---> Bootstrap trees
              |
              v
          Consensus tree
              |
              v
          [tree-io] - Save in various formats
```

## Supported Tree Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| `newick` | .nwk, .tre | Standard format with branch lengths |
| `nexus` | .nex | Rich format (MrBayes, PAUP) |
| `phyloxml` | .xml | XML with metadata support |
| `nexml` | .nexml | Modern XML format |

## Example Prompts

### Tree I/O
- "Read this Newick tree file and show the taxa"
- "Convert my Nexus tree to Newick format"
- "Parse all trees from this MrBayes output"

### Tree Visualization
- "Draw this tree and save as PDF"
- "Show bootstrap values on the tree"
- "Create a tree figure with branch lengths labeled"

### Tree Manipulation
- "Root this tree using Mouse as outgroup"
- "Remove all bacterial sequences from the tree"
- "Ladderize the tree for cleaner visualization"

### Distance Calculations
- "Build a neighbor joining tree from this alignment"
- "Create a bootstrap consensus with 1000 replicates"
- "Calculate pairwise distances between all taxa"

## Installation

```bash
pip install biopython matplotlib numpy
```

## Notes

- **`draw_graphviz()` deprecated** - Removed in Biopython 1.79. Use `Phylo.draw()` for rectangular trees. For radial/circular layouts, use external tools like ETE3 or DendroPy.
- **PhyloXML for metadata** - Convert to PhyloXML format for color and annotation support.
- **Bootstrap vs simple trees** - Always assess support with bootstrap analysis for publication.

## Related Skills

- **alignment** - Prepare MSAs for tree building
- **alignment-io** - Read alignment files as input
- **sequence-io** - Read sequences for alignment and tree building
- **database-access** - Fetch sequences from NCBI for phylogenetic analysis

## References

- [Bio.Phylo documentation](https://biopython.org/docs/latest/api/Bio.Phylo.html)
- [Bio.Phylo.TreeConstruction](https://biopython.org/docs/latest/api/Bio.Phylo.TreeConstruction.html)
- [Phylo cookbook](https://biopython.org/wiki/Phylo_cookbook)
