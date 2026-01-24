# Tree Visualization - Usage Guide

## Overview

This skill creates publication-quality phylogenetic tree figures using Bio.Phylo's matplotlib integration. It supports customizing labels, colors, and exporting to various image formats.

## Prerequisites

```bash
pip install biopython matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Draw this Newick tree and save as PNG"
- "Create a publication-quality tree figure with bootstrap values"
- "Highlight the primate clade in red"

## Example Prompts

### Basic Drawing
> "Plot this tree with taxa labels"

> "Show me an ASCII representation of this tree"

> "Draw the tree and save as PDF"

### Customization
> "Add branch length labels to the tree"

> "Show bootstrap support values at internal nodes"

> "Make the tree figure taller to fit all labels"

### Coloring
> "Color the human and chimp clades in red"

> "Highlight clades with bootstrap > 90 in green"

### Export
> "Save the tree as SVG for my publication"

> "Export high-resolution PNG at 300 DPI"

## What the Agent Will Do

1. Read the tree file using Bio.Phylo
2. Pre-process tree (ladderize, set missing lengths)
3. Configure figure size based on number of taxa
4. Apply requested customizations (colors, labels)
5. Render tree using matplotlib
6. Save to requested format

## Output Formats

| Format | Best For |
|--------|----------|
| PNG | Presentations, web |
| PDF | Publications (vector) |
| SVG | Web, editing in Illustrator |

## Tips

- Always `ladderize()` trees before drawing for cleaner appearance
- Scale figure height with number of taxa (0.3 inches per taxon works well)
- Set `do_show=False` when saving to file to avoid display issues
- Convert to PhyloXML for better color support
- Use `bbox_inches='tight'` when saving to avoid cropped labels
