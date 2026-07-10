# Tree Visualization - Usage Guide

## Overview

This skill draws, styles, annotates, and exports phylogenetic tree figures, and chooses the drawing tool. The central idea is that a tree figure is an argument, not a neutral picture: a Newick or Nexus file fixes only the topology and (optionally) branch lengths and one node label, while the layout geometry, node ordering, root placement, and which support value is shown are all supplied by whoever draws the tree. The same file can be rendered to tell mutually contradictory but technically correct stories. Four choices dominate: cladogram vs phylogram vs chronogram (what the lengths mean), ladderization (which manufactures a false arrow of progress), root placement (which decides what is "basal"), and support display (a bare integer is the quietest lie because bootstrap, posterior, SH-aLRT, and UFBoot are different scales).

A second decisive point is that the drawing tool is decided by the I/O layer, not the surrounding language: BEAST/MrBayes HPD intervals and posteriors plot cleanly only through treeio plus ggtree, because Bio.Phylo flattens the tree to topology plus one label and silently drops the uncertainty that was the result.

## Prerequisites

- Python: `pip install biopython matplotlib` for quick rectangular figures from a Python pipeline.
- R alternative for publication figures: `ggtree`, `treeio`, and `ggtreeExtra` (Bioconductor) for circular layouts, metadata, dual support, and BEAST HPD bars.
- Web/desktop: iTOL v6 (browser) for large trees and template-based annotation; FigTree (desktop) for interactive inspection of a BEAST tree.
- Conceptual: know whether the tree carries branch lengths (phylogram vs cladogram), whether it is rooted, and which support measure the node labels hold. A figure that does not declare these makes silent claims.

## Quick Start

Tell your AI agent what you want to do:
- "Draw this Newick tree as a phylogram and save it as a vector PDF"
- "Show bootstrap support at the internal nodes and label it as bootstrap in the caption"
- "Color the branches for one clade red and ladderize for legibility"
- "My tree has 800 tips and the labels collide - what layout should I use?"
- "I have a BEAST MCC tree with HPD bars - how do I draw the node-age uncertainty?"

## Example Prompts

### Drawing and exporting
> "Read tree.nwk, ladderize it, draw it as a phylogram with a title naming the branch-length unit, and export to tree.pdf as vector"
> "Give me an ASCII diagram of this tree so I can sanity-check the topology in a log"

### Layout choice for size
> "This tree has ~600 tips and rectangular labels are an unreadable band - recommend and produce a circular layout with radial labels"
> "I want metadata rings (host species, sampling year) around a circular tree of 400 genomes"

### Support and honesty
> "Show the SH-aLRT/UFBoot dual support at each node and state which is which in the legend"
> "Collapse every node below 70% bootstrap into a polytomy and say so in the caption"

### Annotated Bayesian trees
> "Draw this BEAST chronogram with the 95% HPD bars on node ages and the posterior probabilities labeled"
> "Bio.Phylo dropped my HPD intervals - route this through treeio and ggtree instead"

## What the Agent Will Do

1. Establish the four argument-level choices: branch-length encoding (cladogram/phylogram/chronogram), node ordering, root, and which support measure to show.
2. Choose the drawing tool from the I/O needs, not the language: Bio.Phylo for a fast Python look; ggtree + treeio for metadata, dual support, or HPD bars; ETE4 for programmatic per-node styling; iTOL for large trees or no-code annotation.
3. Pick a layout from the tip count (rectangular below ~150 tips, circular/fan above, strips/rings beyond ~1000) and lock the branch-length scale.
4. Render with support labeled by its measure, color by group rather than narrating ordering, and add a scale bar to any phylogram.
5. Export vector (SVG/PDF/EPS) for publication, write outputs to a temp/namespaced path, and state in the caption the layout, length unit, ordering disclaimer, and support measure.

## Tips

- The geometry is a claim: never draw untrusted branch lengths as a phylogram, and never draw a chronogram without its age uncertainty bars.
- Ladderize for legibility only, and say in the caption that ordering carries no phylogenetic meaning; color by group so non-monophyly cannot be hidden by rotation.
- A bare support number always flatters the result. Label the measure and order; remember a posterior of 0.95 is weaker than a bootstrap of 95.
- Route annotated Bayesian trees through treeio + ggtree (or FigTree for inspection); Bio.Phylo silently drops BEAST HPD bars and posteriors.
- Export vector for any publication figure; only rasterize at final size with >=600 dpi if a journal forces it, and never trust a 150-dpi PNG of a tree.
- An unrooted/radial tree has no basal taxon; root explicitly before making any "early-diverging" claim and show the root.

## Related Skills

- tree-io - parsing and preserving BEAST/MrBayes/IQ-TREE annotations so they survive into the figure
- tree-manipulation - rooting, ladderizing, pruning, and collapsing low-support nodes before drawing
- data-visualization/ggplot2-fundamentals - the grammar, themes, and ggsave vector export underlying ggtree
- data-visualization/multipanel-figures - composing a tree plus aligned metadata panels into one figure
