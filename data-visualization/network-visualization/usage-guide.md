# Network Visualization - Usage Guide

## Overview

Network visualizations render protein-interaction, gene-regulatory, co-expression, and pathway graphs. Layout choice dominates the result: force-directed (Fruchterman-Reingold, ForceAtlas2) is general-purpose but layout positions are NOT biologically meaningful; hive plots (Krzywinski 2012) and circular layouts preserve categorical structure; hierarchical edge bundling (Holten 2006) reduces hairball clutter for many-to-many. NetworkX + matplotlib for static publication; PyVis for interactive HTML; Cytoscape for journal-grade composition.

## Prerequisites

```bash
pip install networkx pyvis py4cytoscape matplotlib numpy
# Large networks:
pip install datashader
```

```r
install.packages(c('igraph', 'ggraph', 'tidygraph'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Plot PPI network with degree-sized nodes and community-colored modules using NetworkX"
- "ForceAtlas2 layout for a 5000-node scale-free network"
- "Interactive PyVis HTML with hover tooltips"
- "Send the network to Cytoscape with degree-based styling via py4cytoscape"
- "Hive plot anchored by module assignment"

## Example Prompts

### Standard PPI

> "PPI network of 200 proteins from STRING. Spring layout with fixed seed. Color nodes by community via greedy_modularity_communities. Size by degree. Label only hubs (degree >= 5)."

### Large network

> "5000-node co-expression network. Use ForceAtlas2 layout. Bundle edges with hierarchical edge bundling. Rasterize for PDF."

### Comparing conditions

> "Two PPI networks (control vs treatment). Compute layout on the UNION graph; render both with the same `pos` so visual changes reflect biology not layout."

### Interactive supplement

> "PyVis HTML with degree-based node size and module-based color; force-atlas physics; export standalone HTML."

### Cytoscape pipeline

> "Send NetworkX graph to Cytoscape via py4cytoscape; apply degree-mapping style; export PDF."

## What the Agent Will Do

1. Load network from edge list / GraphML / SIF / interaction-database query.
2. Compute basic stats (degree distribution, density, modularity).
3. Choose layout: force-directed for general; hive for comparative; circular for small dense; edge bundling for many-to-many.
4. Set random seed (always); reuse `pos` when comparing conditions.
5. Compute node attributes (degree, betweenness, community assignment).
6. Render with NetworkX layered drawing (`draw_networkx_edges`, `nodes`, `labels`) for fine control.
7. Label hubs only (cap at degree threshold or ~30 labels).
8. Encode edge weights via line width when meaningful.
9. Export PDF for static (with raster for >2000 nodes) or HTML for interactive.

## Tips

- **Layout positions are NOT biology.** Frame conclusions on edge existence and node degree, not inter-node distance.

- **Always set seed.** `seed=42` in NetworkX, `set.seed(42)` for R igraph. Without it, layout varies across runs.

- **Comparing conditions requires shared layout.** Compute on union network; reuse `pos`.

- **ForceAtlas2 for scale-free** (>500 nodes, hub-spoke topology). Spring for general; Kamada-Kawai for small dense.

- **Hive plots (Krzywinski 2012)** anchor nodes to axes by categorical attribute - biology-faithful alternative to force-directed.

- **Hierarchical edge bundling (Holten 2006)** dramatically reduces clutter for many-to-many networks within a hierarchy.

- **Cap labels at hubs only** (degree > 5 or genes of interest). Labeling all nodes destroys readability.

- **Encode edge weight via line width** when interaction confidence varies - uniform widths hide information.

- **PyVis max ~2000 nodes** for browser HTML. Above this use Cytoscape.js or Datashader.

- **Datashader for >50000 nodes** - raster aggregation; only honest display at that scale.

- **Cytoscape (desktop) for publication-grade biological networks.** py4cytoscape exposes script control from Python/R.

## Related Skills

- gene-regulatory-networks/coexpression-networks - Build the network
- database-access/interaction-databases - Fetch PPI data
- data-visualization/multipanel-figures - Combine with other plots
- data-visualization/color-palettes - Community palette
- single-cell/cell-communication - Cell-cell interaction networks
