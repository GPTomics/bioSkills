# Network Visualization - Usage Guide

## Overview

This skill enables AI agents to create static, interactive, and publication-quality visualizations of biological networks. Supports protein interaction networks, gene regulatory networks, co-expression modules, and pathway graphs using NetworkX, PyVis, and Cytoscape automation.

## Prerequisites

```bash
# Core
pip install networkx matplotlib numpy

# Interactive HTML networks
pip install pyvis

# Cytoscape automation (requires Cytoscape desktop running)
pip install py4cytoscape
```

## Quick Start

Tell your AI agent what you want to do:

- "Visualize my protein interaction network from STRING"
- "Create an interactive HTML network from my gene interactions"
- "Color network nodes by community membership"
- "Draw my gene regulatory network with activation and repression edges"
- "Send my NetworkX graph to Cytoscape for styling"

## Example Prompts

### Static Network Plots

> "Draw a network from my interaction DataFrame with node sizes proportional to degree"

> "Create a PPI network plot with edge colors based on confidence scores"

> "Make a side-by-side comparison of treated vs control co-expression networks"

### Community and Clustering

> "Detect communities in my network and color the nodes by group"

> "Highlight my differentially expressed genes on the PPI network"

> "Show the top 10 hub genes in a different color"

### Interactive Networks

> "Create an interactive HTML network I can share with collaborators"

> "Build a PyVis network with hover tooltips showing gene annotations"

### Gene Regulatory Networks

> "Draw my GRN with arrows showing activation in green and repression in red"

> "Visualize the SCENIC regulon network with TFs as diamond nodes"

### Cytoscape

> "Send my NetworkX graph to Cytoscape and apply degree-based styling"

> "Export a publication-quality PDF from Cytoscape"

## What the Agent Will Do

1. Load or construct a NetworkX graph from interaction data
2. Select an appropriate layout algorithm based on network size and topology
3. Apply node sizing (degree), coloring (community or attribute), and edge styling (weight or type)
4. Render the network as a static figure, interactive HTML, or Cytoscape session
5. Export in the requested format (PNG, PDF, SVG, HTML, GraphML)

## Tips

- **Layout seed** - Always set `seed=42` in spring_layout for reproducible figures
- **Spring k parameter** - Increase `k` (default ~1/sqrt(n)) for sparser layouts with better label readability
- **Kamada-Kawai** - Use for small-medium networks (under 200 nodes) where spring layout is too tangled
- **Label top nodes only** - For large networks, label only hub nodes (degree >= threshold) to avoid clutter
- **PyVis for exploration** - Use interactive HTML networks for initial exploration, then switch to matplotlib for publication
- **Cytoscape for publication** - py4cytoscape provides the finest control over styling and exports; requires Cytoscape desktop running
- **Edge bundling** - For dense networks, reduce edge alpha (0.1-0.2) and width to prevent visual clutter
- **Community detection** - greedy_modularity_communities works well for undirected networks; use Louvain (community.best_partition) for larger graphs
- **Color palette** - Use colorblind-friendly palettes; the Nature palette (#E64B35, #4DBBD5, #00A087, #3C5488) works well for up to 8 communities
- **Export resolution** - Use dpi=300 for publication PNGs; use PDF/SVG for scalable vector output

## Related Skills

- data-visualization/multipanel-figures - Combine network panels with other plots
- gene-regulatory-networks/coexpression-networks - Build co-expression networks to visualize
- database-access/interaction-databases - Fetch PPI data for network construction
