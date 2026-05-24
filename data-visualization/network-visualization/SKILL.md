---
name: bio-data-visualization-network-visualization
description: Visualize biological networks (PPI, gene-regulatory, co-expression, pathway) with layout algorithm choice (ForceAtlas2, Fruchterman-Reingold, Kamada-Kawai, hive plots), edge bundling, community-based coloring, and reproducible seeds using NetworkX, PyVis, igraph, and Cytoscape automation. Use when rendering biological networks for static publication, interactive HTML exploration, or Cytoscape-format export.
tool_type: python
primary_tool: NetworkX
---

## Version Compatibility

Reference examples tested with: networkx 3.2+, igraph 0.10+ (Python and R), pyvis 0.3+, py4cytoscape 1.9+, matplotlib 3.8+, datashader 0.16+ (for large-graph rasterization).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Network Visualization

**"Plot a biological network"** -> Select a layout algorithm (force-directed for general; hive plot for comparative; ForceAtlas2 for scale-free; circular for small dense), encode node attributes (size by degree/centrality, color by community/module), and choose rendering tier (matplotlib for static publication; PyVis for interactive HTML; Cytoscape for journal-grade compositing). The dominant pitfall is treating layout as biology — node positions in force-directed plots are NOT biologically meaningful; only connectivity is.

- Python: `networkx`, `pyvis.Network`, `py4cytoscape`, `datashader` (large graphs)
- R: `igraph`, `ggraph` (ggplot2-grammar for networks)
- Desktop: Cytoscape (Shannon 2003), Gephi (ForceAtlas2 native)

## The Single Most Important Modern Insight -- Layout Is an Artifact, Not Biology

A force-directed layout (Fruchterman-Reingold, ForceAtlas2, spring) is the result of an optimization that minimizes edge crossing and balances repulsion. The visual position of a node has no biological meaning — it is determined by the layout algorithm + random initialization + iteration count + repulsion parameters.

Two consequences:
1. **Set `random_state` / `seed` for reproducibility.** Without it, the same network produces different layouts across runs.
2. **Do not read "cluster A is closer to cluster B than C" as biology.** Inter-community distances in force-directed layouts are not preserved. Only EDGE existence and node DEGREE are biological signals from the visual.

For biology-faithful layouts, use **hive plots** (Krzywinski 2012) which anchor nodes to fixed axes by metadata, OR **circular** layouts which preserve symmetry but don't claim distance meaning.

## Decision Tree by Network Type and Question

| Network | Recommended layout | Reason |
|---------|--------------------|--------|
| Generic PPI (<500 nodes) | Fruchterman-Reingold OR Kamada-Kawai | General-purpose; clean separation |
| Scale-free PPI (>500 nodes, hub-spoke) | ForceAtlas2 (Jacomy 2014) | Designed for scale-free networks |
| Gene regulatory (directed) | Hierarchical OR ForceAtlas2 with edge direction | Direction matters; hierarchical for cascade |
| Pathway / signaling | Manual or Cytoscape layout | Curated layouts in WikiPathways/Reactome |
| Co-expression module visualization | Hive plot anchored by module assignment | Comparative; nodes by category |
| Many-to-many (>10k edges) | Hierarchical edge bundling (Holten 2006) | Reduces visual clutter |
| Large network (>50k nodes) | Datashader raster + interactive zoom | matplotlib chokes; raster is the only honest display |
| Connectivity-only (no positions) | Adjacency matrix heatmap | Network as matrix avoids layout artifact |
| Comparing two networks | Side-by-side same layout (`pos` reused) | Otherwise layout differences mask biology |

## Layout Algorithms

```python
import networkx as nx

# Spring / Fruchterman-Reingold (general)
pos = nx.spring_layout(G, k=1/np.sqrt(len(G)), iterations=100, seed=42)

# Kamada-Kawai (better for small dense)
pos = nx.kamada_kawai_layout(G)

# Circular
pos = nx.circular_layout(G)

# Shell (hub at center, periphery outside)
pos = nx.shell_layout(G, nlist=[hub_nodes, periphery_nodes])

# Spectral (reveals clusters)
pos = nx.spectral_layout(G)

# Bipartite (two sets)
pos = nx.bipartite_layout(G, top_nodes)

# Hierarchical (DAG)
pos = nx.nx_pydot.graphviz_layout(G, prog='dot')   # requires graphviz
```

For ForceAtlas2 in Python: `fa2_modified` (newer maintained fork) or use Gephi for the canonical implementation. For ggraph in R:

```r
library(ggraph)
ggraph(g, layout = 'fr') +                          # Fruchterman-Reingold
    geom_edge_link(alpha = 0.3) +
    geom_node_point()

ggraph(g, layout = 'kk') +                          # Kamada-Kawai
ggraph(g, layout = 'circle') +
ggraph(g, layout = 'graphopt') +                    # OpenOrd-style for large
```

## Hive Plots (Krzywinski 2012) — Biology-Faithful

A hive plot anchors nodes to 2-3 fixed axes by a categorical attribute (e.g., node type, module, chromosome); edges drawn as arcs between axes. Removes the "hairball" effect by replacing free 2D layout with structured 1D axes.

```python
# HiveNetX or pyveplot for hive layouts
# Or use d3.js HivePlot for interactive
# R: HivePlotData via igraph + custom rendering
```

Use hive plots when comparing networks across conditions OR when nodes have a categorical structure (e.g., TFs vs targets, chromosomes for 3D-genome interactions).

## Hierarchical Edge Bundling (Holten 2006)

For many-to-many networks within a hierarchical structure (gene hierarchies, taxonomies), edge bundling routes edges along the tree backbone, dramatically reducing clutter.

```r
library(ggraph)
ggraph(graph, layout = 'dendrogram', circular = TRUE) +
    geom_conn_bundle(data = get_con(from = from_idx, to = to_idx),
                     alpha = 0.4, tension = 0.8, edge_colour = 'grey60') +
    geom_node_point() +
    theme_void()
```

## NetworkX + matplotlib — Standard Static

**Goal:** Render a PPI network with node size proportional to degree, color by community, and edge width by interaction confidence.

**Approach:** Compute layout once with fixed seed; compute attributes (degree, community); render in layers via `nx.draw_networkx_*` functions for fine control.

```python
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms.community import greedy_modularity_communities
import numpy as np

# Layout with fixed seed for reproducibility
pos = nx.spring_layout(G, k=1.5, seed=42)

# Compute attributes
degrees = dict(G.degree())
communities = list(greedy_modularity_communities(G))
node_to_community = {n: i for i, c in enumerate(communities) for n in c}

# Sizes scaled to degree
sizes = [100 + degrees[n] * 50 for n in G.nodes()]
colors = [node_to_community[n] for n in G.nodes()]

# Render in layers
fig, ax = plt.subplots(figsize=(10, 8))
nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='grey', width=0.5, ax=ax)
nodes = nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color=colors,
                                cmap='tab20', edgecolors='black', linewidths=0.5, ax=ax)
# Label only high-degree (hub) nodes
hubs = [n for n in G.nodes() if degrees[n] >= 10]
nx.draw_networkx_labels(G, pos, labels={n: n for n in hubs}, font_size=8, ax=ax)
ax.axis('off')
plt.tight_layout()
plt.savefig('network.pdf', bbox_inches='tight', dpi=300)
```

## PyVis — Interactive HTML

```python
from pyvis.network import Network

net = Network(height='700px', width='100%', bgcolor='white', font_color='black')
net.from_nx(G)

# Per-node styling
for node in G.nodes():
    net.get_node(node)['size'] = 10 + degrees[node] * 5
    net.get_node(node)['color'] = palette[node_to_community[node] % len(palette)]
    net.get_node(node)['title'] = f'{node}\nDegree: {degrees[node]}'

net.toggle_physics(True)
net.set_options('{"physics": {"forceAtlas2Based": {"gravitationalConstant": -50}}}')
net.save_graph('network.html')
```

PyVis wraps vis.js; produces standalone HTML. Suitable for supplementary HTML; not for static journal figure.

## Cytoscape Automation (py4cytoscape)

```python
import py4cytoscape as p4c
# Cytoscape desktop must be running

p4c.create_network_from_networkx(G, title='PPI')
p4c.layout_network('force-directed')

# Custom style
style_name = 'DegreeStyle'
p4c.create_visual_style(style_name)
p4c.set_node_size_mapping('degree', [1, 5, 20], [30, 60, 120],
                            mapping_type='c', style_name=style_name)
p4c.set_node_color_mapping('degree', [1, 10, 20], ['#FFFFCC', '#FD8D3C', '#BD0026'],
                             mapping_type='c', style_name=style_name)
p4c.set_visual_style(style_name)

# Export
p4c.export_image('network.pdf', type='PDF')
```

Cytoscape is the desktop reference for publication-grade biological networks; py4cytoscape exposes script control from Python or R (via cyREST).

## Per-Method Failure Modes

### Layout positions interpreted as biology

**Trigger:** "Cluster A is between cluster B and C, so it's transitional."

**Mechanism:** Force-directed positions are optimization artifacts.

**Symptom:** Conclusion contradicts orthogonal evidence; not replicable with different seed.

**Fix:** Frame conclusions in terms of edge existence and node degree only. For trajectory claims, use the relevant time-series tool (RNA velocity, pseudotime), not the network layout.

### Layout differs across runs

**Trigger:** No random seed set.

**Mechanism:** Spring / FA2 are stochastic.

**Symptom:** Rerun produces a visibly different figure.

**Fix:** `seed=42` (NetworkX) or `set.seed(42)` (R igraph) before layout.

### Comparing two networks with different layouts

**Trigger:** `spring_layout` run separately for two conditions.

**Mechanism:** Layouts differ; visual change conflated with biological change.

**Symptom:** Concludes "this protein moved" when only the layout moved.

**Fix:** Compute layout on the union network OR pass the same `pos` to both renders.

### Hairball — too many edges with poor layout

**Trigger:** Dense network with default force-directed; >5k edges.

**Mechanism:** Edge crossings dominate; no structure visible.

**Symptom:** Visual is a uniform dense blob.

**Fix:** Hierarchical edge bundling (Holten 2006), filter to top-confidence edges, use a hive plot, OR raster with Datashader.

### Hub labels obscure non-hub structure

**Trigger:** Labeling every node in a network with >100 nodes.

**Mechanism:** Labels overlap; visual clutter.

**Symptom:** Cannot read any labels; figure too busy.

**Fix:** Label only hubs (degree > threshold) OR genes of interest. Use ggrepel-style repulsion in matplotlib via adjustText.

### Edge widths uniform when weights are meaningful

**Trigger:** Default `width=1` for all edges.

**Mechanism:** Edge attribute (correlation, confidence, weight) not encoded.

**Symptom:** Reader cannot tell strong from weak interactions.

**Fix:** `width = [G[u][v]['weight'] for u, v in G.edges()]` with normalization to visible range.

### PyVis HTML size explodes for large networks

**Trigger:** `net.from_nx(G)` with 10000+ nodes.

**Mechanism:** Embedded JavaScript file balloons; browser hangs.

**Symptom:** HTML file 100+ MB; doesn't render.

**Fix:** For large networks switch to Datashader or Cytoscape with Cytoscape.js for web; PyVis is for <2000 nodes.

## Reconciliation: When Layouts Disagree

| Pattern | Cause | Action |
|---------|-------|--------|
| Two layouts of same network look different | Different algorithm or seed | Standardize; report algorithm + seed |
| Cytoscape and NetworkX disagree | Cytoscape default = grid; NetworkX = spring | Pick one; document |
| Communities don't separate visually | Layout doesn't preserve community structure | Use spectral layout OR color-code communities; do not rely on positional separation |
| Same nodes "move" between conditions | Layout re-computed | Reuse layout from union network |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Max edges for spring layout legibility | ~2000 | Practical |
| Max nodes for PyVis HTML | ~2000 | Browser memory |
| When to bundle edges | >5000 edges or many-to-many | Holten 2006 |
| When to use Datashader | >50000 nodes or edges | Standard |
| Min degree for labeling | depends; 5-10 typical | Practical |
| Random seed | always set (42 is convention) | Reproducibility |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Layout differs across runs | No seed | Always `seed=42` |
| "Distance between clusters" interpreted | Layout artifact | Frame conclusions on edges/degree only |
| Hairball | Dense + force-directed | Bundle / hive / filter / Datashader |
| Two networks' layouts not comparable | Computed separately | Use union network layout |
| Edge widths uniform | Default | Encode weight |
| Label clutter | All nodes labeled | Hubs only |
| PyVis 100MB HTML | Too large for PyVis | Switch to Cytoscape.js / Datashader |

## References

- Csardi G, Nepusz T. 2006. The igraph software package for complex network research. *InterJournal Complex Systems* 1695.
- Fruchterman TMJ, Reingold EM. 1991. Graph drawing by force-directed placement. *Softw Pract Exp* 21(11):1129-1164.
- Hagberg A, Schult D, Swart P. 2008. Exploring network structure, dynamics, and function using NetworkX. *Proc 7th Python in Science Conference (SciPy 2008)*.
- Holten D. 2006. Hierarchical edge bundles: visualization of adjacency relations in hierarchical data. *IEEE TVCG* 12(5):741-748.
- Jacomy M, Venturini T, Heymann S, Bastian M. 2014. ForceAtlas2, a continuous graph layout algorithm for handy network visualization designed for the Gephi software. *PLoS ONE* 9(6):e98679.
- Krzywinski M, Birol I, Jones SJM, Marra MA. 2012. Hive plots—rational approach to visualizing networks. *Brief Bioinform* 13(5):627-644.
- Pedersen T. 2024. ggraph (CRAN). https://ggraph.data-imaginist.com
- Shannon P, et al. 2003. Cytoscape: a software environment for integrated models of biomolecular interaction networks. *Genome Res* 13(11):2498-2504.

## Related Skills

- gene-regulatory-networks/coexpression-networks - Build the network to visualize
- database-access/interaction-databases - Fetch PPI data
- data-visualization/multipanel-figures - Combine network with other plots
- data-visualization/color-palettes - Community / module color schemes
- single-cell/cell-communication - Cell-cell interaction networks
