---
name: bio-data-visualization-network-visualization
description: Visualize biological networks including gene regulatory networks, protein interaction networks, and co-expression modules using NetworkX, PyVis, and Cytoscape automation. Produces interactive and publication-quality network figures. Use when creating network diagrams from interaction data, GRN results, or co-expression modules.
tool_type: python
primary_tool: NetworkX
---

# Network Visualization

Visualize biological networks with static (matplotlib), interactive (PyVis), and publication-quality (Cytoscape) approaches.

## NetworkX + Matplotlib

### Basic Network Plot

```python
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

G = nx.karate_club_graph()

fig, ax = plt.subplots(figsize=(10, 8))
pos = nx.spring_layout(G, seed=42)
nx.draw_networkx(G, pos, node_size=300, node_color='#4DBBD5', edge_color='gray',
                 font_size=8, width=0.8, ax=ax)
ax.set_title('Network')
ax.axis('off')
plt.tight_layout()
plt.savefig('network.png', dpi=300, bbox_inches='tight')
```

### Layout Algorithms

| Layout | Function | Best For |
|--------|----------|----------|
| Spring | `nx.spring_layout(G, k=1.5, seed=42)` | General-purpose, force-directed |
| Kamada-Kawai | `nx.kamada_kawai_layout(G)` | Small-medium networks, clean separation |
| Circular | `nx.circular_layout(G)` | Showing all connections, small networks |
| Shell | `nx.shell_layout(G, nlist=[core, periphery])` | Hub-spoke topology |
| Spectral | `nx.spectral_layout(G)` | Revealing clusters |

The spring layout `k` parameter controls node spacing: increase for sparse layouts, decrease for compact. Always set `seed` for reproducibility.

### Degree-Based Node Sizing

```python
def draw_network(G, pos=None, title='', ax=None, node_cmap='YlOrRd'):
    '''Draw network with degree-proportional node sizes and colored by degree'''
    if pos is None:
        pos = nx.spring_layout(G, seed=42, k=1.5)
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    degrees = dict(G.degree())
    node_sizes = [100 + degrees[n] * 150 for n in G.nodes()]
    node_colors = [degrees[n] for n in G.nodes()]

    nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='gray', width=0.8, ax=ax)
    nodes = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors,
                                   cmap=plt.cm.get_cmap(node_cmap), edgecolors='black',
                                   linewidths=0.5, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)

    plt.colorbar(nodes, ax=ax, label='Degree', shrink=0.8)
    ax.set_title(title)
    ax.axis('off')
    return ax
```

### Edge Weight Visualization

```python
def draw_weighted_network(G, pos=None, weight_key='weight'):
    '''Draw network with edge width proportional to weight'''
    if pos is None:
        pos = nx.spring_layout(G, seed=42, k=1.5)

    fig, ax = plt.subplots(figsize=(10, 8))
    weights = [G[u][v].get(weight_key, 1.0) for u, v in G.edges()]

    # Normalize weights to [0.5, 4] for visible edge widths
    min_w, max_w = min(weights), max(weights)
    if max_w > min_w:
        scaled = [0.5 + 3.5 * (w - min_w) / (max_w - min_w) for w in weights]
    else:
        scaled = [1.0] * len(weights)

    nx.draw_networkx_edges(G, pos, width=scaled, alpha=0.5, edge_color='gray', ax=ax)
    nx.draw_networkx_nodes(G, pos, node_size=400, node_color='#4DBBD5', edgecolors='black', linewidths=0.5, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)
    ax.axis('off')
    plt.tight_layout()
    return fig, ax
```

### Community Detection Coloring

```python
from networkx.algorithms.community import greedy_modularity_communities

def draw_communities(G, pos=None):
    '''Color nodes by community membership'''
    if pos is None:
        pos = nx.spring_layout(G, seed=42, k=1.5)

    communities = list(greedy_modularity_communities(G))
    node_to_community = {}
    for i, comm in enumerate(communities):
        for node in comm:
            node_to_community[node] = i
    node_colors = [node_to_community[n] for n in G.nodes()]

    palette = plt.cm.get_cmap('Set2', len(communities))

    fig, ax = plt.subplots(figsize=(10, 8))
    nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color='gray', ax=ax)
    nx.draw_networkx_nodes(G, pos, node_size=400, node_color=node_colors, cmap=palette,
                           edgecolors='black', linewidths=0.5, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)

    for i, comm in enumerate(communities):
        ax.scatter([], [], c=[palette(i)], label=f'Community {i+1} ({len(comm)} nodes)', s=80)
    ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
    ax.set_title(f'Network Communities (modularity, {len(communities)} groups)')
    ax.axis('off')
    plt.tight_layout()
    return fig, ax
```

### Highlight Subnetwork

```python
def highlight_genes(G, highlight_nodes, pos=None, highlight_color='#E64B35', default_color='#cccccc'):
    '''Highlight specific nodes (e.g., DE genes) in a network'''
    if pos is None:
        pos = nx.spring_layout(G, seed=42, k=1.5)

    colors = [highlight_color if n in highlight_nodes else default_color for n in G.nodes()]
    sizes = [600 if n in highlight_nodes else 200 for n in G.nodes()]
    font_sizes = {n: 9 if n in highlight_nodes else 0 for n in G.nodes()}

    fig, ax = plt.subplots(figsize=(10, 8))
    nx.draw_networkx_edges(G, pos, alpha=0.15, edge_color='gray', ax=ax)
    nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color=colors, edgecolors='black', linewidths=0.5, ax=ax)
    labels = {n: n for n in highlight_nodes if n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=9, font_weight='bold', ax=ax)
    ax.axis('off')
    plt.tight_layout()
    return fig, ax
```

## PyVis Interactive Networks

### Basic Interactive Plot

```python
from pyvis.network import Network

def interactive_network(G, output='network.html', height='700px', width='100%'):
    '''Create interactive HTML network with PyVis'''
    net = Network(height=height, width=width, bgcolor='white', font_color='black')
    net.from_nx(G)
    net.toggle_physics(True)
    net.show_buttons(filter_=['physics'])
    net.save_graph(output)
    return output
```

### Styled Interactive Network

```python
def styled_interactive_network(G, output='styled_network.html', community_colors=None):
    '''Interactive network with degree-based sizing and community colors'''
    net = Network(height='700px', width='100%', bgcolor='white', font_color='black')

    degrees = dict(G.degree())
    palette = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4', '#91D1C2', '#DC0000']

    for node in G.nodes():
        size = 10 + degrees[node] * 5
        if community_colors and node in community_colors:
            color = palette[community_colors[node] % len(palette)]
        else:
            color = '#4DBBD5'
        net.add_node(node, size=size, color=color, title=f'{node}\nDegree: {degrees[node]}')

    for u, v, data in G.edges(data=True):
        weight = data.get('weight', 1.0)
        net.add_edge(u, v, width=weight * 2, title=f'Weight: {weight:.2f}')

    net.toggle_physics(True)
    net.set_options('{"physics": {"forceAtlas2Based": {"gravitationalConstant": -50}}}')
    net.save_graph(output)
    return output
```

### PyVis Configuration Options

| Option | Values | Effect |
|--------|--------|--------|
| `toggle_physics(True)` | True/False | Enable force-directed layout |
| `show_buttons()` | filter list | UI controls for layout tuning |
| `barnes_hut()` | - | Fast layout for large networks |
| `force_atlas_2based()` | - | Good cluster separation |
| `repulsion()` | - | Uniform spacing |

## Cytoscape Automation

### py4cytoscape Basics

```python
import py4cytoscape as p4c

def send_to_cytoscape(G, title='Network', layout='force-directed'):
    '''Send NetworkX graph to Cytoscape (must be running)'''
    p4c.create_network_from_networkx(G, title=title)
    p4c.layout_network(layout)
    return p4c.get_network_suid()

def style_by_degree(network_suid=None):
    '''Apply degree-proportional node sizing in Cytoscape'''
    style_name = 'DegreeStyle'
    p4c.create_visual_style(style_name)

    p4c.set_node_size_mapping('degree', [1, 5, 20], [30, 60, 120],
                              mapping_type='c', style_name=style_name)
    p4c.set_node_color_mapping('degree', [1, 10, 20], ['#FFFFCC', '#FD8D3C', '#BD0026'],
                               mapping_type='c', style_name=style_name)
    p4c.set_visual_style(style_name)
```

### Export from Cytoscape

```python
def export_cytoscape_figure(filename='network.pdf', resolution=300):
    '''Export current Cytoscape view as publication figure'''
    if filename.endswith('.pdf'):
        p4c.export_image(filename, type='PDF')
    else:
        p4c.export_image(filename, type='PNG', resolution=resolution)
```

### Cytoscape Layout Options

| Layout | Use Case |
|--------|----------|
| `force-directed` | General-purpose, good default |
| `circular` | Small networks, all connections visible |
| `hierarchical` | Signaling cascades, directed networks |
| `grid` | Uniform placement for comparison |
| `degree-circle` | Hub nodes in center |

## Bioinformatics Styling Patterns

### PPI Network Styling

```python
def style_ppi_network(G, pos=None, score_key='score'):
    '''Style for protein-protein interaction networks'''
    if pos is None:
        pos = nx.kamada_kawai_layout(G)

    fig, ax = plt.subplots(figsize=(12, 10))
    degrees = dict(G.degree())
    node_sizes = [200 + degrees[n] * 100 for n in G.nodes()]

    edges = G.edges(data=True)
    edge_colors = ['#E64B35' if d.get(score_key, 0) >= 900
                   else '#F39B7F' if d.get(score_key, 0) >= 700
                   else '#cccccc' for _, _, d in edges]

    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=1.2, alpha=0.6, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='#4DBBD5',
                           edgecolors='black', linewidths=0.5, ax=ax)
    labels = {n: n for n in G.nodes() if degrees[n] >= 3}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_weight='bold', ax=ax)
    ax.axis('off')
    plt.tight_layout()
    return fig, ax
```

### GRN Directed Network

```python
def draw_grn(G, pos=None):
    '''Draw gene regulatory network with directed edges and TF highlighting'''
    if pos is None:
        pos = nx.spring_layout(G, seed=42, k=2.0)

    fig, ax = plt.subplots(figsize=(12, 10))

    activating = [(u, v) for u, v, d in G.edges(data=True) if d.get('regulation') == 'activation']
    repressing = [(u, v) for u, v, d in G.edges(data=True) if d.get('regulation') == 'repression']

    nx.draw_networkx_edges(G, pos, edgelist=activating, edge_color='#00A087',
                           arrows=True, arrowsize=15, width=1.5, alpha=0.7, ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=repressing, edge_color='#E64B35',
                           arrows=True, arrowsize=15, width=1.5, alpha=0.7,
                           style='dashed', ax=ax)

    tfs = [n for n in G.nodes() if G.out_degree(n) > 0]
    targets = [n for n in G.nodes() if n not in tfs]

    nx.draw_networkx_nodes(G, pos, nodelist=tfs, node_size=600, node_color='#3C5488',
                           node_shape='D', edgecolors='black', linewidths=0.5, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=targets, node_size=300, node_color='#F39B7F',
                           edgecolors='black', linewidths=0.5, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)

    ax.scatter([], [], c='#00A087', marker='>', s=60, label='Activation')
    ax.scatter([], [], c='#E64B35', marker='>', s=60, label='Repression')
    ax.scatter([], [], c='#3C5488', marker='D', s=60, label='TF')
    ax.scatter([], [], c='#F39B7F', marker='o', s=60, label='Target')
    ax.legend(loc='upper left', fontsize=9, framealpha=0.9)
    ax.axis('off')
    plt.tight_layout()
    return fig, ax
```

### Multi-Panel Network Comparison

```python
def compare_networks(graphs, titles, ncols=2):
    '''Side-by-side network comparison'''
    nrows = (len(graphs) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 6 * nrows))
    axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]

    for ax, G, title in zip(axes, graphs, titles):
        pos = nx.spring_layout(G, seed=42, k=1.5)
        degrees = dict(G.degree())
        sizes = [100 + degrees[n] * 100 for n in G.nodes()]
        nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color='gray', ax=ax)
        nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color='#4DBBD5',
                               edgecolors='black', linewidths=0.5, ax=ax)
        nx.draw_networkx_labels(G, pos, font_size=6, ax=ax)
        ax.set_title(f'{title}\n({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)', fontsize=10)
        ax.axis('off')

    for ax in axes[len(graphs):]:
        ax.set_visible(False)
    plt.tight_layout()
    return fig
```

## Export Formats

| Format | Method | Use Case |
|--------|--------|----------|
| PNG | `plt.savefig('net.png', dpi=300)` | Presentations, web |
| PDF | `plt.savefig('net.pdf')` | Publication figures |
| SVG | `plt.savefig('net.svg')` | Editable vector graphics |
| HTML | `net.save_graph('net.html')` | Interactive sharing (PyVis) |
| GraphML | `nx.write_graphml(G, 'net.graphml')` | Import to Cytoscape |
| GML | `nx.write_gml(G, 'net.gml')` | Import to Gephi |

## Related Skills

- data-visualization/multipanel-figures - Combine network panels with other plots
- gene-regulatory-networks/coexpression-networks - Build co-expression networks to visualize
- database-access/interaction-databases - Fetch PPI data for network construction
