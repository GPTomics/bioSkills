'''Static network visualization with NetworkX and matplotlib'''
# Reference: matplotlib 3.8+, numpy 1.26+ | Verify API if version differs

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- ALTERNATIVE: Load real interaction data ---
# Load from STRING or other PPI database:
#
# interactions = pd.read_csv('string_interactions.csv')
# G = nx.Graph()
# for _, row in interactions.iterrows():
#     G.add_edge(row['gene_a'], row['gene_b'], weight=row['score'] / 1000)

np.random.seed(42)

# Simulate a PPI-like network with hub structure
genes = [f'Gene{i}' for i in range(1, 31)]
G = nx.barabasi_albert_graph(30, 2, seed=42)
G = nx.relabel_nodes(G, {i: genes[i] for i in range(30)})
for u, v in G.edges():
    # Simulated confidence score. STRING scores range 0-1000.
    G[u][v]['score'] = np.random.randint(400, 1000)
    G[u][v]['weight'] = G[u][v]['score'] / 1000

pos = nx.spring_layout(G, seed=42, k=1.5)

# --- Plot 1: Degree-sized nodes with colorbar ---
fig, ax = plt.subplots(figsize=(12, 10))

degrees = dict(G.degree())
node_sizes = [100 + degrees[n] * 150 for n in G.nodes()]
node_colors = [degrees[n] for n in G.nodes()]
edge_weights = [G[u][v]['weight'] * 2 for u, v in G.edges()]

nx.draw_networkx_edges(G, pos, width=edge_weights, alpha=0.3, edge_color='gray', ax=ax)
nodes = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors,
                               cmap=plt.cm.YlOrRd, edgecolors='black', linewidths=0.5, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)

plt.colorbar(nodes, ax=ax, label='Degree', shrink=0.8)
ax.set_title(f'PPI Network ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)')
ax.axis('off')
plt.tight_layout()
plt.savefig('network_degree.png', dpi=300, bbox_inches='tight')
plt.close()

# --- Plot 2: Community detection coloring ---
from networkx.algorithms.community import greedy_modularity_communities

communities = list(greedy_modularity_communities(G))
node_to_community = {}
for i, comm in enumerate(communities):
    for node in comm:
        node_to_community[node] = i

palette = plt.cm.Set2
community_colors = [node_to_community[n] for n in G.nodes()]

fig, ax = plt.subplots(figsize=(12, 10))
nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color='gray', ax=ax)
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=community_colors,
                       cmap=palette, edgecolors='black', linewidths=0.5, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)

for i, comm in enumerate(communities):
    ax.scatter([], [], c=[palette(i)], s=80, label=f'Community {i+1} ({len(comm)} nodes)')
ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
ax.set_title(f'Community Structure ({len(communities)} communities)')
ax.axis('off')
plt.tight_layout()
plt.savefig('network_communities.png', dpi=300, bbox_inches='tight')
plt.close()

# --- Plot 3: Highlight specific genes ---
# Top 5 hub genes by degree
hub_genes = set(sorted(degrees, key=degrees.get, reverse=True)[:5])

fig, ax = plt.subplots(figsize=(12, 10))
colors = ['#E64B35' if n in hub_genes else '#cccccc' for n in G.nodes()]
sizes = [800 if n in hub_genes else 200 for n in G.nodes()]

nx.draw_networkx_edges(G, pos, alpha=0.15, edge_color='gray', ax=ax)
nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color=colors, edgecolors='black', linewidths=0.5, ax=ax)
labels = {n: n for n in hub_genes}
nx.draw_networkx_labels(G, pos, labels=labels, font_size=10, font_weight='bold', ax=ax)
ax.set_title(f'Hub Genes Highlighted (top 5 by degree)')
ax.axis('off')
plt.tight_layout()
plt.savefig('network_hubs.png', dpi=300, bbox_inches='tight')
plt.close()

# --- Plot 4: Edge confidence coloring ---
fig, ax = plt.subplots(figsize=(12, 10))

# Score >= 900: highest confidence (red), >= 700: high (orange), < 700: medium (gray)
edge_colors = ['#E64B35' if G[u][v]['score'] >= 900
               else '#F39B7F' if G[u][v]['score'] >= 700
               else '#cccccc' for u, v in G.edges()]

nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=1.5, alpha=0.7, ax=ax)
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='#4DBBD5',
                       edgecolors='black', linewidths=0.5, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)

ax.scatter([], [], c='#E64B35', s=60, label='Score >= 900')
ax.scatter([], [], c='#F39B7F', s=60, label='Score >= 700')
ax.scatter([], [], c='#cccccc', s=60, label='Score < 700')
ax.legend(loc='upper left', fontsize=9, framealpha=0.9, title='Confidence')
ax.set_title('PPI Network Colored by Interaction Confidence')
ax.axis('off')
plt.tight_layout()
plt.savefig('network_confidence.png', dpi=300, bbox_inches='tight')
plt.close()

# --- Summary statistics ---
print(f'Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')
print(f'Density: {nx.density(G):.3f}')
print(f'Communities: {len(communities)}')
print(f'Avg clustering coefficient: {nx.average_clustering(G):.3f}')

degree_df = pd.DataFrame(sorted(G.degree(), key=lambda x: x[1], reverse=True), columns=['gene', 'degree'])
print('\nTop 10 by degree:')
print(degree_df.head(10).to_string(index=False))

print('\nPlots saved: network_degree.png, network_communities.png, network_hubs.png, network_confidence.png')
