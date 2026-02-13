'''Query STRING REST API for protein interactions and build a NetworkX graph'''
# Reference: biopython 1.83+, pandas 2.2+ | Verify API if version differs

import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from io import StringIO

STRING_API = 'https://version-12-0.string-db.org/api'

# --- ALTERNATIVE: Use real gene lists ---
# Load DE genes or pathway members:
#
# de_results = pd.read_csv('de_results.csv')
# genes = de_results[de_results['padj'] < 0.05]['gene'].head(50).tolist()

genes = ['TP53', 'BRCA1', 'MDM2', 'ATM', 'CHEK2', 'CDK2', 'RB1', 'CDKN1A', 'BAX', 'BCL2']

def get_string_interactions(genes, species=9606, score_threshold=400):
    url = f'{STRING_API}/tsv/network'
    params = {
        'identifiers': '%0d'.join(genes),
        'species': species,
        'required_score': score_threshold,
        'caller_identity': 'bioskills'
    }
    response = requests.get(url, params=params)
    return pd.read_csv(StringIO(response.text), sep='\t')


def resolve_string_ids(genes, species=9606):
    url = f'{STRING_API}/tsv/get_string_ids'
    params = {'identifiers': '%0d'.join(genes), 'species': species, 'caller_identity': 'bioskills'}
    response = requests.get(url, params=params)
    return pd.read_csv(StringIO(response.text), sep='\t')


# Resolve gene names to STRING IDs
resolved = resolve_string_ids(genes)
print('Resolved identifiers:')
print(resolved[['queryItem', 'preferredName', 'stringId']].to_string(index=False))

# Score threshold 700: High confidence. Filters text-mining and prediction noise.
# Use 400 for exploratory analysis, 900 for highest-confidence core interactions.
interactions = get_string_interactions(genes, score_threshold=700)
print(f'\nFound {len(interactions)} high-confidence interactions (score >= 700)')
print(interactions[['preferredName_A', 'preferredName_B', 'score']].head(10).to_string(index=False))

# Build NetworkX graph
G = nx.Graph()
for _, row in interactions.iterrows():
    G.add_edge(row['preferredName_A'], row['preferredName_B'],
               weight=row['score'] / 1000, score=row['score'])

print(f'\nNetwork: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')
print(f'Density: {nx.density(G):.3f}')
print(f'Connected components: {nx.number_connected_components(G)}')

# Degree centrality
degree_df = pd.DataFrame(sorted(G.degree(), key=lambda x: x[1], reverse=True), columns=['gene', 'degree'])
print('\nDegree ranking:')
print(degree_df.to_string(index=False))

# Betweenness centrality (identifies bridge nodes)
betweenness = nx.betweenness_centrality(G)
betweenness_df = pd.DataFrame(sorted(betweenness.items(), key=lambda x: x[1], reverse=True), columns=['gene', 'betweenness'])
print('\nBetweenness centrality:')
print(betweenness_df.to_string(index=False))

# Visualize with matplotlib
fig, ax = plt.subplots(figsize=(10, 8))
pos = nx.spring_layout(G, seed=42, k=1.5)

node_sizes = [300 + G.degree(n) * 200 for n in G.nodes()]
edge_weights = [G[u][v]['weight'] * 2 for u, v in G.edges()]

nx.draw_networkx_edges(G, pos, width=edge_weights, alpha=0.4, edge_color='gray', ax=ax)
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='#4DBBD5', edgecolors='black', linewidths=0.5, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold', ax=ax)

ax.set_title(f'STRING PPI Network (score >= 700, {G.number_of_nodes()} nodes, {G.number_of_edges()} edges)')
ax.axis('off')
plt.tight_layout()
plt.savefig('string_network.png', dpi=300, bbox_inches='tight')
plt.close()

print('\nNetwork saved: string_network.png')
