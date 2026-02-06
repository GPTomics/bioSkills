'''Interactive HTML network visualization with PyVis'''

import networkx as nx
import numpy as np
from pyvis.network import Network
from networkx.algorithms.community import greedy_modularity_communities

# --- ALTERNATIVE: Load real interaction data ---
# import pandas as pd
# interactions = pd.read_csv('string_interactions.csv')
# G = nx.Graph()
# for _, row in interactions.iterrows():
#     G.add_edge(row['gene_a'], row['gene_b'], weight=row['score'] / 1000)

np.random.seed(42)

genes = [f'Gene{i}' for i in range(1, 41)]
G = nx.barabasi_albert_graph(40, 2, seed=42)
G = nx.relabel_nodes(G, {i: genes[i] for i in range(40)})
for u, v in G.edges():
    G[u][v]['weight'] = np.random.uniform(0.4, 1.0)
    G[u][v]['score'] = int(G[u][v]['weight'] * 1000)


def create_interactive_network(G, output='network.html', title='Biological Network'):
    '''Basic interactive network with PyVis'''
    net = Network(height='700px', width='100%', bgcolor='white', font_color='black', heading=title)
    net.from_nx(G)
    net.toggle_physics(True)
    net.show_buttons(filter_=['physics'])
    net.save_graph(output)
    return output


def create_styled_network(G, output='styled_network.html', title='Styled Network'):
    '''Interactive network with degree sizing and community coloring'''
    net = Network(height='700px', width='100%', bgcolor='white', font_color='black', heading=title)

    communities = list(greedy_modularity_communities(G))
    node_to_community = {}
    for i, comm in enumerate(communities):
        for node in comm:
            node_to_community[node] = i

    # Nature-style color palette, colorblind-friendly for up to 8 communities
    palette = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4', '#91D1C2', '#DC0000']
    degrees = dict(G.degree())

    for node in G.nodes():
        comm_idx = node_to_community.get(node, 0)
        color = palette[comm_idx % len(palette)]
        size = 10 + degrees[node] * 5
        tooltip = f'{node}\nDegree: {degrees[node]}\nCommunity: {comm_idx + 1}'
        net.add_node(node, size=size, color=color, title=tooltip, label=node)

    for u, v, data in G.edges(data=True):
        weight = data.get('weight', 0.5)
        score = data.get('score', 500)
        tooltip = f'{u} -- {v}\nScore: {score}'
        net.add_edge(u, v, width=weight * 3, title=tooltip)

    options = {
        'physics': {
            'forceAtlas2Based': {'gravitationalConstant': -50, 'centralGravity': 0.01,
                                 'springLength': 100, 'springConstant': 0.08},
            'solver': 'forceAtlas2Based',
            'stabilization': {'iterations': 150}
        },
        'interaction': {'hover': True, 'navigationButtons': True, 'keyboard': True}
    }
    import json
    net.set_options(json.dumps(options))
    net.save_graph(output)
    return output


# --- Create basic interactive network ---
basic_file = create_interactive_network(G, output='network_basic.html', title='Basic Interactive Network')
print(f'Basic network saved: {basic_file}')

# --- Create styled network with communities ---
styled_file = create_styled_network(G, output='network_styled.html', title='PPI Network with Communities')
print(f'Styled network saved: {styled_file}')

# --- Summary ---
communities = list(greedy_modularity_communities(G))
print(f'\nNetwork: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')
print(f'Communities detected: {len(communities)}')
for i, comm in enumerate(communities):
    print(f'  Community {i+1}: {len(comm)} nodes')
print('\nOpen the HTML files in a browser to interact with the networks.')
