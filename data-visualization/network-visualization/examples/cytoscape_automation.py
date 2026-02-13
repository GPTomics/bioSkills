'''Cytoscape automation via py4cytoscape for publication-quality network figures'''
# Reference: matplotlib 3.8+, numpy 1.26+ | Verify API if version differs

import networkx as nx
import numpy as np

# --- ALTERNATIVE: Load real interaction data ---
# import pandas as pd
# interactions = pd.read_csv('string_interactions.csv')
# G = nx.Graph()
# for _, row in interactions.iterrows():
#     G.add_edge(row['gene_a'], row['gene_b'], weight=row['score'] / 1000)

np.random.seed(42)

genes = ['TP53', 'MDM2', 'BRCA1', 'ATM', 'CHEK2', 'CDK2', 'RB1', 'CDKN1A',
         'BAX', 'BCL2', 'CASP3', 'CASP9', 'BID', 'AKT1', 'MAPK1']
G = nx.barabasi_albert_graph(len(genes), 2, seed=42)
G = nx.relabel_nodes(G, {i: genes[i] for i in range(len(genes))})
for u, v in G.edges():
    G[u][v]['score'] = np.random.randint(400, 1000)
    G[u][v]['weight'] = G[u][v]['score'] / 1000

# Add node attributes for styling
degrees = dict(G.degree())
for node in G.nodes():
    G.nodes[node]['degree'] = degrees[node]
    G.nodes[node]['gene_type'] = 'kinase' if node in ['ATM', 'CHEK2', 'CDK2', 'AKT1', 'MAPK1'] else 'other'

# --- Requires Cytoscape desktop application running ---
# Cytoscape must be open before running py4cytoscape commands.
# Download from https://cytoscape.org/

try:
    import py4cytoscape as p4c
    p4c.cytoscape_ping()
    print('Cytoscape is running.')
except Exception:
    print('Cytoscape is not running. Start Cytoscape desktop and re-run.')
    print('Below is the code that would execute with Cytoscape running.\n')

    print('--- py4cytoscape workflow ---')
    print('1. Create network: p4c.create_network_from_networkx(G)')
    print('2. Apply layout: p4c.layout_network("force-directed")')
    print('3. Create style: p4c.create_visual_style("PPIStyle")')
    print('4. Map node size to degree')
    print('5. Map edge width to score')
    print('6. Export: p4c.export_image("network.pdf", type="PDF")')


def send_network_to_cytoscape(G, title='PPI Network'):
    '''Send a NetworkX graph to Cytoscape and apply layout'''
    suid = p4c.create_network_from_networkx(G, title=title)
    p4c.layout_network('force-directed')
    print(f'Network created in Cytoscape (SUID: {suid})')
    return suid


def apply_ppi_style():
    '''Create and apply a PPI-appropriate visual style'''
    style_name = 'PPIStyle'
    defaults = {
        'NODE_SHAPE': 'ELLIPSE',
        'NODE_FILL_COLOR': '#4DBBD5',
        'NODE_BORDER_WIDTH': 1.0,
        'NODE_BORDER_PAINT': '#000000',
        'NODE_LABEL_FONT_SIZE': 10,
        'EDGE_STROKE_UNSELECTED_PAINT': '#999999',
        'EDGE_TARGET_ARROW_SHAPE': 'NONE',
        'NETWORK_BACKGROUND_PAINT': '#FFFFFF'
    }
    p4c.create_visual_style(style_name, defaults=defaults)

    # Node size mapped to degree: small (1) -> large (15+)
    p4c.set_node_size_mapping('degree', [1, 5, 15], [30, 60, 120],
                              mapping_type='c', style_name=style_name)

    # Node color mapped to degree: light yellow (low) -> dark red (high)
    p4c.set_node_color_mapping('degree', [1, 5, 15],
                               ['#FFFFCC', '#FD8D3C', '#BD0026'],
                               mapping_type='c', style_name=style_name)

    # Edge width mapped to interaction score
    p4c.set_edge_line_width_mapping('score', [400, 700, 1000], [1, 2, 4],
                                    mapping_type='c', style_name=style_name)

    # Node shape mapped to gene type
    p4c.set_node_shape_mapping('gene_type', ['kinase', 'other'], ['DIAMOND', 'ELLIPSE'],
                               mapping_type='d', style_name=style_name)

    p4c.set_visual_style(style_name)
    print(f'Applied style: {style_name}')


def export_figure(filename='ppi_network.pdf', resolution=300):
    '''Export the current network view'''
    if filename.endswith('.pdf'):
        p4c.export_image(filename, type='PDF')
    elif filename.endswith('.svg'):
        p4c.export_image(filename, type='SVG')
    else:
        p4c.export_image(filename, type='PNG', resolution=resolution)
    print(f'Exported: {filename}')


def apply_layout_options():
    '''Demonstrate different layout options'''
    layouts = {
        'force-directed': 'General-purpose force-directed layout',
        'circular': 'Nodes arranged in a circle',
        'hierarchical': 'Top-down layout for directed networks',
        'grid': 'Nodes on a grid for comparison'
    }
    for layout, desc in layouts.items():
        print(f'  {layout}: {desc}')

    # Apply force-directed as default
    p4c.layout_network('force-directed')


# --- Run if Cytoscape is available ---
try:
    p4c.cytoscape_ping()
    send_network_to_cytoscape(G)
    apply_ppi_style()
    export_figure('ppi_network.pdf')
    export_figure('ppi_network.png', resolution=300)
    print('\nNetwork visualization complete.')
except Exception:
    pass
