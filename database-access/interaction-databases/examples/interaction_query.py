'''Aggregate interactions from STRING and OmniPath into a unified network'''

import requests
import pandas as pd
import networkx as nx
from io import StringIO

STRING_API = 'https://version-12-0.string-db.org/api'

# --- ALTERNATIVE: Use real gene lists ---
# de_genes = pd.read_csv('de_results.csv')
# genes = de_genes[de_genes['padj'] < 0.05]['gene'].head(30).tolist()

genes = ['TP53', 'MDM2', 'BRCA1', 'ATM', 'CHEK2', 'CDK2', 'CDKN1A', 'RB1',
         'E2F1', 'MYC', 'JUN', 'FOS', 'AKT1', 'MAPK1', 'EGFR']


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


def get_omnipath_interactions(genes):
    url = 'https://omnipathdb.org/interactions'
    params = {'genesymbols': 1, 'fields': 'sources,references', 'partners': ','.join(genes)}
    response = requests.get(url, params=params)
    return pd.read_csv(StringIO(response.text), sep='\t')


def get_string_enrichment(genes, species=9606):
    url = f'{STRING_API}/tsv/enrichment'
    params = {'identifiers': '%0d'.join(genes), 'species': species, 'caller_identity': 'bioskills'}
    response = requests.get(url, params=params)
    return pd.read_csv(StringIO(response.text), sep='\t')


# --- Query STRING ---
# Score threshold 700: High confidence for reliable network construction.
string_df = get_string_interactions(genes, score_threshold=700)
print(f'STRING: {len(string_df)} interactions (score >= 700)')

# --- Query OmniPath ---
omni_df = get_omnipath_interactions(genes)
omni_filtered = omni_df[omni_df['source_genesymbol'].isin(genes) & omni_df['target_genesymbol'].isin(genes)]
print(f'OmniPath: {len(omni_filtered)} interactions between query genes')

# --- Aggregate into unified network ---
G = nx.Graph()

for _, row in string_df.iterrows():
    a, b = row['preferredName_A'], row['preferredName_B']
    if G.has_edge(a, b):
        G[a][b]['sources'].add('STRING')
        G[a][b]['string_score'] = row['score']
    else:
        G.add_edge(a, b, sources={'STRING'}, string_score=row['score'], omnipath_directed=False)

for _, row in omni_filtered.iterrows():
    a, b = row['source_genesymbol'], row['target_genesymbol']
    if G.has_edge(a, b):
        G[a][b]['sources'].add('OmniPath')
        G[a][b]['omnipath_directed'] = True
    else:
        G.add_edge(a, b, sources={'OmniPath'}, string_score=0, omnipath_directed=True)

print(f'\nUnified network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')

# Edges found in both databases (higher confidence)
multi_source_edges = [(u, v) for u, v, d in G.edges(data=True) if len(d['sources']) > 1]
print(f'Edges in both databases: {len(multi_source_edges)}')

# --- Network statistics ---
degree_df = pd.DataFrame(sorted(G.degree(), key=lambda x: x[1], reverse=True), columns=['gene', 'degree'])
print('\nHub genes (by degree):')
print(degree_df.head(10).to_string(index=False))

betweenness = nx.betweenness_centrality(G)
top_bridges = sorted(betweenness.items(), key=lambda x: x[1], reverse=True)[:5]
print('\nBridge genes (by betweenness):')
for gene, score in top_bridges:
    print(f'  {gene}: {score:.3f}')

# --- STRING Enrichment ---
enrichment = get_string_enrichment(genes)
for category in ['Process', 'KEGG', 'Component']:
    cat_df = enrichment[enrichment['category'] == category].head(5)
    if len(cat_df) > 0:
        print(f'\nTop {category} enrichments:')
        for _, row in cat_df.iterrows():
            print(f'  {row["description"]} (FDR={row["fdr"]:.2e}, genes={row["number_of_genes"]})')

# --- Export ---
edge_list = pd.DataFrame([(u, v, d['sources'], d['string_score'])
                           for u, v, d in G.edges(data=True)],
                          columns=['gene_a', 'gene_b', 'sources', 'string_score'])
edge_list['sources'] = edge_list['sources'].apply(lambda x: ','.join(sorted(x)))
edge_list.to_csv('aggregated_interactions.csv', index=False)
print('\nAggregated interactions saved: aggregated_interactions.csv')

nx.write_graphml(G, 'aggregated_network.graphml')
print('Network saved: aggregated_network.graphml (open in Cytoscape)')
