'''STRING v12 network with per-channel inspection; functional-vs-physical filtering; multi-confidence comparison.'''
# Reference: requests 2.31+, pandas 2.2+, networkx 3.2+ | Verify API if version differs
import requests
import pandas as pd
import networkx as nx
from io import StringIO

STRING = 'https://version-12-0.string-db.org/api'
CALLER = 'bioskills-2026'


def string_network(genes, species=9606, threshold=400):
    url = f'{STRING}/tsv/network'
    params = {
        'identifiers': '%0d'.join(genes),
        'species': species,
        'required_score': threshold,
        'caller_identity': CALLER,
    }
    r = requests.get(url, params=params); r.raise_for_status()
    return pd.read_csv(StringIO(r.text), sep='\t')


def resolve_ids(genes, species=9606):
    url = f'{STRING}/tsv/get_string_ids'
    params = {'identifiers': '%0d'.join(genes), 'species': species, 'caller_identity': CALLER}
    r = requests.get(url, params=params); r.raise_for_status()
    return pd.read_csv(StringIO(r.text), sep='\t')


GENES = ['TP53', 'BRCA1', 'MDM2', 'ATM', 'CHEK2', 'CDK2', 'RB1', 'CDKN1A', 'BAX', 'BCL2']

print('=== Resolve gene symbols to STRING IDs ===')
print(resolve_ids(GENES)[['queryItem', 'preferredName', 'stringId']].to_string(index=False))

print('\n=== Compare confidence tiers (number of edges retained) ===')
for thr in [400, 700, 900]:
    df = string_network(GENES, threshold=thr)
    print(f'  Threshold {thr}: {len(df)} edges')

print('\n=== Per-channel breakdown at threshold 700 ===')
df = string_network(GENES, threshold=700)
# Channels: experiments (escore), database (dscore), textmining (tscore),
# coexpression (ascore), neighborhood (nscore), fusion (fscore), cooccurrence (pscore)
df_show = df[['preferredName_A', 'preferredName_B', 'score',
              'escore', 'dscore', 'tscore', 'ascore']].head(8)
print(df_show.to_string(index=False))

print('\n=== Physical-only filter (experiments channel only) ===')
# escore > 0.4 (in 0-1 scale; STRING reports 0-1 for individual channels in TSV)
physical = df[df['escore'] > 0.4]
print(f'  Combined-score 700+: {len(df)} edges total')
print(f'  Filtered to escore > 0.4 (physical-only): {len(physical)} edges')
print(f'  Use the physical subset for "physically interact" claims.')

print('\n=== Build NetworkX graph and rank by centrality ===')
g = nx.Graph()
for _, row in df.iterrows():
    g.add_edge(row['preferredName_A'], row['preferredName_B'],
               score=row['score'], physical=row['escore'])

print(f'  Nodes: {g.number_of_nodes()}; edges: {g.number_of_edges()}; '
      f'density: {nx.density(g):.3f}; components: {nx.number_connected_components(g)}')

deg = pd.DataFrame(sorted(g.degree(), key=lambda x: -x[1]), columns=['gene', 'degree'])
print('\n  Hubs (top 5 by degree):')
print(deg.head(5).to_string(index=False))
