---
name: bio-interaction-databases
description: Query protein-protein and gene interaction databases including STRING, BioGRID, and IntAct via their REST APIs and Python clients. Retrieve interaction networks, confidence scores, and functional enrichment. Use when building protein interaction networks, contextualizing gene lists with known interactions, or retrieving pathway-level interaction data.
tool_type: python
primary_tool: STRINGdb
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Interaction Databases

Query protein-protein interaction (PPI) and gene interaction databases to build and analyze biological networks.

**"Get protein interactions for a gene list"** → Query PPI databases (STRING, BioGRID, IntAct) via REST APIs to retrieve interaction partners and confidence scores.
- Python: `requests.get()` against STRING/BioGRID REST endpoints + `pandas` for parsing

## STRING API

### Query Interactions

```python
import requests
import pandas as pd
from io import StringIO

STRING_API = 'https://version-12-0.string-db.org/api'

def get_string_interactions(genes, species=9606, score_threshold=400):
    '''Fetch PPI from STRING. species: 9606=human, 10090=mouse, 7227=fly'''
    url = f'{STRING_API}/tsv/network'
    params = {
        'identifiers': '%0d'.join(genes),
        'species': species,
        'required_score': score_threshold,
        'caller_identity': 'bioskills'
    }
    response = requests.get(url, params=params)
    return pd.read_csv(StringIO(response.text), sep='\t')

genes = ['TP53', 'BRCA1', 'MDM2', 'ATM', 'CHEK2', 'CDK2']
interactions = get_string_interactions(genes, score_threshold=700)
interactions[['preferredName_A', 'preferredName_B', 'score']].head()
```

### STRING Confidence Score Tiers

| Threshold | Tier | Guidance |
|-----------|------|----------|
| 150 | Low | Includes transferred and text-mined; noisy |
| 400 | Medium | Default. Balanced sensitivity/specificity |
| 700 | High | Recommended for network analysis. Good confidence |
| 900 | Highest | Experimentally validated core. Very stringent |

Use 700+ for publication-quality networks. Use 400 for exploratory analysis where recall matters.

### Resolve Gene Identifiers

```python
def resolve_string_ids(genes, species=9606):
    '''Map gene names to STRING identifiers'''
    url = f'{STRING_API}/tsv/get_string_ids'
    params = {'identifiers': '%0d'.join(genes), 'species': species, 'caller_identity': 'bioskills'}
    response = requests.get(url, params=params)
    return pd.read_csv(StringIO(response.text), sep='\t')

resolved = resolve_string_ids(['TP53', 'BRCA1', 'MDM2'])
resolved[['queryItem', 'preferredName', 'stringId', 'annotation']]
```

### STRING Functional Enrichment

```python
def get_string_enrichment(genes, species=9606):
    '''Get GO/KEGG enrichment for a gene set via STRING'''
    url = f'{STRING_API}/tsv/enrichment'
    params = {'identifiers': '%0d'.join(genes), 'species': species, 'caller_identity': 'bioskills'}
    response = requests.get(url, params=params)
    df = pd.read_csv(StringIO(response.text), sep='\t')
    return df.sort_values('fdr')

enrichment = get_string_enrichment(genes)
enrichment[enrichment['category'] == 'Process'][['term', 'description', 'fdr', 'number_of_genes']].head(10)
```

### STRING Network Image

```python
def get_string_image(genes, species=9606):
    '''Download STRING network image as PNG'''
    url = f'{STRING_API}/image/network'
    params = {'identifiers': '%0d'.join(genes), 'species': species, 'caller_identity': 'bioskills'}
    response = requests.get(url, params=params)
    with open('string_network.png', 'wb') as f:
        f.write(response.content)

get_string_image(genes)
```

## BioGRID REST API

```python
BIOGRID_API = 'https://webservice.thebiogrid.org/interactions/'

def get_biogrid_interactions(gene, access_key, taxon=9606):
    '''Fetch interactions from BioGRID. Requires free API key from thebiogrid.org'''
    params = {
        'accesskey': access_key,
        'format': 'json',
        'searchNames': True,
        'geneList': gene,
        'taxId': taxon,
        'includeInteractors': True,
        'max': 1000
    }
    response = requests.get(BIOGRID_API, params=params)
    data = response.json()
    rows = [{'gene_a': v['OFFICIAL_SYMBOL_A'], 'gene_b': v['OFFICIAL_SYMBOL_B'],
             'system': v['EXPERIMENTAL_SYSTEM'], 'pubmed': v['PUBMED_ID'],
             'throughput': v['THROUGHPUT']} for v in data.values()]
    return pd.DataFrame(rows)

# biogrid_df = get_biogrid_interactions('TP53', access_key='YOUR_KEY')
```

### Filter BioGRID by Evidence

```python
def filter_biogrid(df, systems=None, low_throughput_only=False):
    '''Filter BioGRID by experimental system or throughput'''
    if low_throughput_only:
        df = df[df['throughput'] == 'Low Throughput']
    if systems:
        df = df[df['system'].isin(systems)]
    return df

# High-confidence physical interactions
# physical_systems = ['Affinity Capture-MS', 'Two-hybrid', 'Co-crystal Structure',
#                     'Reconstituted Complex', 'Affinity Capture-Western']
# filtered = filter_biogrid(biogrid_df, systems=physical_systems, low_throughput_only=True)
```

## IntAct / PSICQUIC

```python
def query_intact(gene, species='human'):
    '''Query IntAct via PSICQUIC REST API'''
    url = f'https://www.ebi.ac.uk/intact/ws/interactors/findInteractor/{gene}'
    params = {'species': species}
    response = requests.get(url, params=params)
    return response.json()

def query_psicquic(query):
    '''Query IntAct PSICQUIC for MITAB-formatted interactions'''
    url = 'https://www.ebi.ac.uk/intact/ws/interactors/list/resolve'
    response = requests.post(url, json={'query': query})
    return response.json()
```

## OmniPath

```python
def get_omnipath_interactions(genes):
    '''Query OmniPath for curated signaling interactions'''
    url = 'https://omnipathdb.org/interactions'
    params = {'genesymbols': 1, 'fields': 'sources,references', 'partners': ','.join(genes)}
    response = requests.get(url, params=params)
    df = pd.read_csv(StringIO(response.text), sep='\t')
    return df

omni = get_omnipath_interactions(['TP53', 'MDM2', 'BRCA1'])
omni[['source_genesymbol', 'target_genesymbol', 'is_directed', 'is_stimulation', 'is_inhibition']].head()
```

### OmniPath Endpoint Reference

| Endpoint | Content |
|----------|---------|
| `/interactions` | Signaling interactions (directed) |
| `/enzsub` | Enzyme-substrate relationships |
| `/complexes` | Protein complexes |
| `/annotations` | Functional annotations |
| `/intercell` | Intercellular communication |

## Convert to NetworkX

**Goal:** Build an analyzable graph object from a STRING interactions table for network metrics and visualization.

**Approach:** Filter interactions by confidence score, then create a NetworkX graph with edge weights derived from STRING scores.

```python
import networkx as nx

def string_to_networkx(interactions_df, score_col='score', min_score=700):
    '''Convert STRING interactions DataFrame to NetworkX graph'''
    filtered = interactions_df[interactions_df[score_col] >= min_score]
    G = nx.Graph()
    for _, row in filtered.iterrows():
        G.add_edge(row['preferredName_A'], row['preferredName_B'],
                   weight=row[score_col] / 1000, score=row[score_col])
    return G

G = string_to_networkx(interactions)
print(f'Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')

degree_df = pd.DataFrame(G.degree(), columns=['gene', 'degree']).sort_values('degree', ascending=False)
degree_df.head()
```

### Multi-Database Aggregation

**Goal:** Combine protein interaction evidence from multiple databases into a single high-confidence network.

**Approach:** Normalize edge pairs from STRING and OmniPath into a common format, merge into one graph, and track which databases support each interaction.

```python
def aggregate_interactions(string_df, omnipath_df):
    '''Combine interactions from multiple databases into a unified network'''
    edges = set()
    for _, row in string_df.iterrows():
        a, b = sorted([row['preferredName_A'], row['preferredName_B']])
        edges.add((a, b, 'STRING', row['score']))
    for _, row in omnipath_df.iterrows():
        a, b = row['source_genesymbol'], row['target_genesymbol']
        edges.add((a, b, 'OmniPath', 1000))

    G = nx.Graph()
    for a, b, source, score in edges:
        if G.has_edge(a, b):
            G[a][b]['sources'].append(source)
            G[a][b]['max_score'] = max(G[a][b]['max_score'], score)
        else:
            G.add_edge(a, b, sources=[source], max_score=score)
    return G
```

### Network Statistics

```python
def network_summary(G):
    '''Compute basic network statistics'''
    stats = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': nx.density(G),
        'components': nx.number_connected_components(G),
        'avg_clustering': nx.average_clustering(G),
    }
    if nx.is_connected(G):
        stats['diameter'] = nx.diameter(G)
        stats['avg_path_length'] = nx.average_shortest_path_length(G)
    return stats

network_summary(G)
```

## Related Skills

- database-access/uniprot-access - Protein annotations and cross-references
- pathway-analysis/go-enrichment - Functional enrichment of network genes
- gene-regulatory-networks/coexpression-networks - Co-expression network construction
- data-visualization/network-visualization - Visualize interaction networks
