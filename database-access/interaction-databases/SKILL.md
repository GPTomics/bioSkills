---
name: bio-interaction-databases
description: Query protein-protein and gene interaction databases (STRING, BioGRID, IntAct, SIGNOR, Reactome, HuRI, HuMAP, OmniPath, ConsensusPathDB, DIP). Use when building PPI networks, choosing between physical vs functional vs genetic interactions, signed/directed vs undirected, high-throughput vs curated, picking confidence thresholds, aggregating across resources, or navigating license constraints. Encodes the database decision matrix, STRING v12 channel semantics, OmniPath as meta-database, SIGNOR for signed signaling, and per-resource rate limits.
tool_type: python
primary_tool: requests
---

## Version Compatibility

Reference examples tested with: requests 2.31+, pandas 2.2+, networkx 3.2+; STRING v12.0, BioGRID 4.4+, IntAct (live), SIGNOR 3.0+, OmniPath (live)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show requests pandas networkx`
- API surface: confirm endpoint URLs match each resource's current docs

STRING URL is version-pinned (`version-12-0` as of 2024); older URLs (`version-11-5`) were deprecated 2023. The Cytoscape/Cytoscape.js ecosystem uses different version semantics; check the docs for the targeted version.

# Interaction Databases

**"Get protein-protein interactions for these genes"** -> The choice of database matters more than the choice of API. Different resources index different evidence (physical binding, functional association, genetic interaction, signed signaling), with different curation pipelines (manually curated vs high-throughput vs text-mined), different species coverage, and different licenses.

The decision matrix below is the postdoc-grade view: **what question is being asked, and which resource answers it best?**

- Python: `requests.get()` against REST endpoints; `pandas` for parsing; `networkx` for graphs
- R: `STRINGdb`, `OmnipathR` (mature Bioconductor clients)
- Web: STRING, BioGRID, IntAct, SIGNOR, OmniPath, ConsensusPathDB browsers

## Required Setup

```python
import requests
import pandas as pd
import networkx as nx
from io import StringIO

CALLER = 'bioskills-2026'  # STRING + OmniPath accept caller_identity for usage attribution
```

API key requirements:
- **BioGRID**: free key required (`https://webservice.thebiogrid.org/`)
- **STRING, IntAct, SIGNOR, OmniPath, Reactome**: no key

## Decision matrix: which resource for which question?

| Question | Best resource | Why |
|---|---|---|
| "Build a network around 10 genes" | STRING (medium confidence ~400) | Comprehensive; channels combinable; good viz integration |
| "Only physically interacting proteins" | IntAct or BioGRID physical | Curated physical interactions; PSI-MI standard |
| "Signed/directed signaling (phospho, ubiq, etc.)" | **SIGNOR** | Only major DB with mechanism types and direction |
| "Functional enrichment based on co-mentioned genes" | STRING functional (default) | Includes textmining channel |
| "Genetic interactions (synthetic lethality)" | BioGRID genetic | Largest curated genetic interaction set |
| "High-throughput Y2H interactome" | **HuRI** | Reference yeast-2-hybrid map of human |
| "Mass-spec-derived protein complexes" | **HuMAP v2** or BioPlex | AP-MS complex maps |
| "Curated pathways with interactions" | Reactome | Pathway-organized; gold standard for signaling |
| "Meta-database aggregating 100+ sources" | **OmniPath** | The modern "one-stop"; pre-aggregated |
| "Cross-species or non-human" | STRING | Species coverage broadest |
| "Bacterial interactome" | STRING bacterial | Limited curated alternatives |
| "Phosphorylation site-specific" | PhosphoSitePlus (commercial license) or SIGNOR | PSP has best PTM coverage but requires license |

**OmniPath** (Türei et al. 2021 *Mol Sys Biol* 17:e9923) deserves special mention: it aggregates >100 sources into a unified API with provenance tracking. For "give me all available interactions for X" workflows, OmniPath is the modern default.

## STRING (v12, channels, confidence)

STRING (Szklarczyk et al. 2023 *Nucleic Acids Res* 51:D638) aggregates evidence into a combined confidence score. The seven evidence **channels**:

| Channel | What it captures |
|---|---|
| `experiments` | Direct experimental evidence (BioGRID, IntAct, etc.) |
| `database` | Curated database (Reactome, KEGG) |
| `textmining` | Co-mention in PubMed |
| `coexpression` | Co-expression across conditions |
| `neighborhood` | Genomic neighborhood (prokaryotes mainly) |
| `fusion` | Gene fusion across species |
| `cooccurrence` | Phylogenetic profile co-occurrence |

The default `score` is the combined-evidence score (0-1000). Confidence tiers:

| Threshold | Tier | Use when |
|---|---|---|
| 150 | Low | Exploration; includes weak textmining |
| 400 | Medium | Default; balanced sensitivity/specificity |
| 700 | High | Publication-quality networks |
| 900 | Highest | Experimentally validated core only |

**Version pinning**: STRING URL is `https://version-12-0.string-db.org/api/...` as of 2024. Older `version-11-5` URLs were deprecated 2023 — code using them silently fails. Use the unversioned `https://string-db.org/api/...` to follow the current release, or pin to a specific version for reproducibility.

**`caller_identity`**: STRING requests all programmatic users pass `caller_identity=<app-name>` for usage attribution. The parameter helps STRING identify and contact heavy automated callers.

## BioGRID (physical + genetic, curated + HT)

BioGRID (Oughtred et al. 2021 *Protein Sci* 30:187) covers physical and genetic interactions across the broadest organism set of any major resource. Requires a **free API key** (`https://webservice.thebiogrid.org/`).

Key concepts:
- `EXPERIMENTAL_SYSTEM`: e.g. "Two-hybrid", "Affinity Capture-MS", "Synthetic Growth Defect"
- `THROUGHPUT`: "Low Throughput" vs "High Throughput" — the most important quality flag
- `EVIDENCE_TYPE`: physical vs genetic

For high-confidence physical interactions, filter to physical systems + Low Throughput:

```python
PHYSICAL_LT_SYSTEMS = {
    'Affinity Capture-MS', 'Affinity Capture-Western', 'Affinity Capture-RNA',
    'Co-fractionation', 'Co-purification', 'Reconstituted Complex',
    'Co-crystal Structure', 'Two-hybrid', 'Far Western', 'FRET', 'PCA',
}
```

## IntAct (PSI-MI curated physical)

IntAct (Del Toro et al. 2022 *Nucleic Acids Res* 50:D648) is the IMEx consortium reference for curated physical interactions in PSI-MI standard format. MINT was folded in ~2014; queries to MINT URLs now redirect to IntAct.

Direct REST API has changed multiple times; the most stable access is via OmniPath (which wraps IntAct) or via the PSICQUIC web services.

## SIGNOR (signed signaling)

SIGNOR (Lo Surdo et al. 2023 *Nucleic Acids Res* 51:D631) is **the only major curated database with signed, directed, mechanism-typed interactions**. Each edge has:

- `direction`: A->B
- `effect`: `up-regulates`, `down-regulates`, `unknown`
- `mechanism`: `phosphorylation`, `dephosphorylation`, `ubiquitination`, `binding`, etc.

Essential for any signaling pathway analysis or dynamic modeling. Coverage smaller than STRING/BioGRID but quality is high.

## Reactome (curated pathways)

Reactome (Milacic et al. 2024 *Nucleic Acids Res* 52:D672) is gold-standard for human pathway curation with full interaction reactions. Species-specific (human is by far the most complete).

Reactome ContentService REST: `https://reactome.org/ContentService/`.

## HuRI / HuMAP (human-specific interactomes)

- **HuRI** (Luck et al. 2020 *Nature* 580:402): yeast-2-hybrid interactome of human; ~53K binary interactions; biased toward binary high-confidence.
- **HuMAP v2** (Drew et al. 2021 *Mol Syst Biol* 17:e10016): AP-MS-derived complex map; integrates >15,000 proteomic experiments.

Both available for download; HuRI also has a web portal (`http://www.interactome-atlas.org/`).

## OmniPath (meta-database)

OmniPath (Türei et al. 2021 *Mol Syst Biol* 17:e9923) aggregates 100+ sources with provenance. Key endpoints:

| Endpoint | Content |
|---|---|
| `/interactions` | Signaling interactions (directed); the most useful for cross-DB consensus |
| `/enzsub` | Enzyme-substrate (kinase-substrate, etc.) |
| `/complexes` | Protein complexes |
| `/annotations` | Functional annotations |
| `/intercell` | Intercellular communication |

Each interaction includes `sources` (list of contributing databases) and `references` (PMIDs). For "give me everything anyone has said about A-B", OmniPath is the answer.

R users: `OmnipathR` (Bioconductor) is more ergonomic than raw HTTP.

## License gotchas (critical for commercial use)

| Resource | License | Notes |
|---|---|---|
| STRING | Free for all use | Permissive |
| BioGRID | Free, but registration required for bulk | Academic and commercial |
| IntAct | CC-BY (PSI-MI) | Permissive |
| SIGNOR | CC-BY-SA | Share-alike |
| Reactome | CC-BY | Permissive |
| HuRI / HuMAP | CC-BY | Permissive |
| OmniPath | Per-source (mostly permissive) | Check individual sources for commercial use |
| ConsensusPathDB | **Academic only** | Cannot use commercially |
| PhosphoSitePlus | **Commercial license required** | Best PTM coverage but costly |
| Pathway Commons | Per-source | Check sources |

For commercial pipelines, **stick to STRING + BioGRID + IntAct + SIGNOR + Reactome + HuRI/HuMAP + OmniPath** with appropriate source attribution. Avoid ConsensusPathDB and PhosphoSitePlus without legal review.

## Code patterns

### STRING network with confidence threshold and channel inspection

**Goal:** Build a network of interactions among a gene list at a stated confidence; surface which channels contribute.

**Approach:** REST `/network` endpoint; parse channel-specific scores from response columns.

**Reference (STRING v12.0):**
```python
import requests
import pandas as pd
from io import StringIO

STRING = 'https://version-12-0.string-db.org/api'


def get_string_network(genes, species=9606, threshold=700):
    '''threshold: 150 (low), 400 (medium), 700 (high), 900 (highest).'''
    url = f'{STRING}/tsv/network'
    params = {
        'identifiers': '%0d'.join(genes),
        'species': species,
        'required_score': threshold,
        'caller_identity': 'bioskills-2026',
    }
    r = requests.get(url, params=params); r.raise_for_status()
    df = pd.read_csv(StringIO(r.text), sep='\t')
    # Columns: stringId_A, stringId_B, preferredName_A, preferredName_B, ncbiTaxonId,
    # score, nscore (neighborhood), fscore (fusion), pscore (cooccurrence),
    # ascore (coexpression), escore (experiments), dscore (database), tscore (textmining)
    return df


genes = ['TP53', 'BRCA1', 'MDM2', 'ATM', 'CHEK2', 'CDK2']
df = get_string_network(genes, threshold=700)
print(f'{len(df)} interactions at score >= 700')
print('Channel composition for first 5 edges:')
print(df[['preferredName_A', 'preferredName_B', 'score',
          'escore', 'dscore', 'tscore', 'ascore']].head())
```

### BioGRID physical interactions, low-throughput only

**Reference (BioGRID 4.4+):**
```python
BIOGRID = 'https://webservice.thebiogrid.org/interactions/'

# Use the PHYSICAL_LT_SYSTEMS set defined above (the full curated set).
# Importing or re-defining is equivalent; using the same constant keeps the filter consistent.


def biogrid_lt_physical(gene, api_key, taxon=9606):
    params = {
        'accesskey': api_key,
        'format': 'json',
        'searchNames': True,
        'geneList': gene,
        'taxId': taxon,
        'includeInteractors': True,
        'max': 10000,
    }
    r = requests.get(BIOGRID, params=params); r.raise_for_status()
    data = r.json()
    rows = []
    for v in data.values():
        if v['THROUGHPUT'] == 'Low Throughput' and v['EXPERIMENTAL_SYSTEM'] in PHYSICAL_LT_SYSTEMS:
            rows.append({
                'gene_a': v['OFFICIAL_SYMBOL_A'],
                'gene_b': v['OFFICIAL_SYMBOL_B'],
                'system': v['EXPERIMENTAL_SYSTEM'],
                'pmid': v['PUBMED_ID'],
            })
    return pd.DataFrame(rows)
```

### SIGNOR signed signaling

```python
SIGNOR = 'https://signor.uniroma2.it/getData.php'


def signor_for_gene(gene_symbol):
    '''Return signed, directed signaling interactions involving the gene.'''
    params = {'organism': 'human', 'entity': gene_symbol}
    r = requests.get(SIGNOR, params=params); r.raise_for_status()
    rows = []
    for line in r.text.strip().split('\n')[1:]:
        cols = line.split('\t')
        if len(cols) >= 8:
            rows.append({
                'source': cols[0], 'target': cols[1], 'effect': cols[2],
                'mechanism': cols[3], 'pmid': cols[7],
            })
    return pd.DataFrame(rows)
```

### OmniPath interactions (meta-database)

```python
OMNI = 'https://omnipathdb.org'


def omnipath_interactions(genes, types='post_translational'):
    '''types: post_translational, transcriptional, mirna_target, lncrna_target.'''
    params = {
        'genesymbols': 1,
        'fields': 'sources,references,curation_effort,n_resources',
        'partners': ','.join(genes),
        'types': types,
        'license': 'academic',  # or 'commercial' for permissive-only sources
    }
    r = requests.get(f'{OMNI}/interactions', params=params); r.raise_for_status()
    df = pd.read_csv(StringIO(r.text), sep='\t')
    return df


df = omnipath_interactions(['TP53', 'MDM2', 'BRCA1'])
# Rich metadata: directionality, sign, sources, references, curation effort
print(df[['source_genesymbol', 'target_genesymbol', 'is_directed', 'is_stimulation',
          'is_inhibition', 'n_resources', 'n_references']].head())
```

### Multi-resource aggregation

**Goal:** Build a union network of interactions from multiple resources; track provenance per edge.

**Approach:** Query each resource; normalize gene symbols; merge into a networkx Graph with `sources` attribute per edge.

**Reference (requests 2.31+, networkx 3.2+):**
```python
import networkx as nx


def aggregate_networks(genes, biogrid_key=None):
    g = nx.Graph()

    # STRING (high confidence)
    string_df = get_string_network(genes, threshold=700)
    for _, row in string_df.iterrows():
        a, b = sorted([row['preferredName_A'], row['preferredName_B']])
        edge = g.get_edge_data(a, b, default={'sources': set(), 'max_score': 0})
        edge['sources'].add('STRING')
        edge['max_score'] = max(edge['max_score'], row['score'] / 1000.0)
        g.add_edge(a, b, **edge)

    # OmniPath
    omni_df = omnipath_interactions(genes)
    for _, row in omni_df.iterrows():
        a, b = sorted([row['source_genesymbol'], row['target_genesymbol']])
        edge = g.get_edge_data(a, b, default={'sources': set(), 'max_score': 0})
        edge['sources'].add('OmniPath')
        g.add_edge(a, b, **edge)

    # BioGRID (if key available)
    if biogrid_key:
        for gene in genes:
            biogrid_df = biogrid_lt_physical(gene, biogrid_key)
            for _, row in biogrid_df.iterrows():
                a, b = sorted([row['gene_a'], row['gene_b']])
                edge = g.get_edge_data(a, b, default={'sources': set(), 'max_score': 0})
                edge['sources'].add('BioGRID-LT-physical')
                g.add_edge(a, b, **edge)

    return g
```

### Network statistics

```python
def summary(g):
    return {
        'nodes': g.number_of_nodes(),
        'edges': g.number_of_edges(),
        'density': nx.density(g),
        'components': nx.number_connected_components(g),
        'mean_degree': sum(dict(g.degree()).values()) / max(g.number_of_nodes(), 1),
    }


# Edges supported by multiple resources are higher-confidence
high_conf = [(a, b, d) for a, b, d in g.edges(data=True) if len(d['sources']) >= 2]
print(f'Multi-source edges: {len(high_conf)}/{g.number_of_edges()}')
```

## Failure modes

### Confusing functional with physical interactions
- **Trigger:** Treating STRING's combined score as "physically interact".
- **Mechanism:** STRING aggregates seven channels; textmining and coexpression are functional, not physical.
- **Symptom:** Network includes co-mentioned but non-interacting proteins.
- **Fix:** Filter STRING to `experiments` channel only (`escore > threshold`); or use BioGRID/IntAct for strictly physical.

### Wrong confidence tier for the use case
- **Trigger:** Default STRING threshold 400 for a publication network.
- **Mechanism:** Includes weak textmining hits.
- **Symptom:** Network has many low-quality edges; downstream stats inflated.
- **Fix:** Use threshold 700 (high) for publication; 900 for experimentally-validated-core.

### High-throughput interactions trusted as low-throughput
- **Trigger:** Treating Y2H or AP-MS bulk screens like curated low-throughput evidence.
- **Mechanism:** HT screens have higher false-positive rates.
- **Symptom:** Spurious interactions; network overfit to specific screens.
- **Fix:** Filter `THROUGHPUT = 'Low Throughput'` in BioGRID; or use IntAct with curated MI scores.

### Symbol drift
- **Trigger:** Using `MARCH1` for a query; renamed to `MARCHF1` in 2020.
- **Mechanism:** HGNC renamed Excel-autocorrect-affected genes; some resources updated, some didn't.
- **Symptom:** Empty results; or matches to wrong gene.
- **Fix:** Resolve symbols to HGNC IDs first via UniProt or Ensembl; query by ID when possible.

### STRING version drift
- **Trigger:** Code using `version-11-5.string-db.org`.
- **Mechanism:** Deprecated 2023; v12 is the current release.
- **Symptom:** 404 or silent return of stale data.
- **Fix:** Use `version-12-0` (pinned for reproducibility) or `string-db.org` (live).

### License surprise
- **Trigger:** Building a commercial product using ConsensusPathDB.
- **Mechanism:** ConsensusPathDB is academic-only.
- **Symptom:** License violation downstream.
- **Fix:** Audit each resource's license before commercial use; OmniPath has a `license=commercial` parameter that filters to commercially-permissive sources.

### Asymmetric / directional confusion
- **Trigger:** Treating SIGNOR or OmniPath directional edges as undirected.
- **Mechanism:** Signaling edges carry direction and sign; collapsing loses information.
- **Symptom:** Wrong network topology in pathway analysis.
- **Fix:** Use DiGraph (nx.DiGraph) for directed resources; preserve `effect` and `mechanism` attributes.

### Rate limit on bulk STRING queries
- **Trigger:** Looping STRING /network for 1000 gene sets.
- **Mechanism:** STRING asks for one request at a time per `caller_identity`.
- **Symptom:** Connection errors; throttling.
- **Fix:** Sleep 1-2 seconds between calls; respect `caller_identity` rules; for very large batches, use the bulk download.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| STRING 404 | Deprecated `version-11-5` URL | Use `version-12-0` or unversioned |
| BioGRID empty result | Missing API key or wrong taxId | Get key; use NCBI taxon ID |
| Symbol mismatch | HGNC renaming | Resolve via UniProt/Ensembl ID |
| HT interactions inflate network | No throughput filter | Filter `THROUGHPUT = 'Low Throughput'` |
| Functional vs physical confusion | Mixed STRING channels | Filter to `escore` for physical |
| Directional edges collapsed | Used Graph for directed source | Use DiGraph for SIGNOR/OmniPath |
| License violation in commercial pipeline | ConsensusPathDB or PhosphoSitePlus | Switch to permissive sources |
| OmniPath returns nothing | `license=commercial` filter too strict | Drop the filter for academic use |

## References

- Szklarczyk D, Kirsch R, Koutrouli M, et al. (2023) The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Res* 51:D638-D646.
- Oughtred R, Rust J, Chang C, et al. (2021) The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions. *Protein Sci* 30:187-200.
- Del Toro N, Shrivastava A, Ragueneau E, et al. (2022) The IntAct database: efficient access to fine-grained molecular interaction data. *Nucleic Acids Res* 50:D648-D653.
- Lo Surdo P, Iannuccelli M, Contino S, et al. (2023) SIGNOR 3.0, the SIGnaling network open resource 3.0: 2022 update. *Nucleic Acids Res* 51:D631-D637.
- Milacic M, Beavers D, Conley P, et al. (2024) The Reactome Pathway Knowledgebase 2024. *Nucleic Acids Res* 52:D672-D678.
- Luck K, Kim DK, Lambourne L, et al. (2020) A reference map of the human binary protein interactome. *Nature* 580:402-408.
- Drew K, Wallingford JB, Marcotte EM. (2021) hu.MAP 2.0: integration of over 15,000 proteomic experiments builds a global compendium of human multiprotein assemblies. *Mol Syst Biol* 17:e10016.
- Türei D, Valdeolivas A, Gül L, et al. (2021) Integrated intra- and intercellular signaling knowledge for multicellular omics analysis. *Mol Syst Biol* 17:e9923.

## Related Skills

- uniprot-access - Resolve symbols to UniProt accessions
- ensembl-rest - Cross-reference Ensembl IDs in network nodes
- gene-regulatory-networks/coexpression-networks - Co-expression as a complement to PPI
- pathway-analysis/go-enrichment - Functional enrichment of network genes
- pathway-analysis/reactome-pathways - Use Reactome pathways alongside Reactome interactions
- data-visualization/network-visualization - Visualize the resulting networks
