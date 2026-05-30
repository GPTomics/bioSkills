---
name: bio-ortholog-inference
description: Pull pre-computed ortholog calls from public databases (OrthoDB, Ensembl Compara, OMA browser, eggNOG, PANTHER, KEGG Orthology, HomoloGene) via their REST APIs. Use when orthologs are already curated upstream, when the question is "what is the X ortholog of Y" rather than "how to infer orthology de novo", when batch-mapping gene IDs across species, or when comparing the resources for consensus calls. Encodes confidence-level semantics, 1:1 vs 1:many vs many:many, HomoloGene deprecation, and when to defect to de novo computation.
tool_type: python
primary_tool: requests
---

## Version Compatibility

Reference examples tested with: requests 2.31+, pandas 2.2+; OrthoDB v12 API, Ensembl REST (Ensembl release 112+), OMA REST API, eggNOG 6.0+, PANTHER v18+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show requests pandas`
- API surface: confirm endpoint URLs and JSON schema match the current API docs

If endpoints return 404 or unexpected JSON, check release notes for the resource; schema migrations happen with each major version (Ensembl release is the biggest moving target).

# Ortholog Inference (Database Access)

**"What is the X ortholog of gene Y?"** -> Many ortholog resources have already done the inference at scale. Pulling their answers is faster and often more reliable than re-computing. This skill is the **database-access** view: how to query the major orthology resources programmatically, what their confidence semantics mean, and when their disagreements matter.

For **de novo orthology inference** (running OrthoFinder, SonicParanoid, OMA standalone on local proteomes), see `comparative-genomics/ortholog-inference` — that's a much deeper treatment of the computational side.

This skill is about pulling answers from:
- **OrthoDB v12** — broadest coverage (1700+ species), levels from species-specific to deep
- **Ensembl Compara** — vertebrate-focused, tree-reconciled, confidence scores
- **OMA browser** — high precision, HOG (Hierarchical Orthologous Group) framework
- **eggNOG 6.0** — pre-computed functional groups, deepest functional annotation
- **PANTHER** — protein family + ortholog calls with experimentally validated curation
- **KEGG Orthology (KO)** — pathway-centric orthologous functional units
- **HomoloGene** — deprecated since 2014, but data still queryable for legacy comparison

- Python: `requests.get()` against REST endpoints; `pandas` for parsing
- CLI: `curl` against the same endpoints; OrthoDB also has bulk downloads

## Required Setup

```python
import requests
import pandas as pd
import time
```

No API keys required for any of these resources (as of 2026), but rate limits apply — see per-resource notes below.

## Decision matrix: which resource for which question?

| Question | Resource | Why |
|---|---|---|
| Ortholog of human gene X in mouse | Ensembl Compara | Best-curated for vertebrates; confidence score per call |
| Ortholog of gene X across all 1700+ species | OrthoDB | Broadest taxonomic coverage |
| Single-copy orthologs for phylogenomics | OrthoDB at species-tree level | Pre-computed; large taxonomic groups |
| Functional annotation transfer | eggNOG-mapper or eggNOG API | OG-based functional categories |
| Pathway-centric orthology (KEGG pathways) | KEGG Orthology (KO) | KO IDs link directly to pathway maps |
| Curated function-aware orthologs | PANTHER | Smaller scope; manually curated; experiment-supported |
| Compare resource consensus | All of them + intersect | Disagreement is itself a signal |
| Plant orthology (Ensembl Plants) | Ensembl Compara (plant division) | Better than Ensembl vertebrate for plants |
| Bacterial orthology | OrthoDB or eggNOG bactNOG | Ensembl Bacteria has limited Compara coverage |
| Custom proteomes not in any database | **De novo computation** | See `comparative-genomics/ortholog-inference` |

## Per-resource API reference

### OrthoDB v12

Base URL: `https://data.orthodb.org/v12/`

Key endpoints (all GET, JSON returned):
- `/search?query=<symbol>&species=<NCBI_taxid>` — find ortholog groups by gene symbol
- `/orthologs?id=<og_id>&species=<taxid>` — get orthologs of a group at a specific level
- `/group?id=<og_id>` — full group info (sequences, evidence)
- `/tab?query=<og_id>` — tab-separated bulk dump

Levels are NCBI taxonomy IDs (e.g. 9606 = human, 40674 = Mammalia, 7742 = Vertebrata).

```python
def orthodb_search(symbol, species_taxid=9606):
    r = requests.get('https://data.orthodb.org/v12/search',
                     params={'query': symbol, 'species': species_taxid})
    r.raise_for_status()
    return r.json()['data']  # list of orthogroup IDs
```

### Ensembl Compara (via Ensembl REST)

Base URL: `https://rest.ensembl.org/`. JSON: `Accept: application/json`. Rate limit: 15 req/sec, 55,000 req/hour. Respect `Retry-After` on 429.

Key endpoints:
- `/homology/symbol/<species>/<symbol>` — orthologs of a gene by symbol
- `/homology/id/<ensembl_gene_id>` — orthologs of a gene by Ensembl ID
- `/lookup/symbol/<species>/<symbol>` — resolve symbol to Ensembl ID first
- Add `?type=orthologues` to filter to orthologs only (drop paralogs)
- Add `?target_species=<species>` to filter to one target species

```python
def ensembl_orthologs(symbol, species='human', target=None):
    url = f'https://rest.ensembl.org/homology/symbol/{species}/{symbol}'
    params = {'type': 'orthologues'}
    if target:
        params['target_species'] = target
    r = requests.get(url, params=params, headers={'Accept': 'application/json'})
    r.raise_for_status()
    homologies = r.json()['data'][0]['homologies']
    return [{
        'target_species': h['target']['species'],
        'target_id': h['target']['id'],
        'type': h['type'],  # ortholog_one2one / one2many / many2many / within_species_paralog
        'confidence': h.get('confidence'),  # 0/1; some calls lack this field
        'identity_target': h['target'].get('perc_id'),
        'identity_query': h['source'].get('perc_id'),
    } for h in homologies]
```

### OMA REST API

Base URL: `https://omabrowser.org/api/`. JSON returned. No rate-limit doc but be polite.

Key endpoints:
- `/protein/<id>/orthologs/` — orthologs of a protein (UniProt or OMA ID)
- `/hog/<hog_id>/` — Hierarchical Orthologous Group info
- `/genome/<species_code>/` — list all genomes; species codes are 5-letter (e.g. HUMAN, MOUSE)

```python
def oma_orthologs(uniprot_acc):
    r = requests.get(f'https://omabrowser.org/api/protein/{uniprot_acc}/orthologs/')
    r.raise_for_status()
    return r.json()  # list of ortholog dicts with omaid, canonicalid, taxonId
```

### eggNOG (5/6)

Base URL: `http://eggnog6.embl.de/api/` (web API; lighter than running eggNOG-mapper). Most heavy lifting still uses **eggNOG-mapper** locally (Cantalapiedra et al. 2021 *Mol Biol Evol* 38:5825) — for batch protein-set annotation, mapper > API.

For ad hoc lookup: search the eggNOG web interface for an orthogroup ID, then download the member set.

### KEGG Orthology (KO)

Base URL: `https://rest.kegg.jp/`. Returns plain text TSV by default (NOT JSON).

```python
def kegg_ko_for_gene(species_code, gene):
    '''KEGG species codes: hsa=human, mmu=mouse, dme=fly, etc.'''
    r = requests.get(f'https://rest.kegg.jp/link/ko/{species_code}:{gene}')
    r.raise_for_status()
    return [line.split('\t')[1].replace('ko:', '') for line in r.text.strip().split('\n') if line]


def kegg_orthologs(ko_id):
    r = requests.get(f'https://rest.kegg.jp/link/genes/{ko_id}')
    return [line.split('\t')[1] for line in r.text.strip().split('\n') if line]
```

KEGG license: commercial use requires a paid license; academic use is free for web/API.

### PANTHER

Base URL: `http://pantherdb.org/services/oai/pantherdb/` (note the unusual base). Has a curated, smaller scope than OrthoDB but with experimental evidence.

### HomoloGene (deprecated)

`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=homologene&id=<id>&rettype=xml` — still works; data frozen since 2014. Useful only for backward-compatibility with old pipelines.

## Confidence-level semantics

Each resource defines "confidence" differently. They are NOT directly comparable.

| Resource | Confidence field | Semantics |
|---|---|---|
| Ensembl Compara | `confidence` 0/1 | Binary; 1 = high confidence based on gene-tree topology |
| OrthoDB | `evolutionary_rate` (not a confidence per se) | Inverse proxy; lower = more conserved |
| OMA | Internal QC; not exposed as a per-call score | All calls passed precision filter |
| eggNOG | Tax-level coverage | Member counts per taxonomic level |
| PANTHER | `evidence` codes | Experimentally validated vs predicted |

**Don't average or compare confidence across resources.** Use within-resource cutoffs; for cross-resource comparison, intersect call sets.

## The orthology conjecture (and why resources disagree)

The orthology conjecture (Tatusov 1997; rigorously evaluated by Studer & Robinson-Rechavi 2009 *Trends Genet* 25:210; Altenhoff et al. 2012 *PLoS Comput Biol* 8:e1002514) — orthologs are more likely than paralogs to share function — is supported but weakly. Sub- and neo-functionalization mean a paralog can become the functional equivalent.

This is also why resources disagree. Different algorithms emphasize different evidence:
- **OMA** is strict (RBH + verification + HOG inference) — higher precision, lower recall.
- **Ensembl Compara** is tree-reconciled — best for vertebrates with deep Compara curation.
- **OrthoDB** uses broader hierarchical clustering — broader coverage, more ambiguous calls.
- **eggNOG** uses pre-computed orthogroups at fixed taxonomic levels — fast but coarser.

For high-stakes calls (publication, drug target choice), **intersect at least two resources** and inspect disagreements.

## Code patterns

### Get the human ortholog of a mouse gene (Ensembl Compara, 1:1 only)

**Goal:** Pull Compara's high-confidence 1:1 ortholog of a single mouse gene in human.

**Approach:** REST query with type filter; assert 1:1; record confidence.

**Reference (Ensembl REST, release 112+):**
```python
import requests
import time


def compara_one2one(symbol, source='mouse', target='human'):
    url = f'https://rest.ensembl.org/homology/symbol/{source}/{symbol}'
    r = requests.get(url, params={'type': 'orthologues', 'target_species': target},
                     headers={'Accept': 'application/json'})
    if r.status_code == 429:
        time.sleep(int(r.headers.get('Retry-After', '5')))
        return compara_one2one(symbol, source, target)
    r.raise_for_status()
    hits = r.json()['data'][0]['homologies']
    one2one = [h for h in hits if h['type'] == 'ortholog_one2one']
    if not one2one:
        return None
    h = one2one[0]
    return {
        'source_id': h['source']['id'],
        'target_id': h['target']['id'],
        'confidence': h.get('confidence'),
        'pid_target': h['target'].get('perc_id'),
    }
```

### Cross-resource agreement (Ensembl + OMA + OrthoDB)

**Goal:** Find orthologs agreed on by multiple resources to flag high-confidence calls.

**Approach:** Query each resource; intersect target IDs after normalizing to a common namespace (UniProt or NCBI Gene).

```python
def cross_resource_orthologs(symbol):
    '''Return target-species ortholog calls from multiple resources for cross-validation.'''
    ensembl = ensembl_orthologs(symbol, species='human')
    # (OMA/OrthoDB lookups omitted for brevity -- need namespace conversion via UniProt ID Mapping)
    return {'ensembl': ensembl}
```

### Batch ortholog table for >100 genes

**Goal:** Build a wide table of orthologs across N species for a gene list.

**Approach:** Loop with rate limit; cache responses; respect `Retry-After`.

**Reference (requests 2.31+):**
```python
def batch_ensembl_orthologs(symbols, source='human', target_species=None, sleep=0.07):
    '''sleep=0.07 keeps under the 15-req/sec ceiling with margin.'''
    rows = []
    for sym in symbols:
        try:
            orthologs = ensembl_orthologs(sym, species=source, target=target_species)
            for o in orthologs:
                rows.append({'source_symbol': sym, **o})
        except requests.HTTPError as e:
            if e.response.status_code == 429:
                wait = int(e.response.headers.get('Retry-After', '10'))
                time.sleep(wait)
            else:
                rows.append({'source_symbol': sym, 'error': str(e)})
        time.sleep(sleep)
    return pd.DataFrame(rows)


df = batch_ensembl_orthologs(['BRCA1', 'TP53', 'MYC'], target_species='mouse')
print(df[df['type'] == 'ortholog_one2one'][['source_symbol', 'target_id', 'confidence']])
```

### Pull all human-mouse 1:1 orthologs as a bulk table

For thousands of genes, prefer Ensembl BioMart bulk export — see `biomart-queries`. The REST API is fine for hundreds; BioMart wins at thousands.

### OMA HOG navigation

```python
def oma_hog_for_protein(oma_or_uniprot_id):
    r = requests.get(f'https://omabrowser.org/api/protein/{oma_or_uniprot_id}/')
    r.raise_for_status()
    return r.json().get('oma_hog_id')


def oma_hog_members(hog_id, level=None):
    url = f'https://omabrowser.org/api/hog/{hog_id}/'
    params = {'level': level} if level else {}
    r = requests.get(url, params=params)
    r.raise_for_status()
    return r.json()
```

### KEGG ortholog lookup

```python
ko_ids = kegg_ko_for_gene('hsa', '7157')  # human TP53
for ko in ko_ids:
    orthologs = kegg_orthologs(ko)
    print(f'{ko}: {len(orthologs)} orthologs across all KEGG species')
```

## Failure modes

### Resource disagreement on 1:1
- **Trigger:** Ensembl Compara says 1:1; OrthoDB says 1:many; OMA says no call.
- **Mechanism:** Different algorithmic emphases; different species coverage.
- **Symptom:** Inconsistent ortholog tables across pipeline stages.
- **Fix:** Define the authoritative resource per project; or take intersection; document the choice.

### Stale resource snapshot
- **Trigger:** Using a 2-year-old OrthoDB download or HomoloGene (frozen 2014).
- **Mechanism:** Species coverage and algorithms have improved; gene model updates.
- **Symptom:** Missing orthologs that the live database has; calling defunct ortholog IDs.
- **Fix:** Pin to a release version with date; refresh annually; for HomoloGene, treat as legacy and verify against a current resource.

### Symbol-based lookup ambiguity
- **Trigger:** `compara_one2one('MARCH1', 'human', 'mouse')` -- but MARCH1 was renamed to MARCHF1 in 2020.
- **Mechanism:** HGNC symbol renames break symbol-based lookups; APIs may return empty or wrong gene.
- **Symptom:** No orthologs found; or orthologs of the wrong gene.
- **Fix:** Resolve symbol to canonical Ensembl/HGNC ID first; use the ID-based endpoint.

### Compara confidence missing
- **Trigger:** Some Compara calls lack `confidence` (older calls; certain species pairs).
- **Mechanism:** Field is not populated for all calls.
- **Symptom:** `KeyError`; or filter drops calls that should pass.
- **Fix:** Use `.get('confidence', None)` and treat missing as unknown (not as low-confidence).

### Rate-limit cascade on bulk queries
- **Trigger:** Loop of 5000 Ensembl REST calls.
- **Mechanism:** 15 req/sec ceiling, 55K/hour; hit gives 429 with Retry-After.
- **Symptom:** Cascading retries; pipeline stalls.
- **Fix:** Sleep 0.07s between calls; check Retry-After on 429; for >5K queries use BioMart bulk export instead.

### Custom proteome not in any DB
- **Trigger:** Querying a newly sequenced species absent from all resources.
- **Mechanism:** All ortholog databases require the species to be in their pre-computed set.
- **Symptom:** No ortholog calls.
- **Fix:** Run de novo (OrthoFinder or SonicParanoid) -- see `comparative-genomics/ortholog-inference`.

### KEGG license confusion
- **Trigger:** Building a commercial product on KEGG REST.
- **Mechanism:** KEGG academic-free, commercial-paid.
- **Symptom:** License violation in a commercial pipeline.
- **Fix:** Confirm license for the use case; eggNOG and OrthoDB have more permissive licenses.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `HTTPError 429` (Ensembl) | Rate limit | Sleep per Retry-After; cap at 15 req/sec |
| Empty homologies | Symbol misspelled or stale | Resolve to Ensembl ID first |
| Missing confidence field | Older calls | `.get()` with default |
| OMA `404` | Wrong namespace (used UniProt where OMA needed OMA ID) | Use the protein lookup endpoint to resolve first |
| KEGG returns HTML | Endpoint wrong (use `rest.kegg.jp`) | Check URL; KEGG is text TSV not JSON |
| Resource disagreement | Different algorithms / coverage | Intersect; document choice |

## References

- Tatusov RL, Koonin EV, Lipman DJ. (1997) A genomic perspective on protein families. *Science* 278:631-637.
- Altenhoff AM, Studer RA, Robinson-Rechavi M, Dessimoz C. (2012) Resolving the ortholog conjecture: orthologs tend to be weakly, but significantly, more similar in function than paralogs. *PLoS Comput Biol* 8:e1002514.
- Studer RA, Robinson-Rechavi M. (2009) How confident can we be that orthologs are similar, but paralogs differ? *Trends Genet* 25:210-216.
- Kuznetsov D, Tegenfeldt F, Manni M, Seppey M, Berkeley M, Kriventseva EV, Zdobnov EM. (2023) OrthoDB v11: annotation of orthologs in the widest sampling of organismal diversity. *Nucleic Acids Res* 51:D445-D451.
- Herrero J, Muffato M, Beal K, et al. (2016) Ensembl comparative genomics resources. *Database* 2016:baw053.
- Altenhoff AM, Vesztrocy AW, Bernard C, et al. (2024) OMA orthology in 2024. *Nucleic Acids Res* 52:D513-D521.
- Hernandez-Plaza A, Szklarczyk D, Botas J, et al. (2023) eggNOG 6.0: enabling comparative genomics across 12,535 organisms. *Nucleic Acids Res* 51:D389-D394.
- Cantalapiedra CP, Hernandez-Plaza A, Letunic I, Bork P, Huerta-Cepas J. (2021) eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Mol Biol Evol* 38:5825-5829.
- Thomas PD, Ebert D, Muruganujan A, Mushayahama T, Albou LP, Mi H. (2022) PANTHER: making genome-scale phylogenetics accessible to all. *Protein Sci* 31:8-22.

## Related Skills

- comparative-genomics/ortholog-inference - De novo orthology computation (OrthoFinder, OMA standalone, SonicParanoid)
- ensembl-rest - Broader Ensembl REST workflows beyond Compara
- biomart-queries - Bulk ortholog table export via Ensembl BioMart
- uniprot-access - Resolve UniProt accessions used by OMA
- pathway-analysis/kegg-pathways - KEGG Orthology and pathway mapping
