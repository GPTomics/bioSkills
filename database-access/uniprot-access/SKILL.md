---
name: bio-uniprot-access
description: Query UniProt's REST API (post-2022 endpoint at rest.uniprot.org) for protein sequences, annotations, GO terms, cross-references, ID mappings, and proteomes. Use when fetching UniProtKB entries, navigating the JSON schema, choosing between UniProtKB/UniRef/UniParc/Proteomes resources, deciding stream vs search endpoint for batch retrieval, running ID-mapping jobs with the async pattern, handling isoform suffixes, or filtering reviewed Swiss-Prot vs auto-annotated TrEMBL. Encodes the legacy URL migration (2022), the new JSON schema layout, and bulk-pull patterns.
tool_type: python
primary_tool: requests
---

## Version Compatibility

Reference examples tested with: requests 2.31+, pandas 2.2+; UniProt REST API as of 2024_06 release

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show requests pandas`
- API surface: confirm endpoint URLs match https://www.uniprot.org/help/api

The REST API JSON schema is stable within a release; major schema changes are documented at https://www.uniprot.org/release-notes. The 2022 migration broke the legacy `https://www.uniprot.org/uniprot/...` endpoints.

# UniProt Access

**"Get protein information from UniProt"** -> Two facts dominate every UniProt workflow in 2026: (1) **the API endpoint migrated in 2022** from `https://www.uniprot.org/uniprot/...` to `https://rest.uniprot.org/uniprotkb/...` with a substantially different JSON schema; pre-2022 code does not work as-is. (2) **`?fields=`** is essential — default JSON returns the full entry (~20-30 KB each); for bulk pulls, request only the fields actually needed.

The major databases under the UniProt umbrella have different scopes:

- **UniProtKB**: the curated knowledgebase — Swiss-Prot (manually reviewed, ~570K entries as of 2024) + TrEMBL (auto-annotated, ~250M). Always specify `reviewed:true` for high-quality reference work.
- **UniRef**: clustered sequences at 100%, 90%, 50% identity. UniRef50 is the standard for redundancy reduction.
- **UniParc**: archival "every unique sequence ever seen" — for provenance and historical lookup.
- **Proteomes**: organism-level groupings; reference proteomes (one per species) are the canonical subset.

- Python: `requests.get('https://rest.uniprot.org/uniprotkb/...')` (REST API)
- Python: `Bio.ExPASy.get_sprot_raw()` (BioPython; legacy SwissProt format)
- CLI: `curl https://rest.uniprot.org/uniprotkb/P04637.json`

## Required Setup

```python
import requests
import pandas as pd
import time
```

No API key required. Rate limit is generous (~200 req/sec tolerated empirically); ID-mapping has its own job queue.

## Endpoint reference

Base: `https://rest.uniprot.org/`

| Resource | Endpoint | Use |
|---|---|---|
| Single entry | `/uniprotkb/{accession}` | One protein record |
| Search | `/uniprotkb/search` | Query with up to 500 results per page |
| Stream | `/uniprotkb/stream` | No 500-result limit; for bulk |
| Batch by accession | `/uniprotkb/accessions` | Multiple specific accessions |
| ID Mapping (run) | `/idmapping/run` | Submit conversion job |
| ID Mapping (status) | `/idmapping/status/{jobId}` | Poll |
| ID Mapping (results) | `/idmapping/results/{jobId}` | Retrieve |
| UniRef entry | `/uniref/{cluster_id}` | One cluster |
| UniRef search | `/uniref/search` | UniRef cluster queries |
| Proteome | `/proteomes/{upid}` | Organism proteome |
| Proteome FASTA | `/proteomes/{upid}.fasta.gz` | Download whole proteome |
| Taxonomy | `/taxonomy/{taxid}` | Taxonomy info |

Append `.json`, `.fasta`, `.tsv`, `.xml`, `.txt`, or `.gff` to single-entry URLs to control format.

## Search query syntax

UniProt search queries use a Lucene-like syntax distinct from Entrez:

| Query | Means |
|---|---|
| `gene:TP53` | Gene name TP53 |
| `gene_exact:TP53` | Exact gene name (no wildcard match) |
| `organism_id:9606` | Human (NCBI taxonomy ID) |
| `organism_name:"Homo sapiens"` | By name (slower than taxid) |
| `reviewed:true` | Swiss-Prot only |
| `reviewed:false` | TrEMBL only |
| `length:[100 TO 500]` | Sequence length range |
| `go:0006915` | GO term (apoptosis) |
| `keyword:KW-0067` | UniProt keyword |
| `ec:2.7.1.1` | Enzyme classification |
| `database:pdb` | Has PDB cross-ref |
| `xref:pdb` | Same as above |
| `existence:1` | Evidence at protein level (1 = strongest) |

Combine: `organism_id:9606 AND reviewed:true AND keyword:KW-0067 AND xref:pdb`.

## `?fields=` for bulk pulls

Default JSON entry is ~20-30 KB. For batch work, restrict fields:

```python
fields = 'accession,id,gene_names,protein_name,length,sequence,xref_pdb,xref_alphafolddb'
url = 'https://rest.uniprot.org/uniprotkb/search'
params = {'query': 'organism_id:9606 AND reviewed:true', 'fields': fields, 'format': 'tsv', 'size': 500}
```

Common field selectors:
| Field | Returns |
|---|---|
| `accession`, `id` | Primary accession (P04637), entry name (P53_HUMAN) |
| `gene_names` | All gene names |
| `gene_primary` | Primary gene name only |
| `protein_name` | Recommended name |
| `organism_name`, `organism_id` | Species |
| `length`, `mass` | Sequence stats |
| `sequence` | The actual sequence |
| `cc_function`, `cc_subcellular_location` | Function and localization comments |
| `ft_domain`, `ft_binding`, `ft_active_site` | Domain/site features |
| `go_p`, `go_c`, `go_f` | GO biological process / cellular component / molecular function |
| `xref_pdb`, `xref_alphafolddb`, `xref_ensembl`, `xref_refseq` | Cross-references |
| `keyword` | UniProt keywords |
| `ec` | Enzyme classification |
| `reviewed` | Swiss-Prot vs TrEMBL flag |
| `cc_alternative_products` | Isoforms |

## Stream vs search vs accessions

| Endpoint | When | Limit |
|---|---|---|
| `/uniprotkb/{acc}` | One accession | 1 entry |
| `/uniprotkb/accessions?accessions=...` | Several known accessions | Up to ~100 per call |
| `/uniprotkb/search?query=...` | Query-driven; need pagination | 500 results per page; `cursor=` for paging |
| `/uniprotkb/stream?query=...` | Bulk query (>500) | No hard limit; one HTTP stream |

For 1000+ results, `/stream` is the right endpoint. Stream returns one HTTP response; iterate over the stream to avoid memory blowup.

## JSON schema navigation (the post-2022 layout)

The new schema is deeply nested. Common access patterns:

```python
entry = requests.get('https://rest.uniprot.org/uniprotkb/P04637.json').json()

acc = entry['primaryAccession']                                                # 'P04637'
entry_name = entry['uniProtkbId']                                              # 'P53_HUMAN'
sequence = entry['sequence']['value']                                          # actual AA sequence
length = entry['sequence']['length']

# Names (nested; defensive .get() because some fields are optional)
recommended = entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value')
primary_gene = entry.get('genes', [{}])[0].get('geneName', {}).get('value')

# Cross-references
xrefs_by_db = {}
for xref in entry.get('uniProtKBCrossReferences', []):
    xrefs_by_db.setdefault(xref['database'], []).append(xref['id'])

# Features (domains, binding sites)
domains = [f for f in entry.get('features', []) if f['type'] == 'Domain']
binding = [f for f in entry.get('features', []) if f['type'] == 'Binding site']

# Isoforms
isoforms = []
for comment in entry.get('comments', []):
    if comment.get('commentType') == 'ALTERNATIVE PRODUCTS':
        isoforms = [iso['name']['value'] for iso in comment.get('isoforms', [])]
```

## Isoform handling

Canonical sequence is returned for the bare accession (e.g. `P04637`). Isoforms have `-2`, `-3`, etc. suffixes (`P04637-2`). To fetch a specific isoform:

```python
iso = requests.get('https://rest.uniprot.org/uniprotkb/P04637-2.fasta').text
```

The canonical entry's `comments[type=ALTERNATIVE PRODUCTS]` lists all isoforms with their differences. For workflows needing all isoforms, iterate the list and fetch separately.

## ID Mapping API (async)

Convert between identifier systems (Ensembl Gene -> UniProt; PDB -> UniProt; UniProt -> RefSeq; etc.). The job pattern:

1. **Submit**: `POST /idmapping/run` with `ids`, `from`, `to`.
2. **Poll**: `GET /idmapping/status/{jobId}` — returns `{'jobStatus': 'RUNNING'}` or `{'results': [...]}`.
3. **Fetch**: `GET /idmapping/results/{jobId}` once status is complete.

Job typically completes in 30s; larger batches take 5-10 min. **Always set a poll timeout** — the API doesn't fail-soft on stuck jobs.

| From | To | Notes |
|---|---|---|
| `UniProtKB_AC-ID` | `UniProtKB` | Resolve obsolete to current accessions |
| `Gene_Name` | `UniProtKB` | Symbol -> accession (lossy; check matches) |
| `Ensembl` | `UniProtKB` | Ensembl Gene/Transcript/Protein |
| `EMBL-GenBank-DDBJ` | `UniProtKB` | INSDC nucleotide accessions |
| `RefSeq_Protein` | `UniProtKB` | NP_/XP_ accessions |
| `PDB` | `UniProtKB` | PDB chain to protein |
| `UniProtKB` | `EMBL-GenBank-DDBJ` | Reverse direction |

Full from/to list at https://rest.uniprot.org/configure/idmapping/fields.

## Code patterns

### Single entry with defensive JSON parsing

**Goal:** Fetch one UniProt entry as JSON and extract canonical name, gene, sequence, PDB cross-refs without KeyErrors.

**Approach:** GET `/uniprotkb/{acc}.json`; navigate with `.get()` chains; handle missing fields gracefully.

**Reference (UniProt REST as of 2024_06):**
```python
import requests


def fetch_uniprot_entry(accession):
    r = requests.get(f'https://rest.uniprot.org/uniprotkb/{accession}.json')
    r.raise_for_status()
    e = r.json()
    return {
        'accession': e['primaryAccession'],
        'entry_name': e.get('uniProtkbId'),
        'reviewed': e.get('entryType') == 'UniProtKB reviewed (Swiss-Prot)',
        'protein_name': e.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value'),
        'gene_primary': (e.get('genes') or [{}])[0].get('geneName', {}).get('value'),
        'sequence': e['sequence']['value'],
        'length': e['sequence']['length'],
        'pdb_ids': [x['id'] for x in e.get('uniProtKBCrossReferences', []) if x['database'] == 'PDB'],
        'alphafold_id': next((x['id'] for x in e.get('uniProtKBCrossReferences', []) if x['database'] == 'AlphaFoldDB'), None),
    }


print(fetch_uniprot_entry('P04637'))
```

### Search via TSV with `fields=` (bulk-friendly)

**Goal:** Get a DataFrame of human reviewed kinases with their PDB and AlphaFold IDs.

**Approach:** /search with format=tsv and explicit fields; paginate via `cursor` if results exceed 500.

**Reference (requests 2.31+):**
```python
import pandas as pd
from io import StringIO


def search_uniprot_tsv(query, fields, size=500):
    url = 'https://rest.uniprot.org/uniprotkb/search'
    params = {'query': query, 'fields': ','.join(fields), 'format': 'tsv', 'size': size}
    r = requests.get(url, params=params)
    r.raise_for_status()
    return pd.read_csv(StringIO(r.text), sep='\t')


df = search_uniprot_tsv(
    'organism_id:9606 AND reviewed:true AND keyword:"Kinase"',
    fields=['accession', 'gene_primary', 'protein_name', 'length', 'xref_pdb', 'xref_alphafolddb'],
)
print(f'{len(df)} reviewed human kinases')
print(df.head())
```

### Stream endpoint for >500 results

```python
import requests
import pandas as pd
from io import StringIO


def stream_uniprot(query, fields):
    url = 'https://rest.uniprot.org/uniprotkb/stream'
    params = {'query': query, 'fields': ','.join(fields), 'format': 'tsv'}
    r = requests.get(url, params=params, stream=True)
    r.raise_for_status()
    return pd.read_csv(StringIO(r.text), sep='\t')


# All human reviewed proteins (~20K)
df = stream_uniprot(
    'organism_id:9606 AND reviewed:true',
    fields=['accession', 'gene_primary', 'protein_name', 'length'],
)
print(f'All human Swiss-Prot: {len(df)}')
```

### ID mapping with proper async polling

**Goal:** Convert Ensembl Gene IDs to UniProt accessions.

**Approach:** Submit job; poll with timeout; retrieve results.

**Reference (UniProt REST 2024_06):**
```python
import time


def map_ids(ids, from_db='Ensembl', to_db='UniProtKB', timeout=600, poll_interval=3):
    submit = requests.post('https://rest.uniprot.org/idmapping/run',
                           data={'ids': ','.join(ids), 'from': from_db, 'to': to_db})
    submit.raise_for_status()
    job_id = submit.json()['jobId']
    print(f'Submitted job {job_id}')

    elapsed = 0
    while elapsed < timeout:
        status = requests.get(f'https://rest.uniprot.org/idmapping/status/{job_id}')
        status.raise_for_status()
        js = status.json()
        if 'jobStatus' in js and js['jobStatus'] == 'RUNNING':
            time.sleep(poll_interval)
            elapsed += poll_interval
            continue
        # Completed (results in status response) or has results endpoint
        break
    else:
        raise TimeoutError(f'ID mapping job {job_id} did not complete in {timeout}s')

    results = requests.get(f'https://rest.uniprot.org/idmapping/results/{job_id}')
    results.raise_for_status()
    return results.json()


mapping = map_ids(['ENSG00000141510', 'ENSG00000171862', 'ENSG00000139618'])
for r in mapping.get('results', []):
    print(f"  {r['from']:<20} -> {r['to']}")
for failed in mapping.get('failedIds', []):
    print(f"  {failed:<20} -> NOT MAPPED")
```

### Resolve obsolete accessions

```python
def resolve_obsolete(accessions):
    '''Use ID mapping to update obsolete accessions to current primary IDs.'''
    return map_ids(accessions, from_db='UniProtKB_AC-ID', to_db='UniProtKB')
```

### Download a reference proteome

```python
import gzip


def download_proteome(upid, out_path):
    '''upid: UniProt Proteome ID, e.g. UP000005640 (human reference).'''
    url = f'https://rest.uniprot.org/proteomes/{upid}.fasta.gz'
    r = requests.get(url, stream=True)
    r.raise_for_status()
    with open(out_path, 'wb') as f:
        for chunk in r.iter_content(8192):
            f.write(chunk)
    return out_path


download_proteome('UP000005640', 'human.fasta.gz')  # human reference proteome
```

### UniRef cluster lookup

```python
def uniref_cluster(uniref_id):
    '''e.g. UniRef50_P04637 -- the UniRef50 cluster centered on P04637.'''
    r = requests.get(f'https://rest.uniprot.org/uniref/{uniref_id}.json')
    r.raise_for_status()
    j = r.json()
    return {
        'id': j['id'],
        'representative': j['representativeMember']['memberId'],
        'member_count': j['memberCount'],
        'identity': j.get('entryType'),
    }
```

## Failure modes

### Legacy URL still in code (post-2022)
- **Trigger:** Old code using `https://www.uniprot.org/uniprot/{acc}.json`.
- **Mechanism:** 2022 migration; old URLs redirect but JSON schema is the new one — old parsers break.
- **Symptom:** Either 404 or `KeyError` from old field paths.
- **Fix:** Use `https://rest.uniprot.org/uniprotkb/{acc}.json`; update field navigation to the new nested layout.

### `?fields=` not specified
- **Trigger:** Bulk pull (1000 accessions) returning full JSON entries.
- **Mechanism:** Default returns ~20-30 KB per entry; 1000 entries = 20-30 MB.
- **Symptom:** Slow; memory blowup; rate-limit triggers.
- **Fix:** Always specify `fields=` for bulk; request only the fields actually needed.

### Search hit 500-record cap
- **Trigger:** Query matches 800 records; iterate first page only.
- **Mechanism:** /search returns 500 per page; need `cursor` for next.
- **Symptom:** Silently dropped tail.
- **Fix:** Use `/stream` for >500 results; or paginate /search with `cursor`.

### ID mapping job poll infinite loop
- **Trigger:** Network glitch during job; status forever "RUNNING".
- **Mechanism:** API doesn't time-out stuck jobs.
- **Symptom:** Pipeline hangs.
- **Fix:** Always set `timeout=` on polling; surface TimeoutError.

### Isoform suffix mishandled
- **Trigger:** Storing `P04637` and assuming that's the only sequence.
- **Mechanism:** TP53 has multiple isoforms; default fetch returns canonical only.
- **Symptom:** Missing alternative-product sequences.
- **Fix:** Read `comments[type=ALTERNATIVE PRODUCTS]`; fetch each isoform with `-N` suffix.

### Swiss-Prot vs TrEMBL confusion
- **Trigger:** Search without `reviewed:true` returning millions of TrEMBL hits.
- **Mechanism:** TrEMBL is automatically annotated, often low-quality.
- **Symptom:** "Why does my analysis include 200M proteins?"
- **Fix:** For reference-quality work, always filter `reviewed:true`.

### Obsolete accessions silently fail
- **Trigger:** Old paper-derived accession that has been merged or demerged.
- **Mechanism:** Direct fetch returns 404 or 301.
- **Symptom:** Missing entries in a batch.
- **Fix:** Use ID mapping (UniProtKB_AC-ID -> UniProtKB) to resolve to current accessions first.

### Gene-symbol disambiguation
- **Trigger:** Search `gene:TP53` returns multiple species or duplicates.
- **Mechanism:** Symbol is shared across species; UniProt indexes all.
- **Symptom:** Mixed-species hits.
- **Fix:** Combine with `organism_id:9606` (or specific taxon); use `gene_exact:` to avoid wildcard matches.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| 404 on legacy URL | Pre-2022 endpoint | Use rest.uniprot.org/uniprotkb/ |
| `KeyError` on old field path | Schema migration 2022 | Update to new nested layout; use `.get()` |
| Bulk fetch very slow | Default JSON entry size | Specify `fields=` for TSV bulk |
| Mid-pagination data missing | 500-record cap | Use /stream or paginate with cursor |
| ID mapping job hangs | API doesn't fail stuck jobs | Set `timeout=` on poll loop |
| Mixed-species search results | Symbol shared across species | Add `organism_id:` filter |
| Million-row search returning TrEMBL | No reviewed filter | Add `reviewed:true` |
| Missing isoform | Default returns canonical only | Fetch with `-N` suffix per isoform |

## References

- The UniProt Consortium. (2024) UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res* 53:D609-D617.
- Bursteinas B, Britto R, Bely B, et al. (2016) Minimizing proteome redundancy in the UniProt Knowledgebase. *Database* 2016:baw139.
- UniProt help: https://www.uniprot.org/help/api
- UniProt REST: https://rest.uniprot.org

## Related Skills

- entrez-fetch - NCBI protein records (RefSeq, GenPept) alternative
- biomart-queries - Alternative ID-mapping path via BioMart (preferred for Ensembl-rooted batches >5K; UniProt /idmapping/run is preferred for obsolete-accession resolution and any UniProt-rooted mapping)
- ortholog-inference - Resolve UniProt accessions used by OMA orthology queries
- structural-biology/structure-io - Download PDB structures referenced from UniProt
- structural-biology/alphafold-predictions - AlphaFoldDB entries cross-referenced in UniProt
- pathway-analysis/go-enrichment - Use GO annotations pulled from UniProt
