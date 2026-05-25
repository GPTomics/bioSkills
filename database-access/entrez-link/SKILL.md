---
name: bio-entrez-link
description: Find cross-database references between NCBI databases using Biopython Bio.Entrez (ELink). Use when navigating gene to protein/structure, sequence to publication, PubMed to GEO, BioProject to SRA runs, or discovering all link relationships for a record. Covers linkname semantics, cmd= variants, asymmetric link warnings, neighbor_history for >200 input IDs, and per-database link tables.
tool_type: python
primary_tool: Bio.Entrez
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, Entrez Direct 21.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Bio.Entrez.elink)` to check signatures
- CLI: `elink -version` then `elink -help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Entrez Link

**"Find records linked to this record in another NCBI database"** -> ELink walks the curated, weekly-maintained link tables between Entrez databases. A link is an asserted relationship (e.g. "this PubMed article describes this nucleotide sequence"), not a similarity hit.

ELink is the navigation layer of Entrez. The decision that matters most is **which `linkname` to use** — not which databases. A single (`dbfrom`, `db`) pair can have a dozen `linkname` variants distinguishing curation level, evidence type, and direction. Picking the wrong one is the difference between 5 high-confidence matches and 500 noisy automated assertions.

- Python: `Entrez.elink(dbfrom=..., db=..., id=..., linkname=...)` (BioPython)
- CLI: `elink -db pubmed -target gene -name pubmed_gene_rif` (Entrez Direct)
- R: `entrez_link(dbfrom=..., db=..., id=...)` (rentrez)

## Required Setup

```python
from Bio import Entrez
Entrez.email = 'researcher@institution.edu'
Entrez.api_key = 'optional_api_key'  # raises rate to 10 req/sec
```

## The `linkname` decision (most important)

For most (`dbfrom`, `db`) pairs NCBI exposes multiple link tables. The qualifiers in the name encode the curation level and the evidence source. Choose deliberately.

### gene -> protein (representative example)

| linkname | Returns | When to use |
|---|---|---|
| `gene_protein` | All linked proteins (curated + automated) | Exploration; expect 10-1000x more hits |
| `gene_protein_refseq` | RefSeq proteins only | Reference-quality analyses; orthology |
| `gene_protein_swissprot` | Reviewed UniProt entries with NCBI cross-ref | Functional annotation; literature support |

### pubmed -> gene

| linkname | Returns |
|---|---|
| `pubmed_gene` | Genes mentioned in this paper (text-mined + curated) |
| `pubmed_gene_rif` | Genes with a Reference Into Function (curated, high-quality) |
| `pubmed_gene_pubmed` | Other PubMed records sharing gene linkage (rare use) |

### nucleotide -> protein

| linkname | Returns |
|---|---|
| `nuccore_protein` | All proteins encoded by this nucleotide record (CDS-linked) |
| `nuccore_protein_refseq` | RefSeq proteins only |

### Discover what link names exist for a pair

```python
h = Entrez.elink(dbfrom='gene', db='protein', id='672', cmd='acheck')
record = Entrez.read(h); h.close()
for ls in record[0]['IdCheckList']['IdLinkSet'][0]['LinkInfo']:
    print(f'{ls["Name"]}  -> {ls["DbTo"]} | {ls["MenuTag"]} ({ls["HtmlTag"]})')
```

`cmd='acheck'` is the only authoritative way to enumerate available linknames — they change with each NCBI release.

## Decision table: which `cmd` for which goal

| Goal | cmd | Returns |
|---|---|---|
| Get linked records | `neighbor` (default) | Linked IDs in target db |
| Get linked + relevance scores | `neighbor_score` | IDs with similarity scores (mostly `pubmed_pubmed`) |
| Get >200 source IDs in one go | `neighbor_history` | WebEnv + QueryKey for downstream EFetch |
| Enumerate available links | `acheck` | List of all linknames for source IDs |
| Check if any link exists | `ncheck` | Boolean per source ID |
| Check specific link exists | `lcheck` | Boolean per source ID + linkname |
| Get NCBI HTML link URLs | `llinks` | URLs to Entrez record pages |
| Get external provider links | `prlinks` | URLs to journal sites, etc. |

The `neighbor_history` cmd is essential when source `id` count exceeds ~200 — past that, the URL-length limit makes the comma-joined form fail. With `neighbor_history` ELink puts results on the history server and returns WebEnv/QueryKey for downstream pickup.

## Asymmetric link warning

ELink relationships are **not guaranteed symmetric**. `pubmed_gene` and `gene_pubmed` may return different sets because:
- Direction-dependent curation: gene-to-PubMed is curated by NCBI staff (GeneRIF); PubMed-to-gene includes text-mining.
- Cutoffs: some link tables truncate at N best links in one direction but not the other.
- Index lag asymmetry: when one db updates faster than the other.

If round-trip consistency matters (e.g. "every gene mentioned in this paper, then every paper mentioning each gene"), expect the round-trip set to be larger than the input — and never assume `A -> B -> A` returns the original ID alone.

## Per-database link catalog (curated subset)

### gene

| Target | Common linknames | Notes |
|---|---|---|
| protein | `gene_protein`, `gene_protein_refseq`, `gene_protein_swissprot` | RefSeq is the safe default |
| nuccore | `gene_nuccore`, `gene_nuccore_refseqrna`, `gene_nuccore_refseqgene` | `refseqrna` for mRNA, `refseqgene` for the curated gene region |
| pubmed | `gene_pubmed`, `gene_pubmed_rif` | RIF is curated and high-quality |
| homologene | `gene_homologene` | Deprecated 2014 but data still queryable |
| snp | `gene_snp` | dbSNP entries in gene region |
| clinvar | `gene_clinvar` | Clinical variants |
| omim | `gene_omim` | Disease associations |

### nuccore / nucleotide

| Target | Common linknames |
|---|---|
| protein | `nuccore_protein`, `nuccore_protein_refseq` |
| gene | `nuccore_gene` |
| taxonomy | `nuccore_taxonomy` |
| biosample | `nuccore_biosample` |
| sra | `nuccore_sra` |
| pubmed | `nuccore_pubmed`, `nuccore_pubmed_refseq` |

### protein

| Target | Common linknames |
|---|---|
| nuccore | `protein_nuccore`, `protein_nuccore_cds`, `protein_nuccore_mrna` |
| gene | `protein_gene` |
| structure | `protein_structure` |
| cdd | `protein_cdd` (conserved domains) |
| pubmed | `protein_pubmed` |

### pubmed

| Target | Common linknames |
|---|---|
| pubmed | `pubmed_pubmed`, `pubmed_pubmed_citedin`, `pubmed_pubmed_refs` |
| gene | `pubmed_gene`, `pubmed_gene_rif` |
| protein | `pubmed_protein` |
| nuccore | `pubmed_nuccore` |
| gds | `pubmed_gds` (GEO datasets cited in paper) |
| sra | `pubmed_sra` |

### bioproject

| Target | Common linknames |
|---|---|
| biosample | `bioproject_biosample` |
| sra | `bioproject_sra` |
| pubmed | `bioproject_pubmed` |

## Code patterns

### Single source -> single target

**Goal:** Get RefSeq proteins for a single gene.

**Approach:** ELink with explicit `linkname` to restrict to curated set.

**Reference (BioPython 1.83+):**
```python
def gene_to_refseq_proteins(gene_id):
    h = Entrez.elink(dbfrom='gene', db='protein', id=gene_id, linkname='gene_protein_refseq')
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    return [link['Id'] for link in r[0]['LinkSetDb'][0]['Link']]

print(gene_to_refseq_proteins('672'))  # BRCA1
```

### Batch source -> target (small batch)

**Goal:** Get linked proteins for a list of <200 gene IDs in one call.

**Approach:** Comma-join IDs; one linkset per input in the response.

**Reference (BioPython 1.83+):**
```python
def batch_gene_protein(gene_ids):
    h = Entrez.elink(dbfrom='gene', db='protein', id=','.join(gene_ids), linkname='gene_protein_refseq')
    r = Entrez.read(h); h.close()
    out = {}
    for linkset in r:
        src = linkset['IdList'][0]
        out[src] = [link['Id'] for link in linkset['LinkSetDb'][0]['Link']] if linkset['LinkSetDb'] else []
    return out
```

### Large batch via history server

**Goal:** Link 5,000 gene IDs to proteins without hitting URL-length limits.

**Approach:** EPost the IDs first (chunked at 200), then ELink with `cmd='neighbor_history'` referencing the WebEnv. Downstream EFetch picks up linked IDs from the history server.

**Reference (BioPython 1.83+):**
```python
def post_then_link(gene_ids, target='protein', linkname='gene_protein_refseq'):
    # EPost in chunks of 200
    webenv = None
    for i in range(0, len(gene_ids), 200):
        chunk = gene_ids[i:i+200]
        kwargs = {'db': 'gene', 'id': ','.join(chunk)}
        if webenv:
            kwargs['WebEnv'] = webenv
        h = Entrez.epost(**kwargs)
        r = Entrez.read(h); h.close()
        webenv = r['WebEnv']
        query_key = r['QueryKey']
        time.sleep(0.1 if Entrez.api_key else 0.34)

    # Link with neighbor_history
    h = Entrez.elink(dbfrom='gene', db=target, linkname=linkname,
                     cmd='neighbor_history', WebEnv=webenv, query_key=query_key)
    r = Entrez.read(h); h.close()
    # WebEnv is at the top level of the response; QueryKey is per-LinkSetDbHistory entry.
    return r[0]['WebEnv'], r[0]['LinkSetDbHistory'][0]['QueryKey']

we, qk = post_then_link(['672', '675', '7157'] * 1000)
# Downstream: Entrez.efetch(db='protein', WebEnv=we, query_key=qk, retstart=..., retmax=500)
```

### Discover all available links

**Goal:** Before writing a pipeline, enumerate what link tables NCBI exposes for a (dbfrom, source-id) pair.

**Approach:** `cmd='acheck'` returns the full LinkInfo list per source.

**Reference (BioPython 1.83+):**
```python
def list_link_names(dbfrom, id):
    h = Entrez.elink(dbfrom=dbfrom, id=id, cmd='acheck')
    r = Entrez.read(h); h.close()
    info = r[0]['IdCheckList']['IdLinkSet'][0]['LinkInfo']
    return [(i['Name'], i['DbTo'], i['MenuTag']) for i in info]

for name, target, label in list_link_names('gene', '672'):
    print(f'{name:<40} -> {target:<15} ({label})')
```

### Chain links (gene -> protein -> structure)

```python
def gene_to_structures(gene_id):
    h = Entrez.elink(dbfrom='gene', db='protein', id=gene_id, linkname='gene_protein_refseq')
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    prot_ids = [l['Id'] for l in r[0]['LinkSetDb'][0]['Link'][:10]]
    time.sleep(0.1 if Entrez.api_key else 0.34)
    h = Entrez.elink(dbfrom='protein', db='structure', id=','.join(prot_ids))
    r = Entrez.read(h); h.close()
    out = []
    for ls in r:
        if ls['LinkSetDb']:
            out.extend(l['Id'] for l in ls['LinkSetDb'][0]['Link'])
    return out
```

### Get neighbor_score for related PubMed articles

```python
def related_pubmed(pmid, top=10):
    h = Entrez.elink(dbfrom='pubmed', db='pubmed', id=pmid,
                     linkname='pubmed_pubmed', cmd='neighbor_score')
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    return [(l['Id'], int(l['Score'])) for l in r[0]['LinkSetDb'][0]['Link'][:top]]
```

### BioProject -> SRA runs

For SRA discovery, `pysradb.SRAweb().sra_metadata(prjna, detailed=True)` (see `sra-data`) is the higher-fidelity path — returns SRR accessions directly with run-level metadata in one call. Use ELink only when staying inside Bio.Entrez:

```python
def bioproject_to_sra(prjna):
    # Convert PRJNA to UID first
    h = Entrez.esearch(db='bioproject', term=f'{prjna}[BioProject]')
    r = Entrez.read(h); h.close()
    if not r['IdList']:
        return []
    bp_uid = r['IdList'][0]
    time.sleep(0.1 if Entrez.api_key else 0.34)
    # Link to SRA
    h = Entrez.elink(dbfrom='bioproject', db='sra', id=bp_uid)
    r = Entrez.read(h); h.close()
    return [l['Id'] for l in r[0]['LinkSetDb'][0]['Link']] if r[0]['LinkSetDb'] else []
```

## Failure modes

### Wrong linkname gives wrong order of magnitude
- **Trigger:** Using `gene_protein` when `gene_protein_refseq` was intended.
- **Mechanism:** `gene_protein` includes all automated and predicted entries (XP_* RefSeq plus all GenBank submissions).
- **Symptom:** 500 proteins returned per gene instead of the expected 1-5 canonical isoforms.
- **Fix:** Pick the curated linkname; verify counts on a known gene.

### Empty LinkSetDb on valid input
- **Trigger:** Gene with no linked records in the requested target.
- **Mechanism:** `record[0]['LinkSetDb']` is an empty list, not raising an error.
- **Symptom:** `KeyError` if code assumes `record[0]['LinkSetDb'][0]` always exists.
- **Fix:** Always guard `if not record[0]['LinkSetDb']: return []`.

### Asymmetric round-trip
- **Trigger:** Pipeline does `genes_for_paper(pmid) -> papers_for_each_gene -> set of PMIDs`.
- **Mechanism:** `pubmed_gene` (text-mined + curated) is larger than `gene_pubmed` (curated only); the round-trip set is not closed.
- **Symptom:** Original PMID may not appear in the round-trip set; new PMIDs do.
- **Fix:** Document the directional asymmetry; use the more-curated linkname (`*_rif` variants) when fidelity matters.

### URL length limit on large batches
- **Trigger:** Comma-joined `id=` with 200+ IDs.
- **Mechanism:** HTTP GET URL exceeds NCBI's parsing limit (~2000 chars).
- **Symptom:** HTTP 414 URI Too Long, or silent truncation.
- **Fix:** EPost the IDs first, then ELink with `cmd='neighbor_history'`.

### One linkset per input ID, indexing confusion
- **Trigger:** Sending 5 IDs, then accessing `record[0]['LinkSetDb'][0]['Link']` expecting the union.
- **Mechanism:** ELink returns one `LinkSet` per input UID, indexed by position.
- **Symptom:** Only the first input's links are processed; rest are dropped.
- **Fix:** Iterate `for linkset in record:` and map by `linkset['IdList'][0]`.

### Mismatched dbfrom and id namespace
- **Trigger:** Passing a PMID into `dbfrom='nucleotide'`.
- **Mechanism:** ELink returns no error — it just looks up the PMID as a nucleotide UID, finds nothing.
- **Symptom:** Empty LinkSetDb on a "valid" ID.
- **Fix:** Validate that the ID matches the source db namespace (PMIDs are db=pubmed, GeneIDs are db=gene).

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `KeyError: 'LinkSetDb'` | Empty result not guarded | `if not record[0]['LinkSetDb']: return []` |
| `HTTPError 414` | Comma-joined id too long | Use EPost + `neighbor_history` |
| `HTTPError 400` | Invalid linkname or wrong db namespace | Use `cmd='acheck'` to enumerate valid links |
| 500 hits instead of 5 | Wrong linkname (e.g. `gene_protein` vs `_refseq`) | Pick curated variant |
| Round-trip set differs from input | Asymmetric link tables | Document; use curated variants |

## References

- Sayers EW et al. (2024) Database resources of the National Center for Biotechnology Information in 2024. *Nucleic Acids Res* 52:D33-D43.
- Kans J. (2024) Entrez Direct: E-utilities on the Unix Command Line. NCBI Bookshelf NBK179288.
- NCBI. ELink help. NBK25499.

## Related Skills

- entrez-search - Resolve UIDs before linking
- entrez-fetch - Retrieve linked records' content
- batch-downloads - History-server retrieval after ELink with `neighbor_history`
- geo-data - Specialized gds <-> pubmed/bioproject links (gds->sra ELink unreliable; use pysradb)
- ncbi-datasets-cli - Modern alternative for gene/genome cross-reference queries
