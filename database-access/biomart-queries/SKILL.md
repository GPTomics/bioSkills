---
name: bio-biomart-queries
description: Bulk-query Ensembl BioMart (and other BioMart instances) for cross-database ID mapping, gene/transcript/exon coordinates, and ortholog tables. Use when batch-converting Ensembl IDs to other namespaces (HGNC, RefSeq, UniProt, Entrez), pulling gene coordinate tables for thousands of genes, building ortholog wide-tables across species, or replacing slow Ensembl REST loops with one-shot bulk export. Encodes BioMart's XML query format, R biomaRt vs Python pybiomart trade-off, mart-vs-dataset hierarchy, and the URL endpoint that's BioMart-specific (separate from rest.ensembl.org).
tool_type: mixed
primary_tool: pybiomart
---

## Version Compatibility

Reference examples tested with: pybiomart 0.9+, R biomaRt 2.58+ (Bioconductor); Ensembl BioMart (release 110+)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show pybiomart`
- R: `packageVersion('biomaRt')`

The BioMart XML query format is stable across Ensembl releases; the underlying mart names and attribute IDs can change between Ensembl releases. For published work, pin the Ensembl release via `useEnsembl(version=110)`.

# BioMart Queries

**"Bulk-convert IDs / pull coordinate tables / extract ortholog wide tables"** -> BioMart is the right answer for any Ensembl-rooted query producing >5,000 rows. It is a separate service from the Ensembl REST API, with separate rate behavior and a different query model (XML-based, batch-oriented). For one-off lookups (<100 records), Ensembl REST is more convenient; for bulk anything, BioMart wins.

The single most important fact: **BioMart returns a flat table from a single query**. There is no per-record loop, no rate-limit cascade, no async polling. One XML query in; one TSV out.

- Python: `pybiomart` (https://github.com/jrderuiter/pybiomart) is the lightest client
- R: `biomaRt` Bioconductor (Durinck et al. 2009 *Nat Protoc* 4:1184) is the canonical client
- CLI: `curl` against the XML endpoint works but is rarely used directly
- Web: `https://www.ensembl.org/biomart/martview` for interactive query design

## Installation

```bash
pip install pybiomart pandas
# R:
# BiocManager::install('biomaRt')
```

## BioMart hierarchy

| Level | Examples |
|---|---|
| Mart | `ENSEMBL_MART_ENSEMBL` (genes), `ENSEMBL_MART_SNP` (variants), `ENSEMBL_MART_MOUSE` (mouse-specific) |
| Dataset | `hsapiens_gene_ensembl`, `mmusculus_gene_ensembl`, etc. (per species) |
| Attribute | Fields to return: `ensembl_gene_id`, `external_gene_name`, `chromosome_name`, etc. |
| Filter | Constraints on the query: `chromosome_name = 17`, `biotype = protein_coding`, etc. |

A query is: pick a mart, pick a dataset, list attributes to return, list filters to constrain. BioMart returns a single TSV.

Discovery:
```python
from pybiomart import Server
server = Server(host='http://www.ensembl.org')
print(server.marts)                                          # list marts
mart = server['ENSEMBL_MART_ENSEMBL']
print(mart.datasets)                                         # list datasets (species)
ds = mart['hsapiens_gene_ensembl']
print(ds.attributes)                                         # list attributes
print(ds.filters)                                            # list filters
```

## Decision matrix: BioMart vs Ensembl REST

| Question | BioMart | Ensembl REST |
|---|---|---|
| Bulk ID mapping (>5000 IDs) | yes (1 query) | rate-limited cascade |
| Single-gene lookup | overkill | yes |
| Coordinate tables for thousands of genes | yes | rate-limited |
| Ortholog wide-table across species | yes (multi-species mart) | per-gene loop |
| VEP variant annotation | no | yes (or local VEP) |
| Sequence retrieval | partial | yes |
| Real-time | no (batch) | yes (per-record) |
| Reproducibility (version pin) | `useEnsembl(version=110)` | archive URL `e110.rest.ensembl.org` |

For >5K rows, BioMart is the right tool. For real-time per-record lookups, REST.

## Common attribute selectors

| Attribute | Returns |
|---|---|
| `ensembl_gene_id` | Stable Ensembl Gene ID |
| `ensembl_gene_id_version` | With `.N` version suffix |
| `external_gene_name` | HGNC symbol (or species-equivalent) |
| `hgnc_id`, `hgnc_symbol` | HGNC permanent ID and symbol |
| `entrezgene_id` | NCBI Gene ID |
| `refseq_mrna`, `refseq_peptide` | RefSeq accessions |
| `uniprotswissprot`, `uniprotsptrembl` | UniProt accessions |
| `chromosome_name`, `start_position`, `end_position`, `strand` | Gene coordinates |
| `transcript_count`, `exon_count` | Counts |
| `biotype` | protein_coding, lncRNA, miRNA, etc. |
| `description` | Free-text gene description |
| `go_id`, `name_1006`, `namespace_1003` | GO term ID, name, namespace |

## Common filter selectors

| Filter | Constraint |
|---|---|
| `ensembl_gene_id` | List of Gene IDs |
| `external_gene_name` | List of symbols |
| `entrezgene_id` | List of NCBI Gene IDs |
| `chromosome_name` | One or more chromosomes |
| `start` / `end` | Coordinate range |
| `biotype` | One or more biotypes |
| `with_<source>` | Boolean: has cross-ref to `<source>` (e.g. `with_hpa` = has Human Protein Atlas) |

## Code patterns

### Bulk ID mapping: Ensembl Gene -> HGNC + RefSeq + UniProt

**Goal:** Convert 5,000 Ensembl Gene IDs to HGNC symbols, RefSeq mRNA accessions, and UniProt accessions in one query.

**Approach:** pybiomart query with three attributes; ID list as a filter; returns one TSV.

**Reference (pybiomart 0.9+, Ensembl release 110+):**
```python
from pybiomart import Server
import pandas as pd

server = Server(host='http://www.ensembl.org')
mart = server['ENSEMBL_MART_ENSEMBL']
ds = mart['hsapiens_gene_ensembl']

ensembl_ids = ['ENSG00000139618', 'ENSG00000141510', 'ENSG00000171862']  # ...up to 5K+

df = ds.query(
    attributes=['ensembl_gene_id', 'external_gene_name', 'hgnc_id',
                'refseq_mrna', 'uniprotswissprot'],
    filters={'ensembl_gene_id': ensembl_ids},
)
print(df.head())
# One row per (gene, cross-ref) pair; genes with multiple RefSeq mRNAs get multiple rows.
```

### Pull gene coordinate table for a chromosome

```python
df = ds.query(
    attributes=['ensembl_gene_id', 'external_gene_name', 'chromosome_name',
                'start_position', 'end_position', 'strand', 'biotype'],
    filters={'chromosome_name': '17', 'biotype': 'protein_coding'},
)
print(f'{len(df)} protein-coding genes on chr17')
```

### Bulk ortholog wide-table (human <-> mouse <-> zebrafish)

**Goal:** One TSV with human Ensembl ID, mouse ortholog Ensembl ID, zebrafish ortholog Ensembl ID per row.

**Approach:** Ortholog attributes from the human mart query both species' orthologs.

```python
df = ds.query(
    attributes=['ensembl_gene_id', 'external_gene_name',
                'mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_orthology_type',
                'drerio_homolog_ensembl_gene', 'drerio_homolog_orthology_type'],
    filters={'chromosome_name': '17'},
)
# pybiomart columns use the mart display names, which can vary across releases.
# Resolve column names defensively rather than hardcoding strings:
mouse_type_col = next(c for c in df.columns if 'Mouse' in c and 'type' in c)
zebra_type_col = next(c for c in df.columns if 'Zebrafish' in c and 'type' in c)
df_one2one = df[(df[mouse_type_col] == 'ortholog_one2one') &
                (df[zebra_type_col] == 'ortholog_one2one')]
print(f'{len(df_one2one)} 1:1 orthologs across all three species on chr17')
```

### GO term annotation for a gene set

```python
df = ds.query(
    attributes=['ensembl_gene_id', 'external_gene_name',
                'go_id', 'name_1006', 'namespace_1003'],
    filters={'external_gene_name': ['TP53', 'BRCA1', 'MYC', 'EGFR']},
)
# Long format: one row per (gene, GO term) pair
```

### Version-pinned query (R biomaRt)

```r
# Reference: Bioconductor biomaRt 2.58+ | Verify API if version differs
library(biomaRt)

# Pin to release 110 for reproducibility
ensembl <- useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl', version=110)

# Or via host URL (for older or specific assemblies)
# ensembl <- useMart('ENSEMBL_MART_ENSEMBL',
#                     dataset='hsapiens_gene_ensembl',
#                     host='https://nov2020.archive.ensembl.org')

df <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id',
                   'uniprotswissprot', 'refseq_mrna'),
    filters = 'ensembl_gene_id',
    values = c('ENSG00000139618', 'ENSG00000141510'),
    mart = ensembl
)
head(df)
```

### Discover attributes / filters programmatically

```python
# What attributes are available?
attrs = ds.attributes
ortho_attrs = [a for a in attrs if 'homolog' in a]
print(f'{len(ortho_attrs)} ortholog attributes; first 5: {ortho_attrs[:5]}')

# What filters?
filts = ds.filters
chrom_filts = [f for f in filts if 'chrom' in f]
```

## Failure modes

### Trying to pull >100K rows in one query
- **Trigger:** Query without any filter (e.g. all attributes for the whole human genome).
- **Mechanism:** BioMart times out or truncates on very large queries.
- **Symptom:** Empty or partial result.
- **Fix:** Chunk by chromosome; combine results client-side.

### No version pinning
- **Trigger:** `useMart('ensembl', ...)` without `version=`.
- **Mechanism:** Defaults to current release; gene model versions change quarterly.
- **Symptom:** Re-running a year later produces different rows.
- **Fix:** Pin with `useEnsembl(version=110)` or archive host URL.

### Multiple cross-refs balloon row count
- **Trigger:** Query for `ensembl_gene_id, refseq_mrna`; a gene with 10 RefSeq mRNAs produces 10 rows.
- **Mechanism:** BioMart joins on cross-refs; many-to-many produces row multiplication.
- **Symptom:** "Why do I have 50K rows for 5K input IDs?"
- **Fix:** Filter to one isoform per gene downstream; or use `ensembl_canonical` filter where available.

### Symbol-based filter misses HGNC renames
- **Trigger:** `filters={'external_gene_name': ['MARCH1']}` post-2020.
- **Mechanism:** HGNC renamed to MARCHF1; BioMart mirrors the new symbol.
- **Symptom:** Empty result for that gene.
- **Fix:** Filter by `ensembl_gene_id` or `hgnc_id`; these are stable.

### Multi-species mart query slow
- **Trigger:** Querying `mmusculus_homolog_ensembl_gene` for 30K human genes.
- **Mechanism:** Ortholog attributes are heavy; large queries take minutes.
- **Symptom:** Timeout or slow.
- **Fix:** Chunk by chromosome; or use Ensembl Compara REST for targeted lookups.

### REST loops where BioMart belongs
- **Trigger:** Loop of 5,000 Ensembl REST `/lookup/symbol` calls.
- **Mechanism:** Rate-limit cascade; 5,000 * 0.07s = 6 minutes just for the rate gate, plus HTTP overhead.
- **Symptom:** Slow; 429 errors.
- **Fix:** Switch to one BioMart query.

### Wrong mart for the question
- **Trigger:** Querying gene info from `ENSEMBL_MART_SNP`.
- **Mechanism:** SNP mart has variant attributes, not gene attributes.
- **Symptom:** Empty result or wrong fields.
- **Fix:** Discover marts with `server.marts`; pick `ENSEMBL_MART_ENSEMBL` for genes.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| Empty result | Wrong attribute / filter name | List with `ds.attributes` and `ds.filters` |
| Timeout on big query | No filter, too many rows | Chunk by chromosome |
| Drift between re-runs | No version pinning | `useEnsembl(version=110)` |
| Row count > expected | Many-to-many cross-ref joins | Filter to canonical isoform |
| Symbol filter returns nothing | HGNC rename | Filter by Ensembl ID or HGNC ID |
| Slow on ortholog wide-table | Multi-species join expensive | Chunk by chromosome |

## References

- Durinck S, Spellman PT, Birney E, Huber W. (2009) Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. *Nat Protoc* 4:1184-1191.
- Kinsella RJ, Kahari A, Haider S, et al. (2011) Ensembl BioMarts: a hub for data retrieval across taxonomic space. *Database* 2011:bar030.
- Smedley D, Haider S, Durinck S, et al. (2015) The BioMart community portal: an innovative alternative to large, centralized data repositories. *Nucleic Acids Res* 43:W589-W598.
- pybiomart documentation: https://github.com/jrderuiter/pybiomart

## Related Skills

- ensembl-rest - Per-record Ensembl queries (BioMart's complement)
- ortholog-inference - Compara ortholog calls with confidence semantics
- uniprot-access - UniProt ID mapping (preferred for UniProt-rooted lookups and obsolete-accession resolution; BioMart is preferred for Ensembl-rooted batches >5K)
- ncbi-datasets-cli - NCBI-side bulk path for genome / gene data
- entrez-search - NCBI alternative for non-Ensembl queries
