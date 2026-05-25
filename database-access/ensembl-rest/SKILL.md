---
name: bio-ensembl-rest
description: Query the Ensembl REST API for gene/transcript/protein lookup, sequence retrieval, comparative genomics (Compara), variant effect prediction (VEP), regulatory features, and cross-species ortholog/paralog calls. Use when pulling Ensembl-native data (Ensembl Gene IDs, version-pinned releases, archive endpoints for reproducibility), gene/transcript/exon structure with stable IDs, or VEP for variant annotation. Encodes the 15 req/sec rate limit, archive (e110.rest.ensembl.org) for reproducibility, Ensembl divisions (vertebrates / plants / fungi / metazoa / bacteria), and the symbol-vs-ID stability problem.
tool_type: python
primary_tool: requests
---

## Version Compatibility

Reference examples tested with: requests 2.31+, Ensembl REST API (release 110+); Ensembl release schedule is roughly quarterly

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show requests`
- API surface: check release notes at https://rest.ensembl.org

Each Ensembl release has an archive REST endpoint (e.g. `https://e110.rest.ensembl.org`) for reproducibility.

# Ensembl REST

**"Pull Ensembl-native gene / transcript / variant data programmatically"** -> Ensembl REST is distinct from NCBI Entrez and BioMart. It is the right answer for: stable Ensembl IDs, transcript / exon structure, VEP (Variant Effect Predictor) annotation, Compara orthologs at vertebrate scale, regulatory feature annotation, and any workflow rooted in Ensembl's coordinate system.

Two facts dominate Ensembl REST work: (1) the **15 req/sec / 55,000 req/hour rate limit** — high enough for hundreds of queries, low enough that bulk work (>5,000) belongs in BioMart instead; (2) **versioned archive endpoints** — `https://e110.rest.ensembl.org` pins to release 110 for reproducibility, while `https://rest.ensembl.org` follows the current release.

- Python: `requests.get('https://rest.ensembl.org/...')`
- Web: https://rest.ensembl.org (interactive doc with try-it-now)
- R: `biomaRt` for bulk (see `biomart-queries`); REST via `httr`

## Required Setup

```python
import requests
import time

BASE = 'https://rest.ensembl.org'
HEADERS = {'Accept': 'application/json'}
SLEEP = 0.07   # 15 req/sec ceiling
```

No API key required. Respect `Retry-After` header on 429.

## Ensembl divisions

Ensembl is divided by clade. Different REST hosts:

| Division | Host | Scope |
|---|---|---|
| Vertebrates | https://rest.ensembl.org | Human, mouse, fish, etc. (the "main" Ensembl) |
| Plants | https://rest.ensembl.org (plants division also accessible) | Arabidopsis, rice, etc. via Ensembl Genomes |
| Fungi | https://rest.ensemblgenomes.org | Yeasts, Aspergillus, etc. |
| Metazoa | https://rest.ensemblgenomes.org | Insects, nematodes, etc. |
| Bacteria | https://rest.ensemblgenomes.org | Limited (most bacteria in NCBI) |

For non-vertebrate work, check `ensemblgenomes.org` mirrors. As of 2024, Ensembl Genomes was being consolidated; check current host.

## Version pinning

| URL | Behavior |
|---|---|
| `https://rest.ensembl.org` | Current release (rolling) |
| `https://e110.rest.ensembl.org` | Pinned to release 110 |
| `https://e111.rest.ensembl.org` | Pinned to release 111 |
| `https://grch37.rest.ensembl.org` | Pinned to GRCh37 (legacy assembly) |

For any published analysis, pin the release. Ensembl releases change gene model versions, exon coordinates, and transcript annotations — re-running a pipeline a year later against the live endpoint may produce different results.

## Major endpoint groups

| Group | Example | Purpose |
|---|---|---|
| Lookup | `/lookup/symbol/human/BRCA1` | Resolve symbol or ID to stable record |
| Sequence | `/sequence/id/{id}` | DNA/protein sequence for ID |
| Cross References | `/xrefs/symbol/human/BRCA1` | Cross-refs to other DBs |
| Homology / Compara | `/homology/symbol/human/BRCA1` | Orthologs and paralogs |
| Gene Tree | `/genetree/id/{tree_id}` | Compara gene tree |
| VEP | `/vep/human/region/{region}/{allele}` | Variant effect prediction |
| Overlap | `/overlap/id/{id}` or `/overlap/region/{region}` | Genes/regulatory in interval |
| Regulatory | `/regulatory/species/{species}/feature/{id}` | Regulatory features |
| Variant | `/variation/{species}/{id}` | dbSNP / 1000G / ClinVar via Ensembl |
| LD | `/ld/{species}/pairwise/{var1}/{var2}` | LD between variants |
| GA4GH | `/ga4gh/...` | GA4GH-compliant subset |

Full reference: https://rest.ensembl.org (interactive).

## The symbol-vs-ID stability problem

Gene symbols are unstable (MARCH1 -> MARCHF1 in 2020 due to Excel autocorrect; SEPT* family also renamed). Ensembl Gene IDs (ENSG...) are stable across releases when the gene model is preserved.

**Best practice:**
1. Resolve symbol -> Ensembl ID once at pipeline start: `/lookup/symbol/{species}/{symbol}`.
2. Persist the Ensembl ID.
3. Run downstream queries by ID, not symbol.

Symbol-based endpoints are convenient for interactive use; ID-based endpoints are for reproducible pipelines.

## Rate-limit math

| Limit | Value |
|---|---|
| Burst | 15 req/sec |
| Hourly | 55,000 req/hour |
| Concurrent | Not enforced; courtesy 1-2 |

Respect `Retry-After` header on HTTP 429. For >5,000 queries, switch to BioMart bulk export (see `biomart-queries`) — BioMart has separate, more permissive limits.

## VEP (Variant Effect Predictor)

VEP via REST is the right call for ad hoc variant annotation. For batch variant annotation (>1000 variants), download VEP and run locally (`variant-calling/variant-annotation` skill).

REST modes:
- `/vep/{species}/region/{region}/{allele}` — single variant by coordinate
- `/vep/{species}/id/{variant_id}` — by dbSNP / Ensembl variant ID
- `/vep/{species}/hgvs/{hgvs_notation}` — by HGVS notation

VEP returns rich annotation: consequence (missense, synonymous, intron), SIFT/PolyPhen scores, gnomAD frequencies (if available), ClinVar significance.

## Compara homology

Compara orthology calls via Ensembl REST are covered in detail in `ortholog-inference` (database-access view). The relevant endpoint:

- `/homology/symbol/{species}/{symbol}` — all orthologs across Ensembl species
- `/homology/id/{ensembl_id}` — same, by ID
- `?target_species=` to restrict to one target
- `?type=orthologues` or `paralogues` to filter

## Code patterns

### Symbol -> stable Ensembl Gene ID

**Goal:** Resolve a gene symbol to its current Ensembl Gene ID once, then use the ID for all downstream queries.

**Approach:** `/lookup/symbol/{species}/{symbol}` returns the canonical record.

**Reference (requests 2.31+):**
```python
import requests
import time

BASE = 'https://rest.ensembl.org'
HEADERS = {'Accept': 'application/json'}


def get_with_retry(url, params=None, max_retries=3):
    for attempt in range(max_retries):
        r = requests.get(url, params=params, headers=HEADERS)
        if r.status_code == 429:
            time.sleep(int(r.headers.get('Retry-After', '5')))
            continue
        r.raise_for_status()
        return r
    raise RuntimeError(f'Failed after {max_retries} retries')


def symbol_to_ensembl(species, symbol):
    r = get_with_retry(f'{BASE}/lookup/symbol/{species}/{symbol}')
    return r.json()


info = symbol_to_ensembl('human', 'BRCA1')
print(f'  Ensembl Gene ID: {info["id"]}')
print(f'  Biotype:         {info["biotype"]}')
print(f'  Chromosome:      {info["seq_region_name"]}:{info["start"]}-{info["end"]}')
print(f'  Strand:          {info["strand"]}')
```

### Sequence retrieval by Ensembl ID

```python
def get_sequence(ensembl_id, seq_type='cdna'):
    '''seq_type: cdna, cds, protein, genomic'''
    r = get_with_retry(f'{BASE}/sequence/id/{ensembl_id}',
                       params={'type': seq_type, 'content-type': 'application/json'})
    return r.json()


prot = get_sequence('ENSG00000139618', seq_type='protein')
print(f'  Length: {len(prot["seq"])} aa')
```

### Overlap: what genes are in this region

```python
def genes_in_region(species, region):
    '''region as "chr:start-end" e.g. "17:43000000-44000000".'''
    r = get_with_retry(f'{BASE}/overlap/region/{species}/{region}',
                       params={'feature': 'gene'})
    return r.json()


for g in genes_in_region('human', '17:43000000-43200000'):
    print(f'  {g["external_name"]:<12} {g["id"]} {g["biotype"]:<20} {g["start"]}-{g["end"]}')
```

### VEP for a single variant

**Goal:** Get full annotation for a variant by coordinate.

**Approach:** `/vep/{species}/region/{region}/{allele}` returns transcript consequences, SIFT/PolyPhen, frequencies.

**Reference (Ensembl REST release 110+):**
```python
def vep_region(species, region, allele):
    r = get_with_retry(f'{BASE}/vep/{species}/region/{region}/{allele}')
    return r.json()


# BRCA1 missense variant in GRCh38 coordinates (rest.ensembl.org defaults to GRCh38);
# for GRCh37 coords use https://grch37.rest.ensembl.org instead.
results = vep_region('human', '17:43044295-43044295:1', 'A')
if results:
    for tc in results[0].get('transcript_consequences', [])[:5]:
        print(f'  {tc["gene_symbol"]:<8} {tc["consequence_terms"]}')
        if 'sift_prediction' in tc:
            print(f'    SIFT: {tc["sift_prediction"]} ({tc.get("sift_score", "?")})')
```

### Compara orthologs (Compara via REST)

```python
def orthologs(species, symbol, target_species=None):
    params = {'type': 'orthologues'}
    if target_species:
        params['target_species'] = target_species
    r = get_with_retry(f'{BASE}/homology/symbol/{species}/{symbol}', params=params)
    return r.json()['data'][0]['homologies']


for o in orthologs('human', 'BRCA1', target_species='mouse'):
    print(f'  {o["target"]["species"]:<15} {o["target"]["id"]}  type={o["type"]}  confidence={o.get("confidence")}')
```

### Batch lookup with rate-limit handling

```python
def batch_symbols(species, symbols):
    out = {}
    for sym in symbols:
        try:
            out[sym] = symbol_to_ensembl(species, sym)
        except requests.HTTPError as e:
            out[sym] = {'error': str(e)}
        time.sleep(0.07)  # 15 req/sec ceiling
    return out
```

### Archive endpoint for reproducibility

```python
# Pin to release 110
ARCHIVE = 'https://e110.rest.ensembl.org'
r = requests.get(f'{ARCHIVE}/lookup/symbol/human/BRCA1', headers={'Accept': 'application/json'})
print(r.json()['id'])
# Re-runs against e110 in 2030 will return the same Gene ID even if the live release has moved on.
```

### LD between two variants

```python
def ld_pairwise(species, var1, var2, population='1000GENOMES:phase_3:CEU'):
    r = get_with_retry(f'{BASE}/ld/{species}/pairwise/{var1}/{var2}',
                       params={'population_name': population})
    return r.json()
```

## Failure modes

### Symbol-based pipeline breaks on HGNC rename
- **Trigger:** `/lookup/symbol/human/MARCH1` after the 2020 rename.
- **Mechanism:** HGNC renamed Excel-autocorrect-affected genes; Ensembl mirrors the rename.
- **Symptom:** 404 or wrong gene returned.
- **Fix:** Resolve symbol to Ensembl Gene ID once at pipeline start; use ID downstream.

### No version pinning, results drift
- **Trigger:** Re-running a pipeline against `https://rest.ensembl.org` a year later.
- **Mechanism:** Live endpoint follows current release; gene models update quarterly.
- **Symptom:** Different transcript coordinates, exon counts, sometimes Gene ID changes.
- **Fix:** Pin to an archive endpoint (`https://e110.rest.ensembl.org`) for reproducibility.

### Rate-limit cascade
- **Trigger:** Loop of 5,000 REST calls without sleep.
- **Mechanism:** 15 req/sec ceiling; 429 with Retry-After.
- **Symptom:** Pipeline stalls; cascade of failures.
- **Fix:** Sleep 0.07s between calls; honor Retry-After; for >5K queries use BioMart bulk export.

### VEP for bulk variants
- **Trigger:** VEP REST for 100K variants.
- **Mechanism:** REST is per-variant; rate-limit makes this infeasible.
- **Symptom:** Days-long runtime; many 429s.
- **Fix:** Download VEP and run locally (`variant-calling/variant-annotation` skill); reserve REST for ad hoc <1K variants.

### Wrong species name
- **Trigger:** Using `Homo_sapiens` instead of `human`, or arbitrary capitalization.
- **Mechanism:** Ensembl species names are lowercase with underscores or common names; case matters.
- **Symptom:** 404 or empty result.
- **Fix:** Use `https://rest.ensembl.org/info/species` to enumerate valid names.

### Non-vertebrate species not in vertebrate host
- **Trigger:** Querying Arabidopsis on `rest.ensembl.org`.
- **Mechanism:** Plants live in Ensembl Genomes (`rest.ensemblgenomes.org`) as of 2024 (consolidation ongoing).
- **Symptom:** 404 or species not recognized.
- **Fix:** Check the right division host; use the lookup endpoint `info/divisions` to confirm.

### Archive endpoint TLS / connectivity
- **Trigger:** Older archive endpoints (e80, e90) gradually decommissioned.
- **Mechanism:** Very old archives are retired ~5 years after release.
- **Symptom:** DNS failure or 503.
- **Fix:** Check the current archive list at https://www.ensembl.org/info/website/archives/index.html.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `404` on symbol lookup | HGNC rename or wrong species | Resolve to Ensembl ID; check species code |
| HTTP 429 | Rate limit | Sleep per Retry-After; cap to 15 req/sec |
| Different results 6 months later | No version pinning | Use archive endpoint (`eXX.rest.ensembl.org`) |
| `404` on non-vertebrate | Wrong host division | Switch to `rest.ensemblgenomes.org` |
| VEP infeasible bulk | REST is per-variant | Local VEP for bulk |
| Old archive 503 | Decommissioned | Use a current archive release |

## References

- Yates AD, Allen J, Amode RM, et al. (2022) Ensembl Genomes 2022: an expanding genome resource for non-vertebrates. *Nucleic Acids Res* 50:D996-D1003.
- Martin FJ, Amode MR, Aneja A, et al. (2023) Ensembl 2023. *Nucleic Acids Res* 51:D933-D941.
- McLaren W, Gil L, Hunt SE, et al. (2016) The Ensembl Variant Effect Predictor. *Genome Biol* 17:122.
- Yates A, Beal K, Keenan S, et al. (2015) The Ensembl REST API: Ensembl data for any language. *Bioinformatics* 31:143-145.

## Related Skills

- biomart-queries - Ensembl BioMart for bulk (>5K) ID mapping
- ortholog-inference - Compara orthologs via Ensembl REST and other resources
- uniprot-access - Cross-reference Ensembl IDs in UniProt entries
- variant-calling/variant-annotation - Local VEP for bulk variant annotation
- ncbi-datasets-cli - NCBI alternative for genome / gene data
- entrez-search - NCBI alternative for non-Ensembl queries
