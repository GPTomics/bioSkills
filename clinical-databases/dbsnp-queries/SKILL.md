---
name: bio-clinical-databases-dbsnp-queries
description: Resolves rsIDs, navigates RsMergeArch/SNPHistory merge chains, and converts between rsID, SPDI, HGVS, and VCF representations using the dbSNP Build 156 JSON architecture. Use when normalizing variant identifiers, joining variant databases by cluster ID, or tracking deprecated rsIDs through historical merges.
tool_type: python
primary_tool: myvariant
---

## Version Compatibility

Reference examples tested with: myvariant 1.0+, requests 2.31+, biopython 1.83+, Entrez Direct 21.0+. dbSNP Build 156 (September 2022) is the current schema; Build 151 (2017) was the last with relational SQL dumps. Builds 152-155 dual-released JSON+SQL; 156+ is JSON-only.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. The Variation Services REST API uses path-based versioning (`/variation/v0/`); E-utilities `db=snp` returns thin legacy summaries missing build-156 schema fields.

# dbSNP Queries and rsID Normalization

**'Look up this rsID / normalize variant representations'** -> Resolve rsIDs through merge chains, compute canonical SPDI, and convert between rsID, HGVS-g, HGVS-c, and VCF allele representations.

- Python (aggregator): `myvariant.MyVariantInfo().getvariant(rsid, fields=['dbsnp', 'clinvar', 'gnomad_exome'])`
- Python (direct): `requests.get(f'https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid_int}')`
- Python (E-utilities, legacy): `Bio.Entrez.esearch(db='snp', term=rsid)`; returns thin summary
- Bulk: `ftp.ncbi.nlm.nih.gov/snp/latest_release/JSON/refsnp-chr{N}.json.bz2`

## rsID Is a Cluster Identifier, Not a Variant Identifier

This is the load-bearing concept. dbSNP cluster definition: ss records (submitted SNPs) are mapped to the genome and clustered into RefSNPs by *position + variant type*, not by allele. A single rsID can point to a *locus* with multiple alleles:

- `rs12345` may resolve to {A>G, A>T, A>C} at one position; the RefSNP JSON `primary_snapshot_data.placements_with_allele[*].alleles` enumerates them.
- ~6-8% of dbSNP rsIDs are multi-allelic (Phan 2025 *NAR* 25-year review).
- PLINK and many older tools historically misuse rsIDs as if they were variant identifiers, which fails for multi-allelic sites and yields wrong genotype assignments.

**Rule:** Use rsID as a *human-facing label only*; use SPDI or ClinGen Allele Registry CA ID for joins.

## Build 156 Schema Overhaul: What Changed

| Aspect | Build 151 (2017) | Build 156 (2022) and current |
|--------|------------------|------------------------------|
| Distribution | Relational SQL dumps + XML | JSON per RefSNP, partitioned by chromosome |
| FTP path | `ftp/snp/organisms/human_9606/` | `ftp.ncbi.nlm.nih.gov/snp/latest_release/JSON/` |
| Primary key | `snp_id`, `ss_id` | `refsnp_id`, with `primary_snapshot_data` block |
| Frequency data | Embedded sparse | ALFA aggregated populations |
| Merge tracking | `RsMergeArch.bcp.gz` | `refsnp-merged.json.bz2` (also `RsMergeArch.bcp.gz` retained for legacy) |
| Withdrawn | `SNPHistory.bcp.gz` | `refsnp-withdrawn.json.bz2` |
| API access | Legacy E-utilities `db=snp` only | Variation Services REST `/v0/refsnp/{id}` returns the full JSON |

E-utilities still works for `db=snp` but returns a thin pre-156 summary missing key fields like `primary_snapshot_data.placements_with_allele`; pipelines reliant on Entrez get out-of-date data.

## RsMergeArch: The Multi-Hop Merge Footgun

When two rsIDs are found to refer to the same allele cluster, the higher (later-assigned) rsID is merged into the lower. `RsMergeArch.bcp.gz` stores `(rsHigh, rsLow, rsCurrent)` tuples.

**The trap:** `rsCurrent` in any given row is the merge target at the time of that merge event, NOT the current dbSNP rsID. A multi-merge chain (rs3 -> rs2 -> rs1, then later rs1 -> rs0) appears as multiple rows. Naive one-hop lookup resolves to a stale ID.

Withdrawn rsIDs (submitter-withdrawn or QC-failed) live in `SNPHistory.bcp.gz`, not RsMergeArch. Both tables must be consulted to resolve any historical rsID.

## SPDI: The Canonical Variant Representation

SPDI (Sequence:Position:Deletion:Insertion) format: `NC_000017.11:43044294:G:A`. Position is **0-based, half-open** (differs from HGVS's 1-based, fully-closed).

The **Contextual Allele** transformation (Variant Overprecision Correction Algorithm) returns the right-aligned, normalized canonical form across left-aligned VCFs and right-aligned HGVS conventions. This is the basis for ClinGen Allele Registry CA ID computation.

| Representation | Build dependency | Transcript dependency | Bijective? | Best for |
|----------------|-----------------|----------------------|------------|----------|
| **VCF (chrom-pos-ref-alt)** | Yes | No | Yes (same build) | Pipelines, bulk |
| **SPDI** | Yes (via RefSeq accession) | No | Yes for SNV/small indel | Canonical normalization |
| **HGVS-g** | Yes (`NC_xxxxx.N`) | No | Yes for SNV/small indel | Human-readable genomic |
| **HGVS-c** | Indirect (via transcript) | Yes | No (one HGVS-c -> many HGVS-g) | Clinical reporting |
| **HGVS-p** | Indirect | Yes | Degenerate (one HGVS-p -> many HGVS-c) | Protein-level annotation |
| **rsID** | None (cluster identifier) | None | NO (multi-allelic) | Human label only |
| **CA ID** | None (canonical) | None | Yes | Cross-database join |

## Decision Tree by Query Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Resolve single rsID to coordinates + alleles | Variation Services `/v0/refsnp/{id}` | Returns full Build 156 JSON, including merge history |
| Resolve historical/deprecated rsID | Variation Services `/v0/refsnp/{id}` -> follow `merged_snapshot_data` chain | Single-hop RsMergeArch lookup misses multi-hop chains |
| Batch query 100-10k rsIDs | myvariant.info `getvariants(rsids)` | Aggregated with ClinVar/gnomAD overlay; rate-limit safe |
| Convert coords <-> rsID | myvariant.info HGVS query or Variation Services `/spdi/{spdi}/rsid` | SPDI is the canonical bridge |
| Normalize variant representations | Variation Services `/hgvs/{hgvs}/contextuals` | Returns canonical SPDI, right-aligned |
| Bulk genomic-wide rsID -> coords | Local download of `refsnp-chr{N}.json.bz2` + parser | No rate limits; weekly snapshots |
| Joining dbSNP with gnomAD by ID | Use SPDI or CA ID, never rsID alone | rsID is a cluster; alleles may not match |
| Get population AF for common variant | ALFA (via Variation Services) for array-genotyped variants; gnomAD for sequencing-derived | Different sample compositions |

## Single rsID Resolution

**Goal:** Resolve an rsID to full Build 156 RefSNP JSON, including coordinates, alleles, gene context, and merge history.

**Approach:** Hit Variation Services `/v0/refsnp/{id_without_rs}`; the response includes `primary_snapshot_data` (current) and `merged_snapshot_data` (if this rsID is itself a merge target).

```python
import requests

VARSVC = 'https://api.ncbi.nlm.nih.gov/variation/v0'

def refsnp(rsid):
    '''Fetch full Build 156 RefSNP JSON. rsid can be 'rs121913529' or 121913529.'''
    rs_int = str(rsid).lstrip('rs')
    r = requests.get(f'{VARSVC}/refsnp/{rs_int}', timeout=30)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()

def summarize_refsnp(payload):
    '''Extract minimal fields. Handles multi-allelic cluster correctly.

    The placement JSON nests assembly metadata; the precise path varies by
    Build / API version. Common variants seen in the wild:
        placement['seq_id_traits_by_assembly'][0]['assembly_name']
        placement['placement_annot']['seq_id_traits_by_assembly'][0]['assembly_name']
    Inspect the actual JSON returned for the current dbSNP Build before
    relying on either path in production.
    '''
    if payload is None or payload.get('is_withdrawn'):
        return None
    primary = payload.get('primary_snapshot_data', {})
    placements = primary.get('placements_with_allele', [])
    def assembly_name(p):
        traits = (p.get('placement_annot') or p).get('seq_id_traits_by_assembly') or []
        return traits[0].get('assembly_name') if traits else ''
    grch38 = next((p for p in placements if 'GRCh38' in (assembly_name(p) or '')), None)
    if grch38 is None:
        return None
    alleles = []
    for allele in grch38.get('alleles', []):
        spdi = allele.get('allele', {}).get('spdi', {})
        alleles.append({
            'ref': spdi.get('deleted_sequence'),
            'alt': spdi.get('inserted_sequence'),
            'seq_id': spdi.get('seq_id'),
            'pos_0based': spdi.get('position')
        })
    return {
        'rsid': payload.get('refsnp_id'),
        'gene': primary.get('allele_annotations', [{}])[0].get('assembly_annotation', [{}])[0].get('genes', [{}])[0].get('locus'),
        'placements_grch38': alleles,
        'is_multiallelic': len(alleles) > 2,
        'merge_history': payload.get('merged_snapshot_data', [])
    }
```

## Multi-Hop Merge Resolution

**Goal:** Resolve a possibly-deprecated rsID to the current canonical rsID, following the full merge chain.

**Approach:** Recursively follow `merged_snapshot_data` until the response has no further merge entries, with cycle detection.

```python
def resolve_merge_chain(rsid, max_hops=10):
    '''Follow multi-hop merge chain. Cycle-safe with max_hops cap.'''
    seen = set()
    current = str(rsid).lstrip('rs')
    for _ in range(max_hops):
        if current in seen:
            return {'error': 'merge cycle detected', 'chain': list(seen)}
        seen.add(current)
        payload = refsnp(current)
        if payload is None:
            return {'error': 'not found', 'final_rsid': current, 'chain': list(seen)}
        if payload.get('is_withdrawn'):
            return {'status': 'withdrawn', 'final_rsid': current, 'chain': list(seen)}
        primary = payload.get('primary_snapshot_data')
        if primary is not None:
            return {'status': 'resolved', 'final_rsid': payload.get('refsnp_id'), 'chain': list(seen)}
        merged = payload.get('merged_snapshot_data', [])
        if not merged:
            return {'status': 'orphan', 'final_rsid': current, 'chain': list(seen)}
        current = str(merged[0].get('merged_into', ''))
    return {'error': 'hop limit', 'chain': list(seen)}
```

## SPDI <-> HGVS <-> VCF Conversion

**Goal:** Move between variant representations using Variation Services as the canonical bridge.

**Approach:** SPDI endpoints handle build resolution and right-alignment; HGVS contextuals applies the Variant Overprecision Correction Algorithm.

```python
def hgvs_to_spdi_canonical(hgvs):
    '''Resolve HGVS to canonical SPDI via the Variant Overprecision Correction Algorithm.'''
    r = requests.get(f'{VARSVC}/hgvs/{hgvs}/contextuals', timeout=30)
    if not r.ok:
        return None
    contextuals = r.json().get('data', {}).get('spdis', [])
    return contextuals[0] if contextuals else None

def spdi_to_rsid(spdi_str):
    '''SPDI 'NC_000017.11:43044294:G:A' -> rsID if a cluster exists.'''
    r = requests.get(f'{VARSVC}/spdi/{spdi_str}/rsids', timeout=30)
    if not r.ok:
        return None
    rsids = r.json().get('data', {}).get('rsids', [])
    return rsids[0] if rsids else None

def vcf_to_canonical_spdi(chrom, pos, ref, alt, assembly='GRCh38'):
    '''VCF (1-based) -> SPDI (0-based, right-aligned).'''
    refseq_map = {('1', 'GRCh38'): 'NC_000001.11', ('17', 'GRCh38'): 'NC_000017.11'}
    refseq = refseq_map.get((str(chrom).lstrip('chr'), assembly))
    if refseq is None:
        return None
    raw_spdi = f'{refseq}:{pos - 1}:{ref}:{alt}'
    r = requests.get(f'{VARSVC}/spdi/{raw_spdi}/canonical_representative', timeout=30)
    return r.json().get('data', {}).get('spdi') if r.ok else None
```

## ALFA Frequencies vs gnomAD

| Source | Sample basis | Variants covered | When to use |
|--------|-------------|------------------|-------------|
| **ALFA** | ~1M dbGaP subjects (array + WGS) across 12 ancestry groups | 447M+ sites; broader (includes array-only) | Common variants, dbGaP-deposited cohorts |
| **gnomAD v4** | 807k WGS+WES individuals across 9 ancestry groups | Sequencing-derived (deeper at rare variants) | Rare-variant FAF95, ACMG BS1/BA1 |

ALFA does NOT provide FAF95-style upper-bound CIs; raw AF only. ALFA captures consent-tier metadata enabling consent-respecting lookups for variants gnomAD doesn't carry.

```python
def alfa_frequency(rsid, ancestry='Total'):
    '''Pull ALFA per-population AF via Variation Services.'''
    payload = refsnp(rsid)
    if payload is None:
        return None
    freq_records = payload.get('primary_snapshot_data', {}).get('allele_annotations', [{}])[0].get('frequency', [])
    alfa_records = [f for f in freq_records if 'ALFA' in f.get('study_name', '')]
    for record in alfa_records:
        if record.get('common_name') == ancestry:
            return {
                'allele': record.get('observation', {}).get('inserted_sequence'),
                'count': record.get('allele_count'),
                'total': record.get('total_count'),
                'freq': record.get('allele_count') / record.get('total_count') if record.get('total_count') else None
            }
    return None
```

## Per-Operation Failure Modes

**1. Treating rsID as a unique variant identifier**
- Trigger: Join two databases by rsID expecting a single variant.
- Mechanism: rsID is a cluster identifier; multi-allelic clusters have 2-4 alleles at one position.
- Symptom: Allele mismatches at low rate (~6-8% of sites); silent merger of unrelated variants.
- Fix: Normalize both sides to SPDI or CA ID before joining.

**2. Single-hop RsMergeArch lookup**
- Trigger: Read one row of `RsMergeArch.bcp.gz` and treat `rsCurrent` as the final answer.
- Mechanism: Multi-hop merges (rs3 -> rs2 -> rs1 -> rs0) span multiple rows; each row records one hop only.
- Symptom: Resolved rsID is itself stale; subsequent queries return outdated annotation.
- Fix: Follow merge chains recursively via Variation Services `merged_snapshot_data` (handles multi-hop in one call).

**3. Confusing withdrawn vs merged**
- Trigger: Query a withdrawn rsID and find no merge target.
- Mechanism: Withdrawn rsIDs live in `SNPHistory.bcp.gz`, not RsMergeArch. They are NOT merged into anything; the cluster was QC-failed or submitter-retracted.
- Symptom: 404 from naive lookup; pipelines proceed with stale rsID.
- Fix: Check `is_withdrawn` field in RefSNP JSON; flag the variant for manual review.

**4. E-utilities thin summary**
- Trigger: Use `Entrez.esummary(db='snp')` and treat output as authoritative.
- Mechanism: E-utilities `db=snp` returns a pre-Build-156 summary missing `primary_snapshot_data.placements_with_allele`, frequency data, and merge history.
- Symptom: Annotations look incomplete or out-of-date.
- Fix: Use Variation Services REST `/v0/refsnp/{id}` for full JSON.

**5. SPDI 0-based vs HGVS 1-based mismatch**
- Trigger: Build SPDI string from VCF position without converting to 0-based.
- Mechanism: SPDI position is 0-based half-open; VCF is 1-based.
- Symptom: Coordinate off-by-one; SPDI does not resolve to expected rsID.
- Fix: SPDI position = (VCF position - 1). For indels, also normalize ref/alt.

**6. Strand/orientation ambiguity for A/T C/G variants**
- Trigger: Merge variants across builds or platforms relying on rsID alone.
- Mechanism: rsID is locus-level; opposite-strand alleles get the same rsID with different ref/alt representation.
- Symptom: Strand-flipped genotypes after merge.
- Fix: Use SPDI (which encodes strand via deleted/inserted sequence) or MAF-match for ambiguous variants.

## Reconciliation: When Sources Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| dbSNP rsID returns 404 in current build | Withdrawn (in `refsnp-withdrawn.json.bz2`) | Check withdrawal reason; consider manual curation |
| rsID resolves to different coords across builds | Genome assembly change (GRCh37 -> GRCh38) | Use SPDI with explicit RefSeq accession; lift over via `pyliftover` |
| ALFA AF and gnomAD AF disagree by >2x | Different sample compositions; ALFA includes array-only sites under-represented in gnomAD | Trust gnomAD for sequencing data; ALFA for array-derived; use the more relevant source |
| Multiple rsIDs map to one SPDI | True duplicates from independent submissions; rare since Build 152 enforced cluster merging | Pick the lowest rsID per RsMergeArch convention |
| One rsID has different ref allele in dbSNP vs gnomAD | dbSNP uses NCBI's ref; gnomAD aligns to its build assembly | Normalize to SPDI before joining |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| Build 156 | Current as of Sep 2022; JSON-only distribution | dbSNP NCBI |
| Multi-allelic rate | ~6-8% of dbSNP rsIDs are multi-allelic | Phan 2025 *NAR* |
| Variation Services rate limit | 10 req/s with API key; 3 req/s without | NCBI E-utilities policy |
| SPDI position | 0-based, half-open | NCBI SPDI specification |
| HGVS position | 1-based, fully-closed | HGVS nomenclature |
| ALFA samples | ~1M individuals across 12 ancestries (2024 release) | NCBI dbGaP aggregation |
| Bulk download chunks | One file per chromosome (`refsnp-chr{N}.json.bz2`) | NCBI FTP |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| 404 from Variation Services on valid rsID | rsID is withdrawn or never assigned | Check `refsnp-withdrawn.json.bz2`; consider strand-flipped equivalent |
| Merge chain resolves but final rsID has different alleles | Multi-allelic cluster; pick allele matching the variant | Filter `placements_with_allele.alleles[*]` by allele match |
| ALFA frequency missing for common variant | Variant not in dbGaP-deposited studies | Fall back to gnomAD (sequencing-derived) |
| `Entrez.esummary` returns old data | Legacy E-utilities, not Build 156 schema | Switch to Variation Services REST `/v0/refsnp/{id}` |
| HGVS-c conversion fails for synonymous variants | Some HGVS-c rely on non-MANE transcripts not in NCBI default | Specify transcript explicitly; use VEP `--mane_select` |
| SPDI for indel does not round-trip | Left/right alignment mismatch | Use `/spdi/{spdi}/canonical_representative` for normalization |
| Bulk JSON parse OOM | `refsnp-chr1.json.bz2` is ~20GB uncompressed | Stream parse with `bz2.BZ2File` + line-by-line JSON; do not load whole file |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why not just use rsID for the join?" | rsID is a cluster identifier; ~6-8% of clusters are multi-allelic, causing silent mismatches. We use SPDI / CA ID. |
| "This annotation says rs12345 but the literature says rs67890" | rsIDs are merged; we resolved through `RsMergeArch` / `merged_snapshot_data` to the current canonical rsID. |
| "dbSNP frequency != gnomAD frequency" | ALFA (dbSNP-embedded) and gnomAD use different sample sets and different ascertainment (array vs sequencing); reconciled per use case. |
| "Why wasn't Entrez used?" | Entrez `db=snp` returns the pre-Build-156 summary missing key fields; we use Variation Services REST for the full JSON. |
| "Coordinate off-by-one in the SPDI" | SPDI is 0-based half-open; VCF is 1-based; intentional conversion applied. |

## References

- Phan L et al. 2025. dbSNP 25-year retrospective. *Nucleic Acids Res* 53:D940.
- Sayers EW et al. 2024. Database resources of the National Center for Biotechnology Information. *Nucleic Acids Res* 52:D33.
- Holmes JB et al. 2020. SPDI: data model for variants and applications at NCBI. *Bioinformatics* 36:1902.
- NCBI Variation Services API: `https://api.ncbi.nlm.nih.gov/variation/v0/`
- dbSNP FTP layout: `ftp.ncbi.nlm.nih.gov/snp/latest_release/JSON/`
- ALFA release notes: `https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/`
- ClinGen Allele Registry: `https://reg.clinicalgenome.org/docs/cg-car/`

## Related Skills

- clinical-databases/myvariant-queries - Aggregated rsID + annotation queries
- clinical-databases/clinvar-lookup - ClinVar VariationID vs rsID linkage
- clinical-databases/gnomad-frequencies - Frequency lookups by canonical SPDI
- clinical-databases/variant-prioritization - Pipeline using normalized variant IDs
- database-access/entrez-search - General Entrez query patterns
