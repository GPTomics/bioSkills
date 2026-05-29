---
name: bio-clinical-databases-myvariant-queries
description: Queries myvariant.info BioThings aggregator for ClinVar, gnomAD, dbSNP, dbNSFP, COSMIC, CADD, and CIViC annotations in batched, version-tracked requests. Use when annotating variant lists from multiple databases simultaneously without managing per-source APIs, and when reproducibility-grade analyses require recording source data versions via _meta.
tool_type: python
primary_tool: myvariant
---

## Version Compatibility

Reference examples tested with: myvariant 1.0.0+, requests 2.31+, pandas 2.2+. myvariant.info aggregates >=21 sources; the operative version of each source is queryable via the `_meta` field and the `/v1/metadata` endpoint.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. dbNSFP version drift is the dominant staleness vector: AlphaMissense was added to dbNSFP v4.4 (~2024); querying `dbnsfp.alphamissense.score` returns whatever version of dbNSFP is currently loaded; check `_meta.src.dbnsfp.version`.

# MyVariant.info Queries; Aggregated Annotation

**'Annotate my variants with ClinVar + gnomAD + CADD + AlphaMissense in one batch'** -> Query the BioThings myvariant.info aggregator with field selection and version tracking, then parse nested responses.

- Python: `myvariant.MyVariantInfo().getvariant(hgvs_or_rsid, fields=['clinvar', 'gnomad_exome', 'dbnsfp'])`
- Python (batch): `mv.getvariants(ids_list, fields=...)`; up to 1000 IDs per request
- Python (search): `mv.query('clinvar.gene.symbol:BRCA1 AND clinvar.clinical_significance:Pathogenic')`
- REST: `GET https://myvariant.info/v1/variant/{hgvs_or_id}?fields=...`
- Bulk: `POST https://myvariant.info/v1/variant` with comma-separated IDs

## BioThings Architecture (Wu 2022 *Bioinformatics*)

myvariant.info is one of three flagship BioThings APIs (with MyGene.info and MyChem.info). All three share the BioThings SDK, which auto-deploys an Elasticsearch index from heterogeneous source files via per-source dataloaders. The 2022 paper formalized the SDK; the architecture itself is older (Xin 2016 *Genome Biol*).

- Elasticsearch-backed: queries use Lucene operators (AND, OR, NOT, range like `dbnsfp.cadd.phred:>20`)
- Dotted-field-name syntax for nested JSON
- The `_id` field is canonical HGVS-g per record (e.g., `chr7:g.117199644G>A`)

## Aggregated Sources: ~21 and Counting

| Source | What | Notes |
|--------|------|-------|
| ClinVar | Pathogenicity | Weekly refresh |
| gnomAD v4 exomes + genomes | Population AF | grpmax_faf95 surfaced |
| dbSNP Build 156 | rsID + alleles | RsMergeArch resolved |
| dbNSFP v4.x | Meta-aggregator of 40+ in silico predictors | Includes AlphaMissense, REVEL, BayesDel |
| CADD | Deleteriousness | Genome-wide |
| CIViC | Cancer interpretation | Per-disease |
| COSMIC | Somatic variants | Catalogue of Somatic Mutations |
| EVS | Exome Variant Server | Legacy (deprecated by gnomAD) |
| ExAC | ExAC frequencies | Legacy (superseded by gnomAD) |
| GRASP | GWAS associations | -- |
| GWAS Catalog | Curated GWAS | -- |
| Wellderly | Disease-resistant elderly cohort | -- |
| EMV | -- | -- |
| DOCM | Database of Curated Mutations | -- |
| ICGC | International cancer | -- |
| MutDB | -- | -- |
| GO | Gene Ontology | -- |
| Snpeff | snpEff annotations | -- |
| GeneReviews | Disease/gene reviews | -- |
| MutPred | Functional impact | -- |

**dbNSFP is itself an aggregator.** Querying `dbnsfp.alphamissense.score` returns the version that dbNSFP loaded, not AlphaMissense direct. The lag from publication (Cheng 2023 *Science*) to integration into myvariant.info is typically 6-18 months via dbNSFP.

## Scopes and Query Forms

| Endpoint | Method | Use |
|----------|--------|-----|
| `/v1/variant/{id}` | GET | Single canonical-ID lookup |
| `/v1/variant` | POST (batched IDs) | Batch lookup, up to 1000 IDs |
| `/v1/query?q={lucene}` | GET | Flexible Elasticsearch search |
| `/v1/metadata` | GET | Per-source versions |

**Scopes** (the `scopes` parameter on `/v1/query` POST) specifies which fields to match an input ID against: `hgvs`, `rsid`, `dbsnp.rsid`, `dbnsfp.genename`, `chrom`, `_id`. The `_id` is canonical HGVS-g.

## Reproducibility: The `_meta` Field

Every record carries `_meta.src` showing per-source version:

```python
mv = myvariant.MyVariantInfo()
record = mv.getvariant('chr7:g.140453136A>T', fields=['_meta', 'clinvar', 'dbnsfp.alphamissense'])
print(record['_meta']['src']['dbnsfp']['version'])  # e.g., '4.7a'
print(record['_meta']['src']['clinvar']['version'])  # e.g., '20250901'
```

For reproducibility, record per-source versions in analysis output alongside results.

## Comparison to Alternatives

| Tool | Approach | When to use |
|------|----------|-------------|
| **myvariant.info** | Cloud aggregator, ES-backed | Quick batch annotation, no local setup |
| **OpenCRAVAT** (Pagel 2020 *Cancer Res*) | Local install, modular annotators | Offline / PHI-sensitive |
| **VarSome** (Kopanos 2019 *Bioinformatics* 35:1978; commercial) | Hosted, 22 sources | 82% ACMG criteria auto-application (highest); clinical labs |
| **Franklin / Genoox** | Commercial hosted | 59 data sources; family/cohort analysis |
| **GeneBe.net** (Stawinski 2024 *Clin Genet*) | Open-source web + API | Free Tavtigian-point-system-based ACMG; comparable to VarSome |
| **ANNOVAR / VEP / snpEff** | Local annotation tools | Pipeline integration, batch annotation, no ACMG |

**myvariant.info does NOT produce ACMG calls**; it is purely an annotation aggregator. Pair with InterVar, GeneBe, or the `acmg-classification` skill for classification.

## Decision Tree by Query Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Single variant batch annotation | `getvariant(hgvs, fields=...)` | One call, all aggregated sources |
| 10-1000 variants | `getvariants(list, fields=...)` | Batch endpoint, up to 1000 |
| > 1000 variants | Chunk to 1000 + sleep | Rate limit + JSON size |
| Search by gene + pathogenicity | `mv.query('clinvar.gene.symbol:BRCA1 AND clinvar.clinical_significance:Pathogenic', size=200)` | Elasticsearch Lucene |
| ACMG-grade pipeline | myvariant for annotation -> InterVar / GeneBe for classification | myvariant does not produce ACMG calls |
| Offline / PHI-sensitive | OpenCRAVAT or VEP locally | myvariant requires HTTP |
| Reproducibility | Always record `_meta.src.<source>.version` | dbNSFP version is the dominant staleness vector |
| Source-specific deep dive | Use the source-specific skill (clinvar-lookup, gnomad-frequencies) | myvariant is aggregator-grade, not source-deep |

## Standard Annotation Workflow

**Goal:** Annotate a list of variants with the canonical clinical fields for downstream prioritization.

**Approach:** Batch `getvariants` with explicit field list; record `_meta` versions; convert to DataFrame.

```python
import myvariant
import pandas as pd

mv = myvariant.MyVariantInfo()

CLINICAL_FIELDS = [
    'clinvar.clinical_significance',
    'clinvar.review_status',
    'clinvar.variant_id',
    'gnomad_exome.faf95',
    'gnomad_exome.af.af',
    'gnomad_exome.an.an',
    'gnomad_genome.faf95',
    'gnomad_genome.af.af',
    'dbsnp.rsid',
    'dbnsfp.alphamissense.score',
    'dbnsfp.alphamissense.pred',
    'dbnsfp.revel.score',
    'dbnsfp.cadd.phred',
    'dbnsfp.spliceai.master_pred',
    'dbnsfp.spliceai.ds_max',
    'cosmic.cosmic_id',
    'civic.openCravatUrl',
    '_meta'
]

def annotate_variant_list(hgvs_list):
    '''Batch-annotate variants with ClinVar / gnomAD / dbNSFP / COSMIC / CIViC fields.'''
    chunked = [hgvs_list[i:i+1000] for i in range(0, len(hgvs_list), 1000)]
    rows = []
    versions = None
    for chunk in chunked:
        results = mv.getvariants(chunk, fields=CLINICAL_FIELDS)
        for r in results:
            if versions is None and r.get('_meta'):
                versions = {src: meta.get('version') for src, meta in r['_meta'].get('src', {}).items()}
            clinvar = r.get('clinvar', {}) or {}
            gnomad_e = r.get('gnomad_exome', {}) or {}
            gnomad_g = r.get('gnomad_genome', {}) or {}
            dbnsfp = r.get('dbnsfp', {}) or {}
            faf95 = (gnomad_e.get('faf95', {}) or gnomad_g.get('faf95', {})) or {}
            rows.append({
                'variant': r.get('query'),
                'clinvar_sig': clinvar.get('clinical_significance'),
                'clinvar_review': clinvar.get('review_status'),
                'gnomad_grpmax_faf95': faf95.get('popmax'),
                'grpmax_ancestry': faf95.get('popmax_population'),
                'gnomad_af': gnomad_e.get('af', {}).get('af') or gnomad_g.get('af', {}).get('af'),
                'rsid': r.get('dbsnp', {}).get('rsid'),
                'alphamissense': dbnsfp.get('alphamissense', {}).get('score'),
                'revel': dbnsfp.get('revel', {}).get('score'),
                'cadd_phred': dbnsfp.get('cadd', {}).get('phred'),
                'spliceai_ds_max': dbnsfp.get('spliceai', {}).get('ds_max')
            })
    return pd.DataFrame(rows), versions
```

## Elasticsearch Query Patterns

**Goal:** Search beyond canonical IDs; e.g., all pathogenic variants in a gene, all variants in a genomic region with CADD > 20.

**Approach:** Lucene syntax in `mv.query()`; support boolean operators, ranges, wildcards.

```python
def find_pathogenic_in_gene(gene_symbol, max_results=500):
    '''Find ClinVar P/LP variants in a gene.'''
    query = f'clinvar.gene.symbol:{gene_symbol} AND '\
            'clinvar.clinical_significance:(Pathogenic OR "Likely pathogenic")'
    hits = mv.query(query, size=max_results, fields=['_id', 'clinvar.clinical_significance',
                                                       'clinvar.review_status'])
    return hits.get('hits', [])


def find_high_cadd_in_region(chrom, start, end, min_cadd=25):
    '''Find variants in region with CADD phred above threshold.'''
    query = f'chrom:{chrom} AND hg19.start:[{start} TO {end}] AND '\
            f'dbnsfp.cadd.phred:>{min_cadd}'
    return mv.query(query, size=500, fields=['_id', 'dbnsfp.cadd.phred', 'clinvar.clinical_significance'])


def find_alphamissense_pathogenic(gene, min_score=0.564):
    '''Find AlphaMissense pathogenic missense in a gene.

    Note: Cheng 2023 dev cutoff is 0.564 BUT this is NOT the Pejaver-style calibrated
    PP3 threshold. ClinGen has not endorsed AlphaMissense thresholds as of May 2026;
    use AlphaMissense as supporting evidence only.
    '''
    query = f'dbnsfp.genename:{gene} AND dbnsfp.alphamissense.score:>{min_score}'
    return mv.query(query, size=500, fields=['_id', 'dbnsfp.alphamissense', 'clinvar.clinical_significance'])
```

## Per-Operation Failure Modes

**1. Stale dbNSFP version**
- Trigger: Query `dbnsfp.alphamissense.score` and trust as current.
- Mechanism: dbNSFP version lags new tools by 6-18 months; AlphaMissense (Cheng 2023) was integrated into dbNSFP 4.4 (~2024).
- Symptom: Predictions appear missing for variants newly scored by AlphaMissense.
- Fix: Check `_meta.src.dbnsfp.version`; for cutting-edge predictions query AlphaMissense API directly.

**2. Treating AlphaMissense dev threshold as PP3-calibrated**
- Trigger: Apply AlphaMissense score >0.564 as PP3 evidence.
- Mechanism: Cheng 2023 developer-recommended threshold is NOT the Pejaver 2022 calibration framework; ClinGen has not endorsed strength-graded thresholds for AlphaMissense as of May 2026.
- Symptom: Over-application of PP3 in ACMG classification.
- Fix: Treat AlphaMissense as supporting evidence; defer to `clinical-databases/acmg-classification` for calibrated thresholds.

**3. Stacking REVEL + BayesDel + AlphaMissense as independent evidence**
- Trigger: Use multiple in silico predictors as additive PP3 evidence.
- Mechanism: REVEL, BayesDel, VEST4 share ClinVar/HGMD training labels; AlphaMissense is partially independent but correlates strongly with conservation.
- Symptom: Inflated PP3 strength; double-counting.
- Fix: Apply ONE predictor per variant (Pejaver 2022 explicit recommendation).

**4. Rate-limit ignorance**
- Trigger: Loop sequentially over 10k variants with `getvariant()`.
- Mechanism: ~1000 req/sec aggregate per source IP; loops trigger throttling or 429s.
- Symptom: Increasing latency, eventual failures.
- Fix: Use `getvariants(chunk, fields=...)` with chunk size 1000; sleep ~0.5s between chunks.

**5. Field-path errors silently return None**
- Trigger: Query `gnomad_exome.faf95.popmax` but typo as `gnomad_exome.faf` or `gnomad.exomes.faf95`.
- Mechanism: Elasticsearch returns null for non-existent paths; no error raised.
- Symptom: All values None; no error.
- Fix: Use `print(mv.getvariant(test_id))` first to inspect actual field structure; check `/v1/metadata/fields`.

**6. Multi-allelic rsID returns one variant only**
- Trigger: Query `rs12345` and treat returned variant as the variant of interest.
- Mechanism: rsID is a cluster identifier; multi-allelic clusters return multiple records.
- Symptom: Wrong allele returned for ~6-8% of rsIDs.
- Fix: Use HGVS-g instead of rsID for unambiguous lookups; for rsID queries, inspect all returned variants and filter by allele.

**7. Sample overlap between sources**
- Trigger: Treat ClinVar + gnomAD as independent corroboration.
- Mechanism: ClinVar variants often *derive from* gnomAD-included individuals; not statistically independent.
- Symptom: False sense of independent evidence in ACMG application.
- Fix: ClinVar P + gnomAD rare is operationally additive evidence per ACMG but understand the populations may overlap.

## Reconciliation: When Sources Inside myvariant Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| dbNSFP REVEL != ClinVar PP3 strength | Different curation cohort | Use Pejaver 2022 calibrated thresholds (see `acmg-classification`) |
| ClinVar P + AlphaMissense benign | NMD-escape region, alternative isoform, ClinVar P stale | Cross-check with conservation, splicing predictions |
| gnomAD AF differs across exome vs genome | Sample sizes differ; exome has 730k, genome 76k | Use exome FAF95 when available; genome as fallback |
| COSMIC + ClinVar overlap | True dual-classification (germline + somatic) | Report both contexts |
| Variant missing from one source | Source-specific coverage gaps | Cross-check directly with primary source skill (clinvar-lookup, gnomad-frequencies) |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| Batch endpoint cap | 1000 IDs per POST | myvariant.info docs |
| Rate limit | ~1000 req/sec aggregate; lower per IP | myvariant.info docs |
| dbNSFP refresh lag | 6-18 months from primary source release | dbNSFP release history |
| `_meta.src` field | Per-source version is always available | BioThings SDK convention |
| Lucene escape | Special chars need `\` (e.g., `chr7\:140453136`) | Elasticsearch convention |
| Multi-allelic rsID | ~6-8% of dbSNP rsIDs are multi-allelic | Phan 2025 *NAR* |
| AlphaMissense PP3 calibration | NOT yet ClinGen-endorsed (as of May 2026) | ClinGen SVI |
| REVEL PP3_Strong calibration | >= 0.932 per Pejaver 2022 | Pejaver 2022 *AJHG* |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| `KeyError: 'gnomad_exome'` | Variant absent from gnomAD exome dataset | Use `.get('gnomad_exome', {})` defensively |
| `None` for AlphaMissense on rare variants | dbNSFP coverage gap; variant in alt-spliced isoform | Query AlphaMissense API directly, or accept None |
| Search returns 0 hits despite known matches | Lucene escape on `:` in chr coords | Quote the chrom-position term or escape `:` |
| Batch returns < input IDs | Some IDs not in any source | Check `notfound` field in response |
| Different AF in myvariant vs gnomAD browser | dbNSFP version != current gnomAD release | Check `_meta.src.gnomad_exome.version` |
| 503 on bulk query | Rate limit | Reduce chunk to 500; sleep 1s between |
| `_id` doesn't match input | myvariant uses canonical HGVS-g; input was rsID or non-canonical | Re-query by `_id` after first resolution |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "myvariant.info is just an aggregator; why not query sources directly?" | Aggregator avoids per-source API setup; sufficient for single-source unique queries we defer to source-specific skills. |
| "This annotation differs from VarSome" | VarSome uses its own ACMG implementation; myvariant.info does NOT produce ACMG calls; we pair with `acmg-classification`. |
| "dbNSFP REVEL differs from REVEL website" | dbNSFP version is on the order of 1 year behind primary; check `_meta.src.dbnsfp.version`. |
| "AlphaMissense calibration thresholds were missed" | AlphaMissense is integrated via dbNSFP; PP3 calibration is in `clinical-databases/acmg-classification` skill. |
| "Why not OpenCRAVAT?" | OpenCRAVAT requires local install; myvariant is faster for batch annotation. Switch to OpenCRAVAT for PHI-sensitive or offline workflows. |

## References

- Wu C et al. 2022. BioThings SDK: a toolkit for building high-performance data APIs in biomedical research. *Bioinformatics* 38:2316.
- Xin J et al. 2016. High-performance web services for querying gene and variant annotation. *Genome Biol* 17:91.
- Cheng J et al. 2023. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 381:eadg7492.
- Pejaver V et al. 2022. Calibration of computational tools for missense variant pathogenicity classification. *Am J Hum Genet* 109:2163.
- Pagel KA et al. 2020. Integrated informatics analysis of cancer-related variants. *Cancer Res* 80:e29. (OpenCRAVAT)
- Kopanos C et al. 2019. VarSome: the human genomic variant search engine. *Bioinformatics* 35:1978.
- Stawinski PM et al. 2024. GeneBe; a tool for automated ACMG/AMP variant interpretation. *Clin Genet* (verify exact volume/pages in the published record before citing).
- myvariant.info docs: `https://docs.myvariant.info/en/latest/`
- BioThings field metadata: `https://myvariant.info/v1/metadata/fields`

## Related Skills

- clinical-databases/clinvar-lookup - Source-level deep ClinVar queries
- clinical-databases/gnomad-frequencies - Source-level deep gnomAD queries
- clinical-databases/dbsnp-queries - Source-level rsID resolution
- clinical-databases/acmg-classification - ACMG framework with Pejaver calibration
- clinical-databases/variant-prioritization - Pipeline using aggregated annotations
