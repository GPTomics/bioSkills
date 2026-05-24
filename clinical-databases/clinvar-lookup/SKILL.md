---
name: bio-clinical-databases-clinvar-lookup
description: Queries ClinVar for variant pathogenicity classifications, ClinGen VCEP curations, and somatic-vs-germline interpretations via REST API, weekly VCF, or bulk XML. Use when determining clinical significance, triangulating conflicting interpretations, or aggregating evidence against the ACMG/AMP framework with ClinGen SVI specifications.
tool_type: python
primary_tool: requests
---

## Version Compatibility

Reference examples tested with: requests 2.31+, cyvcf2 0.30+, pandas 2.2+, bcftools 1.19+, Entrez Direct 21.0+, lxml 5.0+ (for v2 XML schema).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. ClinVar XML schema v2 (rolled out in 2024) replaces `<ClinVarSet>` with `<VariationArchive>` as the top-level anchor; XSLT or parsers targeting the legacy element silently emit zero records.

# ClinVar Lookup and Clinical-Significance Triangulation

**'Look up the clinical significance of this variant'** -> Retrieve ClinVar VCV-level aggregate, SCV-level submissions, ClinGen Variant Curation Expert Panel (VCEP) overrides, and conflict-resolution status.

- Python (REST): `requests.get()` against the E-utilities `clinvar` database
- Python (local VCF): `cyvcf2.VCF('clinvar.vcf.gz')` for batch queries against the weekly snapshot
- CLI: `bcftools annotate -a clinvar.vcf.gz -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN`
- Cross-database: ClinGen Allele Registry CA ID via `https://reg.clinicalgenome.org/`

## The Identifier Hierarchy (VCV / SCV / RCV): Get This Wrong and Everything Downstream Breaks

| Level | Format | What it aggregates | When to use | Fails when |
|-------|--------|--------------------|-------------|-----------|
| **SCV** | `SCVxxxxxxxxx.N` | One submitter, one variant, one condition (atomic submission unit) | Auditing who said what; conflict triangulation | Aggregated reporting (use VCV); cross-condition analysis |
| **RCV** | `RCVxxxxxxxxx.N` | All SCVs for a single (variant, condition) pair | Condition-stratified analysis; legacy aggregation | Variant-level reporting across all conditions (use VCV) |
| **VCV** | `VCVxxxxxxxxx.N` | All RCVs for one variant across all conditions | Canonical anchor since 2017; default API entrypoint | Condition-specific clinical action (use RCV); CLNSIG collapses multi-condition |

**Operational footgun:** the `clinvar.vcf.gz` `CLNSIG` field is the *variant-level* (VCV) aggregate. A variant Pathogenic for disease A but VUS for disease B collapses to "Pathogenic/Conflicting". For condition-stratified analysis, parse RCV-level XML, never `CLNSIG` alone.

**2024 XML schema overhaul:** ClinVar v2 XML separates `GermlineClassification`, `SomaticClinicalImpact`, and `OncogenicityClassification` under one `<VariationArchive>` anchor. The legacy `<ClinicalSignificance>` element is gone. Pipelines built before September 2024 against `<ClinVarSet>` silently emit zero records on new XML. The dual-release period ended December 2024.

## Star Ratings and the Override Hierarchy

| Stars | Review status | What it means operationally |
|-------|--------------|---------------------------|
| 4 | Practice guideline | ACMG/CAP CFTR-level (vanishingly rare) |
| 3 | Expert panel reviewed (ClinGen VCEP) | **FDA-recognized tier**; overrides lower-star records for clinical action |
| 2 | Multiple submitters, criteria provided, no conflicts | Reliable aggregate |
| 1 | Single submitter OR conflicting interpretations (often mis-reported as 2-star) | Use with scrutiny |
| 0 | No assertion criteria provided | Literature-only or legacy submissions |

ClinVar does NOT retract or hide lower-star records when a VCEP publishes; a variant can simultaneously display "Pathogenic (3-star VCEP)" and "Conflicting interpretations (1-star)". Tools handle this differently (VarSeq, Franklin, GenoOx each pick a winner via different rules); this is a major source of inter-tool disagreement.

## ClinGen Variant Curation Expert Panels (VCEPs)

As of 2025, ~80-90 VCEPs are approved or in progress across RASopathies, hereditary cancer (ENIGMA BRCA1/2, InSiGHT MMR), cardiomyopathy (sarcomere genes), hearing loss, RPE65/IRD, inborn errors of metabolism, and FH. The current count is moving; the authoritative directory is the Criteria Specification Registry at `https://cspec.genome.network/cspec/ui/svi/all`.

Each VCEP publishes a **gene-disease-specific CSpec** that re-weights ACMG/AMP criteria. The Hearing Loss VCEP downgrades PM2 to supporting by default and upgrades PS3 thresholds for OTOF. Treating "ACMG/AMP" as a single rubric across all genes is the most common error in non-specialist tooling.

## ACMG/AMP, ClinGen SVI Specifications, and the Bayesian Point System

The Richards 2015 28-criterion framework is the foundation, but **every modern automated classifier (InterVar, GeneBe, Franklin, VarSome) implements the Tavtigian 2018/2020 Bayesian point system**, not the original combining rules. Strengths map to points: Supporting=1, Moderate=2, Strong=4, Very Strong=8 (benign codes negative). Final categories: P >=10, LP 6-9, VUS 0-5, LB -1 to -6, B <=-7.

For variant interpretation framework details, calibrated in-silico thresholds, and PVS1 decision-tree logic, defer to `clinical-databases/acmg-classification`. This skill focuses on querying ClinVar; it intentionally does not re-implement classification.

## Conflicting Interpretations and Conflict Resolution

Harrison 2017 *Genet Med* 19:1096 (PMID 28301460) showed 87% of inter-lab conflicts were resolvable by reassessment plus data sharing. As of 2024, only 3.8% of conflicting BRCA1 missense VUS reached consensus despite years of effort; conflict resolution is slow even in best-curated genes.

**Submission staleness** is non-trivial: ClinVar does not push reclassifications to submitters; a 2017 SCV can persist on an active label in 2026 if the lab has not re-submitted. Genome Alert! (Yauy 2022 *Genet Med*) was built specifically to detect classification drift between weekly releases. The median delta is ~1,247 classification changes per month with potential clinical impact.

## Decision Tree by Query Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Single variant, known gene/condition | E-utilities `esummary` against `clinvar` DB | Lowest latency, returns VCV-level summary |
| Batch (10-1000 variants) by HGVS or rsID | myvariant.info with `fields=clinvar` | Aggregated, includes ClinVar review status |
| Batch (>1000) or coordinate-based | Local `clinvar.vcf.gz` with `bcftools annotate` or `cyvcf2` | No rate limits; weekly snapshot |
| Condition-stratified (variant in disease A vs B) | Bulk XML `VariationArchive` parsing | RCV is the only level that preserves per-condition classification |
| Cross-database join with gnomAD / dbSNP / COSMIC | ClinGen Allele Registry CA ID | Build-agnostic, transcript-agnostic canonical identifier |
| Reproducible analysis with citable date | First-Thursday-of-month archive on FTP | Only monthly snapshots are archived; weekly releases disappear |

## REST API Query (E-utilities)

**Goal:** Retrieve VCV-level ClinVar summary for a single variant by ID, gene, or HGVS.

**Approach:** Hit `esummary.fcgi` or `esearch.fcgi` against `db=clinvar`, parse JSON, then optionally hydrate to full record with `efetch`.

```python
import requests

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'

def clinvar_summary(variation_id):
    '''Retrieve VCV-level summary by ClinVar VariationID (do not confuse with CA ID).

    The germline / somatic / oncogenicity classification nesting shown below
    follows the ClinVar 2024 eSummary v2 schema described in the data-access
    documentation. Field names have changed between API versions -- inspect
    the actual JSON returned by eSummary for the live ClinVar version before
    pinning these key paths in production.
    '''
    r = requests.get(f'{EUTILS}/esummary.fcgi',
                     params={'db': 'clinvar', 'id': variation_id, 'retmode': 'json'},
                     timeout=30)
    r.raise_for_status()
    record = r.json()['result'][str(variation_id)]
    return {
        'vcv': record.get('accession'),
        'name': record.get('title'),
        'germline_class': record.get('germline_classification', {}).get('description'),
        'germline_review_status': record.get('germline_classification', {}).get('review_status'),
        'somatic_clinical': record.get('clinical_impact_classification', {}).get('description'),
        'oncogenicity': record.get('oncogenicity_classification', {}).get('description'),
        'last_evaluated': record.get('germline_classification', {}).get('last_evaluated')
    }

def clinvar_search_gene(gene, pathogenic_only=False, retmax=500):
    term = f'{gene}[gene]'
    if pathogenic_only:
        term += ' AND (clinsig_pathogenic[Properties] OR clinsig_likely_pathogenic[Properties])'
    r = requests.get(f'{EUTILS}/esearch.fcgi',
                     params={'db': 'clinvar', 'term': term, 'retmax': retmax, 'retmode': 'json'},
                     timeout=30)
    return r.json()['esearchresult']['idlist']
```

## Local VCF Query (Weekly Snapshot)

**Goal:** Annotate or look up thousands of variants without rate limits.

**Approach:** Download the weekly `clinvar.vcf.gz` (note: only first-Thursday-of-month is archived; for longitudinal stability pin to monthly archives), query by genomic coordinates with cyvcf2 or annotate VCFs with bcftools.

```bash
mkdir -p clinvar/$(date +%Y%m); cd clinvar/$(date +%Y%m)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

```python
from cyvcf2 import VCF

clinvar = VCF('clinvar.vcf.gz')

def lookup(chrom, pos, ref, alt):
    '''Look up by GRCh38 coords. Returns variant-level (VCV) aggregate; not RCV.'''
    for v in clinvar(f'{chrom}:{pos}-{pos}'):
        if v.REF == ref and alt in v.ALT:
            info = v.INFO
            return {
                'vcv_id': info.get('ALLELEID'),
                'clnsig': info.get('CLNSIG'),
                'clnsig_conf': info.get('CLNSIGCONF'),
                'clnrevstat': info.get('CLNREVSTAT'),
                'clndn': info.get('CLNDN'),
                'clnvc': info.get('CLNVC'),
                'clnhgvs': info.get('CLNHGVS'),
                'clndisdb': info.get('CLNDISDB'),
                'oncdn': info.get('ONCDN'),
                'scidn': info.get('SCIDN')
            }
    return None
```

```bash
bcftools annotate \
    -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVC,INFO/CLNHGVS,INFO/CLNSIGCONF \
    input.vcf.gz -O z -o annotated.vcf.gz
bcftools index -t annotated.vcf.gz
```

## ClinGen Allele Registry (CA IDs): The Real Cross-Database Anchor

ClinGen Allele Registry (`https://reg.clinicalgenome.org/`) computes a build-agnostic, transcript-agnostic CA ID (format `CA######`) for any allele projectable onto NCBI references (GRCh37, GRCh38, T2T-CHM13, any RefSeq transcript). The Registry covers ~700M+ alleles, vastly more than ClinVar. CA ID and ClinVar VariationID are one-to-one *when a variant exists in ClinVar*.

```python
def car_id(hgvs_g):
    '''Resolve HGVS-g to canonical ClinGen Allele Registry CA ID.'''
    r = requests.put(f'https://reg.clinicalgenome.org/allele',
                     headers={'Content-Type': 'text/plain'},
                     data=hgvs_g, timeout=30)
    return r.json().get('@id', '').split('/')[-1] if r.ok else None
```

Use CA ID for any join touching non-ClinVar resources (gnomAD, dbSNP, COSMIC, MAVEdb). VariationID was renumbered during the 2017 ClinVar schema redesign; treating it as a stable cross-build identifier is unsafe.

## Per-Operation Failure Modes

**1. Treating CLNSIG as gospel for condition-specific work**
- Trigger: Pull `CLNSIG=Pathogenic` from `clinvar.vcf.gz` and report variant as pathogenic for the patient's specific phenotype.
- Mechanism: CLNSIG is VCV-level aggregate; a variant can be P for disease A and B for disease B.
- Symptom: Patient phenotype does not match the disease where the variant is actually pathogenic; clinical action is wrong.
- Fix: Parse RCV-level XML (`<RCVAccession>` per condition); cross-check `CLNDN` and report per-condition classifications.

**2. Parsing legacy XML against 2024 schema**
- Trigger: XSLT or parser anchored on `<ClinVarSet>` or `<ClinicalSignificance>`.
- Mechanism: 2024 schema replaces both anchors with `<VariationArchive>` + germline/somatic/oncogenicity tripartite classifications.
- Symptom: Silent zero-record output, no error.
- Fix: Re-target to `<VariationArchive>` and read `GermlineClassification`, `SomaticClinicalImpact`, `OncogenicityClassification` separately.

**3. Counting `variant_summary.txt` rows naively**
- Trigger: `wc -l variant_summary.txt` to estimate variant count.
- Mechanism: One row per assembly per variant (GRCh37 *and* GRCh38); double-counts.
- Symptom: Counts inflated ~2x.
- Fix: `awk -F'\t' '$17=="GRCh38"' variant_summary.txt | wc -l`.

**4. Trusting VariationID as a stable cross-build identifier**
- Trigger: Join gnomAD-v4 records by ClinVar VariationID.
- Mechanism: VariationIDs were renumbered for a subset of variants during the 2017 schema overhaul.
- Symptom: Spurious mismatches at low rate (~1-2%).
- Fix: Use ClinGen Allele Registry CA ID, which is computed deterministically from sequence.

**5. Ignoring star-rating override hierarchy**
- Trigger: Pipeline picks the most-recent SCV regardless of review status.
- Mechanism: A 1-star SCV submitted yesterday outranks a 3-star VCEP curation from 2022 by date.
- Symptom: Clinical reports cite outdated or non-expert assertions over expert-panel decisions.
- Fix: Sort by `review_status` rank (4>3>2>1>0); use the highest-star record. For ties, sort by date.

**6. Aggregating "Conflicting" without inspecting the conflict**
- Trigger: Treat `CLNSIG=Conflicting interpretations` as VUS.
- Mechanism: "Conflicting" can mean (P vs LP), (P vs VUS), or (P vs LB); the clinical meaning is completely different across these.
- Symptom: Patients with high-evidence P variants reported as ambiguous; or true VUS reported as actionable.
- Fix: Parse `CLNSIGCONF` to see exact conflict; weight by submitter star.

**7. Missing somatic interpretations**
- Trigger: Pre-2024 pipeline reads only `CLNSIG`.
- Mechanism: Somatic classifications live in new INFO fields (`ONCDN`, `SCIDN`, `CLNSIGSOMATIC`) since 2024.
- Symptom: Cancer variants appear unclassified.
- Fix: Read `ONCDN` (oncogenicity disease name), `SCIDN` (somatic clinical impact disease name), and the somatic-specific significance fields.

## Reconciliation: When Sources Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| ClinVar P vs gnomAD AF > 1% | Variant is true founder allele in unstratified gnomAD subset, OR ClinVar P is a stale low-star assertion | Check `grpmax_faf95` excluding bottleneck groups; check ClinVar star rating |
| ClinVar P vs AlphaMissense < 0.1 | Variant in NMD-escape region, alternative isoform, or ClinVar P is mis-curated | Check Pejaver 2022 calibration in `acmg-classification` skill; cross-check VCEP |
| VCEP 3-star P vs commercial-lab 1-star B | VCEP supersedes for clinical action | Use VCEP; flag submitter for resubmission |
| ClinVar VCV-level P vs RCV-level VUS for actual condition | VCV averages across conditions | Always report at RCV level for clinical action |
| ClinVar P vs LOVD/HGMD discordant | LOVD/HGMD use different classification systems; HGMD "DM" != ACMG P | Triangulate against published evidence; do not auto-translate labels |
| ClinVar P missing for a known disease variant | Submission lag (~6-12 months typical for new findings) | Check published literature; flag for ClinVar submission |

## Quantitative Thresholds and Operational Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| Star >= 2 | Acceptable confidence for clinical action without further review | ClinGen SVI operational guidance |
| Star = 3 | VCEP-curated; supersedes lower-star records | ClinGen FDA Recognition 2018 |
| `CLNSIG` includes 'Pathogenic' OR 'Likely_pathogenic' | Treat as actionable for ACMG | ClinVar field schema |
| `CLNSIGCONF` present | Multiple SCVs disagree; do NOT auto-action | ClinVar field schema |
| Monthly archive | Use first-Thursday-of-month FTP snapshot for reproducible analyses | NCBI FTP retention policy |
| Submission staleness | Re-check classification annually for active diagnostic variants | Yauy 2022 *Genet Med* (Genome Alert!) |
| AF > 5% in gnomAD | BA1 standalone benign per ClinGen SVI default (VCEP overrides exist) | Richards 2015; SVI specs |
| 1247 changes/month | Median variants with classification change per release | Yauy 2022 |

## ClinVar Somatic vs Germline: 2024 Tripartite

The 2024 schema separates three orthogonal classifications, each with its own `ReviewStatus` and `DateLastEvaluated`:

- **GermlineClassification:** Pathogenic / Likely Pathogenic / VUS / LB / B per ACMG/AMP 2015 + SVI.
- **SomaticClinicalImpact:** Tier I / II / III / IV per AMP/ASCO/CAP 2017 (Li 2017 *J Mol Diagn*).
- **OncogenicityClassification:** Oncogenic / Likely Oncogenic / VUS / Likely Benign / Benign per ClinGen/CGC/VICC 2022 oncogenicity framework.

A single VCV can carry all three with distinct evaluations; the legacy "Pathogenic" label is now ambiguous if not qualified by classification type.

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| Empty result from `efetch db=clinvar` | rsID passed where VariationID expected | Use `esearch` first to resolve rsID to VariationID |
| `CLNSIG` is `_None` or comma-separated mess | Variant has multi-condition RCVs; VCF collapses them | Parse RCV XML for per-condition values |
| Variant present in ClinVar XML but absent from VCF | Variant lacks GRCh38 coordinates (legacy GRCh37-only submission) | Check `<VariationArchive><SequenceLocation>` per assembly |
| 2024-format XML parser silently emits zero records | XML schema v2 incompatibility | Re-target to `<VariationArchive>` |
| Conflicting interpretations with same star rating across two submitters | True scientific disagreement; sometimes resolved by VCEP later | Apply Tavtigian point system to manually reconcile; flag for VCEP review |
| Variant has CA ID but no VariationID | Variant in Allele Registry but never submitted to ClinVar | Use AlleleRegistry as canonical; submit to ClinVar if novel pathogenic |
| `CLNSIG` says Pathogenic but no associated condition `CLNDN` | Orphan classification (older submissions) | Treat as low confidence; cross-check publication |
| Variants pulled by gene return only some isoforms | RefSeq transcript priority differences | Use MANE Select transcript explicitly; cross-check with VEP `--mane_select` |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why is your pathogenic variant 1-star?" | We report star rating per submission; clinical action requires star >=2 OR VCEP curation per ClinGen SVI 2018. |
| "ClinVar says P but gnomAD AF = 2%" | Reconciled via Whiffin FAF95 max-credible-AF framework; bottleneck-group rule applied. |
| "Your VCV count differs from ClinVar.gov" | We pulled from the monthly archive (first-Thursday-of-month) for reproducibility; the live web is post-most-recent-weekly. |
| "Why didn't you flag the somatic variant?" | Pre-2024 XML schema had no separate somatic field; we now read `ONCDN`/`SCIDN`/`SomaticClinicalImpact` per v2 schema. |
| "VarSome says LP but you say VUS" | Tool-specific aggregation rule differences; VarSome auto-applies PP3+PM2 by default per Tavtigian point system; we apply VCEP-specific PP3 calibration per CSpec. |
| "rsID match returned wrong variant" | rsID is a cluster identifier; multi-allelic rsIDs require allele-level resolution; we use SPDI or CA ID. |
| "Why retest a 2022-curated variant?" | Submission staleness median 5-year reclassification cycle in active genes (Harrison 2017); ClinGen recommends annual re-review for active diagnostic variants. |

## References

- Landrum MJ et al. 2025. ClinVar: improvements to accessing data. *Nucleic Acids Res* 53(D1):D1313 (2025 NAR DB issue).
- Harrison SM et al. 2017. Clinical laboratories collaborate to resolve differences in variant interpretations submitted to ClinVar. *Genet Med* 19:1096.
- Yauy K et al. 2022. Genome Alert! a standardized procedure for genomic variant reinterpretation and automated genotype-phenotype reassessment. *Genet Med* 24:S1098.
- Tavtigian SV et al. 2018. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med* 20:1054.
- Tavtigian SV et al. 2020. Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Hum Mutat* 41:1734.
- Richards S et al. 2015. Standards and guidelines for the interpretation of sequence variants. *Genet Med* 17:405. (Original ACMG/AMP 2015)
- Pejaver V et al. 2022. Calibration of computational tools for missense variant pathogenicity classification. *Am J Hum Genet* 109:2163.
- Abou Tayoun AN et al. 2018. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. *Hum Mutat* 39:1517.
- Brnich SE et al. 2020. Recommendations for application of the functional evidence PS3/BS3 criterion using the ACMG/AMP sequence variant interpretation framework. *Genome Med* 12:3.
- Li MM et al. 2017. Standards and guidelines for the interpretation and reporting of sequence variants in cancer. *J Mol Diagn* 19:4. (AMP/ASCO/CAP somatic framework)
- ClinGen Allele Registry: `https://reg.clinicalgenome.org/docs/cg-car/`
- CSpec Registry: `https://cspec.genome.network/cspec/ui/svi/all`

## Related Skills

- clinical-databases/acmg-classification - ACMG/AMP framework, PVS1 decision tree, Pejaver calibrated PP3/BP4 thresholds, Bayesian point system
- clinical-databases/myvariant-queries - Aggregated queries including ClinVar overlay
- clinical-databases/variant-prioritization - Rare-disease filtering pipeline using ClinVar
- clinical-databases/gnomad-frequencies - Population frequency for BS1/BA1 cross-check
- variant-calling/clinical-interpretation - Clinical reporting workflow
