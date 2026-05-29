---
name: bio-clinical-databases-gnomad-frequencies
description: Queries gnomAD v4 (807k samples), v3, v2.1.1, and constraint metrics with grpmax FAF95, bottleneck-group exclusion, LOEUF interpretation, SV/CNV/mtDNA catalogs, and Whiffin max-credible-AF framework. Use when filtering rare variants, applying ACMG BS1/BA1, ranking genes by LoF intolerance, or selecting between v2 (GRCh37 + chrX/Y constraint) and v4 (GRCh38 + 807k samples).
tool_type: python
primary_tool: requests
---

## Version Compatibility

Reference examples tested with: requests 2.31+, hail 0.2.130+, pandas 2.2+, myvariant 1.0+. Current gnomAD release is **v4.1 (May 2024)**; v4.1 fixed the v4.0 AN under-counting issue that inflated rare-variant AF estimates by 5-10%.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- Hail: `hl.version()`; pin to >=0.2.130 for v4 schema

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. The gnomAD browser GraphQL API at `https://gnomad.broadinstitute.org/api` is the supported public endpoint; Hail Tables on Google Cloud Storage at `gs://gcp-public-data--gnomad/` are the supported bulk access.

# gnomAD Frequency Queries and Constraint

**'How rare is this variant in the general population?'** -> Pull allele frequency, grpmax FAF95 (the ACMG-grade frequency), LOEUF gene-level constraint, structural variant catalog, mtDNA frequencies, and the appropriate dataset version per use case.

- Python (single variant): GraphQL via `requests.post('https://gnomad.broadinstitute.org/api', json={'query': ..., 'variables': ...})`
- Python (aggregator): `myvariant.MyVariantInfo().getvariant(hgvs, fields=['gnomad_exome', 'gnomad_genome'])`
- Python (bulk): `hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/exomes/gnomad.exomes.v4.1.sites.ht')`

## v2.1.1 / v3.1.2 / v4.x: When to Use Which

This is the most consequential decision in any gnomAD query. The releases are **not interchangeable**; choice determines what can and cannot be said about a variant.

| Release | Build | Samples | Use when | Fails when |
|---------|-------|---------|----------|-----------|
| **v2.1.1** | GRCh37 | 125,748 exomes + 15,708 genomes | Constraint metrics needed (LOEUF v2 most-validated); chrX/Y constraint required; GRCh37 native non-negotiable | GRCh38 native cohort; modern rare-variant FAF95 (use v4) |
| **v3.1.2** | GRCh38 | 76,156 genomes (NO exomes) | Non-coding region rare variants on GRCh38; mtDNA frequencies | Exome variants needed (no exomes); 76k cohort smaller than v4 |
| **v4.0/v4.1** | GRCh38 | 730,947 exomes + 76,215 genomes = **807,162 total** | Default for everything; rare-variant filtering, FAF95, gene queries | chrX/Y constraint (not released); cancer-cohort analysis (no TCGA in v4) |

**Critical caveats:**
- v4 genomes are the SAME 76,215 v3 samples reprocessed against GRCh38 with updated pipelines; not independent.
- ~81% of v2 genomes are also in v3; joint v2+v3 meta-analysis must dedupe at sample ID.
- v4 includes 416,555 UK Biobank exomes under a specific collaboration agreement; check use terms.
- **v4 does NOT include TCGA**, so the `non_cancer` subset is unnecessary; the v4 subset is `non_ukb` (excludes UKB exomes for ancestry rebalancing).
- Liftover v2 (GRCh37 -> GRCh38) is NOT equivalent to v4 native; variant representation differs at ~0.5-1% of sites due to assembly fixes.

## v4 Ancestry Groups: popmax -> grpmax Terminology

v4 ancestry groups: **AFR, AMR, ASJ, EAS, FIN, MID, NFE, SAS, AMI, REMAINING**. The **MID (Middle Eastern) group was new in v4**; previously absorbed into "OTH". The **REMAINING** group (31,256 v4 samples) is individuals who did not cluster with any reference; they contribute to overall AF but not to grpmax.

**Terminology shift:** gnomAD documentation and ACMG-facing narrative uses **grpmax** (genetic ancestry group max) -- replacing the older **popmax** ("population max") term -- to disambiguate genetic ancestry from self-reported race/ethnicity. The public GraphQL schema still exposes legacy field names containing `popmax` (e.g. `faf95.popmax`, `faf95.popmax_population`); these are the grpmax values under the modern terminology. Always check the schema version when writing queries; new browsers may rename these fields.

`grpmax_faf95` is the operational ACMG field. It computes the maximum 95% lower-CI allele frequency, **excluding bottleneck groups** (AMI, ASJ, FIN, REMAINING) because pathogenic founder variants in those groups would otherwise falsely trigger BS1/BA1. MID is included in grpmax but is the smallest non-bottleneck group with highest per-allele variance.

## Filtering Allele Frequency (FAF95): The ACMG-Grade AF

Whiffin 2017 *Genet Med* 19:1151 introduced FAF95 = Poisson lower bound of 95% CI for AF. By construction, AF > FAF95; FAF95 is the conservative frequency for ACMG application.

**Max-credible-AF formula:** `(prevalence x heterogeneity x allelic-contribution) / (penetrance x 2)`. Plug in disease parameters to get the gene-specific BA1 / BS1 threshold; compare against `grpmax_faf95`.

| Code | Threshold | Notes |
|------|-----------|-------|
| **BA1** | AF > 5% in any non-bottleneck group | ClinGen SVI default; VCEPs may override (Hearing Loss VCEP uses 0.5%) |
| **BS1** | AF > gene-specific max-credible-AF | Computed per gene via Whiffin formula |
| **PM2_Supporting** | Absent or ultra-rare in gnomAD | Downgraded from PM2_Moderate in SVI 2020 |

Use `grpmax_faf95`, not raw AF, for BS1/BA1 application; this is the ClinGen-recommended approach.

## Constraint Metrics: pLI, LOEUF, missense Z

Karczewski 2020 *Nature* 581:434 defined LOEUF as the upper bound of the 90% CI of observed/expected pLoF count per gene. LOEUF is **recommended over pLI** because it is continuous and accounts for gene size more rigorously.

| Metric | What | Interpretation |
|--------|------|----------------|
| **LOEUF** | Upper bound of 90% CI of LoF observed/expected ratio | Lower = more LoF-intolerant; **first decile (LOEUF < 0.35 v2; < 0.6 v4) = strongly intolerant** |
| **pLI** | Probability LoF intolerant | Still used; gnomAD team recommends LOEUF for ranking |
| **Missense Z** | Z-score of observed-vs-expected missense | Z > 3.09 = top 1% missense-constrained |
| **Missense O/E** | Observed/expected missense ratio | Continuous form of missense Z |

**Critical version mismatch:**
- v2.1.1 constraint metrics published 2020; v4 constraint published **March 2024** (4 months after v4 data release).
- **v4 constraint is autosomes only; chrX and chrY constraint metrics in v4 are NOT released**. For X/Y constraint, fall back to v2.1.1.
- **LOEUF first decile shifted v2 to v4**: v2 < 0.35; v4 < 0.6 (larger sample shifted the distribution). Gene rank in deciles is stable across versions but absolute thresholds are NOT interchangeable.

## Subsets: non_cancer, non_neuro, controls

| Release | Subset | Removes | Use when |
|---------|--------|---------|----------|
| v2.1.1 | `non_cancer` | TCGA | Cancer-related variant analysis (avoids circularity) |
| v2.1.1 | `non_neuro` | Psychiatric/neuro cohorts | Neuropsychiatric variant analysis |
| v2.1.1 | `controls` | Cases with known disease (~60k samples) | Disease-association calibration |
| v3.1.2 | `non_v2` | v2 overlapping samples | Independent of v2 |
| v3.1.2 | `controls_and_biobanks` | Disease cases retained, biobanks emphasized | Population-level reference |
| v4 | `non_ukb` | UK Biobank exomes | When EUR-skew of UKB problematic |
| v4 | `non_neuro` | Deprecated | -- |
| v4 | `non_cancer` | Unnecessary (no TCGA in v4) | -- |

## SV Catalog and CNV

| Resource | Release | Samples | Coverage |
|----------|---------|---------|----------|
| gnomAD-SV v2 | Collins 2020 *Nature* 581:444 | 14,891 unrelated WGS | 433k SVs, GRCh37 |
| gnomAD-SV v4 | Nov 2023 | 63,046 unrelated WGS | 1,199,117 high-confidence SVs, GRCh38 |
| gnomAD-CNV v4 | Nov 2023 | 464,297 individuals (exome-derived gCNV) | Rare (AF < 1%) autosomal coding CNVs |

gnomAD-CNV v4 is the resource that democratized exome-derived CNV background frequencies; previously only ExAC-CNV provided this at scale.

## mtDNA (Laricchia 2022 *Genome Res* 32:569)

10,850 unique mtDNA variants across 56,434 individuals (v3.1). Frequencies reported per nuclear-ancestry AND per mitochondrial-haplogroup. Heteroplasmy >=10% threshold; ~1/250 individuals carry pathogenic mtDNA variant at heteroplasmy >=10%. mtDNA inheritance is non-Mendelian; standard ACMG criteria do not apply directly; use MITOMAP and HmtVar in parallel.

## VEP Version Pinning

Each gnomAD release pins to a VEP version:
- **v4 uses VEP 105** with GENCODE 39 / Ensembl 105 transcripts
- v2.1.1 uses VEP 85

A variant's consequence prediction can flip between v2 and v4 due to MANE Select adoption and transcript-set updates. Always pin VEP version when reproducing gnomAD annotations.

## Decision Tree by Query Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Single variant AF lookup | GraphQL API or myvariant.info | Lowest latency; returns full per-ancestry breakdown |
| ACMG BS1/BA1 application | `grpmax_faf95` from v4 | The ClinGen-recommended field |
| Gene-level LoF constraint (autosomes) | LOEUF from v4 March 2024 release | Larger sample, more stable |
| Gene-level LoF constraint (chrX/Y) | LOEUF from v2.1.1 | v4 X/Y constraint NOT released |
| Bulk rare-variant filter (cohort-scale) | Hail Table on GCS | No rate limits; full schema |
| SV frequency | gnomAD-SV v4 (WGS) or gnomAD-CNV v4 (exome) | Choose by data type |
| mtDNA frequency | v3.1 mtDNA release (Laricchia 2022) | Only gnomAD release with mtDNA |
| Cancer-variant analysis | v2.1.1 `non_cancer` subset OR v4 (no TCGA) | Avoid TCGA circularity in v2 |
| Comparison across builds | Use canonical SPDI or CA ID, normalize first | Liftover != native |

## Single Variant Query (GraphQL)

**Goal:** Retrieve exome + genome AF, grpmax, FAF95, and per-ancestry breakdown for one variant.

**Approach:** Hit gnomAD's GraphQL API with explicit dataset version; parse the nested response.

```python
import requests

GNOMAD_API = 'https://gnomad.broadinstitute.org/api'

def query_variant(chrom, pos, ref, alt, dataset='gnomad_r4'):
    '''Query gnomAD GraphQL for variant frequency + grpmax FAF95.

    dataset options: gnomad_r4 (v4.1, default), gnomad_r3, gnomad_r2_1
    '''
    query = '''
    query VariantById($variantId: String!, $dataset: DatasetId!) {
      variant(variantId: $variantId, dataset: $dataset) {
        variant_id
        rsids
        exome {
          ac
          an
          af
          homozygote_count
          filters
          populations { id ac an }
          faf95 { popmax popmax_population }
        }
        genome {
          ac
          an
          af
          homozygote_count
          filters
          populations { id ac an }
          faf95 { popmax popmax_population }
        }
      }
    }
    '''
    variant_id = f'{chrom}-{pos}-{ref}-{alt}'
    r = requests.post(GNOMAD_API,
                      json={'query': query, 'variables': {'variantId': variant_id, 'dataset': dataset}},
                      timeout=30)
    r.raise_for_status()
    return r.json().get('data', {}).get('variant')


def grpmax_faf95(payload):
    '''Extract the grpmax FAF95; the ACMG-grade frequency. Excludes bottleneck groups.'''
    exome = payload.get('exome') if payload else None
    if exome and exome.get('faf95'):
        return {
            'faf95': exome['faf95'].get('popmax'),
            'grpmax_ancestry': exome['faf95'].get('popmax_population'),
            'source': 'exome'
        }
    genome = payload.get('genome') if payload else None
    if genome and genome.get('faf95'):
        return {
            'faf95': genome['faf95'].get('popmax'),
            'grpmax_ancestry': genome['faf95'].get('popmax_population'),
            'source': 'genome'
        }
    return {'faf95': 0.0, 'grpmax_ancestry': None, 'source': 'absent'}
```

## ACMG BS1/BA1 Application

**Goal:** Apply Whiffin max-credible-AF framework to a candidate variant.

**Approach:** Compute the gene-specific BS1 threshold from disease parameters, compare to `grpmax_faf95`.

```python
def max_credible_af(prevalence, max_allelic_contribution=1.0, max_genetic_contribution=1.0,
                    penetrance=1.0):
    '''Whiffin 2017 max-credible-AF formula.

    Args:
        prevalence: disease prevalence (e.g., 1/10000 = 1e-4)
        max_allelic_contribution: max contribution of single allele to disease in any case
        max_genetic_contribution: max contribution of this gene to disease in any case
        penetrance: probability that variant carriers develop disease

    Returns: max-credible per-allele frequency under dominant inheritance (use /2 for AR)
    '''
    return (prevalence * max_genetic_contribution * max_allelic_contribution) / (penetrance * 2)


def apply_bs1_ba1(grpmax_faf95_val, max_credible, ba1_threshold=0.05):
    '''Apply ClinGen SVI BS1/BA1 criteria.

    BA1 default 5% per ClinGen SVI; VCEP-specific overrides exist (Hearing Loss = 0.5%).
    BS1 = max-credible-AF specific to gene+disease.
    '''
    if grpmax_faf95_val is None:
        return 'PM2_Supporting'  # Absent or ultra-rare
    if grpmax_faf95_val > ba1_threshold:
        return 'BA1'
    if grpmax_faf95_val > max_credible:
        return 'BS1'
    return None  # No criterion triggered; variant is consistent with rare-disease causation
```

## Gene-Level Constraint (LOEUF)

**Goal:** Retrieve gene constraint metrics with awareness of version mismatch for chrX/Y.

**Approach:** Use v4 LOEUF for autosomes; fall back to v2.1.1 for chrX/Y. Report LOEUF decile, not raw value, to avoid cross-version comparison errors.

```python
def query_gene_constraint(gene_symbol, dataset='gnomad_r4'):
    '''Pull gene constraint metrics. Note: v4 has no chrX/Y constraint; use v2 fallback.'''
    query = '''
    query GeneById($symbol: String!) {
      gene(gene_symbol: $symbol, reference_genome: GRCh38) {
        gene_id
        symbol
        chrom
        gnomad_constraint {
          oe_lof
          oe_lof_lower
          oe_lof_upper
          oe_mis
          oe_mis_upper
          pli
          mis_z
        }
      }
    }
    '''
    r = requests.post(GNOMAD_API,
                      json={'query': query, 'variables': {'symbol': gene_symbol}},
                      timeout=30)
    r.raise_for_status()
    gene = r.json().get('data', {}).get('gene')
    if gene is None:
        return None
    if gene.get('chrom') in ('X', 'Y'):
        gene['constraint_note'] = ('v4 constraint NOT released for chrX/Y; query v2.1.1 '
                                   'via gnomad_r2_1 dataset on the v2 endpoint')
    return gene
```

## Bulk Query via Hail (cohort-scale)

**Goal:** Filter millions of variants by AF, grpmax, or LOEUF without API rate limits.

**Approach:** Read gnomAD v4 Hail Table from Google Cloud Storage; use `hl.read_table()` + filter operations.

```python
import hail as hl

def init_hail_for_gnomad():
    '''Initialize Hail for gnomAD v4 GCS access. Requires Hail 0.2.130+.'''
    hl.init(default_reference='GRCh38')


def filter_rare_variants_hail(input_vcf, max_grpmax_faf95=0.0001, output_path='filtered.mt'):
    '''Filter input MT to variants below grpmax FAF95 threshold using gnomAD v4 exomes.'''
    ht_v4 = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/exomes/'
                          'gnomad.exomes.v4.1.sites.ht')
    mt = hl.import_vcf(input_vcf, reference_genome='GRCh38')
    mt = mt.annotate_rows(gnomad=ht_v4[mt.locus, mt.alleles])
    mt = mt.filter_rows(
        (hl.is_missing(mt.gnomad.grpmax_faf95)) |
        (mt.gnomad.grpmax_faf95.faf95 < max_grpmax_faf95)
    )
    mt.write(output_path, overwrite=True)
    return mt
```

## Per-Operation Failure Modes

**1. Using popmax/AF where grpmax_faf95 belongs**
- Trigger: Apply BS1 with raw AF instead of FAF95.
- Mechanism: Raw AF inflates for low-N populations; FAF95 is the lower-bound CI; conservative.
- Symptom: Pathogenic variants falsely categorized BS1 in small-N ancestry groups (especially MID with v4's smallest sample size).
- Fix: Use `grpmax_faf95.popmax` field; not `populations[i].af`.

**2. Failing to exclude bottleneck groups**
- Trigger: Compute grpmax including AMI, ASJ, FIN, REMAINING.
- Mechanism: Founder variants in bottleneck groups can reach AF > 5% but are not population-general; would falsely trigger BA1.
- Symptom: Founder-population pathogenic variants reported benign.
- Fix: Use gnomAD's pre-computed `grpmax_faf95` which excludes bottleneck groups by design.

**3. Querying v4 constraint for chrX/Y**
- Trigger: Pull LOEUF for DMD or USP9Y from v4 release.
- Mechanism: v4 March 2024 constraint release excluded sex chromosomes.
- Symptom: Missing or stale constraint metrics for X/Y genes.
- Fix: Query v2.1.1 LOEUF for chrX/Y; use v4 for autosomes; report LOEUF decile rather than raw value.

**4. Comparing LOEUF absolute values across v2/v4**
- Trigger: "v4 LOEUF for GENE-X is 0.45; v2 was 0.30; has it become more tolerant?"
- Mechanism: Larger v4 sample shifts the LOEUF distribution upward; first-decile threshold shifted v2 < 0.35 -> v4 < 0.6.
- Symptom: Genes appear to lose constraint between versions when they have not.
- Fix: Compare deciles, not absolute values; or stay within one version.

**5. v2 -> v4 liftover assumed equivalent**
- Trigger: Project v2 GRCh37 variants onto GRCh38 with CrossMap, treat as v4 native.
- Mechanism: ~0.5-1% of sites have different representations after liftover due to assembly fixes (e.g., gaps closed, contigs joined).
- Symptom: Inconsistent AFs at low rate; failed cross-version reproducibility.
- Fix: Query v4 native by GRCh38 coordinates directly; do not use liftover output as v4-equivalent.

**6. UKB sample contamination of grpmax**
- Trigger: Compute grpmax across v4 default subset; observe inflated NFE/SAS.
- Mechanism: 416,555 UK Biobank exomes dominate the v4 NFE+SAS subsets.
- Symptom: Variants common in UKB but rare globally falsely look common.
- Fix: Use `non_ukb` subset for grpmax when ancestry composition matters.

**7. v3 vs v4 confusion; "I want WGS"**
- Trigger: User says "I want WGS AFs" and pipeline pulls v4 genomes.
- Mechanism: v4 genomes are the SAME 76,215 v3 samples reprocessed against GRCh38; not independent.
- Symptom: WGS AFs appear identical to v3.1.2; not a bug, but worth flagging.
- Fix: Document that v4 genomes = v3 genomes reprocessed; for true independent WGS, no such resource yet exists at scale.

**8. Constraint applied to multi-isoform gene without transcript awareness**
- Trigger: Apply LOEUF "for the gene" when LoF is isoform-specific.
- Mechanism: gnomAD constraint is computed on the canonical transcript; tissue-specific or alternative isoforms may have different LoF tolerance.
- Symptom: Mis-prioritization of variants on minor transcripts.
- Fix: Cross-check with MANE Select; for isoform-specific LoF, use isoform-level constraint where available (rare).

## Reconciliation: When Sources Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| ClinVar P vs gnomAD `grpmax_faf95` > 1% | Founder-population pathogenic; or ClinVar is stale low-star | Apply Whiffin max-credible-AF for the gene; check ClinVar star/freshness |
| v2 LOEUF < 0.35 vs v4 LOEUF = 0.5 | Distribution shifted with v4 sample size, not biology | Use deciles; v4 first decile = < 0.6 |
| v2 AF != v4 AF for same variant | Sample overlap (v3 in v4) + new exomes; expected | Trust v4 default; non-overlapping subsets via `non_v2` or `non_ukb` |
| Variant present in v3 genomes, absent v4 exomes | Variant outside exome capture region (intronic, intergenic) | Use v3.1.2 or v4 genomes for non-coding |
| gnomAD-SV v2 vs v4 different breakpoints | v2 GRCh37, v4 GRCh38; assembly fixes shift coords | Use v4 native; document build |
| Browser shows lower AF than Hail Table | Browser pre-filters with `filters=PASS`; Hail Table includes all | Apply `filters` filter in Hail explicitly |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| BA1 default | grpmax_faf95 > 5% in non-bottleneck group | Richards 2015 + ClinGen SVI |
| BS1 | grpmax_faf95 > gene-specific max-credible-AF | Whiffin 2017 |
| PM2_Supporting | Absent or ultra-rare in gnomAD | SVI 2020 downgrade |
| LOEUF first decile v2 | < 0.35 | Karczewski 2020 |
| LOEUF first decile v4 | < 0.6 | gnomAD constraint release March 2024 |
| Missense Z top 1% | Z > 3.09 | Karczewski 2020 |
| mtDNA heteroplasmy carrier threshold | >=10% heteroplasmy | Laricchia 2022 |
| v4 sample size | 730,947 exomes + 76,215 genomes = 807,162 | gnomAD v4.0 release Nov 2023 |
| Bottleneck groups (excluded from grpmax) | AMI, ASJ, FIN, REMAINING | gnomAD v4 documentation |
| API rate limit | None published; ~10 req/s practical | gnomAD browser GraphQL |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| `Cannot read property 'af' of undefined` | Variant not in dataset; `variant` returned null | Check `if payload is None`; absence is biologically informative |
| FAF95 = 0 for a known common variant | `grpmax_faf95` only computed when AN sufficient | Check AC and AN directly; FAF95 is 0 when N too low to estimate |
| Variant filter status `AC0` or `RF` | Failed gnomAD QC | Variants with non-`PASS` should usually be excluded from analysis |
| Different AFs between gnomAD browser and Hail Table | Browser auto-applies PASS filter; Hail does not | Filter `filters.size() == 0` (i.e., `PASS`) in Hail |
| LOEUF appears worse in v4 vs v2 | Distribution shifted with larger sample | Compare deciles, not absolute values |
| SV not found in v4-SV | v2-SV is GRCh37, v4-SV is GRCh38; or variant not called in WGS | Try v2-SV with liftover; or check gnomAD-CNV for exome-derived |
| mtDNA variant missing | Only v3.1 has mtDNA; not in v4 | Query v3.1 directly |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why FAF95 instead of AF?" | Raw AF is point estimate; FAF95 is Poisson lower-bound 95% CI; ClinGen SVI recommendation for BS1/BA1. |
| "Why exclude FIN and ASJ from grpmax?" | Founder-population pathogenic variants reach high AF locally; including them would trigger false BA1. |
| "This LOEUF differs from the 2020 paper" | We use v4 March 2024 constraint (807k samples); 2020 paper used v2 (141k samples). Decile rank is stable; absolute shifted. |
| "Why not v4 for chrX constraint?" | v4 March 2024 constraint release is autosomes only; chrX/Y not yet released as of 2025. Fall back to v2. |
| "Why v3 if v4 exists?" | v4 genomes = v3 genomes reprocessed; for genome-only analysis they are equivalent. |
| "Variant exists in liftover v2 but not v4" | ~0.5-1% of sites differ post-assembly fixes; use v4 native, not liftover, as ground truth. |
| "Browser AF higher than this value" | Browser includes flagged variants by default; we filter on PASS. |

## References

- Chen S et al. 2024. A genomic mutational constraint map using variation in 76,156 human genomes. *Nature* 625:92.
- Karczewski KJ et al. 2020. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature* 581:434.
- Collins RL et al. 2020. A structural variation reference for medical and population genetics. *Nature* 581:444.
- Laricchia KM et al. 2022. Mitochondrial DNA variation across 56,434 individuals in gnomAD. *Genome Res* 32:569.
- Whiffin N et al. 2017. Using high-resolution variant frequencies to empower clinical genome interpretation. *Genet Med* 19:1151.
- ClinGen guidance on gnomAD v4 (March 2024): `https://clinicalgenome.org/site/assets/files/9445/clingen_guidance_to_vceps_regarding_the_use_of_gnomad_v4_march_2024.pdf`
- gnomAD v4 release notes: `https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/`
- gnomAD v4.1 updates: `https://gnomad.broadinstitute.org/news/2024-05-gnomad-v4-1-updates/`

## Related Skills

- clinical-databases/clinvar-lookup - Pathogenicity classification (gnomAD AF used for BS1/BA1)
- clinical-databases/acmg-classification - Whiffin FAF95 framework applied to ACMG criteria
- clinical-databases/variant-prioritization - Rare-disease pipeline using grpmax_faf95
- clinical-databases/myvariant-queries - Aggregated queries including gnomAD overlay
- population-genetics/population-structure - Population stratification background
