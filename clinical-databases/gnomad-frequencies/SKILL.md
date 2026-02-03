---
name: bio-clinical-databases-gnomad-frequencies
description: Query gnomAD for population allele frequencies to assess variant rarity. Use when filtering variants by population frequency for rare disease analysis or determining if a variant is common in the general population.
tool_type: python
primary_tool: requests
---

# gnomAD Frequency Queries

## gnomAD REST API

### Query Single Variant

```python
import requests

def query_gnomad(chrom, pos, ref, alt, dataset='gnomad_r4'):
    '''Query gnomAD API for variant frequency

    dataset options: gnomad_r4, gnomad_r3, gnomad_r2_1
    '''
    url = 'https://gnomad.broadinstitute.org/api'

    query = '''
    query ($variantId: String!, $dataset: DatasetId!) {
        variant(variantId: $variantId, dataset: $dataset) {
            exome {
                ac
                an
                af
                homozygote_count
            }
            genome {
                ac
                an
                af
                homozygote_count
            }
        }
    }
    '''

    variant_id = f'{chrom}-{pos}-{ref}-{alt}'
    variables = {'variantId': variant_id, 'dataset': dataset}

    response = requests.post(url, json={'query': query, 'variables': variables})
    return response.json()
```

### Parse gnomAD Response

```python
def parse_gnomad_result(result):
    '''Extract allele frequencies from gnomAD response'''
    data = result.get('data', {}).get('variant', {})
    if not data:
        return None

    exome = data.get('exome', {}) or {}
    genome = data.get('genome', {}) or {}

    return {
        'exome_af': exome.get('af'),
        'exome_ac': exome.get('ac'),
        'exome_an': exome.get('an'),
        'exome_hom': exome.get('homozygote_count'),
        'genome_af': genome.get('af'),
        'genome_ac': genome.get('ac'),
        'genome_an': genome.get('an'),
        'genome_hom': genome.get('homozygote_count')
    }
```

## Query via myvariant.info

```python
import myvariant

mv = myvariant.MyVariantInfo()

def get_gnomad_via_myvariant(variant_hgvs):
    '''Get gnomAD frequencies via myvariant.info'''
    result = mv.getvariant(variant_hgvs, fields=['gnomad_exome', 'gnomad_genome'])

    exome = result.get('gnomad_exome', {})
    genome = result.get('gnomad_genome', {})

    return {
        'exome_af': exome.get('af', {}).get('af'),
        'genome_af': genome.get('af', {}).get('af')
    }
```

## Population-Specific Frequencies

```python
def get_population_frequencies(variant_hgvs):
    '''Get gnomAD frequencies by ancestry population'''
    mv = myvariant.MyVariantInfo()
    result = mv.getvariant(variant_hgvs, fields=['gnomad_exome.af'])

    af_data = result.get('gnomad_exome', {}).get('af', {})

    populations = {
        'af': af_data.get('af'),           # Global
        'af_afr': af_data.get('af_afr'),   # African
        'af_amr': af_data.get('af_amr'),   # Admixed American
        'af_asj': af_data.get('af_asj'),   # Ashkenazi Jewish
        'af_eas': af_data.get('af_eas'),   # East Asian
        'af_fin': af_data.get('af_fin'),   # Finnish
        'af_nfe': af_data.get('af_nfe'),   # Non-Finnish European
        'af_sas': af_data.get('af_sas'),   # South Asian
    }
    return populations
```

## Filtering Thresholds

Common frequency cutoffs for variant filtering:

| Threshold | Use Case |
|-----------|----------|
| < 0.01 (1%) | Rare disease, ACMG PM2 |
| < 0.001 (0.1%) | Stringent rare disease |
| < 0.0001 (0.01%) | Ultra-rare |
| Absent | Novel variant |

## Filter Variants by Frequency

```python
def is_rare(gnomad_af, threshold=0.01):
    '''Check if variant is rare based on gnomAD AF

    threshold: Default 0.01 (1%) per ACMG PM2 supporting criterion
    Use 0.001 for more stringent filtering
    '''
    if gnomad_af is None:
        return True  # Absent from gnomAD = rare
    return gnomad_af < threshold

def filter_rare_variants(variants, threshold=0.01):
    '''Filter list of variants to keep only rare ones'''
    rare = []
    for v in variants:
        exome_af = v.get('gnomad_exome_af')
        genome_af = v.get('gnomad_genome_af')
        max_af = max(filter(None, [exome_af, genome_af]), default=None)
        if is_rare(max_af, threshold):
            rare.append(v)
    return rare
```

## Batch Query with Local gnomAD

For large-scale analysis, use local gnomAD VCF/Hail Table:

```python
# Using Hail for gnomAD v4
import hail as hl

ht = hl.read_table('gs://gcp-public-data--gnomad/release/4.0/ht/exomes/gnomad.exomes.v4.0.sites.ht')

# Filter to rare variants
rare_ht = ht.filter(ht.freq[0].AF < 0.01)
```

## Related Skills

- myvariant-queries - Aggregated queries including gnomAD
- variant-prioritization - Filter by frequency thresholds
- population-genetics/population-structure - Population stratification analysis
