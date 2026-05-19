'''gnomAD v4 query patterns: GraphQL, grpmax FAF95, ACMG BS1/BA1, LOEUF constraint.

Reference: requests 2.31+, myvariant 1.0+, hail 0.2.130+ | Verify gnomAD API if v4 schema differs.
v4.1 (May 2024) is current; v4.0 had AN under-counting fixed in v4.1.
'''
import requests
import myvariant
import pandas as pd

GNOMAD_API = 'https://gnomad.broadinstitute.org/api'

BOTTLENECK_GROUPS = {'ami', 'asj', 'fin', 'remaining'}


def query_variant_v4(chrom, pos, ref, alt):
    '''Query a single variant in gnomAD v4 via GraphQL. Returns exome + genome data + grpmax FAF95.

    grpmax_faf95 is the ACMG-grade frequency: 95% lower-CI of AF, excluding bottleneck groups
    (AMI, ASJ, FIN, REMAINING) by design. Use this -- NOT raw AF -- for BS1/BA1 application.
    '''
    query = '''
    query VariantById($variantId: String!) {
      variant(variantId: $variantId, dataset: gnomad_r4) {
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
          filters
          populations { id ac an }
          faf95 { popmax popmax_population }
        }
      }
    }
    '''
    variant_id = f'{chrom}-{pos}-{ref}-{alt}'
    r = requests.post(GNOMAD_API, json={'query': query, 'variables': {'variantId': variant_id}}, timeout=30)
    r.raise_for_status()
    return r.json().get('data', {}).get('variant')


def grpmax_faf95(payload):
    '''Extract grpmax FAF95 -- the ACMG-grade frequency. Prefer exome over genome.

    AMG SVI default: BA1 if grpmax_faf95 > 0.05 in non-bottleneck group.
    Bottleneck groups (AMI, ASJ, FIN, REMAINING) are excluded by design.
    '''
    if payload is None:
        return {'faf95': 0.0, 'grpmax_ancestry': None, 'source': 'absent'}
    exome = payload.get('exome')
    if exome and exome.get('faf95') and exome['faf95'].get('popmax') is not None:
        return {
            'faf95': exome['faf95']['popmax'],
            'grpmax_ancestry': exome['faf95'].get('popmax_population'),
            'source': 'exome',
            'pass': 'PASS' in (exome.get('filters') or []) or not exome.get('filters')
        }
    genome = payload.get('genome')
    if genome and genome.get('faf95') and genome['faf95'].get('popmax') is not None:
        return {
            'faf95': genome['faf95']['popmax'],
            'grpmax_ancestry': genome['faf95'].get('popmax_population'),
            'source': 'genome',
            'pass': 'PASS' in (genome.get('filters') or []) or not genome.get('filters')
        }
    return {'faf95': 0.0, 'grpmax_ancestry': None, 'source': 'absent'}


def max_credible_af(prevalence, max_allelic_contribution=1.0, max_genetic_contribution=1.0,
                    penetrance=1.0):
    '''Whiffin 2017 max-credible-AF formula for ACMG BS1 application.

    Args:
        prevalence: disease prevalence (e.g., HCM = 1/500 = 0.002)
        max_allelic_contribution: max share of cases attributable to a single allele
        max_genetic_contribution: max share of cases attributable to this gene
        penetrance: probability variant carriers develop disease

    Returns: gene-specific max-credible per-allele frequency under dominant inheritance.
             For autosomal recessive, multiply by carrier-frequency-squared appropriately.
    '''
    return (prevalence * max_genetic_contribution * max_allelic_contribution) / (penetrance * 2)


def apply_acmg_freq_codes(faf95_val, max_credible, ba1_threshold=0.05):
    '''Apply ClinGen SVI BS1/BA1/PM2_Supporting criteria from grpmax_faf95.

    BA1 (stand-alone benign): default 5% in non-bottleneck group. VCEPs override
    (Hearing Loss = 0.5% AR). Check the relevant VCEP CSpec.
    BS1: gene-specific max-credible-AF via Whiffin formula.
    PM2_Supporting: absent or ultra-rare (downgraded from PM2_Moderate in SVI 2020).
    '''
    if faf95_val is None or faf95_val == 0.0:
        return 'PM2_Supporting'
    if faf95_val > ba1_threshold:
        return 'BA1'
    if faf95_val > max_credible:
        return 'BS1'
    return None


def query_gene_constraint_v4(gene_symbol):
    '''Pull gene constraint from v4. Note: v4 X/Y constraint NOT released -- fall back to v2.'''
    query = '''
    query GeneConstraint($symbol: String!) {
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
    r = requests.post(GNOMAD_API, json={'query': query, 'variables': {'symbol': gene_symbol}}, timeout=30)
    r.raise_for_status()
    gene = r.json().get('data', {}).get('gene')
    if gene is None:
        return None
    if gene.get('chrom') in ('X', 'Y'):
        gene['note'] = 'v4 constraint NOT released for chrX/Y -- use v2.1.1 LOEUF instead'
    loeuf = gene.get('gnomad_constraint', {}).get('oe_lof_upper')
    gene['loeuf_v4_first_decile'] = loeuf is not None and loeuf < 0.6
    return gene


def annotate_variant_list(variants):
    '''Batch-annotate variants with gnomAD v4 grpmax FAF95 via myvariant.info aggregator.'''
    mv = myvariant.MyVariantInfo()
    fields = ['gnomad_exome.faf95', 'gnomad_exome.af.af', 'gnomad_genome.faf95',
              'gnomad_exome.an.an']
    results = mv.getvariants(variants, fields=fields)
    rows = []
    for r in results:
        exome_faf95 = r.get('gnomad_exome', {}).get('faf95', {})
        rows.append({
            'variant': r.get('query'),
            'grpmax_faf95': exome_faf95.get('popmax') if isinstance(exome_faf95, dict) else None,
            'grpmax_ancestry': exome_faf95.get('popmax_population') if isinstance(exome_faf95, dict) else None,
            'exome_af': r.get('gnomad_exome', {}).get('af', {}).get('af'),
            'exome_an': r.get('gnomad_exome', {}).get('an', {}).get('an')
        })
    return pd.DataFrame(rows)


def hail_bulk_filter_snippet():
    '''Print Hail code for bulk rare-variant filtering against gnomAD v4 exomes.

    Reference: hail 0.2.130+. Requires GCS authentication.
    Use when filtering >10k variants -- API rate limits make GraphQL impractical at scale.
    '''
    return '''
import hail as hl
hl.init(default_reference='GRCh38')

ht_v4 = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/exomes/gnomad.exomes.v4.1.sites.ht')
mt = hl.import_vcf('input.vcf.gz', force_bgz=True, reference_genome='GRCh38')
mt = mt.annotate_rows(gnomad=ht_v4[mt.locus, mt.alleles])

# Filter to variants below grpmax FAF95 threshold; absent variants pass
mt = mt.filter_rows(
    (hl.is_missing(mt.gnomad)) |
    (hl.is_missing(mt.gnomad.grpmax_faf95)) |
    (mt.gnomad.grpmax_faf95.faf95 < 0.0001)
)
mt.write('rare_variants.mt', overwrite=True)
'''


if __name__ == '__main__':
    brca1_variant = query_variant_v4('17', 43094464, 'G', 'A')
    grpmax_info = grpmax_faf95(brca1_variant)
    print(f'BRCA1 c.5096G>A grpmax FAF95: {grpmax_info}')

    hcm_threshold = max_credible_af(prevalence=1/500, max_genetic_contribution=0.30,
                                     max_allelic_contribution=0.10, penetrance=0.80)
    print(f'HCM gene-specific max-credible-AF (MYH7-like): {hcm_threshold:.6f}')
    code = apply_acmg_freq_codes(grpmax_info['faf95'], hcm_threshold)
    print(f'ACMG frequency code: {code}')

    scn2a_constraint = query_gene_constraint_v4('SCN2A')
    if scn2a_constraint:
        print(f'SCN2A LOEUF (v4): {scn2a_constraint.get("gnomad_constraint", {}).get("oe_lof_upper")}; '
              f'first decile: {scn2a_constraint["loeuf_v4_first_decile"]}')
