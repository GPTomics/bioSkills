'''myvariant.info aggregated annotation: batch queries, Lucene search, reproducibility.

Reference: myvariant 1.0+, requests 2.31+ | Verify BioThings SDK signatures if differs.
Records _meta.src per-source versions for reproducibility (dbNSFP version is the
dominant staleness vector for in-silico predictors).
'''
import myvariant
import pandas as pd
import time

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
    'dbnsfp.spliceai.ds_max',
    'dbnsfp.bayesdel.add_af.score',
    'cosmic.cosmic_id',
    'civic.openCravatUrl',
    '_meta'
]


def annotate_variants(hgvs_list, chunk_size=1000, sleep=0.5):
    '''Batch-annotate a list of variants. Returns (DataFrame, versions_dict).

    Records _meta.src per-source versions for reproducibility.
    '''
    rows = []
    versions = None
    for i in range(0, len(hgvs_list), chunk_size):
        chunk = hgvs_list[i:i + chunk_size]
        results = mv.getvariants(chunk, fields=CLINICAL_FIELDS)
        for r in results:
            if versions is None and r.get('_meta'):
                versions = {src: meta.get('version', meta.get('build', '?'))
                            for src, meta in r['_meta'].get('src', {}).items()}
            rows.append(_parse_record(r))
        if i + chunk_size < len(hgvs_list):
            time.sleep(sleep)
    return pd.DataFrame(rows), versions


def _parse_record(r):
    '''Defensively parse a myvariant record into a flat dict.'''
    clinvar = r.get('clinvar', {}) or {}
    gnomad_e = r.get('gnomad_exome', {}) or {}
    gnomad_g = r.get('gnomad_genome', {}) or {}
    dbnsfp = r.get('dbnsfp', {}) or {}
    am = dbnsfp.get('alphamissense', {}) or {}
    revel = dbnsfp.get('revel', {}) or {}
    cadd = dbnsfp.get('cadd', {}) or {}
    spliceai = dbnsfp.get('spliceai', {}) or {}
    bayesdel = dbnsfp.get('bayesdel', {}).get('add_af', {}) or {}
    faf95 = (gnomad_e.get('faf95') or gnomad_g.get('faf95') or {})
    return {
        'variant': r.get('query'),
        'rsid': r.get('dbsnp', {}).get('rsid'),
        'clinvar_sig': clinvar.get('clinical_significance'),
        'clinvar_review': clinvar.get('review_status'),
        'clinvar_var_id': clinvar.get('variant_id'),
        'gnomad_grpmax_faf95': faf95.get('popmax') if isinstance(faf95, dict) else None,
        'grpmax_ancestry': faf95.get('popmax_population') if isinstance(faf95, dict) else None,
        'gnomad_af_exome': gnomad_e.get('af', {}).get('af'),
        'gnomad_af_genome': gnomad_g.get('af', {}).get('af'),
        'alphamissense': am.get('score') if isinstance(am, dict) else None,
        'alphamissense_pred': am.get('pred') if isinstance(am, dict) else None,
        'revel': revel.get('score') if isinstance(revel, dict) else None,
        'cadd_phred': cadd.get('phred'),
        'spliceai_ds_max': spliceai.get('ds_max'),
        'bayesdel_score': bayesdel.get('score'),
        'cosmic_id': r.get('cosmic', {}).get('cosmic_id')
    }


def find_pathogenic_in_gene(gene_symbol, max_results=500, min_star=2):
    '''Lucene search: ClinVar P/LP in a gene with star >= min_star.

    Star convention:
      review_status:"reviewed by expert panel" -> 3 stars (VCEP)
      review_status:"criteria provided, multiple submitters, no conflicts" -> 2 stars
      review_status:"criteria provided, single submitter" -> 1 star
    '''
    review_query = ''
    if min_star >= 3:
        review_query = ' AND clinvar.review_status:"reviewed by expert panel"'
    elif min_star >= 2:
        review_query = (' AND (clinvar.review_status:"reviewed by expert panel" '
                        'OR clinvar.review_status:"criteria provided, multiple submitters")')
    q = (f'clinvar.gene.symbol:{gene_symbol} AND '
         'clinvar.clinical_significance:(Pathogenic OR "Likely pathogenic")' + review_query)
    return mv.query(q, size=max_results, fields=['_id', 'clinvar.clinical_significance',
                                                   'clinvar.review_status', 'clinvar.hgvs'])


def find_high_impact_in_region(chrom, start, end, min_cadd=25, am_threshold=0.564):
    '''Find variants in genomic region with high CADD or high AlphaMissense.

    Note: AlphaMissense developer threshold 0.564 is NOT the Pejaver 2022 calibrated
    PP3 threshold. Treat as supporting evidence only until ClinGen calibrates.
    '''
    q = (f'chrom:{chrom} AND hg19.start:[{start} TO {end}] AND '
         f'(dbnsfp.cadd.phred:>{min_cadd} OR dbnsfp.alphamissense.score:>{am_threshold})')
    return mv.query(q, size=500, fields=['_id', 'dbnsfp.cadd.phred', 'dbnsfp.alphamissense',
                                          'clinvar.clinical_significance'])


def metadata_versions():
    '''Pull current per-source versions across the entire myvariant.info instance.'''
    import requests
    r = requests.get('https://myvariant.info/v1/metadata', timeout=30)
    r.raise_for_status()
    meta = r.json()
    return {name: src.get('version', src.get('build', '?'))
            for name, src in meta.get('src', {}).items()}


if __name__ == '__main__':
    print('Current myvariant.info per-source versions:')
    versions = metadata_versions()
    for src, ver in sorted(versions.items())[:10]:
        print(f'  {src}: {ver}')

    print('\nBatch annotating 4 test variants...')
    test = ['chr7:g.140453136A>T', 'rs121913529', 'rs1800566', 'rs104894155']
    df, src_versions = annotate_variants(test)
    print(df[['variant', 'rsid', 'clinvar_sig', 'gnomad_grpmax_faf95',
              'alphamissense', 'revel', 'cadd_phred']].to_string(index=False))
    print(f'\ndbNSFP version recorded: {src_versions.get("dbnsfp", "n/a")}')
