'''VEP via Ensembl REST: ad hoc variant annotation; for bulk (>1K) use local VEP.'''
# Reference: requests 2.31+, Ensembl REST release 110+ | Verify API if version differs
import requests
import time

BASE = 'https://rest.ensembl.org'
HEADERS = {'Accept': 'application/json'}
SLEEP = 0.07


def get_with_retry(url, params=None, max_retries=3):
    for attempt in range(max_retries):
        r = requests.get(url, params=params, headers=HEADERS)
        if r.status_code == 429:
            time.sleep(int(r.headers.get('Retry-After', '5')))
            continue
        r.raise_for_status()
        return r
    raise RuntimeError(f'{max_retries} retries exhausted')


def vep_region(species, region, allele):
    '''region format: chr:start-end:strand; allele is the alt base(s).'''
    r = get_with_retry(f'{BASE}/vep/{species}/region/{region}/{allele}')
    return r.json()


def vep_hgvs(species, hgvs):
    r = get_with_retry(f'{BASE}/vep/{species}/hgvs/{hgvs}')
    return r.json()


def vep_id(species, variant_id):
    r = get_with_retry(f'{BASE}/vep/{species}/id/{variant_id}')
    return r.json()


def summarize_consequences(vep_result):
    if not vep_result:
        return
    for tc in vep_result[0].get('transcript_consequences', []):
        gene = tc.get('gene_symbol', tc.get('gene_id', '?'))
        consequences = ','.join(tc['consequence_terms'])
        impact = tc.get('impact', '?')
        line = f'  {gene:<10} {tc.get("biotype", "?"):<20} {impact:<10} {consequences}'
        if 'sift_prediction' in tc:
            line += f'  SIFT={tc["sift_prediction"]}({tc.get("sift_score", "?")})'
        if 'polyphen_prediction' in tc:
            line += f'  PolyPhen={tc["polyphen_prediction"]}({tc.get("polyphen_score", "?")})'
        print(line)


print('=== VEP by region (BRCA1 region in GRCh38 coords) ===')
# rest.ensembl.org defaults to GRCh38; for GRCh37 swap base URL to https://grch37.rest.ensembl.org
res = vep_region('human', '17:43044295-43044295:1', 'A')
summarize_consequences(res)
time.sleep(SLEEP)

print('\n=== VEP by dbSNP ID ===')
res = vep_id('human', 'rs55794205')
summarize_consequences(res)
time.sleep(SLEEP)

print('\n=== VEP by HGVS notation ===')
# HGVS coding format: <RefSeq or Ensembl tx>:c.<position><ref>><alt>
res = vep_hgvs('human', 'ENST00000366667:c.803G>A')
summarize_consequences(res)
time.sleep(SLEEP)

print('\n=== Bulk warning ===')
print('  For >1,000 variants, VEP REST is rate-limit-infeasible.')
print('  Install VEP locally and run on a VCF: vep --input_file in.vcf --output_file out.vcf --cache')
print('  See variant-calling/variant-annotation skill.')
