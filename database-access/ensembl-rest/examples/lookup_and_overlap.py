'''Ensembl REST lookup, sequence, overlap; demonstrates symbol -> ID resolution and archive pinning.'''
# Reference: requests 2.31+, Ensembl REST release 110+ | Verify API if version differs
import requests
import time

BASE = 'https://rest.ensembl.org'
ARCHIVE = 'https://e110.rest.ensembl.org'   # pin for reproducibility
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


def symbol_to_id(species, symbol, base=BASE):
    r = get_with_retry(f'{base}/lookup/symbol/{species}/{symbol}')
    return r.json()


def gene_info(ensembl_id, base=BASE):
    r = get_with_retry(f'{base}/lookup/id/{ensembl_id}', params={'expand': 1})
    return r.json()


def sequence_for_id(ensembl_id, seq_type='protein', base=BASE):
    r = get_with_retry(f'{base}/sequence/id/{ensembl_id}', params={'type': seq_type})
    return r.json()


def genes_in_region(species, region, base=BASE):
    r = get_with_retry(f'{base}/overlap/region/{species}/{region}', params={'feature': 'gene'})
    return r.json()


print('=== Symbol resolution (live release) ===')
info = symbol_to_id('human', 'BRCA1')
print(f'  Ensembl Gene ID: {info["id"]}')
print(f'  Biotype:         {info["biotype"]}')
print(f'  Location:        chr{info["seq_region_name"]}:{info["start"]}-{info["end"]} (strand {info["strand"]})')
time.sleep(SLEEP)

print('\n=== Same query against archive (release 110) -- for reproducibility ===')
info_pinned = symbol_to_id('human', 'BRCA1', base=ARCHIVE)
print(f'  Ensembl Gene ID (e110): {info_pinned["id"]}')
if info_pinned["id"] != info["id"]:
    print(f'  WARNING: live and archive Gene IDs differ -- gene model has changed')
time.sleep(SLEEP)

print('\n=== Full gene info (transcripts + exons via ?expand=1) ===')
detail = gene_info(info['id'])
print(f'  Transcripts: {len(detail.get("Transcript", []))}')
for tx in detail.get('Transcript', [])[:3]:
    print(f'    {tx["id"]:<22} {tx["biotype"]:<22} {len(tx.get("Exon", []))} exons')
time.sleep(SLEEP)

print('\n=== Protein sequence ===')
prot = sequence_for_id('ENSG00000139618', seq_type='protein')
print(f'  BRCA2 protein: {len(prot["seq"])} aa  (first 60: {prot["seq"][:60]}...)')
time.sleep(SLEEP)

print('\n=== Genes in interval chr17:43000000-43200000 ===')
for g in genes_in_region('human', '17:43000000-43200000'):
    sym = g.get('external_name', '?')
    print(f'  {sym:<12} {g["id"]} {g["biotype"]:<22} {g["start"]}-{g["end"]}')
