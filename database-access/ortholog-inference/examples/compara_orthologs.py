'''Ensembl Compara REST: single gene + batch, with rate-limit handling and confidence semantics.'''
# Reference: requests 2.31+, Ensembl REST release 112+ | Verify API if version differs
import requests
import time
import pandas as pd

BASE = 'https://rest.ensembl.org'
HEADERS = {'Accept': 'application/json'}
SLEEP = 0.07  # 15 req/sec ceiling


def get_with_retry(url, params=None, max_retries=3):
    for attempt in range(max_retries):
        r = requests.get(url, params=params, headers=HEADERS)
        if r.status_code == 429:
            wait = int(r.headers.get('Retry-After', '5'))
            print(f'  429 rate limit; sleeping {wait}s')
            time.sleep(wait)
            continue
        r.raise_for_status()
        return r
    raise RuntimeError(f'{max_retries} retries exhausted')


def resolve_symbol(species, symbol):
    '''Always resolve symbol to Ensembl Gene ID first -- symbols are unstable (MARCH1->MARCHF1).'''
    r = get_with_retry(f'{BASE}/lookup/symbol/{species}/{symbol}')
    return r.json()['id']


def compara_orthologs(species, symbol, target_species=None):
    url = f'{BASE}/homology/symbol/{species}/{symbol}'
    params = {'type': 'orthologues'}
    if target_species:
        params['target_species'] = target_species
    r = get_with_retry(url, params=params)
    data = r.json()['data']
    if not data:
        return []
    return [{
        'source_id': h['source']['id'],
        'source_pid': h['source'].get('perc_id'),
        'target_species': h['target']['species'],
        'target_id': h['target']['id'],
        'target_pid': h['target'].get('perc_id'),
        'type': h['type'],  # ortholog_one2one / one2many / many2many
        'confidence': h.get('confidence'),  # 0 or 1; may be absent
    } for h in data[0]['homologies']]


def batch_compara(symbols, source='human', target='mouse'):
    rows = []
    for sym in symbols:
        try:
            for o in compara_orthologs(source, sym, target):
                rows.append({'symbol': sym, **o})
        except requests.HTTPError as e:
            rows.append({'symbol': sym, 'error': str(e)})
        time.sleep(SLEEP)
    return pd.DataFrame(rows)


print('=== Single-gene: BRCA1 human -> mouse ===')
for o in compara_orthologs('human', 'BRCA1', target_species='mouse'):
    print(f'  {o["target_id"]}  type={o["type"]}  confidence={o["confidence"]}  '
          f'identity={o["target_pid"]}%')

print('\n=== Symbol resolution example (MARCH1 was renamed in 2020) ===')
try:
    old = resolve_symbol('human', 'MARCH1')
    print(f'  MARCH1 -> {old}')
except requests.HTTPError:
    print('  MARCH1 not found (renamed to MARCHF1 by HGNC 2020)')
new = resolve_symbol('human', 'MARCHF1')
print(f'  MARCHF1 -> {new}')

print('\n=== Batch: 5 genes, human -> zebrafish ===')
df = batch_compara(['TP53', 'BRCA1', 'MYC', 'ATM', 'MDM2'], source='human', target='zebrafish')
one2one = df[df['type'] == 'ortholog_one2one']
print(one2one[['symbol', 'target_id', 'confidence', 'target_pid']].to_string(index=False))
print(f'\n1:1 calls: {len(one2one)}; 1:many or many:many: {len(df) - len(one2one) - df["error"].notna().sum() if "error" in df else 0}')
