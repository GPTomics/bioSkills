'''Ensembl Compara ortholog/paralog queries via REST; includes confidence handling.'''
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


def orthologs_by_symbol(species, symbol, target=None):
    params = {'type': 'orthologues'}
    if target:
        params['target_species'] = target
    r = get_with_retry(f'{BASE}/homology/symbol/{species}/{symbol}', params=params)
    data = r.json().get('data', [])
    if not data:
        return []
    return data[0].get('homologies', [])


def paralogs_by_symbol(species, symbol):
    r = get_with_retry(f'{BASE}/homology/symbol/{species}/{symbol}',
                       params={'type': 'paralogues'})
    data = r.json().get('data', [])
    if not data:
        return []
    return data[0].get('homologies', [])


print('=== BRCA1 human -> all Compara orthologs ===')
all_orth = orthologs_by_symbol('human', 'BRCA1')
type_counts = {}
for h in all_orth:
    type_counts[h['type']] = type_counts.get(h['type'], 0) + 1
for t, c in sorted(type_counts.items(), key=lambda x: -x[1]):
    print(f'  {t:<28} {c}')
time.sleep(SLEEP)

print('\n=== BRCA1 human -> mouse, with confidence ===')
mouse_orth = orthologs_by_symbol('human', 'BRCA1', target='mouse')
for h in mouse_orth:
    print(f'  type={h["type"]}  confidence={h.get("confidence")}  '
          f'target={h["target"]["id"]}  '
          f'identity_target={h["target"].get("perc_id")}  '
          f'identity_query={h["source"].get("perc_id")}')
time.sleep(SLEEP)

print('\n=== BRCA1 within-species paralogs ===')
paralogs = paralogs_by_symbol('human', 'BRCA1')
for p in paralogs[:5]:
    print(f'  paralog: {p["target"]["id"]}  type={p["type"]}  '
          f'taxonomic_level={p.get("taxonomy_level")}')
time.sleep(SLEEP)

print('\n=== Compara confidence semantics ===')
print('  type=ortholog_one2one: 1:1 ortholog (high-confidence single match)')
print('  type=ortholog_one2many / many2one: lineage-specific duplication')
print('  type=ortholog_many2many: ancestral duplication; multiple co-orthologs')
print('  type=within_species_paralog: in-paralog (post-speciation duplication)')
print('  confidence: 0 or 1 (binary; 1 = gene-tree-topology high confidence)')
