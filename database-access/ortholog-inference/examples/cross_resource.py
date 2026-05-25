'''Cross-resource ortholog consensus: query Ensembl Compara, OMA REST, and OrthoDB; surface disagreement.'''
# Reference: requests 2.31+ | Verify API if version differs
import requests
import time

ENSEMBL = 'https://rest.ensembl.org'
OMA = 'https://omabrowser.org/api'
ORTHODB = 'https://data.orthodb.org/v12'
SLEEP = 0.1


def ensembl_orthologs(symbol, source='human', target='mouse'):
    r = requests.get(f'{ENSEMBL}/homology/symbol/{source}/{symbol}',
                     params={'type': 'orthologues', 'target_species': target},
                     headers={'Accept': 'application/json'})
    r.raise_for_status()
    data = r.json()['data']
    if not data:
        return []
    return [h['target']['id'] for h in data[0]['homologies'] if h['type'].startswith('ortholog')]


def oma_orthologs(uniprot_acc):
    r = requests.get(f'{OMA}/protein/{uniprot_acc}/orthologs/')
    if r.status_code == 404:
        return []
    r.raise_for_status()
    return [o.get('canonicalid') or o.get('omaid') for o in r.json()]


def orthodb_groups(symbol, species_taxid=9606):
    r = requests.get(f'{ORTHODB}/search', params={'query': symbol, 'species': species_taxid})
    r.raise_for_status()
    return r.json().get('data', [])


def orthodb_orthologs(og_id, species_taxid=10090):
    '''Members of an OG at a target species level (10090 = mouse).'''
    r = requests.get(f'{ORTHODB}/orthologs', params={'id': og_id, 'species': species_taxid})
    if r.status_code != 200:
        return []
    return [m.get('gene_id', {}).get('id') for m in r.json().get('data', [])]


print('=== Compara: BRCA1 human -> mouse ===')
ensembl_hits = ensembl_orthologs('BRCA1', 'human', 'mouse')
print(f'  {len(ensembl_hits)} hits: {ensembl_hits}')
time.sleep(SLEEP)

print('\n=== OMA: BRCA1 (UniProt P38398) human -> mouse ===')
oma_hits = oma_orthologs('P38398')
mouse_oma = [o for o in oma_hits if o and 'MOUSE' in str(o)]
print(f'  Total OMA orthologs: {len(oma_hits)}; mouse subset: {mouse_oma[:5]}')
time.sleep(SLEEP)

print('\n=== OrthoDB: BRCA1 at human level, mouse members ===')
groups = orthodb_groups('BRCA1', species_taxid=9606)
print(f'  OrthoDB groups containing BRCA1: {groups[:3]}')
if groups:
    mouse_members = orthodb_orthologs(groups[0], species_taxid=10090)
    print(f'  Group {groups[0]} mouse members: {mouse_members[:5]}')

print('\n=== Interpretation ===')
print('Disagreement across resources is informative -- not error.')
print('Compara is tree-reconciled (best for vertebrates).')
print('OMA is strict (high precision, lower recall).')
print('OrthoDB is broad (broad coverage, more ambiguous calls).')
print('For publication-grade calls, intersect 2+ resources; inspect disagreements case by case.')
