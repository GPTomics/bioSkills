'''UniProt REST workflows: single entry parsing, bulk TSV search, /stream for >500, async ID mapping.'''
# Reference: requests 2.31+, pandas 2.2+, UniProt REST 2024_06 | Verify API if version differs
import requests
import pandas as pd
import time
from io import StringIO


def fetch_entry_json(accession):
    r = requests.get(f'https://rest.uniprot.org/uniprotkb/{accession}.json')
    r.raise_for_status()
    e = r.json()
    pdb = [x['id'] for x in e.get('uniProtKBCrossReferences', []) if x['database'] == 'PDB']
    alphafold = next((x['id'] for x in e.get('uniProtKBCrossReferences', []) if x['database'] == 'AlphaFoldDB'), None)
    return {
        'accession': e['primaryAccession'],
        'entry_name': e.get('uniProtkbId'),
        'reviewed': e.get('entryType', '').startswith('UniProtKB reviewed'),
        'protein_name': e.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value'),
        'gene_primary': (e.get('genes') or [{}])[0].get('geneName', {}).get('value'),
        'length': e['sequence']['length'],
        'pdb_count': len(pdb),
        'pdb_ids': pdb[:5],
        'alphafold': alphafold,
    }


def search_tsv(query, fields, size=500):
    url = 'https://rest.uniprot.org/uniprotkb/search'
    params = {'query': query, 'format': 'tsv', 'fields': ','.join(fields), 'size': size}
    r = requests.get(url, params=params); r.raise_for_status()
    return pd.read_csv(StringIO(r.text), sep='\t')


def stream_tsv(query, fields):
    '''No 500-cap; right endpoint for bulk pulls.'''
    url = 'https://rest.uniprot.org/uniprotkb/stream'
    params = {'query': query, 'format': 'tsv', 'fields': ','.join(fields)}
    r = requests.get(url, params=params, stream=True); r.raise_for_status()
    return pd.read_csv(StringIO(r.text), sep='\t')


def map_ids(ids, from_db='Ensembl', to_db='UniProtKB', timeout=600, poll=3):
    submit = requests.post('https://rest.uniprot.org/idmapping/run',
                           data={'ids': ','.join(ids), 'from': from_db, 'to': to_db})
    submit.raise_for_status()
    job_id = submit.json()['jobId']
    elapsed = 0
    while elapsed < timeout:
        s = requests.get(f'https://rest.uniprot.org/idmapping/status/{job_id}')
        s.raise_for_status()
        js = s.json()
        if js.get('jobStatus') == 'RUNNING':
            time.sleep(poll); elapsed += poll
            continue
        break
    else:
        raise TimeoutError(f'job {job_id} did not complete in {timeout}s')
    r = requests.get(f'https://rest.uniprot.org/idmapping/results/{job_id}')
    r.raise_for_status()
    return r.json()


print('=== Single entry parse (P04637 / TP53) ===')
for k, v in fetch_entry_json('P04637').items():
    print(f'  {k:<14} {v}')

print('\n=== Bulk search with fields= (human reviewed kinases) ===')
df = search_tsv(
    'organism_id:9606 AND reviewed:true AND keyword:"Kinase"',
    fields=['accession', 'gene_primary', 'protein_name', 'length', 'xref_pdb'],
    size=500,
)
print(f'  {len(df)} reviewed human kinases')
print(df.head(5).to_string(index=False))

print('\n=== Stream for >500 results (all human Swiss-Prot reviewed) ===')
df_all = stream_tsv(
    'organism_id:9606 AND reviewed:true',
    fields=['accession', 'gene_primary', 'length'],
)
print(f'  All human Swiss-Prot: {len(df_all)}')

print('\n=== ID Mapping: Ensembl Gene -> UniProt ===')
mapping = map_ids(['ENSG00000141510', 'ENSG00000171862', 'ENSG00000139618'])
for r in mapping.get('results', []):
    print(f'  {r["from"]:<22} -> {r["to"]}')
for failed in mapping.get('failedIds', []):
    print(f'  {failed:<22} -> NOT MAPPED')
