'''Chain ELink calls and demonstrate EPost + neighbor_history for batches >200 IDs.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def link_one(dbfrom, db, source_id, linkname=None):
    kwargs = {'dbfrom': dbfrom, 'db': db, 'id': source_id}
    if linkname:
        kwargs['linkname'] = linkname
    h = Entrez.elink(**kwargs)
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    return [link['Id'] for link in r[0]['LinkSetDb'][0]['Link']]


def link_batch_via_history(dbfrom, db, source_ids, linkname=None, chunk=200):
    '''EPost in chunks of 200 (NCBI hard limit), then ELink with neighbor_history.'''
    webenv = None
    query_key = None
    for i in range(0, len(source_ids), chunk):
        kwargs = {'db': dbfrom, 'id': ','.join(source_ids[i:i+chunk])}
        if webenv:
            kwargs['WebEnv'] = webenv
        h = Entrez.epost(**kwargs)
        r = Entrez.read(h); h.close()
        webenv = r['WebEnv']
        query_key = r['QueryKey']
        time.sleep(DELAY)

    kwargs = {'dbfrom': dbfrom, 'db': db, 'cmd': 'neighbor_history',
              'WebEnv': webenv, 'query_key': query_key}
    if linkname:
        kwargs['linkname'] = linkname
    h = Entrez.elink(**kwargs)
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDbHistory']:
        return None, None
    # WebEnv is top-level; QueryKey is per-LinkSetDbHistory entry.
    return r[0]['WebEnv'], r[0]['LinkSetDbHistory'][0]['QueryKey']


print('=== Chain: Gene -> Protein (RefSeq) -> Structure ===')
tp53 = '7157'
proteins = link_one('gene', 'protein', tp53, linkname='gene_protein_refseq')
print(f'TP53 -> {len(proteins)} RefSeq proteins')
time.sleep(DELAY)

if proteins:
    protein_batch = ','.join(proteins[:10])
    h = Entrez.elink(dbfrom='protein', db='structure', id=protein_batch)
    r = Entrez.read(h); h.close()
    structures = []
    for ls in r:
        if ls['LinkSetDb']:
            structures.extend([l['Id'] for l in ls['LinkSetDb'][0]['Link']])
    print(f'  -> {len(structures)} structure UIDs in PDB-MMDB')
time.sleep(DELAY)

print('\n=== Large batch via history server (simulated) ===')
# Simulate 250 gene IDs (over the 200-id-per-EPost ceiling)
fake_genes = ['672', '7157', '1956', '4609', '983'] * 50
print(f'Input: {len(fake_genes)} gene UIDs (with duplicates -- realistic batch)')
we, qk = link_batch_via_history('gene', 'protein', fake_genes, linkname='gene_protein_refseq')
print(f'Got WebEnv (truncated): {we[:30] if we else "<none>"}...  QueryKey: {qk}')
print('Downstream: Entrez.efetch(db=protein, WebEnv=we, query_key=qk, retstart=..., retmax=500) in batches')
