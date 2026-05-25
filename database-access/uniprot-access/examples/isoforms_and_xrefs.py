'''Iterate isoforms and cross-references for a UniProt entry; download reference proteome.'''
# Reference: requests 2.31+, UniProt REST 2024_06 | Verify API if version differs
import requests
import time

DELAY = 0.05  # 200 req/sec tolerated


def list_isoforms(accession):
    '''Read comments[type=ALTERNATIVE PRODUCTS] for isoform IDs and names.'''
    r = requests.get(f'https://rest.uniprot.org/uniprotkb/{accession}.json')
    r.raise_for_status()
    entry = r.json()
    for comment in entry.get('comments', []):
        if comment.get('commentType') == 'ALTERNATIVE PRODUCTS':
            return [{
                'name': iso.get('name', {}).get('value'),
                'ids': iso.get('isoformIds', []),
                'canonical': iso.get('isoformSequenceStatus') == 'Displayed',
            } for iso in comment.get('isoforms', [])]
    return []


def fetch_isoform_fasta(accession_with_suffix):
    r = requests.get(f'https://rest.uniprot.org/uniprotkb/{accession_with_suffix}.fasta')
    r.raise_for_status()
    return r.text


def xref_summary(accession):
    r = requests.get(f'https://rest.uniprot.org/uniprotkb/{accession}.json')
    r.raise_for_status()
    e = r.json()
    by_db = {}
    for x in e.get('uniProtKBCrossReferences', []):
        by_db.setdefault(x['database'], []).append(x['id'])
    return by_db


def download_proteome(upid, out_path):
    '''upid: e.g. UP000005640 (human reference proteome).'''
    url = f'https://rest.uniprot.org/proteomes/{upid}.fasta.gz'
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(out_path, 'wb') as f:
            for chunk in r.iter_content(8192):
                f.write(chunk)
    return out_path


print('=== TP53 isoforms (P04637) ===')
isos = list_isoforms('P04637')
print(f'Found {len(isos)} isoforms:')
for iso in isos[:5]:
    canonical = '[canonical]' if iso['canonical'] else ''
    print(f'  {iso["ids"]:<20} {iso["name"]:<30} {canonical}')

time.sleep(DELAY)

# Fetch a non-canonical isoform
if len(isos) > 1:
    second = isos[1]['ids'][0] if isos[1]['ids'] else None
    if second:
        fasta = fetch_isoform_fasta(second)
        print(f'\n=== Fetched isoform {second} FASTA (first 200 chars) ===')
        print(fasta[:200])
        time.sleep(DELAY)

print('\n=== Cross-references summary (P04637) ===')
xrefs = xref_summary('P04637')
for db in sorted(xrefs, key=lambda k: -len(xrefs[k]))[:10]:
    print(f'  {db:<20} {len(xrefs[db])} entries')

print('\n=== Download human reference proteome (commented out - large file) ===')
print('  download_proteome("UP000005640", "human_reference_proteome.fasta.gz")')
print('  # ~20 MB compressed; ~80 MB unpacked; ~20K proteins')
