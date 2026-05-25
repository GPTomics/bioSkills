'''Download by known ID list: direct EFetch for <200 IDs, chained EPost+history for larger lists, with integrity verification.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez, SeqIO
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34  # 0.10 with API key
EPOST_LIMIT = 200  # NCBI hard limit per EPost call


def direct_efetch(db, ids, out_path, rettype='fasta'):
    '''For <200 IDs: comma-joined EFetch in one call.'''
    assert len(ids) <= EPOST_LIMIT, f'Use chained EPost for >{EPOST_LIMIT} IDs'
    h = Entrez.efetch(db=db, id=','.join(ids), rettype=rettype, retmode='text')
    with open(out_path, 'w') as out:
        out.write(h.read())
    h.close()


def chained_epost_fetch(db, ids, out_path, rettype='fasta', batch_size=500):
    '''For larger ID lists: EPost in chunks of 200 sharing a WebEnv, then EFetch by WebEnv/QueryKey.'''
    delay = 0.1 if Entrez.api_key else DELAY
    webenv = None
    posts = []  # list of (query_key, chunk_size)

    for i in range(0, len(ids), EPOST_LIMIT):
        chunk = ids[i:i+EPOST_LIMIT]
        kwargs = {'db': db, 'id': ','.join(chunk)}
        if webenv:
            kwargs['WebEnv'] = webenv
        h = Entrez.epost(**kwargs)
        r = Entrez.read(h); h.close()
        webenv = r['WebEnv']
        posts.append((r['QueryKey'], len(chunk)))
        time.sleep(delay)
        print(f'  Posted {min(i + EPOST_LIMIT, len(ids))}/{len(ids)} IDs')

    with open(out_path, 'w') as out:
        for query_key, chunk_total in posts:
            for start in range(0, chunk_total, batch_size):
                h = Entrez.efetch(db=db, rettype=rettype, retmode='text',
                                  retstart=start, retmax=batch_size,
                                  webenv=webenv, query_key=query_key)
                out.write(h.read()); h.close()
                time.sleep(delay)


def verify_count(out_path, expected, fmt='fasta'):
    observed = sum(1 for _ in SeqIO.parse(out_path, fmt))
    if observed != expected:
        print(f'  WARNING: expected {expected}, got {observed}')
    else:
        print(f'  Integrity OK: {observed} records')
    return observed == expected


small_list = [
    'NM_007294.4', 'NM_000059.4', 'NM_000546.6',
    'NM_001126112.3', 'NM_004985.5', 'NM_000492.4',
]

print(f'=== Direct EFetch for {len(small_list)} known accessions ===')
direct_efetch('nucleotide', small_list, 'small_list.fasta')
verify_count('small_list.fasta', len(small_list))

# Realistic larger list: 250 IDs (would fail with comma-joined EFetch)
large_list = small_list * 42 + small_list[:3]  # 255 IDs (with duplicates -- ok for demo)
print(f'\n=== Chained EPost + history fetch for {len(large_list)} IDs ===')
chained_epost_fetch('nucleotide', large_list, 'large_list.fasta')
# Dedup expected count: 6 unique * (42 + 3/6 -> 252+3 = 255 entries, fetched as 255)
print(f'  Verification skipped (duplicates in input)')
