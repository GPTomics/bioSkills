'''ESummary for bulk metadata at 5-10x lower cost than EFetch; demonstrates chunked retrieval.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34
CHUNK = 500  # ESummary supports up to 10K; 500 stays well under URL length cap


def bulk_summaries(db, ids, chunk=CHUNK):
    '''Yield docsums in chunks; respects rate limit.'''
    for i in range(0, len(ids), chunk):
        h = Entrez.esummary(db=db, id=','.join(ids[i:i+chunk]))
        for r in Entrez.read(h):
            yield r
        h.close()
        time.sleep(DELAY)


print('=== Nucleotide docsum: organism, length, AccessionVersion ===')
ids = ['NM_007294.4', 'NM_000059.4', 'NM_000546.6', 'NM_001126112.3']
for s in bulk_summaries('nucleotide', ids):
    print(f'  {s["AccessionVersion"]:<18} {s["Length"]:>8} nt   {s["Organism"]}')

print('\n=== PubMed docsum: title, journal, date ===')
pmids = ['35412348', '34502548', '36045532']
for s in bulk_summaries('pubmed', pmids):
    title = s.get('Title', '?')
    authors = s.get('AuthorList', [])
    journal = s.get('Source', '?')
    date = s.get('PubDate', '?')
    print(f'  PMID {s["Id"]}')
    print(f'    {title[:80]}')
    print(f'    {", ".join(authors[:3])}{"..." if len(authors) > 3 else ""}')
    print(f'    {journal}, {date}')

print('\n=== Compare ESummary vs EFetch payload size ===')
import sys
h = Entrez.esummary(db='nucleotide', id=','.join(ids))
summary_bytes = sys.getsizeof(h.read()); h.close()
time.sleep(DELAY)
h = Entrez.efetch(db='nucleotide', id=','.join(ids), rettype='gb', retmode='text')
efetch_bytes = sys.getsizeof(h.read()); h.close()
print(f'ESummary payload: {summary_bytes:>10,} bytes')
print(f'EFetch (gb) payload: {efetch_bytes:>10,} bytes')
print(f'EFetch / ESummary: {efetch_bytes / summary_bytes:.1f}x')
