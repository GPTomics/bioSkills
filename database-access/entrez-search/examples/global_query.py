'''EGQuery for cross-database counts, with a fallback that loops ESearch for authoritative numbers.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'

DELAY = 0.34
CURATED_DBS = ['pubmed', 'pmc', 'nucleotide', 'protein', 'gene', 'sra', 'gds', 'bioproject', 'biosample', 'clinvar']


def egquery_counts(term):
    '''Cross-db counts via EGQuery. Note: index can lag per-database counts by 1-2 days.'''
    handle = Entrez.egquery(term=term)
    record = Entrez.read(handle); handle.close()
    return {r['DbName']: int(r['Count']) for r in record['eGQueryResult']}


def loop_esearch_counts(term, dbs):
    '''Authoritative per-database counts via ESearch. Slower but never lags the indexes.'''
    counts = {}
    for db in dbs:
        h = Entrez.esearch(db=db, term=term, retmax=0)
        r = Entrez.read(h); h.close()
        counts[db] = int(r['Count'])
        time.sleep(DELAY)
    return counts


print('=== EGQuery (fast, may lag) ===')
egq = egquery_counts('CRISPR')
for db in sorted(egq, key=egq.get, reverse=True):
    if egq[db] > 0:
        print(f'  {db:<15} {egq[db]:>10,}')
time.sleep(DELAY)

print('\n=== Loop ESearch (authoritative) for the same term ===')
esq = loop_esearch_counts('CRISPR', CURATED_DBS)
for db in sorted(esq, key=esq.get, reverse=True):
    if esq[db] > 0:
        delta = esq[db] - egq.get(db, 0)
        marker = f'(EGQuery off by {delta:+,})' if delta else ''
        print(f'  {db:<15} {esq[db]:>10,}  {marker}')
