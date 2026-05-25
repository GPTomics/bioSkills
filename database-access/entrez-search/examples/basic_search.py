'''Field-qualified Entrez searches with explicit retmax handling, history-server intersection, and query-translation inspection.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'

DELAY = 0.34  # 3 req/sec ceiling without API key


def safe_count(db, term):
    '''Return count and translated query without retrieving any UIDs.'''
    handle = Entrez.esearch(db=db, term=term, retmax=0)
    record = Entrez.read(handle); handle.close()
    return int(record['Count']), record['QueryTranslation']


print('=== Field-qualified RefSeq mRNA search ===')
term = 'Homo sapiens[ORGN] AND srcdb_refseq[PROP] AND biomol_mrna[PROP] AND BRCA1[Gene Name]'
count, translation = safe_count('nucleotide', term)
print(f'Count: {count}')
print(f'Translation: {translation[:200]}')
time.sleep(DELAY)

print('\n=== retmax cap diagnostic ===')
term = 'Homo sapiens[ORGN] AND srcdb_refseq[PROP] AND biomol_mrna[PROP]'
handle = Entrez.esearch(db='nucleotide', term=term, retmax=20)
record = Entrez.read(handle); handle.close()
total = int(record['Count'])
returned = len(record['IdList'])
if total > returned:
    print(f'WARNING: {total} matched, only {returned} returned. Increase retmax or use history server.')
else:
    print(f'All {total} returned.')
time.sleep(DELAY)

print('\n=== History-server chaining ===')
h1 = Entrez.esearch(db='pubmed', term='CRISPR[Title]', usehistory='y', retmax=0)
r1 = Entrez.read(h1); h1.close()
time.sleep(DELAY)

h2 = Entrez.esearch(db='pubmed', term='2024[PDAT]', usehistory='y', WebEnv=r1['WebEnv'], retmax=0)
r2 = Entrez.read(h2); h2.close()
time.sleep(DELAY)

intersect_term = f'#{r1["QueryKey"]} AND #{r2["QueryKey"]}'
h3 = Entrez.esearch(db='pubmed', term=intersect_term, usehistory='y', WebEnv=r1['WebEnv'], retmax=0)
r3 = Entrez.read(h3); h3.close()
print(f'CRISPR in title: {r1["Count"]}')
print(f'2024 PDAT: {r2["Count"]}')
print(f'Intersection: {r3["Count"]}')
time.sleep(DELAY)

print('\n=== Spell-check before searching ===')
h = Entrez.espell(db='pubmed', term='breast canser')
r = Entrez.read(h); h.close()
print(f'Original: breast canser')
print(f'Corrected: {r["CorrectedQuery"]}')
time.sleep(DELAY)

print('\n=== Diagnose taxonomy-walk over-expansion ===')
count_default, _ = safe_count('nucleotide', 'Mammalia[ORGN] AND insulin[Gene Name]')
time.sleep(DELAY)
count_exact, _ = safe_count('nucleotide', 'Mammalia[Organism:exp] AND insulin[Gene Name]')
print(f'Mammalia (walked):     {count_default:>10}')
print(f'Mammalia (no expand):  {count_exact:>10}')
print(f'Ratio (subtree vs node-only): {count_default / max(count_exact, 1):.1f}x')
