'''Enumerate every linkname available for a (dbfrom, source-id) pair via cmd=acheck.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def discover_links(dbfrom, source_id):
    '''Return (LinkName, DbTo, MenuTag) tuples for the (dbfrom, source-id) pair.'''
    h = Entrez.elink(dbfrom=dbfrom, id=source_id, cmd='acheck')
    r = Entrez.read(h); h.close()
    info = r[0]['IdCheckList']['IdLinkSet'][0]['LinkInfo']
    return [(i['Name'], i['DbTo'], i['MenuTag']) for i in info]


print('=== Gene 672 (BRCA1): all available link tables ===')
for name, target, label in discover_links('gene', '672'):
    print(f'  {name:<42} -> {target:<14} ({label})')
time.sleep(DELAY)

print('\n=== Nucleotide UID 31322957 (NM_007294.4): link tables ===')
for name, target, label in discover_links('nucleotide', '31322957'):
    print(f'  {name:<42} -> {target:<14} ({label})')
time.sleep(DELAY)

print('\n=== PubMed 35412348: link tables ===')
for name, target, label in discover_links('pubmed', '35412348'):
    print(f'  {name:<42} -> {target:<14} ({label})')
