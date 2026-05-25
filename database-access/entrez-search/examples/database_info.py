'''Discover Entrez database structure: fields, links, last update timestamp, sortable terms.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'

DELAY = 0.34


def all_databases():
    h = Entrez.einfo()
    r = Entrez.read(h); h.close()
    return r['DbList']


def db_info(db):
    h = Entrez.einfo(db=db)
    r = Entrez.read(h); h.close()
    return r['DbInfo']


print('=== All Entrez databases ===')
dbs = all_databases()
print(f'Count: {len(dbs)}')
print(', '.join(dbs))
time.sleep(DELAY)

for db in ['nucleotide', 'pubmed', 'sra', 'gds']:
    info = db_info(db)
    print(f'\n=== {info["DbName"]} ===')
    print(f'Description: {info["Description"]}')
    print(f'Count:       {info["Count"]}')
    print(f'LastUpdate:  {info["LastUpdate"]}  # records submitted after this not yet searchable')
    print(f'Field count: {len(info["FieldList"])}')
    for field in info['FieldList'][:6]:
        print(f'    {field["Name"]:<12} ({field["FullName"]:<25}) {field["Description"][:60]}')
    print(f'Link count:  {len(info["LinkList"])}')
    for link in info['LinkList'][:4]:
        print(f'    {link["Name"]:<32} -> {link["DbTo"]}')
    time.sleep(DELAY)
