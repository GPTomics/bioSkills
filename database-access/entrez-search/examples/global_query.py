'''Search across all NCBI databases with Entrez.egquery()'''
from Bio import Entrez

Entrez.email = 'your.email@example.com'

# Global query - search all databases
print('=== Global Query: CRISPR ===')
handle = Entrez.egquery(term='CRISPR')
record = Entrez.read(handle)
handle.close()

print('Databases with matches:')
for result in record['eGQueryResult']:
    count = int(result['Count'])
    if count > 0:
        print(f"  {result['DbName']:20} {count:>10,} records")

# Compare terms across databases
print('\n=== Comparing Terms ===')
terms = ['COVID-19', 'SARS-CoV-2', 'coronavirus']
for term in terms:
    handle = Entrez.egquery(term=term)
    record = Entrez.read(handle)
    handle.close()

    for result in record['eGQueryResult']:
        if result['DbName'] == 'pubmed':
            print(f"PubMed '{term}': {result['Count']} articles")
            break
