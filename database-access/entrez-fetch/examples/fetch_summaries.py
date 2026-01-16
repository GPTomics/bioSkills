'''Get document summaries without downloading full records'''
from Bio import Entrez

Entrez.email = 'your.email@example.com'

# Get nucleotide summary
print('=== Nucleotide Summary ===')
handle = Entrez.esummary(db='nucleotide', id='NM_007294')
records = Entrez.read(handle)
handle.close()

summary = records[0]
print(f"Caption: {summary['Caption']}")
print(f"Title: {summary['Title']}")
print(f"Length: {summary['Length']}")
print(f"Organism: {summary['Organism']}")
print(f"TaxId: {summary['TaxId']}")

# Get multiple summaries
print('\n=== Multiple Summaries ===')
ids = '224589800,224589802,224589804'
handle = Entrez.esummary(db='nucleotide', id=ids)
records = Entrez.read(handle)
handle.close()

for summary in records:
    print(f"{summary['Caption']}: {summary['Length']} bp - {summary['Organism']}")

# PubMed summary
print('\n=== PubMed Summary ===')
handle = Entrez.esummary(db='pubmed', id='35412348')
records = Entrez.read(handle)
handle.close()

summary = records[0]
print(f"Title: {summary['Title'][:80]}...")
print(f"Authors: {', '.join(summary['AuthorList'][:3])}...")
print(f"Journal: {summary['Source']}")
print(f"Date: {summary['PubDate']}")
