'''Explore NCBI database structure with Entrez.einfo()'''
from Bio import Entrez

Entrez.email = 'your.email@example.com'

# List all databases
print('=== Available NCBI Databases ===')
handle = Entrez.einfo()
record = Entrez.read(handle)
handle.close()
print(f"Total databases: {len(record['DbList'])}")
print(f"Databases: {', '.join(record['DbList'][:10])}...")

# Get info about specific database
print('\n=== Nucleotide Database Info ===')
handle = Entrez.einfo(db='nucleotide')
record = Entrez.read(handle)
handle.close()
info = record['DbInfo']
print(f"Name: {info['DbName']}")
print(f"Description: {info['Description']}")
print(f"Total records: {info['Count']}")
print(f"Last updated: {info['LastUpdate']}")

# List searchable fields
print('\n=== Searchable Fields ===')
for field in info['FieldList'][:10]:
    print(f"  {field['Name']:15} - {field['Description']}")

# List available links to other databases
print('\n=== Links to Other Databases ===')
for link in info['LinkList'][:5]:
    print(f"  {link['Name']:30} -> {link['DbTo']}")
