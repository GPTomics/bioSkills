'''Chain multiple database links together'''
from Bio import Entrez

Entrez.email = 'your.email@example.com'

def gene_to_proteins(gene_id):
    handle = Entrez.elink(dbfrom='gene', db='protein', id=gene_id, linkname='gene_protein_refseq')
    record = Entrez.read(handle)
    handle.close()
    if not record[0]['LinkSetDb']:
        return []
    return [link['Id'] for link in record[0]['LinkSetDb'][0]['Link']]

def proteins_to_structures(protein_ids):
    handle = Entrez.elink(dbfrom='protein', db='structure', id=','.join(protein_ids[:10]))
    record = Entrez.read(handle)
    handle.close()
    structures = []
    for linkset in record:
        if linkset['LinkSetDb']:
            structures.extend([link['Id'] for link in linkset['LinkSetDb'][0]['Link']])
    return structures

# Navigate Gene -> Protein -> Structure
print('=== Gene -> Protein -> Structure ===')
gene_id = '7157'  # TP53
print(f"Starting with gene: {gene_id} (TP53)")

protein_ids = gene_to_proteins(gene_id)
print(f"Found {len(protein_ids)} proteins")

if protein_ids:
    structure_ids = proteins_to_structures(protein_ids[:5])
    print(f"Found {len(structure_ids)} structures")
    if structure_ids:
        print(f"Structure IDs: {structure_ids[:10]}")

# Gene to publications
print('\n=== Gene -> Publications ===')
handle = Entrez.elink(dbfrom='gene', db='pubmed', id='672', linkname='gene_pubmed')
record = Entrez.read(handle)
handle.close()

if record[0]['LinkSetDb']:
    pubmed_ids = [link['Id'] for link in record[0]['LinkSetDb'][0]['Link'][:10]]
    print(f"BRCA1 publications (first 10): {pubmed_ids}")
