'''Find GEO datasets from a publication'''
from Bio import Entrez

Entrez.email = 'your.email@example.com'

def find_geo_from_pubmed(pmid):
    '''Find GEO datasets associated with a PubMed article'''
    # Link PubMed to GDS
    handle = Entrez.elink(dbfrom='pubmed', db='gds', id=pmid)
    links = Entrez.read(handle)
    handle.close()

    if not links[0]['LinkSetDb']:
        print(f"No GEO datasets linked to PMID {pmid}")
        return []

    gds_ids = [link['Id'] for link in links[0]['LinkSetDb'][0]['Link']]
    print(f"Found {len(gds_ids)} linked GEO records")

    # Get summaries
    handle = Entrez.esummary(db='gds', id=','.join(gds_ids))
    summaries = Entrez.read(handle)
    handle.close()

    return summaries

# Get article info
print('=== Publication Info ===')
pmid = '32228226'  # Example PMID
handle = Entrez.esummary(db='pubmed', id=pmid)
pubmed_info = Entrez.read(handle)[0]
handle.close()
print(f"PMID: {pmid}")
print(f"Title: {pubmed_info['Title'][:80]}...")

# Find linked GEO
print('\n=== Linked GEO Datasets ===')
geo_datasets = find_geo_from_pubmed(pmid)

for ds in geo_datasets:
    print(f"\n{ds['Accession']}")
    print(f"  Title: {ds['title'][:60]}...")
    print(f"  Organism: {ds['taxon']}")
    print(f"  Samples: {ds['n_samples']}")
    print(f"  Platform: {ds['GPL']}")
    print(f"  Type: {ds['gdsType']}")
