'''PubMed retrieval via MEDLINE parser (stable) and XML parser (rich but schema-drifts).'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez, Medline
import time
from io import StringIO

Entrez.email = 'your.email@example.com'
DELAY = 0.34

pmids = ['35412348', '34502548']


def safe_get(d, *path, default=None):
    '''Defensive nested-dict navigation for drifting XML schemas.'''
    cur = d
    for key in path:
        if isinstance(cur, list):
            if not cur:
                return default
            cur = cur[0]
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


print('=== MEDLINE format (schema-stable, preferred for long-lived parsers) ===')
h = Entrez.efetch(db='pubmed', id=','.join(pmids), rettype='medline', retmode='text')
records = list(Medline.parse(StringIO(h.read()))); h.close()
for r in records:
    print(f'  PMID {r.get("PMID", "?")}')
    print(f'    {r.get("TI", "(no title)")[:80]}')
    print(f'    {r.get("FAU", ["(no authors)"])[0]}')
    print(f'    {r.get("SO", "(no source)")}')
    mesh = r.get('MH', [])
    print(f'    MeSH: {len(mesh)} terms; first 3: {mesh[:3]}')
time.sleep(DELAY)

print('\n=== XML format (richer; check Bio version on schema drift) ===')
h = Entrez.efetch(db='pubmed', id=pmids[0], retmode='xml')
records = Entrez.read(h); h.close()
article = records['PubmedArticle'][0]
citation = article['MedlineCitation']
title = safe_get(citation, 'Article', 'ArticleTitle', default='?')
journal = safe_get(citation, 'Article', 'Journal', 'Title', default='?')
mesh_descriptors = [m['DescriptorName'] for m in citation.get('MeshHeadingList', [])]
grants = citation.get('Article', {}).get('GrantList', [])
pmc_id = next((id['#text'] for id in article.get('PubmedData', {}).get('ArticleIdList', [])
               if hasattr(id, 'attributes') and id.attributes.get('IdType') == 'pmc'), None)
print(f'  Title: {title}')
print(f'  Journal: {journal}')
print(f'  MeSH terms: {len(mesh_descriptors)} (first 3: {[str(m) for m in mesh_descriptors[:3]]})')
print(f'  Grants: {len(grants)}')
print(f'  PMC ID: {pmc_id or "not in PMC"}')
