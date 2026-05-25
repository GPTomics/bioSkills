'''Find GEO datasets cited in a paper via pubmed -> gds ELink, with SuperSeries flagging.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def find_geo_for_pubmed(pmid):
    h = Entrez.elink(dbfrom='pubmed', db='gds', id=pmid)
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    gds_uids = [l['Id'] for l in r[0]['LinkSetDb'][0]['Link']]
    time.sleep(DELAY)
    h = Entrez.esummary(db='gds', id=','.join(gds_uids))
    return Entrez.read(h)


def article_info(pmid):
    h = Entrez.esummary(db='pubmed', id=pmid)
    r = Entrez.read(h)[0]; h.close()
    return r


PMID = '32228226'  # Blanco-Melo et al. 2020 *Cell* (COVID-19 transcriptional response)

print('=== Article ===')
art = article_info(PMID)
print(f'  PMID:    {PMID}')
print(f'  Title:   {art.get("Title", "?")[:90]}')
print(f'  Journal: {art.get("Source", "?")}, {art.get("PubDate", "?")}')
time.sleep(DELAY)

print(f'\n=== GEO datasets cited in PMID {PMID} ===')
for ds in find_geo_for_pubmed(PMID):
    super_marker = ''
    summary = ds.get('summary', '')
    if 'SuperSeries' in str(summary):
        super_marker = '[SuperSeries -- check SubSeries before processing]'
    print(f'  {ds["Accession"]:12}  {ds["n_samples"]:>4} samples  {ds["GPL"]:<8}  {super_marker}')
    print(f'      {ds["title"][:90]}')
