'''Search GEO with field-qualified queries and surface SuperSeries before download.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
import gzip
import urllib.request
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def search_geo(term, study_type='gse', organism=None, gds_type=None, max_results=20):
    parts = [term, f'{study_type}[Entry Type]']
    if organism:
        parts.append(f'{organism}[Organism]')
    if gds_type:
        parts.append(f'{gds_type}[GDS Type]')
    full_term = ' AND '.join(parts)
    h = Entrez.esearch(db='gds', term=full_term, retmax=max_results)
    s = Entrez.read(h); h.close()
    if not s['IdList']:
        return []
    h = Entrez.esummary(db='gds', id=','.join(s['IdList']))
    return Entrez.read(h)


def detect_super_series(gse):
    '''Check the SOFT family file for !Series_relation. Returns (super_of, sub_of).'''
    prefix = gse[:-3] + 'nnn'
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/soft/{gse}_family.soft.gz'
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = resp.read()
    except Exception as e:
        return {'super_of': [], 'sub_of': None, 'error': str(e)}
    super_of, sub_of = [], None
    with gzip.GzipFile(fileobj=__import__('io').BytesIO(data), mode='rt') as f:
        for line in f:
            if line.startswith('^SAMPLE'):
                break
            if line.startswith('!Series_relation'):
                if 'SuperSeries of' in line:
                    super_of.append(line.split('SuperSeries of: ')[1].strip())
                elif 'SubSeries of' in line:
                    sub_of = line.split('SubSeries of: ')[1].strip()
    return {'super_of': super_of, 'sub_of': sub_of}


print('=== Search: breast cancer RNA-seq, human ===')
for s in search_geo(
    'breast cancer',
    organism='Homo sapiens',
    gds_type='expression profiling by high throughput sequencing',
    max_results=10,
):
    print(f'  {s["Accession"]:12} {s["n_samples"]:>4} samples  {s["GPL"]:<10}  {s["title"][:60]}')
    time.sleep(DELAY)

print('\n=== SuperSeries detection (live FTP check) ===')
for gse in ['GSE122288', 'GSE123456']:
    info = detect_super_series(gse)
    if info['super_of']:
        print(f'  {gse}: SuperSeries of {info["super_of"][:3]}{"..." if len(info["super_of"]) > 3 else ""}')
    elif info['sub_of']:
        print(f'  {gse}: SubSeries of {info["sub_of"]}')
    elif 'error' in info:
        print(f'  {gse}: error checking ({info["error"]})')
    else:
        print(f'  {gse}: standalone Series (safe to process as one experiment)')
    time.sleep(DELAY)
