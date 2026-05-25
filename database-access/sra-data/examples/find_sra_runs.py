'''Resolve hierarchical accessions (PRJNA, GSE, SRX, SRP) to SRR runs via pysradb (preferred) or Entrez runinfo.'''
# Reference: pysradb 2.2+, biopython 1.83+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def runs_via_entrez(term, max_results=10000):
    '''Fallback path: Entrez runinfo CSV; less ergonomic than pysradb but no extra dep.'''
    h = Entrez.esearch(db='sra', term=term, retmax=0, usehistory='y')
    s = Entrez.read(h); h.close()
    total = int(s['Count'])
    if total == 0:
        return []
    webenv, qk = s['WebEnv'], s['QueryKey']
    runs = []
    for start in range(0, min(total, max_results), 500):
        h = Entrez.efetch(db='sra', rettype='runinfo', retmode='text',
                          retstart=start, retmax=500,
                          webenv=webenv, query_key=qk)
        text = h.read(); h.close()
        for line in text.strip().split('\n')[1:]:
            if line:
                runs.append(line.split(',')[0])
        time.sleep(DELAY)
    return runs


def runs_via_pysradb(identifier):
    '''Preferred path: pysradb handles GSE/PRJNA/SRX/SRP -> SRR cleanly.'''
    try:
        from pysradb import SRAweb
    except ImportError:
        print('pysradb not installed; pip install pysradb')
        return []
    db = SRAweb()
    meta = db.sra_metadata(identifier, detailed=True)
    if meta.empty:
        return []
    return meta['run_accession'].tolist()


print('=== Entrez runinfo path: BioProject PRJNA398962 ===')
runs = runs_via_entrez('PRJNA398962[BioProject]', max_results=20)
print(f'  {len(runs)} SRR runs (first 5: {runs[:5]})')

print('\n=== Entrez runinfo path: human RNA-Seq sample ===')
runs = runs_via_entrez('Homo sapiens[ORGN] AND RNA-Seq[Strategy] AND transcriptomic[Source]', max_results=10)
print(f'  {len(runs)} runs')

print('\n=== pysradb path: GSE -> SRR ===')
gse_runs = runs_via_pysradb('GSE110009')
print(f'  GSE110009 -> {len(gse_runs)} runs (first 5: {gse_runs[:5]})')

print('\n=== Write accessions file for batch download ===')
with open('accessions.txt', 'w') as f:
    for r in runs:
        f.write(f'{r}\n')
print(f'  Wrote {len(runs)} accessions to accessions.txt')
print('  Now: bash download_batch.sh accessions.txt ./fastq')
