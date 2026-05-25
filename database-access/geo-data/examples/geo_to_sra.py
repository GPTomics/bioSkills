'''Resolve GSE -> SRR via pysradb (preferred) or gds -> bioproject -> sra ELink chain (fallback).'''
# Reference: biopython 1.83+, pysradb 2.2+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def gse_to_srr_pysradb(gse):
    try:
        from pysradb import SRAweb
    except ImportError:
        return None
    db = SRAweb()
    srp_df = db.gse_to_srp(gse)
    if srp_df.empty:
        return []
    srp = srp_df['study_accession'].iloc[0]
    srr_df = db.srp_to_srr(srp)
    return srr_df['run_accession'].tolist()


def gse_to_srr_entrez(gse):
    '''Fallback: gds -> bioproject -> sra ELink chain (direct gds -> sra is not supported).'''
    h = Entrez.esearch(db='gds', term=f'{gse}[Accession]')
    s = Entrez.read(h); h.close()
    if not s['IdList']:
        return []
    gds_uid = s['IdList'][0]
    time.sleep(DELAY)

    # gds -> bioproject (acheck on gds shows: bioproject, gds, pmc, pubmed, taxonomy; NOT sra)
    h = Entrez.elink(dbfrom='gds', db='bioproject', id=gds_uid)
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    bp_uids = [l['Id'] for l in r[0]['LinkSetDb'][0]['Link']]
    time.sleep(DELAY)

    # bioproject -> sra
    h = Entrez.elink(dbfrom='bioproject', db='sra', id=','.join(bp_uids))
    r = Entrez.read(h); h.close()
    sra_uids = []
    for ls in r:
        if ls['LinkSetDb']:
            sra_uids.extend(l['Id'] for l in ls['LinkSetDb'][0]['Link'])
    if not sra_uids:
        return []
    time.sleep(DELAY)

    h = Entrez.efetch(db='sra', id=','.join(sra_uids), rettype='runinfo', retmode='text')
    text = h.read(); h.close()
    runs = []
    for line in text.strip().split('\n')[1:]:
        if line:
            runs.append(line.split(',')[0])
    return runs


GSE = 'GSE147507'  # COVID-19 RNA-seq dataset (multi-cohort)

print(f'=== {GSE} -> SRR via pysradb (preferred) ===')
runs = gse_to_srr_pysradb(GSE)
if runs is None:
    print('  pysradb not installed; skip')
else:
    print(f'  {len(runs)} SRR runs (first 5: {runs[:5]})')

print(f'\n=== {GSE} -> SRR via Entrez chain (fallback: gds -> bioproject -> sra) ===')
runs = gse_to_srr_entrez(GSE)
print(f'  {len(runs)} SRR runs (first 5: {runs[:5]})')

print('\n=== Write run accessions for download ===')
with open(f'{GSE}_sra_runs.txt', 'w') as f:
    for r in runs:
        f.write(f'{r}\n')
print(f'  Wrote {len(runs)} accessions to {GSE}_sra_runs.txt')
print(f'  Hand off to sra-data: bash download_batch.sh {GSE}_sra_runs.txt ./fastq')
