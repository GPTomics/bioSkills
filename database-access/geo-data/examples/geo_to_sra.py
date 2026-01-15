'''Link GEO series to SRA for raw data download'''
from Bio import Entrez

Entrez.email = 'your.email@example.com'

def get_sra_from_geo(gse_accession):
    '''Find SRA experiments linked to a GEO series'''
    # Find GEO record
    handle = Entrez.esearch(db='gds', term=f'{gse_accession}[accn]')
    search = Entrez.read(handle)
    handle.close()

    if not search['IdList']:
        print(f"No GEO record found for {gse_accession}")
        return []

    gds_id = search['IdList'][0]
    print(f"Found GEO record ID: {gds_id}")

    # Link to SRA
    handle = Entrez.elink(dbfrom='gds', db='sra', id=gds_id)
    links = Entrez.read(handle)
    handle.close()

    if not links[0]['LinkSetDb']:
        print("No SRA data linked to this GEO series")
        return []

    sra_ids = [link['Id'] for link in links[0]['LinkSetDb'][0]['Link']]
    print(f"Found {len(sra_ids)} SRA records")

    return sra_ids

def get_run_info(sra_ids):
    '''Get run information from SRA IDs'''
    if len(sra_ids) > 50:
        sra_ids = sra_ids[:50]

    handle = Entrez.efetch(db='sra', id=','.join(sra_ids), rettype='runinfo', retmode='text')
    runinfo = handle.read()
    handle.close()

    runs = []
    lines = runinfo.strip().split('\n')
    if len(lines) > 1:
        header = lines[0].split(',')
        for line in lines[1:]:
            if line:
                values = line.split(',')
                run_dict = dict(zip(header, values))
                if run_dict.get('Run'):
                    runs.append(run_dict)
    return runs

# Example: Find SRA data for a GEO series
print('=== GEO to SRA ===')
gse = 'GSE147507'  # COVID-19 RNA-seq dataset
sra_ids = get_sra_from_geo(gse)

if sra_ids:
    runs = get_run_info(sra_ids[:20])
    print(f"\nSRA Runs for {gse}:")
    for run in runs[:10]:
        print(f"  {run.get('Run', 'N/A')}: {run.get('LibraryStrategy', 'N/A')} - {run.get('spots', 'N/A')} spots")

    # Write accessions to file
    with open(f'{gse}_sra_runs.txt', 'w') as f:
        for run in runs:
            if run.get('Run'):
                f.write(run['Run'] + '\n')
    print(f"\nWrote {len(runs)} accessions to {gse}_sra_runs.txt")
