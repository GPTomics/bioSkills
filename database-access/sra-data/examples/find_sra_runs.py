'''Find SRA run accessions using Entrez'''
from Bio import Entrez

Entrez.email = 'your.email@example.com'

def search_sra(term, max_results=100):
    '''Search SRA database and return run accessions'''
    handle = Entrez.esearch(db='sra', term=term, retmax=max_results)
    search = Entrez.read(handle)
    handle.close()

    print(f"Found {search['Count']} total results")
    if not search['IdList']:
        return []

    handle = Entrez.efetch(db='sra', id=','.join(search['IdList']), rettype='runinfo', retmode='text')
    runinfo = handle.read()
    handle.close()

    runs = []
    lines = runinfo.strip().split('\n')
    if len(lines) > 1:
        for line in lines[1:]:
            if line:
                fields = line.split(',')
                if fields[0]:
                    runs.append(fields[0])
    return runs

# Search by BioProject
print('=== Search by BioProject ===')
runs = search_sra('PRJNA398962[bioproject]', max_results=10)
print(f"Runs: {runs}")

# Search by organism and experiment type
print('\n=== Search by Criteria ===')
runs = search_sra('human[orgn] AND RNA-Seq[strategy] AND transcriptomic[source]', max_results=10)
print(f"Runs: {runs}")

# Write accessions to file
print('\n=== Writing to File ===')
with open('accessions.txt', 'w') as f:
    for run in runs:
        f.write(run + '\n')
print(f"Wrote {len(runs)} accessions to accessions.txt")
