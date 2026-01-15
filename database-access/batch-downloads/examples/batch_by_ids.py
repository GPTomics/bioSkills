'''Download records by a list of IDs'''
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'

def download_by_ids(db, ids, output_file, rettype='fasta', batch_size=200):
    total = len(ids)
    delay = 0.1 if Entrez.api_key else 0.34

    print(f"Downloading {total} records...")

    with open(output_file, 'w') as out:
        for start in range(0, total, batch_size):
            batch = ids[start:start+batch_size]
            end = start + len(batch)
            print(f"  Batch {start+1}-{end} of {total}")

            handle = Entrez.efetch(db=db, id=','.join(batch), rettype=rettype, retmode='text')
            out.write(handle.read())
            handle.close()

            time.sleep(delay)

    print(f"Saved to {output_file}")

# Example: list of accessions
accessions = [
    'NM_007294', 'NM_000059', 'NM_000546',
    'NM_001126112', 'NM_004985', 'NM_000492',
    'NM_000314', 'NM_001197', 'NM_002467', 'NM_000518'
]
download_by_ids('nucleotide', accessions, 'selected_genes.fasta')

# For large ID lists, use epost first
def download_large_id_list(db, ids, output_file, rettype='fasta', batch_size=500):
    print(f"Posting {len(ids)} IDs to history server...")
    handle = Entrez.epost(db=db, id=','.join(ids))
    result = Entrez.read(handle)
    handle.close()

    webenv = result['WebEnv']
    query_key = result['QueryKey']
    total = len(ids)
    delay = 0.1 if Entrez.api_key else 0.34

    with open(output_file, 'w') as out:
        for start in range(0, total, batch_size):
            end = min(start + batch_size, total)
            print(f"  Fetching {start+1}-{end} of {total}")

            handle = Entrez.efetch(db=db, rettype=rettype, retmode='text',
                                   retstart=start, retmax=batch_size,
                                   webenv=webenv, query_key=query_key)
            out.write(handle.read())
            handle.close()

            time.sleep(delay)

    print(f"Saved to {output_file}")
