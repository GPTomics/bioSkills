'''Robust batch download with retry logic'''
from Bio import Entrez
from urllib.error import HTTPError
import time

Entrez.email = 'your.email@example.com'
# Entrez.api_key = 'your_api_key'

def robust_batch_download(db, term, output_file, rettype='fasta', batch_size=500, max_retries=3):
    # Search with history
    handle = Entrez.esearch(db=db, term=term, usehistory='y')
    search = Entrez.read(handle)
    handle.close()

    webenv = search['WebEnv']
    query_key = search['QueryKey']
    total = int(search['Count'])

    print(f"Found {total} records to download")
    if total == 0:
        return

    delay = 0.1 if Entrez.api_key else 0.34
    downloaded = 0

    with open(output_file, 'w') as out:
        for start in range(0, total, batch_size):
            end = min(start + batch_size, total)

            for retry in range(max_retries):
                try:
                    handle = Entrez.efetch(db=db, rettype=rettype, retmode='text',
                                           retstart=start, retmax=batch_size,
                                           webenv=webenv, query_key=query_key)
                    data = handle.read()
                    handle.close()

                    if data.strip():
                        out.write(data)
                        downloaded = end
                        print(f"Progress: {downloaded}/{total} ({100*downloaded/total:.1f}%)")
                    break

                except HTTPError as e:
                    if e.code == 429:
                        wait_time = 10 * (retry + 1)
                        print(f"Rate limited. Waiting {wait_time}s (attempt {retry+1}/{max_retries})")
                        time.sleep(wait_time)
                    elif e.code == 400 and retry < max_retries - 1:
                        print(f"Bad request, retrying... (attempt {retry+1}/{max_retries})")
                        time.sleep(5)
                    else:
                        print(f"Error: {e}")
                        raise

                except Exception as e:
                    if retry < max_retries - 1:
                        print(f"Error: {e}, retrying...")
                        time.sleep(5)
                    else:
                        raise

            time.sleep(delay)

    print(f"Download complete: {downloaded} records saved to {output_file}")

# Example
robust_batch_download('nucleotide', 'mouse[orgn] AND hemoglobin[gene] AND mRNA[fkey]', 'mouse_hemoglobin.fasta')
