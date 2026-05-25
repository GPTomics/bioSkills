'''Production-grade batch download: history server + disk checkpoint + WebEnv-expiry detection + jittered backoff.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
from urllib.error import HTTPError
import json
import random
import time
from pathlib import Path

Entrez.email = 'your.email@example.com'
# Entrez.api_key = 'your_api_key'


def truncate_to_last_newline(path):
    '''If a previous run crashed mid-chunk, truncate the trailing partial record.'''
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return
    with open(p, 'rb+') as f:
        f.seek(-min(8192, p.stat().st_size), 2)
        tail = f.read()
        last_nl = tail.rfind(b'\n')
        if last_nl >= 0:
            f.seek(-len(tail) + last_nl + 1, 2)
            f.truncate()


def checkpointed_download(db, term, out_path, ckpt_path, rettype='fasta',
                           batch_size=500, max_retries=5):
    delay = 0.1 if Entrez.api_key else 0.34
    ckpt = Path(ckpt_path)
    start = json.loads(ckpt.read_text())['start'] if ckpt.exists() else 0

    def refresh_session():
        h = Entrez.esearch(db=db, term=term, usehistory='y', retmax=0)
        s = Entrez.read(h); h.close()
        return s['WebEnv'], s['QueryKey'], int(s['Count'])

    webenv, query_key, total = refresh_session()
    print(f'{total:,} records matched; resuming at {start:,}')
    if start >= total:
        print('Already complete.')
        return

    if start == 0:
        Path(out_path).write_text('')
    else:
        truncate_to_last_newline(out_path)

    with open(out_path, 'a') as out:
        while start < total:
            success = False
            for attempt in range(max_retries):
                try:
                    h = Entrez.efetch(db=db, rettype=rettype, retmode='text',
                                      retstart=start, retmax=batch_size,
                                      webenv=webenv, query_key=query_key)
                    body = h.read(); h.close()
                    if not body.strip():
                        raise RuntimeError('Empty body')
                    if '<ERROR>' in body[:500]:
                        raise RuntimeError(f'WebEnv expired or server error: {body[:200]}')
                    out.write(body); out.flush()
                    success = True
                    break
                except HTTPError as e:
                    backoff = min(120, (2 ** attempt) + random.uniform(0, 1))
                    if e.code == 429:
                        print(f'  HTTP 429; backing off {backoff:.1f}s')
                        time.sleep(backoff)
                    elif e.code in (500, 502, 503, 504):
                        print(f'  HTTP {e.code}; backing off {backoff:.1f}s')
                        time.sleep(backoff)
                    else:
                        raise
                except RuntimeError as e:
                    print(f'  {e}; refreshing WebEnv')
                    webenv, query_key, total = refresh_session()
                    time.sleep(1)

            if not success:
                raise RuntimeError(f'Failed after {max_retries} retries at start={start}')

            start += batch_size
            ckpt.write_text(json.dumps({'start': start, 'total': total}))
            time.sleep(delay)
            if (start // batch_size) % 10 == 0 or start >= total:
                print(f'  {min(start, total):,}/{total:,}')

    ckpt.unlink(missing_ok=True)
    print(f'Done: {total:,} records -> {out_path}')


if __name__ == '__main__':
    checkpointed_download(
        db='nucleotide',
        term='Mus musculus[ORGN] AND hemoglobin[Gene Name] AND srcdb_refseq[PROP] AND biomol_mrna[PROP]',
        out_path='mouse_hemoglobin.fasta',
        ckpt_path='mouse_hemoglobin.ckpt.json',
    )
