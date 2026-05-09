'''Search a query structure against AlphaFoldDB or a custom Foldseek database.

Returns best hits with TM-score, e-value, and structural alignment metadata.
--alignment-type 2 (3Di+AA Gotoh local) is Foldseek's default and is the right
choice for the same-fold regime; --alignment-type 1 (TMalign) refines top hits
with a full global TM-score at higher cost.
'''
# Reference: foldseek 8+ | Verify CLI flags if version differs

import subprocess
import csv
import math

def foldseek_search(query_pdb, database, output_m8, tmp_dir='tmp/', alignment_type=2, max_seqs=200):
    cmd = [
        'foldseek', 'easy-search', query_pdb, database, output_m8, tmp_dir,
        '--alignment-type', str(alignment_type),
        '--max-seqs', str(max_seqs),
        '--format-output', 'query,target,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,alnlen,pident',
    ]
    subprocess.run(cmd, check=True)

def parse_results(output_m8):
    columns = ['query', 'target', 'evalue', 'bits', 'alntmscore', 'qtmscore', 'ttmscore', 'lddt', 'alnlen', 'pident']
    hits = []
    with open(output_m8) as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) != len(columns):
                print(f'Warning: skipping row with {len(row)} fields (expected {len(columns)})')
                continue
            hit = dict(zip(columns, row))
            for col in ['evalue', 'bits', 'alntmscore', 'qtmscore', 'ttmscore', 'lddt', 'pident']:
                value = hit[col]
                hit[col] = float(value) if value not in ('', 'NA', 'nan') else math.nan
            hit['alnlen'] = int(hit['alnlen']) if hit['alnlen'].isdigit() else 0
            hits.append(hit)
    return hits

if __name__ == '__main__':
    foldseek_search('query.pdb', '/path/to/afdb', 'result.m8')
    hits = parse_results('result.m8')

    print('Top 20 Foldseek hits:')
    print(f'{"target":<25} {"evalue":>10} {"alnTM":>7} {"qTM":>7} {"lddt":>6} {"%id":>5}')
    for hit in hits[:20]:
        print(f'{hit["target"][:25]:<25} {hit["evalue"]:>10.2e} {hit["alntmscore"]:>7.3f} '
              f'{hit["qtmscore"]:>7.3f} {hit["lddt"]:>6.3f} {hit["pident"]:>5.1f}')

    confident_homologs = [h for h in hits if not math.isnan(h['alntmscore']) and h['alntmscore'] > 0.5 and h['evalue'] < 1e-3]
    print(f'\nConfident structural homologs (alnTM > 0.5, e-value < 1e-3): {len(confident_homologs)}')
