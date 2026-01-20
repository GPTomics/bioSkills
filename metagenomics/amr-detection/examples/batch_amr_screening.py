#!/usr/bin/env python3
'''Batch AMR screening for multiple assemblies'''

import subprocess
import pandas as pd
from pathlib import Path
import sys

def run_amrfinder(fasta_path, output_path):
    '''Run AMRFinderPlus on a single assembly'''
    cmd = [
        'amrfinder',
        '-n', str(fasta_path),
        '-o', str(output_path),
        '--plus',
        '--threads', '4'
    ]
    subprocess.run(cmd, check=True, capture_output=True)

def batch_amr_screen(input_dir, output_dir):
    '''Screen all FASTA files in directory'''
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    fasta_files = list(input_dir.glob('*.fasta')) + list(input_dir.glob('*.fa'))

    all_results = []
    for fasta in fasta_files:
        sample = fasta.stem
        output_file = output_dir / f'{sample}_amr.tsv'

        print(f'Processing {sample}...')
        try:
            run_amrfinder(fasta, output_file)
            df = pd.read_csv(output_file, sep='\t')
            df['sample'] = sample
            all_results.append(df)
        except Exception as e:
            print(f'  Error: {e}')

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(output_dir / 'combined_amr.tsv', sep='\t', index=False)

        print('\n=== Summary ===')
        print(f'Samples processed: {len(all_results)}')
        print(f'Total AMR genes: {len(combined)}')
        print('\nTop drug classes:')
        print(combined['Class'].value_counts().head(10))

        pivot = combined.pivot_table(index='sample', columns='Class', aggfunc='size', fill_value=0)
        pivot.to_csv(output_dir / 'amr_matrix.tsv', sep='\t')
        print(f'\nAMR matrix saved to: {output_dir}/amr_matrix.tsv')

if __name__ == '__main__':
    input_dir = sys.argv[1] if len(sys.argv) > 1 else 'assemblies'
    output_dir = sys.argv[2] if len(sys.argv) > 2 else 'amr_results'
    batch_amr_screen(input_dir, output_dir)
