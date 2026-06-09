"""Batch ARG presence calling across assembled MAGs with AMRFinderPlus.

Reports ARG presence per MAG and a sample-by-drug-class matrix. These are presence calls against a
curated database, NOT resistance phenotypes - a metagenomic ARG has no host link and no expression.
Run --organism only on a single resolved species, never on mixed community contigs.
"""
# Reference: AMRFinderPlus 3.12+, pandas 2.2+ | Verify API if version differs
import subprocess
import sys
from pathlib import Path
import pandas as pd


def run_amrfinder(fasta_path, output_path):
    cmd = ['amrfinder', '-n', str(fasta_path), '-o', str(output_path), '--plus', '--threads', '4']
    subprocess.run(cmd, check=True, capture_output=True)


def batch_amr_screen(input_dir, output_dir):
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    fasta_files = list(input_dir.glob('*.fasta')) + list(input_dir.glob('*.fa'))

    all_results = []
    for fasta in fasta_files:
        sample = fasta.stem
        output_file = output_dir / f'{sample}_amr.tsv'
        print(f'Processing {sample}...')
        run_amrfinder(fasta, output_file)
        df = pd.read_csv(output_file, sep='\t')
        df['sample'] = sample
        all_results.append(df)

    if not all_results:
        return
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(output_dir / 'combined_amr.tsv', sep='\t', index=False)
    print(f'\nMAGs processed: {len(all_results)}, ARG calls: {len(combined)}')
    print('Top drug classes (presence, not phenotype):')
    print(combined['Class'].value_counts().head(10))
    pivot = combined.pivot_table(index='sample', columns='Class', aggfunc='size', fill_value=0)
    pivot.to_csv(output_dir / 'amr_matrix.tsv', sep='\t')
    print(f'Presence matrix saved to {output_dir}/amr_matrix.tsv')


if __name__ == '__main__':
    in_dir = sys.argv[1] if len(sys.argv) > 1 else 'mags'
    out_dir = sys.argv[2] if len(sys.argv) > 2 else 'amr_results'
    batch_amr_screen(in_dir, out_dir)
