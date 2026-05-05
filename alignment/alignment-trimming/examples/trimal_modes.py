'''Compare trimAl modes on the same input alignment.

trimAl modes (Capella-Gutierrez et al 2009 Bioinformatics):
- automated1: heuristically picks gappyout, strict, or strictplus
- gappyout: aggressive gap removal (recommended for HMM profile prep)
- strict: combined gap + similarity (recommended for phylogenetics)
- strictplus: ~50% more aggressive than strict
- gt N: manual gap threshold (e.g. -gt 0.5)
'''
# Reference: trimAl 1.4+ | Verify CLI flags if version differs

import subprocess
from Bio import AlignIO

def run_trimal(input_fasta, output_fasta, mode_flags, kept_columns_file=None):
    cmd = ['trimal', '-in', input_fasta, '-out', output_fasta] + mode_flags
    if kept_columns_file:
        cmd += ['-colnumbering']
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    if kept_columns_file:
        with open(kept_columns_file, 'w') as f:
            f.write(result.stdout)

def get_alignment_length(fasta_path):
    return AlignIO.read(fasta_path, 'fasta').get_alignment_length()

if __name__ == '__main__':
    modes = {
        'automated1': ['-automated1'],
        'gappyout': ['-gappyout'],
        'strict': ['-strict'],
        'strictplus': ['-strictplus'],
        'gt0.5': ['-gt', '0.5'],
    }

    original_length = get_alignment_length('input.fasta')
    print(f'Original alignment: {original_length} columns\n')
    print(f'{"Mode":<12} {"Columns":>8} {"Retained":>10} {"Removed":>9}')
    print('-' * 45)
    for mode_name, flags in modes.items():
        output_path = f'trimmed_{mode_name}.fasta'
        run_trimal('input.fasta', output_path, flags, f'cols_{mode_name}.txt')
        length = get_alignment_length(output_path)
        retention = length / original_length
        print(f'{mode_name:<12} {length:>8} {retention*100:>9.1f}% {original_length - length:>9}')
