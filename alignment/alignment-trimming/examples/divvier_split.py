'''Split ambiguous alignment columns with Divvier instead of removing them.

Divvier (Ali, Bogusz & Whelan 2019 MBE) preserves phylogenetic signal that pure
column removal would discard, by splitting ambiguous columns into multiple
columns each containing only confidently-assigned residues.
'''
# Reference: Divvier 1.01+ | Verify CLI flags if version differs

import subprocess
from Bio import AlignIO

def run_divvier(input_fasta, partial=False, mincol=None):
    cmd = ['divvier']
    if partial:
        cmd += ['-partial']
    else:
        cmd += ['-divvy']
    if mincol is not None:
        cmd += ['-mincol', str(mincol)]
    cmd += [input_fasta]
    subprocess.run(cmd, check=True)
    suffix = '.divvy.fas' if not partial else '.partial.fas'
    return f'{input_fasta}{suffix}'

if __name__ == '__main__':
    original = AlignIO.read('input.fasta', 'fasta')
    print(f'Original: {len(original)} sequences, {original.get_alignment_length()} columns')

    divvy_path = run_divvier('input.fasta')
    divvy_alignment = AlignIO.read(divvy_path, 'fasta')
    print('\nFull divvying (-divvy):')
    print(f'  {len(divvy_alignment)} sequences, {divvy_alignment.get_alignment_length()} columns')
    print('  (Note: column count may INCREASE due to splitting, unlike traditional trimming)')

    partial_path = run_divvier('input.fasta', partial=True, mincol=4)
    partial_alignment = AlignIO.read(partial_path, 'fasta')
    print('\nPartial divvying (-partial -mincol 4):')
    print(f'  {len(partial_alignment)} sequences, {partial_alignment.get_alignment_length()} columns')
