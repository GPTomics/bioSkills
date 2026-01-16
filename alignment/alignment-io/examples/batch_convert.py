'''Batch convert alignment files'''

from pathlib import Path
from Bio import AlignIO

input_dir = Path('alignments/')
output_dir = Path('converted/')
output_dir.mkdir(exist_ok=True)

input_format = 'clustal'
output_format = 'fasta'

for input_file in input_dir.glob('*.aln'):
    alignment = AlignIO.read(input_file, input_format)
    output_file = output_dir / f'{input_file.stem}.fasta'
    AlignIO.write(alignment, output_file, output_format)
    print(f'{input_file.name} -> {output_file.name} ({len(alignment)} seqs)')

print(f'\nConverted all files to {output_format} format')
