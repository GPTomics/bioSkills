'''Convert alignment between different formats'''

from Bio import AlignIO

input_file = 'alignment.aln'
input_format = 'clustal'

conversions = [
    ('output.fasta', 'fasta'),
    ('output.phy', 'phylip-relaxed'),
    ('output.nex', 'nexus'),
]

alignment = AlignIO.read(input_file, input_format)
print(f'Read alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns')

for output_file, output_format in conversions:
    AlignIO.write(alignment, output_file, output_format)
    print(f'Wrote: {output_file} ({output_format})')
