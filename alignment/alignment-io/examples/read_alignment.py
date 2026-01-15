'''Read and inspect a multiple sequence alignment'''

from Bio import AlignIO

alignment = AlignIO.read('sample_alignment.aln', 'clustal')

print(f'Alignment length: {alignment.get_alignment_length()} columns')
print(f'Number of sequences: {len(alignment)}')
print(f'\nSequence IDs:')
for record in alignment:
    print(f'  {record.id}: {len(record.seq)} bp')

print(f'\nFirst 60 columns of first sequence:')
print(alignment[0].seq[:60])
