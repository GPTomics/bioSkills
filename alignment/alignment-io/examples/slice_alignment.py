'''Slice and subset alignments'''

from Bio import AlignIO

alignment = AlignIO.read('alignment.aln', 'clustal')
print(f'Original: {len(alignment)} sequences, {alignment.get_alignment_length()} columns')

subset_seqs = alignment[0:5]
print(f'First 5 sequences: {len(subset_seqs)} sequences')

trimmed = alignment[:, 50:150]
print(f'Columns 50-150: {trimmed.get_alignment_length()} columns')

region = alignment[0:5, 50:150]
print(f'Combined slice: {len(region)} sequences, {region.get_alignment_length()} columns')

AlignIO.write(region, 'trimmed_subset.fasta', 'fasta')
print('Wrote trimmed subset to trimmed_subset.fasta')
