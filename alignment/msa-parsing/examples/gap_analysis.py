'''Analyze gap distribution in alignment'''

from Bio import AlignIO

alignment = AlignIO.read('alignment.fasta', 'fasta')
print(f'Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns\n')

print('Gaps per sequence:')
for record in alignment:
    gaps = str(record.seq).count('-')
    gap_pct = gaps / len(record.seq) * 100
    print(f'  {record.id}: {gaps} gaps ({gap_pct:.1f}%)')

print('\nGaps per column (showing columns with gaps):')
for col_idx in range(alignment.get_alignment_length()):
    column = alignment[:, col_idx]
    gaps = column.count('-')
    if gaps > 0:
        gap_pct = gaps / len(alignment) * 100
        print(f'  Column {col_idx}: {gaps} gaps ({gap_pct:.1f}%)')
