'''Analyze gap distribution in alignment'''

from Bio import AlignIO

alignment = AlignIO.read('alignment.fasta', 'fasta')
num_seqs = len(alignment)
num_cols = alignment.get_alignment_length()

print(f'Alignment: {num_seqs} sequences, {num_cols} columns\n')

gaps_per_seq = [(r.id, str(r.seq).count('-')) for r in alignment]
gaps_per_col = [alignment[:, i].count('-') for i in range(num_cols)]

total_gaps = sum(g for _, g in gaps_per_seq)
total_positions = num_seqs * num_cols

print('Gap Statistics:')
print(f'  Total gaps: {total_gaps}')
print(f'  Total positions: {total_positions}')
print(f'  Overall gap fraction: {total_gaps/total_positions*100:.1f}%')

print('\nGaps per sequence:')
for seq_id, gaps in sorted(gaps_per_seq, key=lambda x: -x[1]):
    pct = gaps / num_cols * 100
    print(f'  {seq_id}: {gaps} ({pct:.1f}%)')

gap_free = sum(1 for g in gaps_per_col if g == 0)
all_gap = sum(1 for g in gaps_per_col if g == num_seqs)
print(f'\nColumn statistics:')
print(f'  Gap-free columns: {gap_free} ({gap_free/num_cols*100:.1f}%)')
print(f'  All-gap columns: {all_gap}')

if any(g > 0 for g in gaps_per_col):
    # 0.5 (50%) threshold: columns with majority gaps often indicate alignment uncertainty
    # or indel events. Lower threshold (0.3) for stringent analysis; higher (0.7) if gaps expected.
    gappy_cols = [(i, g) for i, g in enumerate(gaps_per_col) if g > num_seqs * 0.5]
    if gappy_cols:
        print(f'\nColumns with >50% gaps: {len(gappy_cols)}')
        for i, g in gappy_cols[:10]:
            print(f'  Column {i}: {g}/{num_seqs} gaps')
