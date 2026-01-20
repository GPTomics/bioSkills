'''Design PCR primers with primer3-py'''

import primer3
from Bio import SeqIO

record = SeqIO.read('template.fasta', 'fasta')
sequence = str(record.seq)
print(f'Template: {record.id}, {len(sequence)}bp')

result = primer3.design_primers(
    seq_args={'SEQUENCE_TEMPLATE': sequence, 'SEQUENCE_ID': record.id},
    global_args={
        'PRIMER_PRODUCT_SIZE_RANGE': [[150, 300]],  # Standard PCR range; adjust for qPCR (70-150) or cloning (500+).
        'PRIMER_NUM_RETURN': 5,
        'PRIMER_OPT_TM': 60.0,  # Industry standard optimal Tm per primer3/IDT recommendations.
        'PRIMER_MIN_TM': 57.0,  # Allows +/-3C from optimal; tighter range improves specificity.
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 40.0,  # 40-60% GC ensures balanced primer stability and specificity.
        'PRIMER_MAX_GC': 60.0,
    }
)

num_returned = result['PRIMER_PAIR_NUM_RETURNED']
print(f'\nFound {num_returned} primer pairs:\n')

for i in range(num_returned):
    left_seq = result[f'PRIMER_LEFT_{i}_SEQUENCE']
    right_seq = result[f'PRIMER_RIGHT_{i}_SEQUENCE']
    left_tm = result[f'PRIMER_LEFT_{i}_TM']
    right_tm = result[f'PRIMER_RIGHT_{i}_TM']
    left_pos = result[f'PRIMER_LEFT_{i}']
    right_pos = result[f'PRIMER_RIGHT_{i}']
    product_size = result[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
    penalty = result[f'PRIMER_PAIR_{i}_PENALTY']

    print(f'Pair {i+1} (penalty: {penalty:.2f})')
    print(f'  Forward: {left_seq}')
    print(f'    Position: {left_pos[0]}, Length: {left_pos[1]}, Tm: {left_tm:.1f}C')
    print(f'  Reverse: {right_seq}')
    print(f'    Position: {right_pos[0]}, Length: {right_pos[1]}, Tm: {right_tm:.1f}C')
    print(f'  Product size: {product_size}bp\n')
