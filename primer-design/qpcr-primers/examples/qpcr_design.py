'''Design qPCR primers and TaqMan probe'''

import primer3
from Bio import SeqIO

record = SeqIO.read('target_gene.fasta', 'fasta')
sequence = str(record.seq)
print(f'Template: {record.id}, {len(sequence)}bp')

result = primer3.design_primers(
    seq_args={'SEQUENCE_TEMPLATE': sequence, 'SEQUENCE_ID': record.id},
    global_args={
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_PRODUCT_SIZE_RANGE': [[70, 150]],
        'PRIMER_NUM_RETURN': 3,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 58.0,
        'PRIMER_MAX_TM': 62.0,
        'PRIMER_INTERNAL_OPT_TM': 70.0,
        'PRIMER_INTERNAL_MIN_TM': 68.0,
        'PRIMER_INTERNAL_MAX_TM': 72.0,
        'PRIMER_INTERNAL_OPT_SIZE': 25,
    }
)

num_returned = result['PRIMER_PAIR_NUM_RETURNED']
print(f'\nFound {num_returned} primer/probe sets:\n')

for i in range(num_returned):
    print(f'Set {i+1}:')
    print(f'  Forward: {result[f"PRIMER_LEFT_{i}_SEQUENCE"]}')
    print(f'    Tm: {result[f"PRIMER_LEFT_{i}_TM"]:.1f}C')
    print(f'  Reverse: {result[f"PRIMER_RIGHT_{i}_SEQUENCE"]}')
    print(f'    Tm: {result[f"PRIMER_RIGHT_{i}_TM"]:.1f}C')
    print(f'  Probe: {result[f"PRIMER_INTERNAL_{i}_SEQUENCE"]}')
    print(f'    Tm: {result[f"PRIMER_INTERNAL_{i}_TM"]:.1f}C')
    print(f'  Product: {result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]}bp\n')
