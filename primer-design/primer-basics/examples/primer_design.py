'''Design and rank PCR primer pairs with primer3-py, then read the penalty ranking and the
explain tallies. primer3 scores only this template -- route the chosen pair to a genome
specificity check (primer-specificity) before ordering. Self-contained: no input file, no output file.'''
# Reference: primer3-py 2.3+ | Verify API if version differs

import primer3

TM_OPT, TM_MIN, TM_MAX = 60.0, 58.0, 62.0   # primer Tm window (C); narrow band keeps the pair matched
PAIR_DIFF_TM_MAX = 2.0                       # max Tm gap between the two primers; >2-3 C unbalances amplification
GC_MIN, GC_MAX = 40.0, 60.0                  # percent; primer3 default 20-80 is far too wide
PRODUCT_RANGE = [[150, 400]]                 # bp; standard PCR amplicon
MV, DV, DNTP, DNA = 50.0, 1.5, 0.6, 50.0     # mM monovalent, mM Mg2+, mM dNTP, nM oligo; free Mg2+ ~= DV - DNTP

template = ('GCAATTCGGCATTAGCCTGAGCGTACGTACGTTAGCCAGTCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCATCGG'
            'ATCCAGTCAGTCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCATGCTAGCGCGATCGTAGCATCGATCGATG'
            'CATCGTAGCTAGCTAGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCATCGATCGATCGTAGCTGATCGATCGTAGCA'
            'TCGATCGATCGATCGTAGCTAGCTAGCATCGATCGTAGCATCGTTACGGCATTAGCCTGAGCGTACGTACGTTAGCCAGT')

result = primer3.design_primers(
    seq_args={'SEQUENCE_ID': 'demo', 'SEQUENCE_TEMPLATE': template, 'SEQUENCE_TARGET': [150, 50]},
    global_args={
        'PRIMER_PICK_LEFT_PRIMER': 1, 'PRIMER_PICK_RIGHT_PRIMER': 1, 'PRIMER_NUM_RETURN': 5,
        'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': TM_OPT, 'PRIMER_MIN_TM': TM_MIN, 'PRIMER_MAX_TM': TM_MAX,
        'PRIMER_PAIR_MAX_DIFF_TM': PAIR_DIFF_TM_MAX,
        'PRIMER_MIN_GC': GC_MIN, 'PRIMER_MAX_GC': GC_MAX,
        'PRIMER_GC_CLAMP': 1,
        'PRIMER_PRODUCT_SIZE_RANGE': PRODUCT_RANGE,
        'PRIMER_SALT_MONOVALENT': MV, 'PRIMER_SALT_DIVALENT': DV,
        'PRIMER_DNTP_CONC': DNTP, 'PRIMER_DNA_CONC': DNA,
        'PRIMER_EXPLAIN_FLAG': 1})

n = result['PRIMER_PAIR_NUM_RETURNED']
print(f'Returned {n} pairs (ranked by penalty; index 0 = lowest)')

if n == 0:
    print('Zero pairs: read the explain tallies and loosen the dominant bucket, one at a time.')
    print('  left :', result['PRIMER_LEFT_EXPLAIN'])
    print('  right:', result['PRIMER_RIGHT_EXPLAIN'])
    print('  pair :', result['PRIMER_PAIR_EXPLAIN'])
else:
    for i in range(n):
        fwd, rev = result[f'PRIMER_LEFT_{i}_SEQUENCE'], result[f'PRIMER_RIGHT_{i}_SEQUENCE']
        fwd_tm, rev_tm = result[f'PRIMER_LEFT_{i}_TM'], result[f'PRIMER_RIGHT_{i}_TM']
        size, penalty = result[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'], result[f'PRIMER_PAIR_{i}_PENALTY']
        print(f'pair {i} penalty={penalty:.2f} product={size}bp dTm={abs(fwd_tm - rev_tm):.1f}C')
        print(f'  F {fwd} Tm={fwd_tm:.1f}')
        print(f'  R {rev} Tm={rev_tm:.1f}')
    print('\nNext step: BLAST / in-silico PCR this pair for genome specificity (primer-specificity).')
