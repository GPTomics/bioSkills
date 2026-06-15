'''Co-design qPCR primers and a TaqMan probe: short amplicon, Tm-matched primers, probe Tm 8-10 C
above the primers (raised explicitly -- primer3 internal defaults equal the primer Tm), and no 5'-G
on the probe (enforced with PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME='HNNNN', IUPAC H = not G).
A passing design still needs a genome specificity check (primer-specificity) and a standard curve.
Self-contained, no file IO.'''
# Reference: primer3-py 2.3+ | Verify API if version differs

import primer3

PRIMER_TM = (60.0, 58.0, 62.0)        # opt/min/max C; tight window, pair matched within 2 C
PROBE_TM = (70.0, 68.0, 72.0)         # opt/min/max C; 8-10 C above primers so the probe is bound when Taq cleaves it
AMPLICON = [[80, 150]]                 # bp; short amplicon protects ~100% efficiency

# AT-balanced primer flanks with a GC-rich central window so a high-Tm probe site exists
template = ('ATGGCACCGTCAAGGCTGAGAACGGGAAGCTTGTCATCAATGGAAATCCCATCACCATCTTCCAGGAGCGAGATCCCTCCAAA'
            'GCCGGCTGCCGGAGGCCGCCGGCAGCCGGCCGCAGGCCG'
            'ATCAAGTGGGGCGATGCTAGCTACGATCAGTACCATGGAGAAGGCTGGGTCATCATCTCTAATGCCCATGTTCGTCATGGTGT')

result = primer3.design_primers(
    seq_args={'SEQUENCE_ID': 'assay1', 'SEQUENCE_TEMPLATE': template},
    global_args={
        'PRIMER_PICK_LEFT_PRIMER': 1, 'PRIMER_PICK_RIGHT_PRIMER': 1, 'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_PRODUCT_SIZE_RANGE': AMPLICON, 'PRIMER_NUM_RETURN': 3,
        'PRIMER_OPT_TM': PRIMER_TM[0], 'PRIMER_MIN_TM': PRIMER_TM[1], 'PRIMER_MAX_TM': PRIMER_TM[2],
        'PRIMER_PAIR_MAX_DIFF_TM': 2.0,
        'PRIMER_INTERNAL_OPT_TM': PROBE_TM[0], 'PRIMER_INTERNAL_MIN_TM': PROBE_TM[1], 'PRIMER_INTERNAL_MAX_TM': PROBE_TM[2],
        'PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME': 'HNNNN',
        'PRIMER_EXPLAIN_FLAG': 1})

n = result['PRIMER_PAIR_NUM_RETURNED']
print(f'Returned {n} primer/probe sets')

if n == 0:
    print('No set: probe Tm window may be unreachable on this template.')
    print('  internal:', result['PRIMER_INTERNAL_EXPLAIN'])
else:
    for i in range(n):
        fwd, rev = result[f'PRIMER_LEFT_{i}_SEQUENCE'], result[f'PRIMER_RIGHT_{i}_SEQUENCE']
        probe = result[f'PRIMER_INTERNAL_{i}_SEQUENCE']
        ptm = result[f'PRIMER_INTERNAL_{i}_TM']
        print(f'set {i} amplicon={result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]}bp')
        print(f'  F {fwd} Tm={result[f"PRIMER_LEFT_{i}_TM"]:.1f}')
        print(f'  R {rev} Tm={result[f"PRIMER_RIGHT_{i}_TM"]:.1f}')
        print(f'  P {probe} 5prime={probe[0]} Tm={ptm:.1f} (offset +{ptm - result[f"PRIMER_LEFT_{i}_TM"]:.1f}C)')
    print('\nNext: genome specificity (primer-specificity) and a standard curve for efficiency (MIQE).')
