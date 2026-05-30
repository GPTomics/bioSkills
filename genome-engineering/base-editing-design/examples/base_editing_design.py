'''Scan base-editor guides: place the target base in the window, read bystanders + context.

The job is product PURITY, not efficiency: report which bases in the window will edit and the
sequence context of each, not a single efficiency number. The position-efficiency values are
COARSE ILLUSTRATIVE GRADIENTS (not measurements) to show the peak-at-5-7 shape -- for a real
genotype spectrum use BE-Hive/DeepBE, then validate by amplicon NGS. Numbering is PAM-distal:
position 1 = 5' end of the spacer, position 20 = adjacent to the PAM.
'''
# Reference: biopython 1.83+ | Verify API if version differs

import re

CBE_WINDOW = (4, 8)   # activity gradient, peak ~5-7
ABE_WINDOW = (4, 8)   # ABE7.10 tighter (~4-7); ABE8e WIDER (~3-11) -- editor-specific, re-check on a switch

# Illustrative activity gradient (shape only; NOT measured efficiency). Use BE-Hive for real outcomes.
POSITION_GRADIENT = {3: 0.2, 4: 0.7, 5: 0.95, 6: 1.0, 7: 0.85, 8: 0.5, 9: 0.2}

# APOBEC1 (CBE) 5'-neighbor context preference: TC favored, GC strongly disfavored (Komor 2016).
CBE_CONTEXT = {'T': 1.0, 'C': 0.8, 'A': 0.6, 'G': 0.4}

STOP_CODONS_VIA_CBE = {'CAA': 'TAA', 'CAG': 'TAG', 'CGA': 'TGA'}   # C->T on the sense strand -> stop


def base_pair_change(ref, alt):
    '''Write the edit as a base-pair change to pick the editor family (resolves strand confusion).'''
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    pair = f'{ref}*{comp[ref]} -> {alt}*{comp[alt]}'
    if (ref, alt) in (('C', 'T'), ('G', 'A')):
        return pair, 'CBE'
    if (ref, alt) in (('A', 'G'), ('T', 'C')):
        return pair, 'ABE'
    if (ref, alt) == ('C', 'G'):
        return pair, 'CGBE'
    return pair, 'prime-editing (transversion outside base-editor chemistry)'


def find_targets(sequence, target_pos, editor='CBE'):
    '''Forward-strand NGG guides that place the base at target_pos in the editing window.

    Also scan the reverse strand (CCN PAMs) in practice: a forward G is a C on the displaced
    strand for a reverse guide. Returns guides with bystanders and per-base context.
    '''
    seq = sequence.upper()
    window = CBE_WINDOW if editor == 'CBE' else ABE_WINDOW
    base = 'C' if editor == 'CBE' else 'A'
    guides = []
    for m in re.finditer(r'(?=([ACGT]GG))', seq):
        start = m.start() - 20
        if start < 0:
            continue
        pos_in_spacer = target_pos - start + 1   # PAM-distal numbering (1 = 5' end)
        if not (window[0] <= pos_in_spacer <= window[1] and seq[target_pos] == base):
            continue
        spacer = seq[start:m.start()]
        bystanders = [i + 1 for i in range(window[0] - 1, window[1])
                      if i < len(spacer) and spacer[i] == base and (start + i) != target_pos]
        context = CBE_CONTEXT.get(spacer[pos_in_spacer - 2], 0.5) if (editor == 'CBE' and pos_in_spacer > 1) else None
        guides.append({'spacer': spacer, 'pam': seq[m.start():m.start() + 3], 'target_position': pos_in_spacer,
                       'gradient': POSITION_GRADIENT.get(pos_in_spacer, 0.05), 'bystanders': bystanders,
                       'context_pref': context})
    return sorted(guides, key=lambda g: (len(g['bystanders']), -g['gradient']))


def find_cbe_stops(cds):
    '''Find in-frame CAA/CAG/CGA codons a CBE could convert to a stop (CRISPR-STOP/iSTOP).

    Prefer EARLY stops (NMD-competent, before functional domains); a base-editor KO is only
    as complete as the editing, so verify protein.
    '''
    return [{'codon_index': i // 3, 'codon': cds[i:i + 3], 'stop': STOP_CODONS_VIA_CBE[cds[i:i + 3]]}
            for i in range(0, len(cds) - 2, 3) if cds[i:i + 3] in STOP_CODONS_VIA_CBE]


if __name__ == '__main__':
    print('Editor triage (write the edit as a base-pair change first):')
    for ref, alt in [('C', 'T'), ('G', 'A'), ('A', 'G'), ('C', 'G'), ('A', 'T')]:
        pair, editor = base_pair_change(ref, alt)
        print(f'  {ref}->{alt}:  {pair}  -> {editor}')

    seq = 'GGCATCGATCGATCGATCGATCAGCTAGCATCGATCGATCGATCGATCGATCGGCATGCTAGCTAGCATCGTAGCTAGCTGACT'
    target = next(i for i in range(20, len(seq)) if seq[i] == 'C' and find_targets(seq, i, 'CBE'))
    print(f'\nCBE guides placing the C at sequence index {target} inside the editing window:')
    for g in find_targets(seq, target, 'CBE')[:3]:
        ctx = f", 5'-context pref {g['context_pref']:.1f}" if g['context_pref'] is not None else ''
        flag = 'PURE (no bystanders)' if not g['bystanders'] else f"bystander C at {g['bystanders']}"
        print(f"  {g['spacer']} PAM={g['pam']} target@pos{g['target_position']} "
              f"(illustrative activity {g['gradient']:.2f}{ctx}) -> {flag}")

    cds = 'ATGGCACAAGCTCGAGGCTAGCAGAAATAG'
    print(f'\nCRISPR-STOP candidates (early CAA/CAG/CGA -> stop via CBE): {find_cbe_stops(cds)}')
    print('Report the genotype spectrum (BE-Hive/DeepBE), then validate by NGS; rank by purity, not efficiency.')
