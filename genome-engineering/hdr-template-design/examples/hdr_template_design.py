'''Design an HDR ssODN with a MANDATORY codon-checked blocking mutation, plus junction primers.

The single most important rule: a donor with no blocking mutation is self-destructing -- the
corrected allele keeps an intact PAM, Cas9 re-cuts it, and the indel reads out as "HDR failed."
This script disrupts the PAM synonymously (verified against a codon table) so the edit survives,
emits both ssODN strands (the Richardson non-target-strand/asymmetric rule is a prior to TEST,
not a law), and designs genotyping primers with primer3-py. primer3 designs PRIMERS, not arms.
'''
# Reference: biopython 1.83+, primer3-py 2.0+ | Verify API if version differs

from Bio.Seq import Seq
from Bio.Data import CodonTable
import primer3

FWD = CodonTable.standard_dna_table.forward_table   # codon -> amino acid (stop codons absent)
SSODN_ARM = 50          # per-arm nt; 30-60 workable, total stays under the ~200 nt synthesis ceiling
EDIT_TO_CUT_MAX = 10    # HDR incorporation falls sharply beyond ~10 bp (Paquet 2016)


def aa(codon):
    return FWD.get(codon.upper(), '*')   # '*' for stop / unknown


def synonymous_swap(codon, pos_in_codon):
    '''Return a synonymous codon differing at pos_in_codon (to break a PAM), or None if impossible.'''
    for base in 'ACGT':
        if base == codon[pos_in_codon]:
            continue
        cand = codon[:pos_in_codon] + base + codon[pos_in_codon + 1:]
        if aa(cand) == aa(codon) and aa(codon) != '*':
            return cand
    return None


def block_pam_synonymously(donor, pam_start, cds_offset):
    '''Disrupt an NGG PAM with a synonymous change so the edited allele cannot be re-cut.

    pam_start: index of the PAM's N in `donor`; cds_offset: index where the reading frame begins.
    Targets the 2nd or 3rd PAM base (the GG); returns (blocked_donor, description) or (donor, reason).
    '''
    for g in (pam_start + 2, pam_start + 1):   # prefer the 3rd PAM base, then the 2nd
        if donor[g] != 'G':
            continue
        codon_i = cds_offset + 3 * ((g - cds_offset) // 3)
        codon = donor[codon_i:codon_i + 3]
        if len(codon) < 3:
            continue
        swapped = synonymous_swap(codon, g - codon_i)
        if swapped:
            return donor[:codon_i] + swapped + donor[codon_i + 3:], f'silent PAM block at codon {codon}->{swapped}'
    return donor, 'no synonymous PAM block found -- introduce silent SEED-region mutations instead'


def design_ssodn(target, cut, edit_pos, new_base, pam_start, cds_offset, arm=SSODN_ARM):
    '''Build an ssODN donor with the edit and a codon-checked blocking mutation; return both strands.'''
    if abs(edit_pos - cut) > EDIT_TO_CUT_MAX:
        return {'error': f'edit is {abs(edit_pos - cut)} bp from the cut (>{EDIT_TO_CUT_MAX}); '
                         'choose a closer-cutting guide or switch to base/prime editing'}
    if cut < arm or cut + arm > len(target):
        return {'error': f'need >={arm} bp of flank on both sides of the cut'}
    edited = target[:edit_pos] + new_base + target[edit_pos + 1:]
    blocked, block_note = block_pam_synonymously(edited, pam_start, cds_offset)
    sense = blocked[cut - arm:cut + arm]
    return {'sense': sense, 'antisense': str(Seq(sense).reverse_complement()), 'length': len(sense),
            'blocking': block_note,
            'note': 'TEST both strands (Richardson rule is a prior, not a law); add 2-3 phosphorothioate '
                    'caps per end; for >50 bp inserts use lssDNA, not a longer oligo'}


def junction_primers(donor, product_min=90, product_max=150):
    '''Design genotyping primers across the insertion junction with primer3-py.'''
    res = primer3.bindings.design_primers(
        {'SEQUENCE_ID': 'hdr_junction', 'SEQUENCE_TEMPLATE': donor},
        {'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 24,
         'PRIMER_OPT_TM': 60.0, 'PRIMER_PRODUCT_SIZE_RANGE': [[product_min, product_max]]})
    return {'left': res.get('PRIMER_LEFT_0_SEQUENCE'), 'right': res.get('PRIMER_RIGHT_0_SEQUENCE'),
            'pairs_returned': res.get('PRIMER_PAIR_NUM_RETURNED', 0)}


if __name__ == '__main__':
    # In-frame coding target (cds_offset=0). Search for an NGG PAM that (a) sits where a
    # synonymous block is possible and (b) leaves >=arm flank on both sides of its cut.
    target = ('ATGGCCGAGACCAAGTTCGACGAGCTGAAGGCCTACATCGAGAAGCTGGGCGAGAAGGGCTTCGAC'
              'GAGGTGAAGCTGTACGGCGACGAGCTGAAGTTCTACGGCGAGAAGCTGTTCGAGAAGCTGTTCTAA')
    cds_offset = 0
    chosen = None
    for i in range(SSODN_ARM, len(target) - SSODN_ARM - 3):
        if target[i + 1:i + 3] == 'GG':                 # NGG PAM at i
            blocked, note = block_pam_synonymously(target, i, cds_offset)
            if blocked != target:
                chosen = (i, i - 3, note); break          # SpCas9 cut ~3 bp 5' of the PAM
    pam_start, cut, _ = chosen
    edit_pos = cut                                        # install the edit at the cut (within 10 bp)
    new_base = 'A' if target[edit_pos] != 'A' else 'C'

    print(f'Blockable NGG PAM at {pam_start}, cut at {cut}, edit {target[edit_pos]}->{new_base} at {edit_pos}')
    print('HDR ssODN (codon-checked blocking mutation is mandatory):')
    donor = design_ssodn(target, cut, edit_pos, new_base, pam_start, cds_offset)
    for k, v in donor.items():
        print(f'  {k}: {v}')

    full = target[:edit_pos] + new_base + target[edit_pos + 1:]
    print('\nJunction genotyping primers (primer3-py):')
    print(' ', junction_primers(full))
    print('\nIf edit-to-cut > 10 bp, or the cell is post-mitotic, reconsider: base/prime editing, or HITI.')
