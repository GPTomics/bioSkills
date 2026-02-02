'''Design pegRNAs for prime editing experiments'''

from Bio.Seq import Seq

# Standard Cas9 scaffold sequence
CAS9_SCAFFOLD = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'

# Optimal parameter ranges
# PBS length: 13-17 nt (shorter = less stable, longer = more secondary structure)
# RT template: 10-20 nt for substitutions, longer for insertions
# Nick distance (PE3): 40-100 bp (too close = indels, too far = low efficiency)
PBS_LENGTH_RANGE = (13, 17)
RT_TEMPLATE_MIN = 10
NICK_DISTANCE_OPTIMAL = (40, 100)


def calculate_gc_content(seq):
    return sum(1 for nt in seq.upper() if nt in 'GC') / len(seq)


def estimate_melting_temp(seq):
    '''Estimate Tm using simple 2+4 rule for short oligos'''
    seq = seq.upper()
    at = sum(1 for nt in seq if nt in 'AT')
    gc = len(seq) - at
    return 2 * at + 4 * gc


def design_pegrna_substitution(target_seq, edit_pos, new_base, pbs_length=13, rt_length=15):
    '''Design pegRNA for a single nucleotide substitution

    Args:
        target_seq: ~100bp sequence centered on edit site
        edit_pos: Position of nucleotide to change (0-indexed)
        new_base: New nucleotide (A, C, G, or T)
        pbs_length: Primer binding site length (13-17nt optimal)
        rt_length: RT template length (10-20nt for substitutions)

    Note: This assumes a PAM (NGG) is available near the edit site.
    In practice, search for PAMs within ~30bp of edit for optimal efficiency.
    '''
    target_seq = target_seq.upper()

    # Find nick position (typically 3bp upstream of PAM)
    # For this example, assume nick is 3bp after edit position
    nick_pos = edit_pos + 3

    # Design spacer (20nt upstream of nick + 3bp for PAM)
    spacer_start = nick_pos - 17
    spacer = target_seq[spacer_start:spacer_start + 20]

    # Design PBS (reverse complement of sequence upstream of nick)
    pbs_region = target_seq[nick_pos - pbs_length:nick_pos]
    pbs = str(Seq(pbs_region).reverse_complement())

    # Design RT template with edit
    rt_region = list(target_seq[nick_pos:nick_pos + rt_length])
    edit_offset = edit_pos - nick_pos
    if 0 <= edit_offset < len(rt_region):
        rt_region[edit_offset] = new_base
    rt_template = str(Seq(''.join(rt_region)).reverse_complement())

    return {
        'spacer': spacer,
        'pbs': pbs,
        'rt_template': rt_template,
        'scaffold': CAS9_SCAFFOLD,
        'pbs_length': pbs_length,
        'pbs_gc': calculate_gc_content(pbs),
        'pbs_tm': estimate_melting_temp(pbs),
        'rt_length': rt_length,
        'edit_type': 'substitution',
        'original_base': target_seq[edit_pos],
        'new_base': new_base
    }


def optimize_pbs(nick_region, min_len=10, max_len=17):
    '''Find optimal PBS length based on GC content and Tm

    Optimal PBS characteristics:
    - Length: 13-17nt
    - GC content: 40-60%
    - Melting temperature: 45-65C
    '''
    options = []
    for length in range(min_len, max_len + 1):
        pbs_region = nick_region[-length:]
        pbs = str(Seq(pbs_region).reverse_complement())
        gc = calculate_gc_content(pbs)
        tm = estimate_melting_temp(pbs)

        score = 1.0
        if gc < 0.4 or gc > 0.6:
            score -= 0.2
        if tm < 45 or tm > 65:
            score -= 0.2
        if length < 13:
            score -= 0.1

        options.append({'length': length, 'sequence': pbs, 'gc': gc, 'tm': tm, 'score': score})

    return sorted(options, key=lambda x: x['score'], reverse=True)


def assemble_pegrna(spacer, rt_template, pbs, scaffold=CAS9_SCAFFOLD):
    '''Assemble full pegRNA sequence

    Order: 5'-[spacer]-[scaffold]-[RT template]-[PBS]-3'

    For U6 promoter: If spacer doesn't start with G, consider adding one
    (though this may affect targeting)
    '''
    full_seq = spacer + scaffold + rt_template + pbs
    return {
        'full_sequence': full_seq,
        'length': len(full_seq),
        'components': {
            'spacer': spacer,
            'scaffold': scaffold,
            'rt_template': rt_template,
            'pbs': pbs
        }
    }


def design_pe3_nick(target_seq, pegrna_nick_pos, min_dist=40, max_dist=100):
    '''Design PE3 second nicking guide

    PE3 strategy: Second nick on non-edited strand improves efficiency
    by biasing DNA repair toward the edited strand.

    Optimal distance: 40-100bp from pegRNA nick site
    - <40bp: Higher indel frequency
    - >100bp: Reduced efficiency benefit
    '''
    candidates = []

    # Search for PAM sites on opposite strand
    for offset in range(min_dist, max_dist + 1):
        # Downstream position
        pos = pegrna_nick_pos + offset
        if pos + 23 <= len(target_seq):
            region = target_seq[pos:pos + 23]
            if region[21:23] == 'GG':
                candidates.append({
                    'spacer': region[:20],
                    'position': pos,
                    'distance': offset,
                    'strand': '+',
                    'score': 1.0 - abs(offset - 70) / 100  # Prefer ~70bp
                })

    return sorted(candidates, key=lambda x: x['score'], reverse=True)[:5]


if __name__ == '__main__':
    # Example: Design pegRNA for a point mutation
    # Target: Correct a hypothetical A>G mutation

    target = 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAGG'
    edit_position = 40  # Position of the base to change
    original = target[edit_position]
    new_base = 'G' if original != 'G' else 'A'

    print(f'Designing pegRNA to change {original} -> {new_base} at position {edit_position}')
    print(f'Target region: ...{target[edit_position-10:edit_position+10]}...')
    print()

    # Design pegRNA
    pegrna = design_pegrna_substitution(target, edit_position, new_base)

    print('pegRNA Components:')
    print(f"  Spacer (20nt):     {pegrna['spacer']}")
    print(f"  PBS ({pegrna['pbs_length']}nt):        {pegrna['pbs']} (GC: {pegrna['pbs_gc']:.0%}, Tm: {pegrna['pbs_tm']}C)")
    print(f"  RT template:       {pegrna['rt_template']}")
    print()

    # Optimize PBS
    nick_region = target[edit_position - 17:edit_position + 3]
    pbs_options = optimize_pbs(nick_region)
    print('PBS length optimization:')
    for opt in pbs_options[:3]:
        print(f"  {opt['length']}nt: {opt['sequence']} (GC: {opt['gc']:.0%}, Tm: {opt['tm']}C, Score: {opt['score']:.2f})")
    print()

    # Assemble full pegRNA
    assembled = assemble_pegrna(pegrna['spacer'], pegrna['rt_template'], pegrna['pbs'])
    print(f"Full pegRNA length: {assembled['length']}nt")
