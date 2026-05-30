'''Build a PBS x RTT pegRNA panel for one nick and search for a PE3/PE3b nicking guide.

Prime editing has no universal PBS/RTT optimum, so the deliverable is a PANEL to test or
rank with a trained model (PRIDICT2.0 / DeepPrime), not a single pegRNA. This script
enforces the hard rules a hand design gets wrong: prepend (never replace) the U6 5' G, and
drop extensions whose first templated base is C (PrimeDesign's --filter_c1_extension). For
real designs, prefer PrimeDesign with its exact inline notation; this illustrates the geometry.
'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio.Seq import Seq

PBS_LADDER = (10, 12, 14, 16)   # scan PBS lengths; optimum is locus-specific (GC/Tm), not fixed
RTT_HOMOLOGY = (11, 16, 22)     # RTT = nick-to-edit + edit + 3' homology tail (~10-16 nt typical)


def prepend_u6_g(spacer):
    '''U6 (Pol III) initiates best with a 5' G: prepend one, never replace the first base.'''
    return spacer if spacer.startswith('G') else 'G' + spacer


def wallace_tm(seq):
    '''Illustrative short-oligo Tm (2*AT + 4*GC). Real PBS tuning matches Tm across a panel.'''
    at = sum(c in 'AT' for c in seq.upper())
    return 2 * at + 4 * (len(seq) - at)


def build_pegrna_panel(target, nick, edit_pos, new_base):
    '''Return a PBS x RTT panel of pegRNA 3' extensions for a substitution at edit_pos.

    target: forward strand; nick: forward-strand index of the nick (3 bp 5' of the PAM on
    the protospacer strand); edit_pos: forward index to change; new_base: replacement.
    The 3' extension is [RTT][PBS] (5'->3'): RTT = revcomp of the edited region 3' of the
    nick; PBS = revcomp of the region 5' of the nick.
    '''
    target = target.upper()
    panel = []
    for rtt_homology in RTT_HOMOLOGY:
        rtt_len = (edit_pos - nick) + 1 + rtt_homology
        flap = list(target[nick:nick + rtt_len])
        off = edit_pos - nick
        if not (0 <= off < len(flap)):
            continue
        flap[off] = new_base
        rtt = str(Seq(''.join(flap)).reverse_complement())
        starts_c = rtt.startswith('C')   # C at the +1 extension position lowers efficiency
        for pbs_len in PBS_LADDER:
            pbs = str(Seq(target[nick - pbs_len:nick]).reverse_complement())
            extension = rtt + pbs
            panel.append({'pbs_len': pbs_len, 'rtt_len': rtt_len, 'extension': extension,
                          'pbs': pbs, 'pbs_tm': wallace_tm(pbs), 'flap_starts_with_c': starts_c})
    return [p for p in panel if not p['flap_starts_with_c']]   # apply the c1 filter


def find_pe3_nick(target, peg_nick, edited=False, edit_pos=None):
    '''Search the non-edited (complementary) strand for an ngRNA nick ~40-100 bp from the pegRNA nick.

    The pegRNA nicks the protospacer/PAM strand (forward here); a PE3/PE3b ngRNA must nick the
    OPPOSITE strand, so its NGG protospacer lies on the reverse-complement strand -- searched here
    by scanning the reverse complement and mapping each candidate nick back to forward coordinates.
    edited=False -> PE3 (any nearby complementary-strand protospacer). edited=True -> PE3b (the ngRNA
    protospacer must span the edit so it is created only after editing); only possible when the edit
    makes/breaks a protospacer.
    '''
    target = target.upper()
    n = len(target)
    rc = str(Seq(target).reverse_complement())
    hits = []
    for j in range(n - 22):
        if rc[j + 21:j + 23] != 'GG':
            continue
        nick_fwd = n - 1 - (j + 17)                       # ngRNA nick (3 bp 5' of PAM) in forward coords
        distance = abs(nick_fwd - peg_nick)
        if not (40 <= distance <= 100):
            continue
        proto_lo, proto_hi = n - 1 - (j + 22), n - 1 - j  # forward span of the ngRNA protospacer+PAM
        if edited and edit_pos is not None and not (proto_lo <= edit_pos <= proto_hi):
            continue
        hits.append({'spacer': rc[j:j + 20], 'position': nick_fwd, 'distance': distance,
                     'mode': 'PE3b' if edited else 'PE3'})
    hits.sort(key=lambda h: abs(h['distance'] - 70))   # sweet spot ~50-90 bp
    return hits[:5]


if __name__ == '__main__':
    target = ('CTGACCTGTAGCAATTCGGCAGTCAGGTACCATGGCTAGCTAGGGCCTAGACTTCGATCCAGGTACGT'
              'TACGGCATTCGATCGGATCCAAGTTCCGATCGATCGTAGCTAGCTAGCTAGGCTAGCATCGATCGTAG')
    nick = 30          # forward-strand nick position (3 bp 5' of the PAM on the protospacer strand)
    edit_pos = 36      # base to change
    new_base = 'A'

    print(f'Edit: {target[edit_pos]}->{new_base} at {edit_pos}; pegRNA nick at {nick}')
    panel = build_pegrna_panel(target, nick, edit_pos, new_base)
    print(f'\nPBS x RTT panel ({len(panel)} candidates after the c1 filter; rank with PRIDICT2.0/DeepPrime, then test):')
    for p in panel[:6]:
        print(f"  PBS={p['pbs_len']} (Tm~{p['pbs_tm']}) RTT={p['rtt_len']} ext={p['extension']}")

    print('\nPE3 nicking-guide candidates (non-edited strand, ~40-100 bp away):')
    for h in find_pe3_nick(target, nick):
        print(f"  {h['mode']} spacer={h['spacer']} dist={h['distance']} pos={h['position']}")
    print('\nReminders: prepend (not replace) the 5-prime G; add a PAM-disrupting + MMR-evading '
          'silent edit; append a tevopreQ1 3-prime motif for expressed pegRNAs; report edit:indel purity.')
