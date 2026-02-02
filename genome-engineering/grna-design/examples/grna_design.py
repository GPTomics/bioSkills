'''Design CRISPR guide RNAs with activity scoring'''

from Bio.Seq import Seq
import re

# Rule Set 2 position-weight matrix (simplified from Doench et al. 2016)
# Position 0 = PAM-distal, Position 19 = PAM-proximal
# Full model available at: https://github.com/MicrosoftResearch/Azimuth
RULE_SET_2_WEIGHTS = {
    0: {'A': 0, 'C': 0, 'G': 0.08, 'T': -0.08},
    1: {'A': 0.02, 'C': -0.06, 'G': 0.06, 'T': -0.02},
    2: {'A': 0.01, 'C': 0.02, 'G': -0.02, 'T': -0.01},
    17: {'A': -0.04, 'C': 0.02, 'G': 0.03, 'T': -0.01},
    18: {'A': -0.07, 'C': 0.13, 'G': -0.01, 'T': -0.05},
    19: {'A': -0.07, 'C': 0.03, 'G': 0.11, 'T': -0.07},
}


def calculate_gc_content(sequence):
    gc = sum(1 for nt in sequence.upper() if nt in 'GC')
    return gc / len(sequence)


def find_pam_sites(sequence, pam='NGG', guide_length=20):
    '''Find all PAM sites and extract guide sequences

    PAM patterns:
    - NGG: SpCas9 (most common, high activity)
    - TTTN: Cas12a/Cpf1 (lower off-targets, staggered cut)
    - NNGRRT: SaCas9 (smaller, fits in AAV)
    '''
    sequence = sequence.upper()
    guides = []

    if pam == 'NGG':
        for match in re.finditer(r'(?=(.GG))', sequence):
            pos = match.start()
            if pos >= guide_length:
                guide = sequence[pos - guide_length:pos]
                guides.append({
                    'sequence': guide,
                    'pam': sequence[pos:pos + 3],
                    'position': pos - guide_length,
                    'strand': '+'
                })

        rc_seq = str(Seq(sequence).reverse_complement())
        for match in re.finditer(r'(?=(.GG))', rc_seq):
            pos = match.start()
            if pos >= guide_length:
                guide = rc_seq[pos - guide_length:pos]
                original_pos = len(sequence) - pos
                guides.append({
                    'sequence': guide,
                    'pam': rc_seq[pos:pos + 3],
                    'position': original_pos,
                    'strand': '-'
                })

    return guides


def score_guide_activity(guide_seq):
    '''Score guide on-target activity (0-1 scale)

    Scoring criteria:
    - GC content 40-70%: optimal range for Cas9 activity
    - No poly-T stretches: terminates U6/H1 Pol III transcription
    - Position-specific nucleotide preferences from Rule Set 2

    Interpretation:
    - >0.6: High activity expected (use these first)
    - 0.4-0.6: Moderate activity (acceptable if no better options)
    - <0.4: Low activity (avoid if possible)
    '''
    guide_seq = guide_seq.upper()
    score = 0.5

    # GC content penalty (optimal 40-70%)
    gc = calculate_gc_content(guide_seq)
    if gc < 0.4 or gc > 0.7:
        score -= 0.15

    # Poly-T penalty (terminates transcription)
    if 'TTTT' in guide_seq:
        score -= 0.3

    # Position-specific scoring
    for pos, weights in RULE_SET_2_WEIGHTS.items():
        if pos < len(guide_seq):
            nt = guide_seq[pos]
            score += weights.get(nt, 0)

    return max(0, min(1, score))


def design_guides_for_gene(gene_sequence, n_guides=5):
    '''Design top N guides for a gene sequence'''
    all_guides = find_pam_sites(gene_sequence)

    for guide in all_guides:
        guide['activity_score'] = score_guide_activity(guide['sequence'])
        guide['gc_content'] = calculate_gc_content(guide['sequence'])

    all_guides.sort(key=lambda x: x['activity_score'], reverse=True)
    return all_guides[:n_guides]


if __name__ == '__main__':
    # Example: Design guides for a target region
    # Use BRCA1 exon 2 as example (first 200bp)
    target_seq = 'ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGA'

    guides = design_guides_for_gene(target_seq, n_guides=5)

    print('Top 5 guides for target region:')
    for i, g in enumerate(guides, 1):
        print(f"{i}. {g['sequence']} (PAM: {g['pam']}, Score: {g['activity_score']:.2f}, GC: {g['gc_content']:.0%})")
