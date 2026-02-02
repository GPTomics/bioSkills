'''Off-target prediction for CRISPR guides'''

import pandas as pd


# CFD score position-specific mismatch penalties (simplified)
# Full matrix from Doench et al. 2016 Nature Biotechnology
# Position 1 = PAM-proximal, Position 20 = PAM-distal
CFD_POSITION_WEIGHTS = {
    1: 0.3,   # PAM-proximal: mismatches severely reduce cutting
    2: 0.35,
    3: 0.4,
    4: 0.45,
    5: 0.5,
    10: 0.6,
    15: 0.75,
    18: 0.85,
    19: 0.9,
    20: 0.95,  # PAM-distal: mismatches have less effect
}

# PAM variant penalties (SpCas9 prefers NGG)
CFD_PAM_SCORES = {
    'AGG': 0.26,  # NAG PAM has ~26% activity vs NGG
    'CGG': 0.11,
    'TGG': 0.02,
    'GAG': 0.07,
    'GCG': 0.03,
    'GTG': 0.02,
    'GGA': 0.01,
    'GGC': 0.01,
    'GGT': 0.01,
}


def calculate_cfd_score(guide, off_target_seq):
    '''Calculate CFD score for an off-target site

    CFD (Cutting Frequency Determination) predicts cleavage probability
    relative to the on-target site.

    Interpretation:
    - 1.0: Perfect match (on-target)
    - >0.5: High cleavage probability (concerning off-target)
    - 0.1-0.5: Moderate probability
    - <0.1: Low probability (likely acceptable)
    '''
    guide = guide.upper()
    off_target = off_target_seq.upper()

    score = 1.0

    # Apply position-specific mismatch penalties
    for i in range(min(20, len(guide), len(off_target))):
        if guide[i] != off_target[i]:
            pos = i + 1  # 1-based position
            penalty = CFD_POSITION_WEIGHTS.get(pos, 0.5)
            score *= penalty

    # Apply PAM penalty if sequence includes PAM
    if len(off_target) >= 23:
        pam = off_target[20:23]
        if pam != 'NGG' and pam in CFD_PAM_SCORES:
            score *= CFD_PAM_SCORES[pam]

    return score


def count_mismatches(seq1, seq2):
    '''Count mismatches between two sequences'''
    return sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a != b)


def prepare_cas_offinder_input(guides, genome_dir, max_mismatches=4):
    '''Prepare input file for Cas-OFFinder

    Args:
        guides: List of 20nt guide sequences
        genome_dir: Path to directory with genome .2bit or indexed fasta
        max_mismatches: Search tolerance (4 is good balance)
                       0-2: Only find very similar sites
                       3-4: Standard search (recommended)
                       5-6: Comprehensive but slow
    '''
    lines = [genome_dir, 'N' * 20 + 'NGG']
    for guide in guides:
        lines.append(f'{guide.upper()}NNN {max_mismatches}')
    return '\n'.join(lines)


def parse_cas_offinder_output(filepath):
    '''Parse Cas-OFFinder output file'''
    columns = ['guide', 'chrom', 'position', 'sequence', 'strand', 'mismatches']
    return pd.read_csv(filepath, sep='\t', header=None, names=columns)


def analyze_off_targets(guide, off_target_df):
    '''Analyze off-targets for a single guide

    Returns:
        dict with specificity metrics and flagged sites
    '''
    guide_ots = off_target_df[off_target_df['guide'] == guide].copy()

    # Calculate CFD for each off-target
    guide_ots['cfd_score'] = guide_ots['sequence'].apply(
        lambda seq: calculate_cfd_score(guide, seq)
    )

    # Aggregate specificity (higher = better)
    # Formula: 1 / (1 + sum of CFD scores)
    cfd_sum = guide_ots['cfd_score'].sum()
    specificity = 1 / (1 + cfd_sum)

    # Count by mismatch level
    mm_counts = guide_ots.groupby('mismatches').size().to_dict()

    # Flag high-risk off-targets (CFD > 0.5)
    high_risk = guide_ots[guide_ots['cfd_score'] > 0.5]

    return {
        'guide': guide,
        'specificity_score': specificity,
        'total_off_targets': len(guide_ots),
        'mismatch_counts': mm_counts,
        'high_risk_count': len(high_risk),
        'high_risk_sites': high_risk.to_dict('records') if len(high_risk) > 0 else []
    }


if __name__ == '__main__':
    # Example: Analyze off-targets for a guide
    guide = 'ATCGATCGATCGATCGATCG'

    # Simulated off-target data (normally from Cas-OFFinder)
    off_targets = [
        {'sequence': 'ATCGATCGATCGATCGATCGAGG', 'mismatches': 0},  # On-target
        {'sequence': 'ATCGATCGATCGATCGATCCAGG', 'mismatches': 1},  # 1 mismatch
        {'sequence': 'ATCGATCGATCGATCGATTCAGG', 'mismatches': 2},  # 2 mismatches
        {'sequence': 'ATCGATCGATCGAACGATCGAGG', 'mismatches': 2},  # 2 mismatches
        {'sequence': 'ATCGATCGATCGTTCGATCGAGG', 'mismatches': 2},  # 2 mismatches
    ]

    print(f'Guide: {guide}')
    print('\nOff-target CFD scores:')
    for ot in off_targets:
        cfd = calculate_cfd_score(guide, ot['sequence'])
        print(f"  {ot['sequence']} ({ot['mismatches']} mm): CFD = {cfd:.3f}")

    # Calculate aggregate specificity
    cfd_sum = sum(calculate_cfd_score(guide, ot['sequence']) for ot in off_targets[1:])
    specificity = 1 / (1 + cfd_sum)
    print(f'\nAggregate specificity score: {specificity:.3f}')
    print('(>0.9 = excellent, 0.7-0.9 = good, 0.5-0.7 = moderate, <0.5 = poor)')
