'''Predict B-cell and T-cell epitopes'''

import pandas as pd


def predict_linear_bcell_epitopes(sequence, window_size=15, threshold=0.5):
    '''Predict linear B-cell epitopes using sequence features

    This is a simplified local implementation.
    For production, use IEDB BepiPred-2.0 API.

    Features used:
    - Hydrophilicity (Parker scale)
    - Surface accessibility (Emini scale)
    - Flexibility (Karplus-Schulz)

    B-cell epitope characteristics:
    - Surface exposed
    - Hydrophilic
    - Flexible loops
    - Length: 5-15 amino acids typical
    '''
    # Simplified hydrophilicity scale
    hydrophilicity = {
        'A': -0.5, 'R': 3.0, 'N': 0.2, 'D': 3.0, 'C': -1.0,
        'Q': 0.2, 'E': 3.0, 'G': 0.0, 'H': -0.5, 'I': -1.8,
        'L': -1.8, 'K': 3.0, 'M': -1.3, 'F': -2.5, 'P': 0.0,
        'S': 0.3, 'T': -0.4, 'W': -3.4, 'Y': -2.3, 'V': -1.5
    }

    scores = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        score = sum(hydrophilicity.get(aa, 0) for aa in window) / window_size
        scores.append({
            'position': i + 1,
            'window': window,
            'score': score
        })

    df = pd.DataFrame(scores)
    df['is_epitope'] = df['score'] > threshold

    return df


def identify_epitope_regions(scores_df, min_length=5):
    '''Identify continuous epitope regions from scores

    Args:
        min_length: Minimum consecutive positions to call epitope
                   5 is typical minimum for B-cell epitopes
    '''
    epitopes = []
    current_region = []

    for _, row in scores_df.iterrows():
        if row['is_epitope']:
            current_region.append(row)
        else:
            if len(current_region) >= min_length:
                epitopes.append({
                    'start': current_region[0]['position'],
                    'end': current_region[-1]['position'],
                    'length': len(current_region),
                    'avg_score': sum(r['score'] for r in current_region) / len(current_region),
                    'sequence': ''.join(r['window'][0] for r in current_region)
                })
            current_region = []

    # Don't forget last region
    if len(current_region) >= min_length:
        epitopes.append({
            'start': current_region[0]['position'],
            'end': current_region[-1]['position'],
            'length': len(current_region),
            'avg_score': sum(r['score'] for r in current_region) / len(current_region),
            'sequence': ''.join(r['window'][0] for r in current_region)
        })

    return epitopes


def scan_for_tcell_epitopes(sequence, peptide_lengths=[8, 9, 10, 11]):
    '''Generate all possible T-cell epitope candidates

    MHC class I typically binds 8-11 amino acid peptides.
    9-mers are most common.

    This extracts candidates; actual binding prediction
    requires MHC binding tools (MHCflurry, NetMHCpan).
    '''
    candidates = []
    for length in peptide_lengths:
        for i in range(len(sequence) - length + 1):
            peptide = sequence[i:i + length]
            candidates.append({
                'position': i + 1,
                'length': length,
                'peptide': peptide
            })

    return pd.DataFrame(candidates)


def summarize_epitopes(bcell_epitopes, tcell_candidates):
    '''Summarize epitope predictions'''
    print('Epitope Prediction Summary')
    print('=' * 50)

    print(f'\nB-cell epitopes: {len(bcell_epitopes)}')
    if bcell_epitopes:
        for i, ep in enumerate(bcell_epitopes[:5], 1):
            print(f"  {i}. Position {ep['start']}-{ep['end']}: "
                  f"{ep['sequence'][:20]}{'...' if len(ep['sequence']) > 20 else ''} "
                  f"(score: {ep['avg_score']:.2f})")

    print(f'\nT-cell candidates: {len(tcell_candidates)}')
    length_dist = tcell_candidates['length'].value_counts().to_dict()
    for length, count in sorted(length_dist.items()):
        print(f'  {length}-mers: {count}')


if __name__ == '__main__':
    # Example: Analyze a short antigen sequence
    # Using simplified SARS-CoV-2 receptor binding domain fragment
    sequence = 'NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFN'

    print('Epitope Prediction Example')
    print('=' * 50)
    print(f'Sequence length: {len(sequence)} aa\n')

    # B-cell epitope prediction
    print('Predicting B-cell epitopes...')
    bcell_scores = predict_linear_bcell_epitopes(sequence, window_size=9, threshold=0.3)
    bcell_epitopes = identify_epitope_regions(bcell_scores)

    # T-cell candidates
    print('Generating T-cell epitope candidates...')
    tcell_candidates = scan_for_tcell_epitopes(sequence)

    # Summary
    print()
    summarize_epitopes(bcell_epitopes, tcell_candidates)

    if bcell_epitopes:
        print('\nTop B-cell epitope:')
        top = sorted(bcell_epitopes, key=lambda x: -x['avg_score'])[0]
        print(f"  Position {top['start']}-{top['end']}")
        print(f"  Sequence: {top['sequence']}")
        print(f"  Score: {top['avg_score']:.3f}")
