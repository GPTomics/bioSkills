'''TCR-epitope specificity prediction'''

import pandas as pd
from difflib import SequenceMatcher


def parse_tcr_sequences(tcr_data):
    '''Parse and validate TCR sequences

    CDR3 (Complementarity Determining Region 3) is the primary
    determinant of TCR antigen specificity.

    Expected format:
    - cdr3_beta: Required, most informative
    - cdr3_alpha: Optional, improves specificity
    - v_beta, j_beta: V/J gene usage (optional)
    '''
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')

    def validate_cdr3(seq):
        if pd.isna(seq) or not seq:
            return False
        return all(aa in valid_aa for aa in seq.upper())

    df = pd.DataFrame(tcr_data)
    df['valid'] = df['cdr3_beta'].apply(validate_cdr3)

    invalid_count = (~df['valid']).sum()
    if invalid_count > 0:
        print(f'Warning: {invalid_count} invalid CDR3 sequences removed')

    return df[df['valid']]


def levenshtein_distance(s1, s2):
    '''Calculate edit distance between two sequences'''
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]


def match_tcr_to_database(query_cdr3, database_df, max_distance=1):
    '''Match TCR to known epitopes in database

    VDJdb contains curated TCR-epitope pairs from:
    - Published literature
    - IEDB
    - 10x Genomics data

    Args:
        max_distance: Maximum edit distance for fuzzy matching
                     0 = exact only
                     1 = allow 1 mismatch (recommended)
                     2 = more permissive
    '''
    matches = []

    for _, row in database_df.iterrows():
        db_cdr3 = row['cdr3']
        distance = levenshtein_distance(query_cdr3, db_cdr3)

        if distance <= max_distance:
            matches.append({
                'query': query_cdr3,
                'db_cdr3': db_cdr3,
                'distance': distance,
                'epitope': row.get('epitope', 'Unknown'),
                'antigen': row.get('antigen', 'Unknown'),
                'species': row.get('species', 'Unknown'),
                'mhc': row.get('mhc', 'Unknown')
            })

    return pd.DataFrame(matches)


def cluster_tcrs(cdr3_sequences, max_distance=3):
    '''Cluster TCRs by CDR3 similarity

    TCRs within 1-3 edit distance often recognize the same epitope.
    This allows grouping for specificity prediction.
    '''
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    import numpy as np

    n = len(cdr3_sequences)
    distances = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            d = levenshtein_distance(cdr3_sequences[i], cdr3_sequences[j])
            distances[i, j] = distances[j, i] = d

    if n > 1:
        condensed = squareform(distances)
        Z = linkage(condensed, method='average')
        clusters = fcluster(Z, t=max_distance, criterion='distance')
    else:
        clusters = [1]

    return dict(zip(cdr3_sequences, clusters))


def simple_binding_score(cdr3, epitope):
    '''Simple heuristic binding score

    For accurate predictions, use ERGO-II or similar deep learning models.

    This simple score considers:
    - Length compatibility
    - Charge complementarity
    '''
    # Length score (similar CDR3 lengths within specificity groups)
    length_score = max(0, 1 - abs(len(cdr3) - 13) / 5)  # Optimal ~13aa

    # Charge complementarity
    positive = set('RKH')
    negative = set('DE')

    cdr3_charge = sum(1 if aa in positive else -1 if aa in negative else 0 for aa in cdr3)
    epitope_charge = sum(1 if aa in positive else -1 if aa in negative else 0 for aa in epitope)

    charge_score = 0.5 - (cdr3_charge * epitope_charge) / 10

    combined = (length_score + max(0, min(1, charge_score))) / 2

    return {
        'cdr3': cdr3,
        'epitope': epitope,
        'length_score': length_score,
        'charge_score': charge_score,
        'combined_score': combined
    }


if __name__ == '__main__':
    print('TCR-Epitope Binding Analysis')
    print('=' * 50)

    # Example TCR sequences
    tcrs = [
        {'cdr3_beta': 'CASSLAPGTTNEKLFF'},
        {'cdr3_beta': 'CASSLGQAYEQYF'},
        {'cdr3_beta': 'CASSIRSSYEQYF'},
        {'cdr3_beta': 'CASSLGQAYEQYF'},  # Duplicate
        {'cdr3_beta': 'CASSLAPGTTNEKLFF'},  # Similar to first
    ]

    # Simulated VDJdb entries
    vdjdb = pd.DataFrame({
        'cdr3': ['CASSLAPGTTNEKLFF', 'CASSLGQAYEQYF', 'CASSIRSSYEQYF'],
        'epitope': ['GILGFVFTL', 'NLVPMVATV', 'GLCTLVAML'],
        'antigen': ['Influenza M1', 'CMV pp65', 'EBV BMLF1'],
        'species': ['Influenza A', 'CMV', 'EBV'],
        'mhc': ['HLA-A*02:01', 'HLA-A*02:01', 'HLA-A*02:01']
    })

    # Parse TCRs
    tcr_df = parse_tcr_sequences(tcrs)
    print(f'\nValid TCRs: {len(tcr_df)}')

    # Match to database
    print('\nDatabase matching:')
    unique_cdr3 = tcr_df['cdr3_beta'].unique()
    for cdr3 in unique_cdr3:
        matches = match_tcr_to_database(cdr3, vdjdb, max_distance=1)
        if len(matches) > 0:
            for _, m in matches.iterrows():
                print(f"  {cdr3} -> {m['epitope']} ({m['antigen']}) [distance={m['distance']}]")

    # Cluster TCRs
    print('\nTCR clusters:')
    clusters = cluster_tcrs(list(unique_cdr3), max_distance=3)
    cluster_groups = {}
    for cdr3, cluster_id in clusters.items():
        if cluster_id not in cluster_groups:
            cluster_groups[cluster_id] = []
        cluster_groups[cluster_id].append(cdr3)

    for cluster_id, members in cluster_groups.items():
        print(f'  Cluster {cluster_id}: {len(members)} TCRs')
        for m in members[:3]:
            print(f'    {m}')

    # Simple binding score example
    print('\nBinding score (heuristic):')
    score = simple_binding_score('CASSLAPGTTNEKLFF', 'GILGFVFTL')
    print(f"  CDR3: {score['cdr3']}")
    print(f"  Epitope: {score['epitope']}")
    print(f"  Combined score: {score['combined_score']:.3f}")
