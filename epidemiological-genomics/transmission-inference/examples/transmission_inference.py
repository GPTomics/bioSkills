'''Infer pathogen transmission networks'''

import pandas as pd
import numpy as np
from collections import Counter


def infer_transmission_pairs(snp_matrix, collection_dates, max_snps=5, max_days=30):
    '''Infer likely transmission pairs from genomic and temporal data

    Criteria for transmission (A -> B):
    1. A collected before B (temporal directionality)
    2. SNP distance consistent with direct transmission
    3. Time difference compatible with generation time

    Args:
        snp_matrix: Pairwise SNP distances (DataFrame or array)
        collection_dates: Dict mapping sample -> decimal year date
        max_snps: Maximum SNPs for direct transmission (pathogen-specific)
                 - SARS-CoV-2: 2-3 SNPs
                 - TB: 5-12 SNPs
                 - E. coli outbreak: 0-5 SNPs
        max_days: Maximum days between collection for direct link
                 Related to serial interval of pathogen

    Returns:
        DataFrame of transmission pairs with confidence scores
    '''
    samples = list(collection_dates.keys())
    n = len(samples)

    pairs = []
    for i, infector in enumerate(samples):
        for j, infectee in enumerate(samples):
            if i == j:
                continue

            # Temporal check: infector must be earlier
            days_diff = (collection_dates[infectee] - collection_dates[infector]) * 365
            if days_diff <= 0 or days_diff > max_days:
                continue

            # Genomic check
            snp_dist = snp_matrix[i, j] if isinstance(snp_matrix, np.ndarray) else snp_matrix.iloc[i, j]
            if snp_dist > max_snps:
                continue

            # Score confidence
            # Lower SNPs + appropriate timing = higher confidence
            if snp_dist <= 1 and days_diff <= 14:
                confidence = 'high'
            elif snp_dist <= 3:
                confidence = 'medium'
            else:
                confidence = 'low'

            pairs.append({
                'infector': infector,
                'infectee': infectee,
                'snp_distance': snp_dist,
                'days_between': days_diff,
                'confidence': confidence
            })

    return pd.DataFrame(pairs)


def identify_superspreaders(transmission_df, threshold=3):
    '''Identify superspreading individuals

    Superspreader definition:
    - 80/20 rule: ~20% of cases cause ~80% of transmission
    - Common threshold: >3 secondary cases

    Superspreaders are important for:
    - Outbreak control (targeted interventions)
    - Understanding transmission dynamics
    - R0 estimation (overdispersion parameter k)
    '''
    infector_counts = Counter(transmission_df['infector'])

    superspreaders = []
    for infector, count in infector_counts.most_common():
        if count >= threshold:
            superspreaders.append({
                'case': infector,
                'secondary_infections': count
            })

    if superspreaders:
        total = sum(infector_counts.values())
        ss_total = sum(s['secondary_infections'] for s in superspreaders)
        print(f'Superspreaders account for {ss_total/total:.1%} of transmissions')

    return superspreaders


def build_transmission_chain(pairs_df, start_case=None):
    '''Reconstruct transmission chain from pairs

    Returns ordered list from source to terminal cases
    '''
    # Build graph
    graph = {}
    all_cases = set(pairs_df['infector']) | set(pairs_df['infectee'])

    for case in all_cases:
        graph[case] = {
            'infected_by': None,
            'infected': []
        }

    for _, row in pairs_df.iterrows():
        if graph[row['infectee']]['infected_by'] is None:
            graph[row['infectee']]['infected_by'] = row['infector']
            graph[row['infector']]['infected'].append(row['infectee'])

    # Find root (index case)
    if start_case is None:
        roots = [c for c, data in graph.items() if data['infected_by'] is None]
        if len(roots) == 1:
            start_case = roots[0]
        else:
            print(f'Multiple potential index cases: {roots}')
            start_case = roots[0]

    return graph, start_case


def print_transmission_tree(graph, case, indent=0):
    '''Print transmission tree in text format'''
    prefix = '  ' * indent + ('└─ ' if indent > 0 else '')
    print(f'{prefix}{case}')
    for child in graph[case]['infected']:
        print_transmission_tree(graph, child, indent + 1)


if __name__ == '__main__':
    print('Transmission Inference Example')
    print('=' * 50)

    # Simulated outbreak data
    np.random.seed(42)
    samples = [f'Case_{i}' for i in range(10)]

    # Collection dates (decimal year)
    dates = {
        'Case_0': 2020.10,  # Index case
        'Case_1': 2020.12,
        'Case_2': 2020.13,
        'Case_3': 2020.15,
        'Case_4': 2020.16,
        'Case_5': 2020.18,
        'Case_6': 2020.20,
        'Case_7': 2020.22,
        'Case_8': 2020.25,
        'Case_9': 2020.28,
    }

    # SNP distance matrix (simulated transmission chain)
    # Index -> 1,2 -> 3,4 -> 5,6 -> etc.
    snp_matrix = np.random.randint(0, 3, size=(10, 10))
    np.fill_diagonal(snp_matrix, 0)
    snp_matrix = (snp_matrix + snp_matrix.T) // 2  # Make symmetric

    # Infer transmission
    pairs = infer_transmission_pairs(snp_matrix, dates, max_snps=3, max_days=50)

    print(f'\nInferred {len(pairs)} potential transmission pairs:')
    print(pairs[['infector', 'infectee', 'snp_distance', 'confidence']].head(10).to_string(index=False))

    # Identify superspreaders
    print('\nSuperspreader Analysis:')
    ss = identify_superspreaders(pairs, threshold=2)
    for s in ss:
        print(f"  {s['case']}: {s['secondary_infections']} secondary infections")

    # Build and display tree
    print('\nTransmission Chain:')
    graph, root = build_transmission_chain(pairs)
    print_transmission_tree(graph, root)
