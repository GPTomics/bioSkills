'''Multi-locus sequence typing for bacterial isolates'''

import subprocess
import pandas as pd
import numpy as np
from pathlib import Path


def run_mlst(fasta_files, scheme=None):
    '''Run MLST typing on genome assemblies

    MLST identifies isolates using ~7 housekeeping gene alleles.
    Returns sequence type (ST) which can be used for:
    - Outbreak investigation
    - Surveillance tracking
    - Epidemiological studies

    Args:
        fasta_files: List of FASTA file paths
        scheme: MLST scheme (auto-detected if None)

    Returns:
        DataFrame with typing results
    '''
    cmd = ['mlst'] + [str(f) for f in fasta_files]
    if scheme:
        cmd.extend(['--scheme', scheme])

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f'Error: {result.stderr}')
        return None

    # Parse output
    lines = result.stdout.strip().split('\n')
    rows = []
    for line in lines:
        parts = line.split('\t')
        rows.append({
            'file': parts[0],
            'scheme': parts[1],
            'ST': parts[2],
            'allele_profile': '\t'.join(parts[3:])
        })

    return pd.DataFrame(rows)


def calculate_allelic_distance(profile1, profile2):
    '''Calculate number of allele differences between two profiles

    Distance interpretation (general guidelines):
    - 0 differences: Identical ST
    - 1-2 differences: Single locus variant (SLV), closely related
    - 3-4 differences: Double locus variant (DLV)
    - 5+ differences: Different clonal complexes
    '''
    alleles1 = profile1.split('\t')
    alleles2 = profile2.split('\t')

    differences = 0
    for a1, a2 in zip(alleles1, alleles2):
        # Skip if either is missing/novel
        if '?' in a1 or '?' in a2 or '~' in a1 or '~' in a2:
            continue
        # Extract allele number
        num1 = a1.split('(')[1].rstrip(')') if '(' in a1 else a1
        num2 = a2.split('(')[1].rstrip(')') if '(' in a2 else a2
        if num1 != num2:
            differences += 1

    return differences


def identify_clusters(typing_results, max_distance=2):
    '''Identify clonal clusters from MLST results

    Args:
        max_distance: Maximum allele differences for same cluster
                     2 = SLVs grouped together (common for outbreaks)
    '''
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    n = len(typing_results)
    distance_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i+1, n):
            dist = calculate_allelic_distance(
                typing_results.iloc[i]['allele_profile'],
                typing_results.iloc[j]['allele_profile']
            )
            distance_matrix[i, j] = distance_matrix[j, i] = dist

    # Cluster
    condensed = squareform(distance_matrix)
    Z = linkage(condensed, method='single')
    clusters = fcluster(Z, t=max_distance, criterion='distance')

    typing_results['cluster'] = clusters
    return typing_results, distance_matrix


def summarize_typing(results):
    '''Summarize MLST typing results'''
    print('MLST Typing Summary')
    print('=' * 50)

    # Count STs
    st_counts = results['ST'].value_counts()
    print(f"\nTotal isolates: {len(results)}")
    print(f"Unique STs: {len(st_counts)}")

    print('\nST distribution:')
    for st, count in st_counts.head(10).items():
        print(f'  ST{st}: {count} isolates')

    # Check for novel types
    novel = results[results['ST'].str.contains('-|~', na=False)]
    if len(novel) > 0:
        print(f'\nNovel/incomplete profiles: {len(novel)}')


if __name__ == '__main__':
    print('MLST Typing Example')
    print('=' * 50)

    # Example: Would run on actual FASTA files
    # results = run_mlst(['isolate1.fasta', 'isolate2.fasta'])

    # Simulated results for demonstration
    demo_results = pd.DataFrame({
        'file': ['isolate1.fasta', 'isolate2.fasta', 'isolate3.fasta',
                 'isolate4.fasta', 'isolate5.fasta'],
        'scheme': ['ecoli'] * 5,
        'ST': ['131', '131', '131', '73', '73'],
        'allele_profile': [
            'adk(53)\tfumC(40)\tgyrB(47)\ticd(13)\tmdh(36)\tpurA(28)\trecA(29)',
            'adk(53)\tfumC(40)\tgyrB(47)\ticd(13)\tmdh(36)\tpurA(28)\trecA(29)',
            'adk(53)\tfumC(40)\tgyrB(47)\ticd(13)\tmdh(36)\tpurA(28)\trecA(7)',  # SLV
            'adk(6)\tfumC(4)\tgyrB(4)\ticd(16)\tmdh(24)\tpurA(8)\trecA(14)',
            'adk(6)\tfumC(4)\tgyrB(4)\ticd(16)\tmdh(24)\tpurA(8)\trecA(14)'
        ]
    })

    summarize_typing(demo_results)

    # Identify clusters
    print('\nCluster Analysis (max 2 allele differences):')
    clustered, distances = identify_clusters(demo_results.copy(), max_distance=2)
    for cluster_id in sorted(clustered['cluster'].unique()):
        members = clustered[clustered['cluster'] == cluster_id]['file'].tolist()
        print(f'  Cluster {cluster_id}: {members}')
