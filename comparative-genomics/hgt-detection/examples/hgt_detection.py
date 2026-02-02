'''Horizontal gene transfer detection'''

import numpy as np
import pandas as pd
from collections import Counter


def calculate_gc_content(sequence):
    '''Calculate GC content as a fraction'''
    seq = sequence.upper()
    gc = sum(1 for nt in seq if nt in 'GC')
    total = sum(1 for nt in seq if nt in 'ACGT')
    return gc / total if total > 0 else 0


def genome_gc_stats(cds_sequences):
    '''Calculate genome-wide GC statistics

    Returns mean and std for identifying anomalous genes
    '''
    gc_values = [calculate_gc_content(seq) for seq in cds_sequences.values()]

    return {
        'mean': np.mean(gc_values),
        'std': np.std(gc_values),
        'median': np.median(gc_values),
        'n_genes': len(gc_values)
    }


def identify_gc_anomalies(cds_sequences, z_threshold=2):
    '''Identify genes with anomalous GC content

    Z-score thresholds:
    - |Z| > 2: Moderate anomaly (5% expected by chance)
    - |Z| > 3: Strong anomaly (0.3% expected by chance)

    Note: Ancient HGT may have ameliorated (adapted to host GC)
    so absence of anomaly doesn't rule out HGT
    '''
    stats = genome_gc_stats(cds_sequences)

    anomalies = []
    for gene_id, seq in cds_sequences.items():
        gc = calculate_gc_content(seq)
        z_score = (gc - stats['mean']) / stats['std'] if stats['std'] > 0 else 0

        if abs(z_score) > z_threshold:
            anomalies.append({
                'gene': gene_id,
                'gc': gc,
                'z_score': z_score,
                'direction': 'high' if z_score > 0 else 'low'
            })

    return anomalies, stats


def calculate_codon_usage(cds_sequence):
    '''Calculate relative codon frequencies'''
    if len(cds_sequence) % 3 != 0:
        cds_sequence = cds_sequence[:len(cds_sequence) - len(cds_sequence) % 3]

    codons = [cds_sequence[i:i+3].upper() for i in range(0, len(cds_sequence), 3)]
    codons = [c for c in codons if len(c) == 3 and all(nt in 'ACGT' for nt in c)]

    counts = Counter(codons)
    total = sum(counts.values())

    return {codon: count / total for codon, count in counts.items()} if total > 0 else {}


def genome_codon_usage(cds_sequences):
    '''Calculate genome-wide reference codon usage'''
    all_codons = Counter()

    for seq in cds_sequences.values():
        if len(seq) % 3 != 0:
            seq = seq[:len(seq) - len(seq) % 3]
        codons = [seq[i:i+3].upper() for i in range(0, len(seq), 3)]
        codons = [c for c in codons if len(c) == 3 and all(nt in 'ACGT' for nt in c)]
        all_codons.update(codons)

    total = sum(all_codons.values())
    return {codon: count / total for codon, count in all_codons.items()}


def codon_usage_deviation(gene_codon_usage, genome_codon_usage):
    '''Calculate deviation from genome codon usage

    Uses chi-squared-like metric
    Higher values indicate more deviation from host
    '''
    deviation = 0
    n_codons = 0

    for codon, gene_freq in gene_codon_usage.items():
        genome_freq = genome_codon_usage.get(codon, 0)
        if genome_freq > 0:
            deviation += ((gene_freq - genome_freq) ** 2) / genome_freq
            n_codons += 1

    return deviation / n_codons if n_codons > 0 else 0


def identify_codon_anomalies(cds_sequences, percentile_threshold=95):
    '''Identify genes with anomalous codon usage

    Genes in top 5% of deviation scores are flagged
    '''
    genome_usage = genome_codon_usage(cds_sequences)

    deviations = []
    for gene_id, seq in cds_sequences.items():
        gene_usage = calculate_codon_usage(seq)
        dev = codon_usage_deviation(gene_usage, genome_usage)
        deviations.append({'gene': gene_id, 'codon_deviation': dev})

    df = pd.DataFrame(deviations)
    threshold = np.percentile(df['codon_deviation'], percentile_threshold)
    df['anomalous'] = df['codon_deviation'] > threshold

    return df[df['anomalous']].to_dict('records'), threshold


def identify_genomic_islands(gene_annotations, max_gap=10000, min_genes=3):
    '''Cluster anomalous genes into genomic islands

    Parameters:
    - max_gap: Maximum distance between genes to cluster
    - min_genes: Minimum genes to call an island

    Islands often contain:
    - Mobile elements (integrases, transposases)
    - Pathogenicity/fitness genes
    - Antibiotic resistance genes
    '''
    # Sort by position
    sorted_genes = sorted(gene_annotations, key=lambda x: x.get('start', 0))

    islands = []
    current_island = []

    for gene in sorted_genes:
        if not gene.get('anomalous', False):
            if current_island and len(current_island) >= min_genes:
                islands.append(current_island)
            current_island = []
            continue

        if not current_island:
            current_island = [gene]
        elif gene['start'] - current_island[-1].get('end', 0) <= max_gap:
            current_island.append(gene)
        else:
            if len(current_island) >= min_genes:
                islands.append(current_island)
            current_island = [gene]

    if current_island and len(current_island) >= min_genes:
        islands.append(current_island)

    return islands


def summarize_hgt_analysis(gc_anomalies, gc_stats, codon_anomalies, islands):
    '''Print HGT analysis summary'''
    print('Horizontal Gene Transfer Analysis')
    print('=' * 50)

    print(f'\nGenome GC statistics:')
    print(f"  Mean GC: {gc_stats['mean']*100:.1f}%")
    print(f"  Std GC: {gc_stats['std']*100:.2f}%")
    print(f"  Total genes: {gc_stats['n_genes']}")

    print(f'\nGC content anomalies:')
    print(f'  Genes with |Z| > 2: {len(gc_anomalies)}')
    high_gc = sum(1 for g in gc_anomalies if g['direction'] == 'high')
    low_gc = len(gc_anomalies) - high_gc
    print(f'    High GC: {high_gc}')
    print(f'    Low GC: {low_gc}')

    print(f'\nCodon usage anomalies:')
    print(f'  Genes with unusual codon usage: {len(codon_anomalies)}')

    print(f'\nGenomic islands detected: {len(islands)}')
    for i, island in enumerate(islands, 1):
        n_genes = len(island)
        start = island[0].get('start', 'N/A')
        end = island[-1].get('end', 'N/A')
        print(f'  Island {i}: {n_genes} genes, position {start}-{end}')


if __name__ == '__main__':
    print('HGT Detection Analysis')
    print('=' * 50)

    # Example: Simulated gene data
    np.random.seed(42)

    # Simulate a genome with ~50% GC
    n_genes = 100
    example_genes = {}
    gene_annotations = []

    for i in range(n_genes):
        gene_id = f'gene_{i:04d}'

        # Most genes: normal GC (~50%)
        # Some genes (10%): anomalous GC (simulated HGT)
        if i < 10:  # Simulate HGT cluster
            gc_target = 0.35  # Low GC foreign DNA
        elif i > 90:  # Another HGT cluster
            gc_target = 0.65  # High GC foreign DNA
        else:
            gc_target = 0.50

        # Generate sequence with target GC
        seq = ''.join(
            np.random.choice(['G', 'C'] if np.random.random() < gc_target else ['A', 'T'])
            for _ in range(900)
        )
        example_genes[gene_id] = seq

        gene_annotations.append({
            'gene': gene_id,
            'start': i * 1000,
            'end': i * 1000 + 900,
        })

    # Run analysis
    gc_anomalies, gc_stats = identify_gc_anomalies(example_genes)
    codon_anomalies, codon_threshold = identify_codon_anomalies(example_genes)

    # Mark anomalous genes in annotations
    anomalous_genes = set(g['gene'] for g in gc_anomalies)
    for ann in gene_annotations:
        ann['anomalous'] = ann['gene'] in anomalous_genes

    islands = identify_genomic_islands(gene_annotations)

    summarize_hgt_analysis(gc_anomalies, gc_stats, codon_anomalies, islands)

    print('\n\nTo run on real data:')
    print('1. Load CDS sequences from FASTA')
    print('2. gc_anomalies, gc_stats = identify_gc_anomalies(cds_seqs)')
    print('3. codon_anomalies, _ = identify_codon_anomalies(cds_seqs)')
    print('4. islands = identify_genomic_islands(gene_annotations)')
