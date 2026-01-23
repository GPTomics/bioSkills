'''Detect ribosome stalling sites from Ribo-seq data'''

import numpy as np
import pandas as pd
from collections import defaultdict

def find_pause_sites(occupancy_per_codon, zscore_threshold=3.0):
    '''Find positions with elevated ribosome occupancy

    Pause sites are identified by z-score normalization within each transcript.
    Positions with z-score > threshold are considered pauses.

    Args:
        occupancy_per_codon: Dict of {transcript: [codon occupancies]}
        zscore_threshold: Cutoff for calling pause (default 3.0 = top ~0.1%)

    Returns:
        List of pause site dicts
    '''
    pause_sites = []

    for tx, occupancy in occupancy_per_codon.items():
        occ_array = np.array(occupancy)

        # Need sufficient data
        if len(occ_array) < 20 or occ_array.sum() < 100:
            continue

        # Z-score within transcript
        mean_occ = occ_array.mean()
        std_occ = occ_array.std()

        if std_occ == 0:
            continue

        zscores = (occ_array - mean_occ) / std_occ

        # Find pause positions
        for pos, (occ, zscore) in enumerate(zip(occ_array, zscores)):
            if zscore > zscore_threshold:
                pause_sites.append({
                    'transcript': tx,
                    'codon_position': pos,
                    'occupancy': occ,
                    'zscore': zscore,
                    'mean_transcript': mean_occ
                })

    return pause_sites

def calculate_codon_occupancy(codon_counts):
    '''Calculate average occupancy per codon type

    Higher occupancy for a codon indicates slower translation.
    Often correlates with:
    - Low tRNA abundance (rare codons)
    - Amino acid properties (proline)
    - Context effects
    '''
    codon_means = {}
    for codon, counts in codon_counts.items():
        if len(counts) >= 100:  # Need sufficient observations
            codon_means[codon] = np.mean(counts)

    # Normalize to mean = 1
    overall_mean = np.mean(list(codon_means.values()))
    codon_normalized = {c: v / overall_mean for c, v in codon_means.items()}

    return codon_normalized

def analyze_pause_context(pause_sites, cds_sequences, window=5):
    '''Analyze amino acid context around pause sites

    Known pause motifs:
    - PPP: Polyproline stalls due to ribosome tunnel geometry
    - D/E runs: Negatively charged nascent chains
    - Stop codon context: Nucleotides around stop affect termination
    '''
    from Bio.Seq import Seq

    contexts = []
    for site in pause_sites:
        tx = site['transcript']
        pos = site['codon_position']

        if tx not in cds_sequences:
            continue

        cds = cds_sequences[tx]

        # Extract window around pause
        start_nt = max(0, (pos - window) * 3)
        end_nt = min(len(cds), (pos + window + 1) * 3)

        seq_window = cds[start_nt:end_nt]
        if len(seq_window) % 3 == 0:
            aa_seq = str(Seq(seq_window).translate())
            contexts.append({
                'transcript': tx,
                'position': pos,
                'aa_context': aa_seq,
                'zscore': site['zscore']
            })

    return contexts

def count_motifs(contexts, motif='PP'):
    '''Count occurrences of motif near pause sites

    Enrichment of motif at pause sites vs random indicates
    causal relationship with stalling.
    '''
    pause_with_motif = sum(1 for c in contexts if motif in c['aa_context'])
    total = len(contexts)

    return pause_with_motif, total, pause_with_motif / total if total > 0 else 0

# Example usage
if __name__ == '__main__':
    print('Ribosome stalling analysis')
    print('')
    print('Z-score thresholds for pause calling:')
    print('  > 2.0: Moderate pause (top 2.5%)')
    print('  > 3.0: Strong pause (top 0.1%) [recommended]')
    print('  > 4.0: Very strong pause')
    print('')
    print('Common pause motifs:')
    print('  PPP/PP: Polyproline - ribosome tunnel stalling')
    print('  DDD/EEE: Poly-acidic - electrostatic effects')
    print('  Rare codons: tRNA limitation')
    print('')
    print('Load your Ribo-seq data to detect pause sites')
