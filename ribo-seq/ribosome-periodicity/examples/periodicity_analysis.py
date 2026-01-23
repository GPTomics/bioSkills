'''Analyze 3-nucleotide periodicity in Ribo-seq data'''

import pysam
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

def get_reads_around_start(bam_path, gtf_path, upstream=50, downstream=100):
    '''Extract read positions relative to start codons'''
    from plastid import GTF2_TranscriptAssembler

    # Load transcripts
    transcripts = {}
    for tx in GTF2_TranscriptAssembler(gtf_path):
        if tx.cds_start is not None:
            transcripts[tx.get_name()] = tx

    # Count reads by position
    position_counts = defaultdict(int)

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            if read.is_unmapped:
                continue

            # Get 5' position (where ribosome sits)
            if read.is_reverse:
                pos_5p = read.reference_end
            else:
                pos_5p = read.reference_start

            # Map to transcript and get relative position
            # (simplified - full implementation needs transcript mapping)

    return position_counts

def calculate_periodicity_score(frame_counts):
    '''Calculate periodicity score from frame counts

    Score = (frame0 - mean(frame1, frame2)) / total
    Higher score = stronger periodicity

    Good libraries: score > 0.5
    Marginal: 0.3-0.5
    Poor: < 0.3
    '''
    total = sum(frame_counts.values())
    if total == 0:
        return 0

    frame0 = frame_counts.get(0, 0)
    frame1 = frame_counts.get(1, 0)
    frame2 = frame_counts.get(2, 0)

    # Frame 0 should dominate in good Ribo-seq
    frame0_fraction = frame0 / total
    other_mean = (frame1 + frame2) / 2 / total

    return frame0_fraction - other_mean

def analyze_periodicity_by_length(bam_path, gtf_path):
    '''Analyze periodicity for each read length

    Returns optimal P-site offset per length.
    '''
    # Group reads by length
    reads_by_length = defaultdict(list)

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            if not read.is_unmapped:
                reads_by_length[read.query_length].append(read)

    # Analyze each length
    results = []
    for length, reads in sorted(reads_by_length.items()):
        if len(reads) < 1000:  # Need sufficient reads
            continue

        # Count frames at CDS starts
        frame_counts = {0: 0, 1: 0, 2: 0}
        # (simplified frame counting)

        score = calculate_periodicity_score(frame_counts)

        # Typical P-site offset from 5' end
        # Ribosome covers ~28 nt, P-site at position 12-15
        offset = 12 if length <= 29 else 13 if length <= 31 else 14

        results.append({
            'length': length,
            'count': len(reads),
            'periodicity_score': score,
            'psite_offset': offset
        })
        print(f'Length {length}: {len(reads)} reads, periodicity={score:.3f}, offset={offset}')

    return results

def plot_metagene_profile(position_counts, output='metagene.pdf'):
    '''Plot metagene profile around start codons'''
    positions = sorted(position_counts.keys())
    counts = [position_counts[p] for p in positions]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Position plot colored by frame
    colors = ['#e41a1c', '#377eb8', '#4daf4a']  # Red, blue, green for frames 0,1,2
    for frame in range(3):
        frame_pos = [p for p in positions if p % 3 == frame]
        frame_counts = [position_counts[p] for p in frame_pos]
        axes[0].bar(frame_pos, frame_counts, alpha=0.7, color=colors[frame],
                   label=f'Frame {frame}')

    axes[0].axvline(0, color='black', linestyle='--', linewidth=2, label='Start codon')
    axes[0].set_xlabel('Position relative to start codon (nt)')
    axes[0].set_ylabel('Normalized read count')
    axes[0].set_title('Metagene Profile')
    axes[0].legend()

    # Periodicity spectrum (FFT)
    if len(counts) > 10:
        from scipy.fft import fft, fftfreq
        fft_result = np.abs(fft(counts))
        freq = fftfreq(len(counts))
        # Only positive frequencies
        pos_freq = freq[1:len(freq)//2]
        pos_fft = fft_result[1:len(fft_result)//2]
        period = 1 / pos_freq

        axes[1].plot(period, pos_fft, 'b-', linewidth=1.5)
        axes[1].axvline(3, color='red', linestyle='--', linewidth=2, label='3-nt period')
        axes[1].set_xlabel('Period (nucleotides)')
        axes[1].set_ylabel('Power')
        axes[1].set_title('Periodicity Spectrum')
        axes[1].set_xlim(1, 10)
        axes[1].legend()

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    print(f'Saved metagene plot to {output}')

# Example usage
if __name__ == '__main__':
    print('Ribo-seq periodicity analysis')
    print('Load your BAM and GTF files to analyze periodicity')
    print('')
    print('Quality thresholds:')
    print('  Periodicity score > 0.5: Good library')
    print('  Periodicity score 0.3-0.5: Marginal')
    print('  Periodicity score < 0.3: Poor - check digestion')
