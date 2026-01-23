'''miRge3 miRNA quantification and isomiR analysis'''

import subprocess
import pandas as pd
from pathlib import Path

def run_mirge3(samples, output_dir, organism='human', adapter='TGGAATTCTCGGGTGCCAAGG', threads=8):
    '''Run miRge3 annotation pipeline

    Args:
        samples: List of FASTQ file paths
        output_dir: Output directory
        organism: Organism name (human, mouse, etc.)
        adapter: 3' adapter sequence (Illumina TruSeq default)
        threads: Number of threads
    '''
    cmd = [
        'miRge3.0', 'annotate',
        '-s', ','.join(samples),
        '-lib', 'miRge3_libs',
        '-on', organism,
        '-db', 'mirbase',
        '-a', adapter,
        '-o', output_dir,
        '--threads', str(threads),
        '--isomir'  # Enable isomiR detection
    ]

    print(f'Running: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)
    return output_dir

def load_mirge3_results(output_dir):
    '''Load miRge3 count matrices'''
    output_dir = Path(output_dir)

    results = {}

    # Raw counts
    counts_file = output_dir / 'miR.Counts.csv'
    if counts_file.exists():
        results['counts'] = pd.read_csv(counts_file, index_col=0)
        print(f"Loaded {len(results['counts'])} miRNAs")

    # RPM normalized
    rpm_file = output_dir / 'miR.RPM.csv'
    if rpm_file.exists():
        results['rpm'] = pd.read_csv(rpm_file, index_col=0)

    # IsomiR counts
    isomir_file = output_dir / 'isomiR.Counts.csv'
    if isomir_file.exists():
        results['isomirs'] = pd.read_csv(isomir_file, index_col=0)
        print(f"Loaded {len(results['isomirs'])} isomiRs")

    return results

def filter_and_normalize(counts, min_total=10):
    '''Filter low-expressed miRNAs and normalize

    Args:
        counts: Raw count DataFrame
        min_total: Minimum total reads across samples
                   Default 10 - removes very low/noise miRNAs
    '''
    import numpy as np

    # Filter by total expression
    counts_filtered = counts[counts.sum(axis=1) >= min_total]
    print(f'Kept {len(counts_filtered)}/{len(counts)} miRNAs with >= {min_total} total reads')

    # RPM normalization
    total_per_sample = counts_filtered.sum(axis=0)
    rpm = counts_filtered / total_per_sample * 1e6

    # Log2 transform (pseudocount=1 to handle zeros)
    log2_rpm = np.log2(rpm + 1)

    return counts_filtered, rpm, log2_rpm

def summarize_isomirs(isomir_counts):
    '''Summarize isomiR diversity per canonical miRNA

    Returns number of isomiRs and dominant variant per miRNA.
    High isomiR diversity may indicate active modification.
    '''
    # Extract canonical miRNA name from isomiR ID
    # Format: hsa-miR-21-5p_variant
    canonical = isomir_counts.index.str.extract(r'(hsa-\w+-\d+[a-z]*-[35]p)')[0]
    isomir_counts = isomir_counts.copy()
    isomir_counts['canonical'] = canonical

    summary = []
    for mirna, group in isomir_counts.groupby('canonical'):
        total = group.drop('canonical', axis=1).sum().sum()
        n_variants = len(group)
        summary.append({
            'miRNA': mirna,
            'total_reads': total,
            'n_isomirs': n_variants,
            'isomir_diversity': n_variants / total if total > 0 else 0
        })

    return pd.DataFrame(summary).sort_values('total_reads', ascending=False)

# Example usage
if __name__ == '__main__':
    # Load example results
    output_dir = 'mirge3_output'

    # If running fresh analysis:
    # samples = ['sample1.fastq.gz', 'sample2.fastq.gz']
    # run_mirge3(samples, output_dir)

    # Load results
    results = load_mirge3_results(output_dir)

    if 'counts' in results:
        counts, rpm, log2_rpm = filter_and_normalize(results['counts'])
        print('\nTop 10 expressed miRNAs:')
        print(rpm.mean(axis=1).sort_values(ascending=False).head(10))

    if 'isomirs' in results:
        isomir_summary = summarize_isomirs(results['isomirs'])
        print('\nIsomiR diversity summary:')
        print(isomir_summary.head(10))
