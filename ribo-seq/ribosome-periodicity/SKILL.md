---
name: bio-ribo-seq-ribosome-periodicity
description: Validate Ribo-seq data quality by checking 3-nucleotide periodicity and calculating P-site offsets. Use when assessing library quality or determining read offsets for downstream analysis.
tool_type: python
primary_tool: Plastid
---

# Ribosome Periodicity Analysis

## 3-Nucleotide Periodicity

Ribosomes move 3 nucleotides per codon. Good Ribo-seq data shows strong periodicity:

```python
from plastid import BAMGenomeArray, FivePrimeMapFactory, GenomicSegment
import numpy as np
import matplotlib.pyplot as plt

# Load aligned reads
alignments = BAMGenomeArray('riboseq.bam', mapping=FivePrimeMapFactory())

# Get metagene around start codons
# Expect strong 3-nt periodicity
```

## Calculate P-site Offset

```python
from plastid import metagene_analysis

# The P-site offset varies by read length
# Typically 12-15 nt from 5' end for 28-30 nt reads

def determine_psite_offset(bam_path, annotation_file):
    '''Determine optimal P-site offset from metagene analysis'''
    from plastid import GTF2_TranscriptAssembler, BAMGenomeArray

    # Load annotations
    transcripts = list(GTF2_TranscriptAssembler(annotation_file))

    # Load reads
    alignments = BAMGenomeArray(bam_path, mapping=FivePrimeMapFactory())

    # Metagene around start codons
    # Peak should align with start codon position
    metagene_data = metagene_analysis(
        transcripts,
        alignments,
        upstream=50,
        downstream=100
    )

    return metagene_data
```

## Metagene Plots

```python
def plot_metagene(metagene_data, offset=12):
    '''Plot metagene profile around start codon'''
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Frame 0, 1, 2 around start codon
    positions = np.arange(-50, 100)

    # Plot by frame
    for frame in range(3):
        frame_positions = positions[positions % 3 == frame]
        counts = metagene_data[positions % 3 == frame]
        axes[0].bar(frame_positions, counts, alpha=0.7, label=f'Frame {frame}')

    axes[0].set_xlabel('Position relative to start codon')
    axes[0].set_ylabel('Normalized counts')
    axes[0].legend()
    axes[0].axvline(0, color='red', linestyle='--', label='Start')

    # Periodicity
    from scipy.fft import fft
    fft_result = np.abs(fft(metagene_data))
    freq = np.fft.fftfreq(len(metagene_data))

    axes[1].plot(1/freq[1:len(freq)//2], fft_result[1:len(freq)//2])
    axes[1].set_xlabel('Period (nt)')
    axes[1].set_ylabel('Power')
    axes[1].axvline(3, color='red', linestyle='--')

    plt.tight_layout()
    plt.savefig('periodicity.pdf')
```

## Assess by Read Length

```python
def periodicity_by_length(bam_path, annotation_file):
    '''Calculate periodicity score for each read length'''
    import pysam

    # Group reads by length
    reads_by_length = {}
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            if not read.is_unmapped:
                length = read.query_length
                if length not in reads_by_length:
                    reads_by_length[length] = []
                reads_by_length[length].append(read)

    # Calculate periodicity for each length
    # Good lengths show strong 3-nt periodicity
    results = {}
    for length, reads in reads_by_length.items():
        if len(reads) > 1000:  # Need sufficient reads
            periodicity = calculate_periodicity(reads, annotation_file)
            results[length] = periodicity

    return results
```

## P-site Offset Table

Common P-site offsets by read length (5' end mapping):

| Read Length | P-site Offset |
|-------------|---------------|
| 28 nt | 12 |
| 29 nt | 12 |
| 30 nt | 13 |
| 31 nt | 13 |
| 32 nt | 14 |

## Validate with RiboCode

```bash
# RiboCode includes periodicity analysis
RiboCode_onestep \
    -g annotation.gtf \
    -r riboseq.bam \
    -f genome.fa \
    -o output_dir

# Check output for periodicity plots
```

## Related Skills

- riboseq-preprocessing - Generate aligned BAM
- orf-detection - Uses P-site offsets
- translation-efficiency - Requires proper positioning
