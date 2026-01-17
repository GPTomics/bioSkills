# ATAC-seq Quality Control Usage Guide

## Overview

ATAC-seq QC assesses both technical quality (library complexity, contamination) and biological signal (chromatin accessibility patterns). Good QC is essential before differential analysis.

## Key QC Metrics

| Metric | What It Measures | Target |
|--------|------------------|--------|
| TSS Enrichment | Signal at promoters | >7 (ENCODE standard) |
| FRiP | Signal in peaks vs background | >0.2 |
| Fragment Size | Nucleosome periodicity | Clear pattern |
| MT Fraction | Contamination | <20% |
| NRF | Library complexity | >0.9 |

## Complete QC Pipeline

```bash
#!/bin/bash
set -euo pipefail

BAM=$1
PEAKS=$2
SAMPLE=$(basename $BAM .bam)
OUTDIR=qc_results/$SAMPLE

mkdir -p $OUTDIR

echo "=== ATAC-seq QC Report for $SAMPLE ===" > $OUTDIR/qc_report.txt
echo "" >> $OUTDIR/qc_report.txt

# 1. Basic alignment stats
echo "1. Alignment Statistics" >> $OUTDIR/qc_report.txt
samtools flagstat $BAM > $OUTDIR/flagstat.txt
cat $OUTDIR/flagstat.txt >> $OUTDIR/qc_report.txt
echo "" >> $OUTDIR/qc_report.txt

# 2. Mitochondrial fraction
echo "2. Mitochondrial Contamination" >> $OUTDIR/qc_report.txt
MT=$(samtools view -c $BAM chrM)
TOTAL=$(samtools view -c -F 4 $BAM)
MT_FRAC=$(echo "scale=4; $MT / $TOTAL" | bc)
echo "Mitochondrial reads: $MT" >> $OUTDIR/qc_report.txt
echo "Total mapped reads: $TOTAL" >> $OUTDIR/qc_report.txt
echo "MT Fraction: $MT_FRAC" >> $OUTDIR/qc_report.txt
echo "" >> $OUTDIR/qc_report.txt

# 3. FRiP
echo "3. FRiP (Fraction of Reads in Peaks)" >> $OUTDIR/qc_report.txt
READS_IN_PEAKS=$(bedtools intersect -a $BAM -b $PEAKS -u | samtools view -c)
FRIP=$(echo "scale=4; $READS_IN_PEAKS / $TOTAL" | bc)
echo "Reads in peaks: $READS_IN_PEAKS" >> $OUTDIR/qc_report.txt
echo "FRiP: $FRIP" >> $OUTDIR/qc_report.txt
echo "" >> $OUTDIR/qc_report.txt

# 4. Peak statistics
echo "4. Peak Statistics" >> $OUTDIR/qc_report.txt
PEAK_COUNT=$(wc -l < $PEAKS)
echo "Total peaks: $PEAK_COUNT" >> $OUTDIR/qc_report.txt
MEDIAN_WIDTH=$(awk '{print $3-$2}' $PEAKS | sort -n | awk '{a[NR]=$1} END {print a[int(NR/2)]}')
echo "Median peak width: $MEDIAN_WIDTH" >> $OUTDIR/qc_report.txt
echo "" >> $OUTDIR/qc_report.txt

# 5. Fragment size distribution
echo "Generating fragment size distribution..."
samtools view -f 66 $BAM | awk '{print sqrt($9^2)}' | sort -n | uniq -c | \
    awk '{print $2"\t"$1}' > $OUTDIR/fragment_sizes.txt

# 6. Insert size histogram with Picard
java -jar $PICARD CollectInsertSizeMetrics \
    I=$BAM \
    O=$OUTDIR/insert_metrics.txt \
    H=$OUTDIR/insert_histogram.pdf \
    M=0.5 2>/dev/null

echo "QC complete. Report: $OUTDIR/qc_report.txt"
```

## Interpreting Fragment Size Distribution

Expected ATAC-seq pattern shows nucleosome periodicity:

```
Fragment Size    Peak Type
<100 bp         Nucleosome-free regions (NFR)
180-247 bp      Mono-nucleosome
315-473 bp      Di-nucleosome
558-615 bp      Tri-nucleosome
```

Plot and check:
```python
import pandas as pd
import matplotlib.pyplot as plt

frag_sizes = pd.read_csv('fragment_sizes.txt', sep='\t', names=['size', 'count'])
frag_sizes = frag_sizes[frag_sizes['size'] <= 1000]

plt.figure(figsize=(10, 4))
plt.fill_between(frag_sizes['size'], frag_sizes['count'], alpha=0.5)
plt.axvline(100, color='red', linestyle='--', label='NFR cutoff')
plt.axvline(180, color='green', linestyle='--', label='Mono-nuc')
plt.xlabel('Fragment Size (bp)')
plt.ylabel('Count')
plt.title('ATAC-seq Fragment Size Distribution')
plt.legend()
plt.savefig('fragment_distribution.png', dpi=150)
```

Good: Clear peaks at ~50bp (NFR), ~200bp (mono), ~400bp (di)
Bad: Single peak or no periodicity

## TSS Enrichment Calculation

```python
import numpy as np
import pyBigWig
import pandas as pd

def tss_enrichment_score(bigwig_file, tss_bed, window=2000):
    '''Calculate ENCODE-style TSS enrichment score.'''
    bw = pyBigWig.open(bigwig_file)

    profiles = []
    for line in open(tss_bed):
        fields = line.strip().split('\t')
        chrom = fields[0]
        tss = int(fields[1])
        strand = fields[5] if len(fields) > 5 else '+'

        start = max(0, tss - window)
        end = tss + window

        try:
            vals = bw.values(chrom, start, end)
            if vals is None:
                continue
            vals = np.array(vals)
            if strand == '-':
                vals = vals[::-1]
            profiles.append(vals)
        except:
            continue

    if not profiles:
        return None, None

    avg_profile = np.nanmean(profiles, axis=0)

    # Background: edges of profile
    background = np.nanmean(np.concatenate([avg_profile[:100], avg_profile[-100:]]))

    # Signal: center
    center_start = window - 50
    center_end = window + 50
    signal = np.nanmean(avg_profile[center_start:center_end])

    tsse = signal / background if background > 0 else 0

    return tsse, avg_profile

# Generate TSS BED from GTF
# awk '$3=="gene" {print $1"\t"$4"\t"$4+1"\t.\t.\t"$7}' genes.gtf > tss.bed

score, profile = tss_enrichment_score('sample.bw', 'tss.bed')
print(f'TSS Enrichment Score: {score:.1f}')
```

ENCODE threshold: TSSE > 7

## Library Complexity

```python
from collections import Counter
import pysam

def library_complexity(bam_file, sample_size=1000000):
    '''Calculate NRF, PBC1, PBC2.'''
    bam = pysam.AlignmentFile(bam_file, 'rb')

    positions = Counter()
    n_reads = 0

    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.is_read2:  # Only count read1 for pairs
            continue

        n_reads += 1
        pos = (read.reference_name, read.reference_start, read.is_reverse)
        positions[pos] += 1

        if n_reads >= sample_size:
            break

    distinct = len(positions)
    m1 = sum(1 for v in positions.values() if v == 1)
    m2 = sum(1 for v in positions.values() if v == 2)

    nrf = distinct / n_reads if n_reads > 0 else 0
    pbc1 = m1 / distinct if distinct > 0 else 0
    pbc2 = m1 / m2 if m2 > 0 else float('inf')

    return {
        'NRF': nrf,      # Non-redundant fraction
        'PBC1': pbc1,    # PCR bottleneck coefficient 1
        'PBC2': pbc2,    # PCR bottleneck coefficient 2
        'distinct_positions': distinct,
        'total_reads': n_reads
    }

complexity = library_complexity('sample.bam')
print(f"NRF: {complexity['NRF']:.3f}")
print(f"PBC1: {complexity['PBC1']:.3f}")
print(f"PBC2: {complexity['PBC2']:.1f}")
```

## MultiQC Integration

```bash
# Run various QC tools, then aggregate with MultiQC
fastqc sample.fastq.gz
picard CollectInsertSizeMetrics ...
picard CollectAlignmentSummaryMetrics ...

# Aggregate all QC
multiqc . -o multiqc_report
```

## Troubleshooting

### Low TSS Enrichment
- Check alignment parameters
- Verify correct genome annotation
- May indicate poor chromatin preparation

### High MT Fraction
- Filter MT reads before analysis
- Consider re-extracting nuclei

### Poor Fragment Distribution
- May indicate over-transposition
- Could be degraded chromatin
