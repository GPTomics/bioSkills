# Alignment Validation - Usage Guide

## Overview
Validate alignment quality with insert size distribution, proper pairing rates, GC bias, strand balance, and other post-alignment metrics before downstream analysis.

## Prerequisites
```bash
# samtools and Picard
conda install -c bioconda samtools picard

# Python packages
pip install pysam numpy matplotlib
```

## Quick Start
Tell your AI agent what you want to do:
- "Validate the quality of my aligned BAM file"
- "Check insert size distribution for sample.bam"
- "Generate a comprehensive QC report for my alignment"
- "Check for GC bias in my sequencing data"

## Example Prompts

### Basic Validation
> "Run alignment validation on sample.bam and tell me if it passes QC"

> "Check the mapping rate and proper pairing percentage for my BAM file"

> "Generate a flagstat report and interpret the results"

### Insert Size Analysis
> "Calculate insert size distribution for sample.bam and plot a histogram"

> "Check if the insert size matches my library prep protocol (expected 350bp)"

> "Compare insert sizes across multiple samples to check for consistency"

### GC Bias Detection
> "Check for GC bias in my WGS data using Picard"

> "Generate a GC bias plot and tell me if correction is needed"

> "Compare GC coverage across my samples to identify problematic ones"

### Strand Balance
> "Calculate the forward/reverse strand ratio for sample.bam"

> "Check strand balance per chromosome to identify bias"

> "Validate strand balance meets expected 0.5 ratio"

### Comprehensive QC
> "Run a complete alignment validation pipeline with all metrics"

> "Generate QC report with mapping rate, pairing, insert size, and strand balance"

> "Identify any quality issues in my alignment before variant calling"

## What the Agent Will Do

1. Load the BAM file and verify it is indexed
2. Calculate mapping statistics using samtools flagstat
3. Extract insert size distribution from properly paired reads
4. Compute strand balance (forward/reverse ratio)
5. Check mapping quality distribution
6. Generate per-chromosome coverage statistics
7. Optionally run Picard metrics for GC bias and alignment summary
8. Compare metrics against quality thresholds
9. Generate a summary report with pass/fail status
10. Recommend corrective actions for any failing metrics

## Key Metrics and Thresholds

### Mapping Rate
Percentage of reads successfully aligned to the reference.

| Status | WGS | Exome | RNA-seq |
|--------|-----|-------|---------|
| Good | > 95% | > 95% | > 80% |
| Warning | 90-95% | 90-95% | 70-80% |
| Fail | < 90% | < 90% | < 70% |

### Proper Pairing
Percentage of paired reads aligned with correct orientation and insert size.

| Status | Threshold |
|--------|-----------|
| Good | > 90% |
| Warning | 80-90% |
| Fail | < 80% |

### Insert Size
Fragment size distribution should match library preparation protocol.

| Library Type | Expected Size |
|--------------|---------------|
| Standard WGS | 300-500 bp |
| PCR-free WGS | 350-550 bp |
| RNA-seq | 150-300 bp |
| ChIP-seq | 150-300 bp |
| ATAC-seq | Multimodal (nucleosome pattern) |

### Strand Balance
Forward/reverse strand ratio should be approximately 0.5 (balanced).

| Status | F/R Ratio |
|--------|-----------|
| Good | 0.48-0.52 |
| Warning | 0.45-0.55 |
| Fail | Outside 0.45-0.55 |

### Mapping Quality
MAPQ indicates confidence in alignment position.

| MAPQ | Meaning |
|------|---------|
| 0 | Multi-mapper or uncertain |
| 1-10 | Low confidence |
| 20-30 | Moderate confidence |
| 40+ | High confidence |
| 60 | Unique mapping (BWA) |

## Common Validation Commands

### Quick Flagstat Check
```bash
samtools flagstat input.bam
```

### Calculate Metrics
```bash
# Mapping rate
mapped=$(samtools view -c -F 4 input.bam)
total=$(samtools view -c input.bam)
rate=$(echo "scale=2; $mapped / $total * 100" | bc)
echo "Mapping rate: ${rate}%"

# Proper pairing rate
proper=$(samtools view -c -f 2 input.bam)
pair_rate=$(echo "scale=2; $proper / $mapped * 100" | bc)
echo "Proper pairing: ${pair_rate}%"

# Strand balance
forward=$(samtools view -c -F 16 input.bam)
reverse=$(samtools view -c -f 16 input.bam)
ratio=$(echo "scale=4; $forward / $reverse" | bc)
echo "F/R ratio: $ratio"
```

### Insert Size with Picard
```bash
java -jar picard.jar CollectInsertSizeMetrics \
    I=input.bam \
    O=insert_metrics.txt \
    H=insert_histogram.pdf
```

### GC Bias with Picard
```bash
java -jar picard.jar CollectGcBiasMetrics \
    I=input.bam \
    O=gc_bias_metrics.txt \
    CHART=gc_bias_chart.pdf \
    S=gc_summary.txt \
    R=reference.fa
```

### Comprehensive Stats
```bash
samtools stats input.bam > stats.txt
plot-bamstats -p stats_plots/ stats.txt
```

## Python Validation

### Insert Size Analysis
```python
import pysam
import numpy as np
import matplotlib.pyplot as plt

def get_insert_sizes(bam_file, max_reads=100000):
    sizes = []
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for i, read in enumerate(bam.fetch()):
            if i >= max_reads:
                break
            if read.is_proper_pair and not read.is_secondary and read.template_length > 0:
                sizes.append(read.template_length)
    return sizes

sizes = get_insert_sizes('sample.bam')
print(f'Median: {np.median(sizes):.0f} bp')
print(f'Mean: {np.mean(sizes):.0f} bp')
print(f'Std: {np.std(sizes):.0f} bp')
```

### Comprehensive Validator
```python
import pysam
from collections import Counter

class AlignmentValidator:
    def __init__(self, bam_file):
        self.bam = pysam.AlignmentFile(bam_file, 'rb')

    def mapping_rate(self):
        stats = self.bam.get_index_statistics()
        mapped = sum(s.mapped for s in stats)
        unmapped = sum(s.unmapped for s in stats)
        return mapped / (mapped + unmapped) * 100

    def proper_pair_rate(self, sample_size=100000):
        proper, paired = 0, 0
        for i, read in enumerate(self.bam.fetch()):
            if i >= sample_size:
                break
            if read.is_paired:
                paired += 1
                if read.is_proper_pair:
                    proper += 1
        return proper / paired * 100 if paired > 0 else 0

    def strand_balance(self, sample_size=100000):
        forward, reverse = 0, 0
        for i, read in enumerate(self.bam.fetch()):
            if i >= sample_size:
                break
            if not read.is_unmapped:
                if read.is_reverse:
                    reverse += 1
                else:
                    forward += 1
        total = forward + reverse
        return forward / total if total > 0 else 0.5

    def report(self):
        map_rate = self.mapping_rate()
        pair_rate = self.proper_pair_rate()
        strand = self.strand_balance()

        print(f'Mapping rate: {map_rate:.1f}% {"PASS" if map_rate > 95 else "WARN" if map_rate > 90 else "FAIL"}')
        print(f'Proper pairing: {pair_rate:.1f}% {"PASS" if pair_rate > 90 else "WARN" if pair_rate > 80 else "FAIL"}')
        print(f'Strand balance: {strand:.3f} {"PASS" if 0.48 <= strand <= 0.52 else "WARN" if 0.45 <= strand <= 0.55 else "FAIL"}')

    def close(self):
        self.bam.close()
```

## Troubleshooting

### Low Mapping Rate (< 90%)
- Verify correct reference genome version
- Check sample species matches reference
- Look for contamination with other species
- Check read quality before alignment

### Poor Proper Pairing (< 80%)
- Insert size may not match aligner expectations
- Check for chimeric reads from structural variants
- Verify library preparation protocol
- Consider re-aligning with adjusted parameters

### Abnormal Insert Size
- Compare to library prep protocol expectations
- Bimodal distribution may indicate mixed libraries
- Very wide distribution suggests degraded DNA
- Check for adapter contamination

### Strand Imbalance
- May indicate capture bias in targeted sequencing
- Check for PCR amplification artifacts
- Examine specific chromosomes for local bias

### GC Bias
- Strong correlation indicates PCR amplification issues
- Consider library complexity problems
- May need GC correction before analysis
- Check for low-input DNA preparation issues

### Low MAPQ
- High proportion of repetitive regions
- Multi-mapping reads in duplicated genomic regions
- Consider longer reads or paired-end sequencing
- May need to filter low MAPQ reads for variant calling

## Tips
- Always validate alignments before variant calling or quantification
- Run validation on a subset first to check for major issues quickly
- Compare metrics across samples to identify outliers
- Keep validation reports for reproducibility and troubleshooting
- Consider experiment type when interpreting thresholds (RNA-seq vs WGS differ)
- Use Picard metrics for publication-quality QC reports
- Automate validation in pipelines with pass/fail thresholds
