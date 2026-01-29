# VCF Statistics Usage Guide

## Overview

This guide covers generating variant statistics and quality metrics.

## Prerequisites

- bcftools installed (`conda install -c bioconda bcftools`)
- plot-vcfstats (comes with bcftools, requires python and matplotlib)

## Quick Start

Tell your AI agent what you want to do:
- "Generate comprehensive statistics for my VCF including Ti/Tv ratio"
- "Compare variant statistics before and after filtering"
- "Create QC plots showing depth and quality distributions"
- "Check for sample swaps or contamination in my cohort VCF"

## bcftools stats Overview

`bcftools stats` generates comprehensive statistics about VCF files:

- Variant counts (SNPs, indels, multiallelic)
- Transition/transversion ratio
- Quality score distribution
- Depth distribution
- Allele frequency spectrum
- Indel length distribution
- Per-sample statistics

## Basic Usage

### Generate Statistics

```bash
bcftools stats input.vcf.gz > stats.txt
```

### View Summary Numbers

```bash
bcftools stats input.vcf.gz | grep "^SN"
```

Output:
```
SN      0       number of samples:      3
SN      0       number of records:      150000
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 130000
SN      0       number of MNPs: 100
SN      0       number of indels:       19900
SN      0       number of others:       0
SN      0       number of multiallelic sites:   5000
```

## Understanding Output Sections

### Section Codes

| Code | Description |
|------|-------------|
| `SN` | Summary numbers |
| `TSTV` | Transition/transversion stats |
| `SiS` | Singleton stats |
| `AF` | Allele frequency distribution |
| `QUAL` | Quality score distribution |
| `IDD` | Indel length distribution |
| `ST` | Substitution types |
| `DP` | Depth distribution |
| `PSC` | Per-sample counts |
| `PSI` | Per-sample indels |

### Transition/Transversion Ratio

```bash
bcftools stats input.vcf.gz | grep "^TSTV"
```

Output:
```
TSTV    0       100000  48000   2.08    99500   47800   2.08
```

Columns: transitions, transversions, Ti/Tv ratio

Expected ratios:
- Whole genome: ~2.0-2.1
- Exome: ~2.8-3.3
- Low ratio may indicate sequencing artifacts

### Quality Distribution

```bash
bcftools stats input.vcf.gz | grep "^QUAL"
```

Shows distribution of QUAL scores.

### Depth Distribution

```bash
bcftools stats input.vcf.gz | grep "^DP"
```

Shows read depth distribution.

## Per-Sample Statistics

### Enable Per-Sample Mode

```bash
bcftools stats -s - input.vcf.gz > stats.txt
```

The `-s -` means "all samples".

### View Per-Sample Counts

```bash
bcftools stats -s - input.vcf.gz | grep "^PSC"
```

Output columns:
- Sample name
- Homozygous reference sites
- Heterozygous sites
- Homozygous alternate sites
- Transitions
- Transversions
- Missing sites

### View Per-Sample Indels

```bash
bcftools stats -s - input.vcf.gz | grep "^PSI"
```

## Comparing VCF Files

### Compare Two Files

```bash
bcftools stats file1.vcf.gz file2.vcf.gz > comparison.txt
```

Each file gets an index (0, 1) in the output.

### View Comparison

```bash
# File 1 stats
grep "^SN.*0" comparison.txt

# File 2 stats
grep "^SN.*1" comparison.txt
```

## Region-Based Statistics

### Specific Region

```bash
bcftools stats -r chr1:1000000-2000000 input.vcf.gz > region_stats.txt
```

### From BED File

```bash
bcftools stats -R targets.bed input.vcf.gz > target_stats.txt
```

### Exome vs Whole Genome

```bash
# Exome regions only
bcftools stats -R exome.bed input.vcf.gz > exome_stats.txt

# Compare with whole genome
bcftools stats input.vcf.gz > genome_stats.txt
```

## Generating Plots

### Create All Plots

```bash
bcftools stats input.vcf.gz > stats.txt
plot-vcfstats -p output_directory stats.txt
```

Creates:
- `summary.pdf` - Multi-page PDF with all plots
- Individual PNG files for each plot

### Plot Types Generated

- Substitution types
- Indel length distribution
- Ts/Tv by quality
- Indels per sample
- SNPs per sample
- Quality score distribution
- Depth distribution

### Compare Multiple Files with Plots

```bash
bcftools stats file1.vcf.gz file2.vcf.gz > comparison.txt
plot-vcfstats -p comparison_plots comparison.txt
```

## Sample Identity Checking

### Check Sample Concordance

```bash
bcftools gtcheck -g reference.vcf.gz query.vcf.gz
```

Reports how well samples match between files.

### Detect Sample Swaps

```bash
bcftools gtcheck -G 1 cohort.vcf.gz
```

Compares all samples pairwise to detect:
- Sample swaps
- Related individuals
- Duplicate samples

### Interpret Results

Output columns:
- Sample pair
- Discordance rate
- Sites compared
- Discordant sites

Low discordance (<0.1) suggests same individual or close relatives.

## Quick Statistics with Query

### Variant Counts

```bash
# Total variants
bcftools view -H input.vcf.gz | wc -l

# SNPs only
bcftools view -v snps -H input.vcf.gz | wc -l

# Indels only
bcftools view -v indels -H input.vcf.gz | wc -l

# PASS variants
bcftools view -f PASS -H input.vcf.gz | wc -l
```

### Mean Quality Score

```bash
bcftools query -f '%QUAL\n' input.vcf.gz | \
    awk '{sum+=$1; n++} END {print "Mean QUAL:", sum/n}'
```

### Mean Depth

```bash
bcftools query -f '%INFO/DP\n' input.vcf.gz | \
    awk '{sum+=$1; n++} END {print "Mean DP:", sum/n}'
```

### Variants per Chromosome

```bash
bcftools view -H input.vcf.gz | cut -f1 | sort | uniq -c
```

## Complete QC Workflow

### Generate QC Report

```bash
#!/bin/bash
VCF="$1"
OUTDIR="qc_report"
mkdir -p "$OUTDIR"

# Generate stats
bcftools stats -s - "$VCF" > "$OUTDIR/stats.txt"

# Extract key metrics
echo "=== Summary ===" > "$OUTDIR/summary.txt"
grep "^SN" "$OUTDIR/stats.txt" >> "$OUTDIR/summary.txt"

echo "" >> "$OUTDIR/summary.txt"
echo "=== Ti/Tv Ratio ===" >> "$OUTDIR/summary.txt"
grep "^TSTV" "$OUTDIR/stats.txt" | head -1 | awk '{print $5}' >> "$OUTDIR/summary.txt"

# Generate plots
plot-vcfstats -p "$OUTDIR/plots" "$OUTDIR/stats.txt"

echo "QC report in $OUTDIR"
```

### Filter Assessment

```bash
# Stats before filtering
bcftools stats raw.vcf.gz > before.txt

# Filter
bcftools filter -e 'QUAL<30 || INFO/DP<10' raw.vcf.gz -Oz -o filtered.vcf.gz

# Stats after filtering
bcftools stats filtered.vcf.gz > after.txt

# Compare
echo "Before filtering:"
grep "^SN.*number of SNPs" before.txt
echo "After filtering:"
grep "^SN.*number of SNPs" after.txt
```

## Python Statistics with cyvcf2

### Comprehensive Stats

```python
from cyvcf2 import VCF
import numpy as np

vcf = VCF('input.vcf.gz')

stats = {
    'total': 0,
    'snps': 0,
    'indels': 0,
    'pass': 0,
    'quals': [],
    'depths': []
}

for variant in vcf:
    stats['total'] += 1
    if variant.is_snp:
        stats['snps'] += 1
    elif variant.is_indel:
        stats['indels'] += 1
    if variant.FILTER is None:
        stats['pass'] += 1
    if variant.QUAL:
        stats['quals'].append(variant.QUAL)
    dp = variant.INFO.get('DP')
    if dp:
        stats['depths'].append(dp)

print(f'Total variants: {stats["total"]}')
print(f'SNPs: {stats["snps"]}')
print(f'Indels: {stats["indels"]}')
print(f'PASS: {stats["pass"]}')
print(f'Mean QUAL: {np.mean(stats["quals"]):.1f}')
print(f'Mean DP: {np.mean(stats["depths"]):.1f}')
```

## Troubleshooting

### "No data found"

VCF may be empty or have wrong format:

```bash
bcftools view -H input.vcf.gz | head
```

### "plot-vcfstats not found"

Install matplotlib and ensure bcftools scripts are in PATH:

```bash
pip install matplotlib
which plot-vcfstats
```

### Very Low Ti/Tv Ratio

May indicate:
- Sequencing artifacts
- Contamination
- Wrong reference genome

Check:
```bash
# Verify reference matches
bcftools norm -c w -f reference.fa input.vcf.gz > /dev/null
```

## Example Prompts

> "Generate comprehensive statistics for my VCF file including Ti/Tv ratio"

> "Compare variant statistics before and after filtering"

> "Create QC plots showing depth and quality distributions"

> "Check for sample swaps or contamination in my cohort VCF"

## See Also

- [bcftools stats documentation](http://www.htslib.org/doc/bcftools.html#stats)
- [bcftools gtcheck documentation](http://www.htslib.org/doc/bcftools.html#gtcheck)
- **filtering-best-practices** - Filter based on statistics
- **vcf-basics** - Query VCF data
