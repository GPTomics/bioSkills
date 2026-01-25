---
name: bio-copy-number-cnvkit-analysis
description: Detect copy number variants from targeted/exome sequencing using CNVkit. Supports tumor-normal pairs, tumor-only, and germline CNV calling. Use when detecting CNVs from WES or targeted panel sequencing data.
tool_type: cli
primary_tool: cnvkit
---

# CNVkit CNV Analysis

## Basic Workflow

```bash
# Complete pipeline for tumor-normal pair
cnvkit.py batch tumor.bam \
    --normal normal.bam \
    --targets targets.bed \
    --fasta reference.fa \
    --output-reference my_reference.cnn \
    --output-dir results/
```

## Build Reference from Normal Samples

```bash
# Step 1: Build reference from multiple normals (recommended)
cnvkit.py batch \
    --normal normal1.bam normal2.bam normal3.bam \
    --targets targets.bed \
    --fasta reference.fa \
    --output-reference pooled_reference.cnn

# Step 2: Run on tumor samples using pre-built reference
cnvkit.py batch tumor1.bam tumor2.bam \
    --reference pooled_reference.cnn \
    --output-dir results/
```

## Flat Reference (No Matched Normal)

```bash
# When no matched normal is available
cnvkit.py batch tumor.bam \
    --targets targets.bed \
    --fasta reference.fa \
    --output-reference flat_reference.cnn \
    --output-dir results/
```

## WGS Mode

```bash
# For whole genome sequencing (no targets file)
cnvkit.py batch tumor.bam \
    --normal normal.bam \
    --fasta reference.fa \
    --method wgs \
    --output-dir results/
```

## Step-by-Step Pipeline

```bash
# 1. Generate target and antitarget regions
cnvkit.py target targets.bed --annotate refFlat.txt -o targets.target.bed
cnvkit.py antitarget targets.bed -o targets.antitarget.bed

# 2. Calculate coverage
cnvkit.py coverage tumor.bam targets.target.bed -o tumor.targetcoverage.cnn
cnvkit.py coverage tumor.bam targets.antitarget.bed -o tumor.antitargetcoverage.cnn
cnvkit.py coverage normal.bam targets.target.bed -o normal.targetcoverage.cnn
cnvkit.py coverage normal.bam targets.antitarget.bed -o normal.antitargetcoverage.cnn

# 3. Build reference
cnvkit.py reference normal.targetcoverage.cnn normal.antitargetcoverage.cnn \
    --fasta reference.fa -o reference.cnn

# 4. Fix and call
cnvkit.py fix tumor.targetcoverage.cnn tumor.antitargetcoverage.cnn reference.cnn -o tumor.cnr
cnvkit.py segment tumor.cnr -o tumor.cns
cnvkit.py call tumor.cns -o tumor.call.cns
```

## Segmentation Options

```bash
# Default CBS (Circular Binary Segmentation)
cnvkit.py segment sample.cnr -o sample.cns

# Use HMM for better performance
cnvkit.py segment sample.cnr --method hmm -o sample.cns

# Adjust smoothing
cnvkit.py segment sample.cnr --smooth-cbs -o sample.cns
```

## CNV Calling with Ploidy/Purity

```bash
# Specify tumor purity and ploidy
cnvkit.py call sample.cns \
    --purity 0.7 \
    --ploidy 2 \
    -o sample.call.cns

# With B-allele frequencies (from VCF)
cnvkit.py call sample.cns \
    --vcf sample.vcf \
    --purity 0.7 \
    -o sample.call.cns
```

## Export Results

```bash
# Export to BED format
cnvkit.py export bed sample.call.cns -o sample.cnv.bed

# Export to VCF
cnvkit.py export vcf sample.call.cns -o sample.cnv.vcf

# Export segments for GISTIC2
cnvkit.py export seg *.cns -o samples.seg

# Export for Nexus
cnvkit.py export nexus-basic sample.cnr -o sample.nexus.txt
```

## Visualization

```bash
# Scatter plot with segments
cnvkit.py scatter sample.cnr -s sample.cns -o sample_scatter.png

# Single chromosome
cnvkit.py scatter sample.cnr -s sample.cns -c chr17 -o sample_chr17.png

# Diagram (ideogram style)
cnvkit.py diagram sample.cnr -s sample.cns -o sample_diagram.pdf

# Heatmap across samples
cnvkit.py heatmap *.cns -o heatmap.pdf
```

## Key Output Files

| Extension | Description |
|-----------|-------------|
| .cnn | Reference or coverage file |
| .cnr | Copy ratios (log2) per bin |
| .cns | Segmented copy ratios |
| .call.cns | Called copy number states |

## Python API

```python
import cnvlib

# Load data
cnr = cnvlib.read('sample.cnr')
cns = cnvlib.read('sample.cns')

# Filter by chromosome
chr17 = cnr[cnr.chromosome == 'chr17']

# Get amplifications
amps = cns[cns['log2'] > 0.5]

# Get deletions
dels = cns[cns['log2'] < -0.5]

# Export
cnr.to_csv('sample.cnr.tsv', sep='\t', index=False)
```

## Quality Control

```bash
# Check reference quality
cnvkit.py metrics *.cnr -s *.cns

# Check for gender mismatches
cnvkit.py sex *.cnr *.cnn

# Median absolute deviation (lower is better)
# Biweight midvariance (sample heterogeneity)
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| --method | hybrid | hybrid, wgs, amplicon |
| --segment-method | cbs | cbs, hmm, haar, none |
| --drop-low-coverage | on | Drop low-coverage bins |
| --purity | 1.0 | Tumor purity (0-1) |
| --ploidy | 2 | Sample ploidy |
| --thresholds | -1.1,-0.25,0.2,0.7 | CN state thresholds |

## Related Skills

- alignment-files/bam-statistics - QC of input BAMs
- copy-number/cnv-visualization - Advanced plotting
- copy-number/cnv-annotation - Gene-level annotation
- long-read-sequencing/structural-variants - Complementary SV calling
