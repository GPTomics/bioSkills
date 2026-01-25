# GATK CNV Workflow Usage Guide

## Overview

GATK provides a robust, best-practices CNV calling workflow. It uses principal component analysis to denoise samples against a panel of normals and integrates allelic information for improved accuracy.

## Installation

```bash
conda install -c bioconda gatk4
# or download from Broad Institute
```

## When to Use GATK vs CNVkit

| Factor | GATK | CNVkit |
|--------|------|--------|
| WGS | Excellent | Good |
| WES | Good | Excellent |
| Panel | Limited | Good |
| Germline | Yes (cohort mode) | Limited |
| Speed | Slower | Faster |
| PON required | Recommended | Recommended |

## Prerequisites

1. **Reference files**:
   - Reference FASTA + index + dict
   - Interval list (targets for WES)
   - Common SNP sites VCF (for allelic counts)

2. **Sample requirements**:
   - Aligned, sorted, indexed BAM files
   - 5-10+ normal samples for PON

## Quick Start

Tell your AI agent what you want to do:
- "Run the GATK CNV workflow on my tumor-normal WGS pair"
- "Create a panel of normals from my cohort for GATK CNV calling"
- "Compare GATK and CNVkit results to identify high-confidence CNVs"
- "Generate CNV plots from GATK ModelSegments output"

## Creating Interval Files

```bash
# From BED file
gatk BedToIntervalList \
    -I targets.bed \
    -SD reference.dict \
    -O targets.interval_list

# Preprocess for CNV calling
gatk PreprocessIntervals \
    -R reference.fa \
    -L targets.interval_list \
    --bin-length 0 \
    -O preprocessed.interval_list
```

## Building a Panel of Normals

Critical for accurate denoising:

```bash
# Step 1: Collect counts from each normal
for bam in normal*.bam; do
    sample=$(basename $bam .bam)
    gatk CollectReadCounts \
        -R reference.fa \
        -I $bam \
        -L preprocessed.interval_list \
        -O ${sample}.counts.hdf5
done

# Step 2: Create PON
gatk CreateReadCountPanelOfNormals \
    -I normal1.counts.hdf5 \
    -I normal2.counts.hdf5 \
    -I normal3.counts.hdf5 \
    -I normal4.counts.hdf5 \
    -I normal5.counts.hdf5 \
    --minimum-interval-median-percentile 5.0 \
    --number-of-eigensamples 20 \
    -O cnv_pon.hdf5
```

## Complete Somatic Workflow

```bash
#!/bin/bash
set -euo pipefail

TUMOR_BAM=$1
NORMAL_BAM=$2
SAMPLE_NAME=$3
OUTDIR=results/${SAMPLE_NAME}
mkdir -p $OUTDIR

REF=reference.fa
INTERVALS=preprocessed.interval_list
PON=cnv_pon.hdf5
SNPS=common_biallelic_snps.vcf.gz

# 1. Collect read counts
echo "Collecting read counts..."
gatk CollectReadCounts -R $REF -I $TUMOR_BAM -L $INTERVALS \
    -O $OUTDIR/tumor.counts.hdf5

gatk CollectReadCounts -R $REF -I $NORMAL_BAM -L $INTERVALS \
    -O $OUTDIR/normal.counts.hdf5

# 2. Denoise tumor
echo "Denoising..."
gatk DenoiseReadCounts \
    -I $OUTDIR/tumor.counts.hdf5 \
    --count-panel-of-normals $PON \
    --standardized-copy-ratios $OUTDIR/tumor.standardized.tsv \
    --denoised-copy-ratios $OUTDIR/tumor.denoised.tsv

# 3. Collect allelic counts
echo "Collecting allelic counts..."
gatk CollectAllelicCounts -R $REF -I $TUMOR_BAM -L $SNPS \
    -O $OUTDIR/tumor.allelicCounts.tsv

gatk CollectAllelicCounts -R $REF -I $NORMAL_BAM -L $SNPS \
    -O $OUTDIR/normal.allelicCounts.tsv

# 4. Model segments
echo "Modeling segments..."
gatk ModelSegments \
    --denoised-copy-ratios $OUTDIR/tumor.denoised.tsv \
    --allelic-counts $OUTDIR/tumor.allelicCounts.tsv \
    --normal-allelic-counts $OUTDIR/normal.allelicCounts.tsv \
    --output-prefix tumor \
    -O $OUTDIR/

# 5. Call segments
echo "Calling segments..."
gatk CallCopyRatioSegments \
    -I $OUTDIR/tumor.cr.seg \
    -O $OUTDIR/tumor.called.seg

# 6. Plot results
echo "Generating plots..."
gatk PlotModeledSegments \
    --denoised-copy-ratios $OUTDIR/tumor.denoised.tsv \
    --allelic-counts $OUTDIR/tumor.hets.tsv \
    --segments $OUTDIR/tumor.modelFinal.seg \
    --sequence-dictionary ${REF%.fa}.dict \
    --output-prefix tumor \
    -O $OUTDIR/plots/

echo "Done! Results in $OUTDIR"
```

## Interpreting Results

### .called.seg file columns:
- CONTIG, START, END: Genomic coordinates
- NUM_POINTS_COPY_RATIO: Number of intervals
- CALL: +, -, 0 (amplification, deletion, neutral)

### Copy ratio interpretation:
| LOG2_COPY_RATIO | Interpretation |
|-----------------|----------------|
| < -0.5 | Deletion |
| -0.5 to 0.5 | Neutral |
| > 0.5 | Amplification |

## Quality Control

```bash
# Check denoising quality
# Good: MAD < 0.3, smooth profiles
# Bad: High variance, many outliers

# Check number of het sites
wc -l tumor.hets.tsv  # Should be >10,000 for WGS

# Compare tumor vs normal profiles
```

## Troubleshooting

### Few segments called
- Check PON quality (need sufficient normals)
- Verify interval preprocessing
- Check sample quality (coverage, contamination)

### Noisy profiles
- Increase PON size
- Filter problematic intervals
- Check for batch effects

### No het sites
- Verify SNP VCF has sites in captured regions
- Check allelic count thresholds

## Comparison with CNVkit

```python
import pandas as pd

# Load GATK results
gatk_segs = pd.read_csv('tumor.called.seg', sep='\t', comment='@')

# Load CNVkit results
cnvkit_segs = pd.read_csv('tumor.cns', sep='\t')

# Compare at gene level
# Merge on overlapping coordinates
```

Both tools should be run and results compared for high-confidence calls.

## Example Prompts

> "Run the GATK CNV workflow on my tumor-normal WGS pair"

> "Create a panel of normals from my cohort for GATK CNV calling"

> "Compare GATK and CNVkit results to identify high-confidence CNV calls"

> "Generate CNV plots from GATK ModelSegments output"
