# CNVkit Analysis Usage Guide

## Overview

CNVkit is the standard tool for detecting copy number variants from targeted sequencing (exome, gene panels). It uses both on-target and off-target reads to infer copy number across the genome.

## Prerequisites

```bash
# Requires Python >= 3.10, NumPy 2.x, Pandas 3.x compatible
conda install -c bioconda cnvkit
# or
pip install cnvkit
# For HMM segmentation: pomegranate >= 1.0 (installed automatically)
```

Dependencies: R with DNAcopy package (for CBS segmentation).

## Quick Start

Tell your AI agent what you want to do:
- "Run CNVkit on my tumor-normal exome pair to call copy number variants"
- "Build a panel of normals from my control samples for CNV calling"
- "Export my CNVkit segments to VCF format for downstream analysis"
- "Check quality metrics for my CNVkit run and identify noisy samples"
- "Generate bedGraph coverage files for sharing without exposing raw sequences"

## Input Requirements

1. **BAM files**: Aligned, sorted, indexed
2. **Target BED**: Capture regions (exome targets)
3. **Reference FASTA**: Same as used for alignment
4. **Gene annotations** (optional): refFlat.txt for gene names

## Workflow Selection

| Scenario | Command |
|----------|---------|
| Tumor-normal pair | `batch tumor.bam --normal normal.bam ...` |
| Multiple normals → reference | `batch --normal *.bam ...` then `batch tumor.bam --reference ref.cnn` |
| Tumor-only (flat reference) | `batch tumor.bam --targets ... --fasta ...` |
| WGS | Add `--method wgs` |
| Germline | Use `--normal` samples only |

## Complete Tumor-Normal Pipeline

```bash
# Variables
TUMOR=tumor.bam
NORMAL=normal.bam
TARGETS=exome_targets.bed
FASTA=reference.fa
REFFLAT=refFlat.txt
OUTDIR=cnvkit_results

mkdir -p $OUTDIR

# Run batch command
cnvkit.py batch $TUMOR \
    --normal $NORMAL \
    --targets $TARGETS \
    --annotate $REFFLAT \
    --fasta $FASTA \
    --output-reference $OUTDIR/reference.cnn \
    --output-dir $OUTDIR

# Call with purity adjustment (if known)
cnvkit.py call $OUTDIR/tumor.cns \
    --purity 0.7 \
    -o $OUTDIR/tumor.call.cns

# Generate plots
cnvkit.py scatter $OUTDIR/tumor.cnr -s $OUTDIR/tumor.cns -o $OUTDIR/tumor_scatter.png
cnvkit.py diagram $OUTDIR/tumor.cnr -s $OUTDIR/tumor.cns -o $OUTDIR/tumor_diagram.pdf
```

## Building a Panel of Normals

For best results, use 5-10+ matched normal samples:

```bash
# Step 1: Generate reference from normals
cnvkit.py batch \
    --normal normal1.bam normal2.bam normal3.bam normal4.bam normal5.bam \
    --targets targets.bed \
    --annotate refFlat.txt \
    --fasta reference.fa \
    --output-reference panel_of_normals.cnn

# Step 2: Process tumors with PON
for tumor in tumor*.bam; do
    sample=$(basename $tumor .bam)
    cnvkit.py batch $tumor \
        --reference panel_of_normals.cnn \
        --output-dir results/$sample
done
```

## Privacy-Preserving Workflow (bedGraph)

CNVkit accepts pre-computed bedGraph files in place of BAMs, enabling coverage sharing without exposing raw sequences:

```bash
# Step 1: Generate bedGraph from BAM
bedtools genomecov -ibam tumor.bam -bg | bgzip > tumor.bed.gz
tabix -p bed tumor.bed.gz

# Step 2: Use bedGraph as input to CNVkit coverage
cnvkit.py coverage tumor.bed.gz targets.target.bed -o tumor.targetcoverage.cnn
cnvkit.py coverage tumor.bed.gz antitargets.bed -o tumor.antitargetcoverage.cnn

# Step 3: Continue with fix/segment/call as normal
```

## Interpreting Results

### CNR file (copy ratios)
- Each row is a bin (target or off-target region)
- `log2` column: log2 ratio vs reference (0 = diploid)
- `weight`: confidence weight

### CNS file (segments)
- Merged adjacent bins with similar log2
- `cn`: absolute copy number (after calling)
- `probes`: number of bins in segment

### Thresholds for calling:
| log2 | CN State | Interpretation |
|------|----------|----------------|
| < -1.1 | 0 | Homozygous deletion |
| -1.1 to -0.25 | 1 | Heterozygous deletion |
| -0.25 to 0.2 | 2 | Diploid |
| 0.2 to 0.7 | 3 | Single copy gain |
| > 0.7 | 4+ | Amplification |

## Quality Control

```bash
# Check metrics for all samples
cnvkit.py metrics *.cnr -s *.cns

# Good quality indicators:
# - Median absolute deviation (MAD) < 0.5
# - Biweight midvariance < 0.5
# - Few bins with extreme values

# Check inferred sex
cnvkit.py sex *.cnr *.cnn
```

## Detailed Metrics

Both `segmetrics` and `genemetrics` share most statistics flags: `--ci`, `--pi`, `--iqr`, `--alpha`, `--bootstrap`, `--mean`, `--median`, `--stdev`, `--sem`, `--mad`.

Two key differences: (1) t-test flag naming -- `segmetrics` uses `--t-test` (hyphenated), `genemetrics` uses `--ttest` (no hyphen); (2) `--smooth-bootstrap` exists only in `segmetrics` (smoothed CI).

```bash
# Per-segment confidence intervals and prediction intervals
cnvkit.py segmetrics sample.cnr -s sample.cns --ci --pi --bootstrap 100 -o sample.segmetrics.cns

# Gene-level CNV report with confidence intervals
# 0.2: CNVkit default gain/loss log2 threshold (2^0.2 ~ 15% copy number change)
cnvkit.py genemetrics sample.cnr -s sample.cns --threshold 0.2 --ci --bootstrap 100 -o sample.genemetrics.tsv
```

## Common Issues

### Noisy data
```bash
# Increase smoothing during segmentation
cnvkit.py segment sample.cnr --smooth-cbs -o sample.cns

# Or use HMM which handles tumor heterogeneity better than CBS
cnvkit.py segment sample.cnr --method hmm-tumor -o sample.cns
```

### Low coverage samples
```bash
# Lower minimum coverage thresholds
cnvkit.py batch tumor.bam --normal normal.bam \
    --targets targets.bed --fasta reference.fa \
    --target-avg-size 200
```

### Wrong ploidy/purity
```bash
# Estimate from data
cnvkit.py call sample.cns --center median -o sample.call.cns

# Or specify if known
cnvkit.py call sample.cns --purity 0.7 --ploidy 2 -o sample.call.cns
```

## Integration with Other Tools

### Export for downstream analysis
```bash
# For GISTIC2 (identifying recurrent CNVs)
cnvkit.py export seg *.cns -o cohort.seg

# For integration with SVs
cnvkit.py export vcf sample.cns -o sample.cnv.vcf
```

### Python analysis
```python
import cnvlib
import pandas as pd

# Load and filter
cns = cnvlib.read('sample.cns')
amplified = cns[(cns['log2'] > 0.5) & (cns['cn'] >= 4)]
deleted = cns[(cns['log2'] < -0.5) & (cns['cn'] <= 1)]

# Genes in amplified regions
print(amplified[['chromosome', 'start', 'end', 'gene', 'log2', 'cn']])
```

## Example Prompts

> "Run CNVkit on my tumor-normal exome pair and call copy number variants"

> "Build a panel of normals from my 10 control samples for CNV calling"

> "Export my CNVkit segments to VCF format for integration with SNV calls"

> "Check the quality metrics for my CNVkit run and identify noisy samples"

> "Generate per-segment confidence intervals and identify genes with significant copy number changes"

## What the Agent Will Do

1. Prepare target and antitarget regions from capture BED
2. Calculate coverage from BAM or bedGraph input
3. Build or apply a reference from normal samples
4. Segment and call copy number states
5. Generate plots and export results
6. Run quality metrics and gene-level reports

## Tips

- **Panel of normals** - Use 5-10+ normal samples for best results
- **HMM segmentation** - hmm-tumor handles noisy tumor data better than CBS
- **MAD threshold** - Check MAD < 0.5 in metrics output as a quality indicator
- **Purity/ploidy** - Use --purity and --ploidy when values are known for more accurate calls
- **bedGraph input** - Enables sharing coverage data without exposing raw sequences

## Related Skills

- alignment-files/bam-statistics - QC of input BAMs
- copy-number/cnv-visualization - Advanced CNV plotting
- copy-number/cnv-annotation - Gene-level annotation of CNV calls
- copy-number/gatk-cnv - Alternative CNV caller
- long-read-sequencing/structural-variants - Complementary SV calling
