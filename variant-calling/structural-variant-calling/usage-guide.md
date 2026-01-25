# Structural Variant Calling (Short Reads) Usage Guide

## Overview

Structural variants (SVs) are genomic alterations typically >50bp that include deletions, insertions, inversions, duplications, and translocations. Short-read SV calling uses paired-end and split-read information to detect these events.

## Prerequisites

```bash
conda install -c bioconda manta delly smoove survivor annotsv bcftools pysam
```

## Quick Start

Tell your AI agent what you want to do:
- "Call structural variants from my BAM file using multiple callers"
- "Run somatic SV calling on my tumor-normal pair with Manta and Delly"
- "Filter SV calls to keep only high-confidence variants with PE and SR support"
- "Annotate structural variants with gene overlaps and pathogenicity predictions"

## SV Detection Signatures

| Signature | SV Type | Evidence |
|-----------|---------|----------|
| Read depth | DEL, DUP | Coverage changes |
| Paired-end | DEL, INS, INV, DUP | Discordant insert sizes |
| Split-read | All types | Reads split across breakpoints |
| Mate-pair | Large SVs | Long-range linking |

## Tool Selection Guide

| Scenario | Recommended Tool |
|----------|------------------|
| General germline | Manta |
| Somatic SVs | Manta (tumor-normal) |
| High sensitivity | Multiple callers + merge |
| WES data | Manta --exome |
| Large cohort | Delly joint calling |
| Simple pipeline | Smoove (LUMPY wrapper) |

## Complete Pipeline

```bash
#!/bin/bash
set -euo pipefail

BAM=$1
REFERENCE=$2
SAMPLE=$(basename $BAM .bam)
OUTDIR=sv_calls/$SAMPLE

mkdir -p $OUTDIR

echo "=== Running SV calling for $SAMPLE ==="

# 1. Manta
echo "Running Manta..."
configManta.py \
    --bam $BAM \
    --referenceFasta $REFERENCE \
    --runDir $OUTDIR/manta

$OUTDIR/manta/runWorkflow.py -j 8
cp $OUTDIR/manta/results/variants/diploidSV.vcf.gz $OUTDIR/${SAMPLE}_manta.vcf.gz

# 2. Delly
echo "Running Delly..."
delly call -g $REFERENCE -o $OUTDIR/${SAMPLE}_delly.bcf $BAM
bcftools view $OUTDIR/${SAMPLE}_delly.bcf -Oz -o $OUTDIR/${SAMPLE}_delly.vcf.gz
bcftools index $OUTDIR/${SAMPLE}_delly.vcf.gz

# 3. Smoove (LUMPY)
echo "Running Smoove..."
smoove call --name $SAMPLE --fasta $REFERENCE --outdir $OUTDIR/smoove -p 8 $BAM
cp $OUTDIR/smoove/${SAMPLE}-smoove.genotyped.vcf.gz $OUTDIR/${SAMPLE}_lumpy.vcf.gz

# 4. Merge callers
echo "Merging callsets..."
gunzip -c $OUTDIR/${SAMPLE}_manta.vcf.gz > $OUTDIR/manta.vcf
gunzip -c $OUTDIR/${SAMPLE}_delly.vcf.gz > $OUTDIR/delly.vcf
gunzip -c $OUTDIR/${SAMPLE}_lumpy.vcf.gz > $OUTDIR/lumpy.vcf

echo "$OUTDIR/manta.vcf" > $OUTDIR/vcf_list.txt
echo "$OUTDIR/delly.vcf" >> $OUTDIR/vcf_list.txt
echo "$OUTDIR/lumpy.vcf" >> $OUTDIR/vcf_list.txt

SURVIVOR merge $OUTDIR/vcf_list.txt 500 2 1 1 0 50 $OUTDIR/${SAMPLE}_merged.vcf

# 5. Filter
echo "Filtering..."
bcftools view -i 'QUAL >= 20' $OUTDIR/${SAMPLE}_merged.vcf -Oz -o $OUTDIR/${SAMPLE}_final.vcf.gz
bcftools index $OUTDIR/${SAMPLE}_final.vcf.gz

# 6. Annotate
echo "Annotating..."
AnnotSV -SVinputFile $OUTDIR/${SAMPLE}_final.vcf.gz \
    -genomeBuild GRCh38 \
    -outputFile $OUTDIR/${SAMPLE}_annotated

echo "Done! Final SVs: $OUTDIR/${SAMPLE}_final.vcf.gz"
```

## Somatic SV Calling

```bash
#!/bin/bash
TUMOR=$1
NORMAL=$2
REFERENCE=$3
OUTDIR=somatic_svs

mkdir -p $OUTDIR

# Manta somatic
configManta.py \
    --tumorBam $TUMOR \
    --normalBam $NORMAL \
    --referenceFasta $REFERENCE \
    --runDir $OUTDIR/manta

$OUTDIR/manta/runWorkflow.py -j 8

# Delly somatic
delly call -g $REFERENCE -o $OUTDIR/delly.bcf $TUMOR $NORMAL

echo -e "tumor\ttumor\nnormal\tcontrol" > $OUTDIR/samples.tsv

delly filter -f somatic -o $OUTDIR/delly_somatic.bcf \
    -s $OUTDIR/samples.tsv $OUTDIR/delly.bcf

bcftools view $OUTDIR/delly_somatic.bcf -Oz -o $OUTDIR/delly_somatic.vcf.gz
```

## Quality Filtering

```python
import pandas as pd
import pysam

def filter_sv_calls(vcf_in, vcf_out, min_qual=20, min_size=50, min_support=3):
    '''Filter SV calls by quality metrics.'''
    vcf = pysam.VariantFile(vcf_in)
    out = pysam.VariantFile(vcf_out, 'w', header=vcf.header)

    kept, filtered = 0, 0
    for rec in vcf:
        # Quality filter
        if rec.qual and rec.qual < min_qual:
            filtered += 1
            continue

        # Size filter
        svlen = abs(rec.info.get('SVLEN', [0])[0])
        if svlen < min_size:
            filtered += 1
            continue

        # Support filter (if available)
        support = rec.info.get('SUPPORT', rec.info.get('SR', [0]))
        if isinstance(support, tuple):
            support = sum(support)
        if support < min_support:
            filtered += 1
            continue

        out.write(rec)
        kept += 1

    out.close()
    print(f'Kept: {kept}, Filtered: {filtered}')
```

## Interpreting SV VCF

### Key INFO Fields
| Field | Description |
|-------|-------------|
| SVTYPE | SV type (DEL, INS, DUP, INV, BND) |
| SVLEN | SV length (negative for DEL) |
| END | End position |
| CIPOS | Confidence interval around POS |
| CIEND | Confidence interval around END |
| PE | Paired-end support |
| SR | Split-read support |

### Format Fields
| Field | Description |
|-------|-------------|
| GT | Genotype (0/0, 0/1, 1/1) |
| PR | Paired-read count (ref, alt) |
| SR | Split-read count (ref, alt) |

## Benchmarking

Compare calls to truth set:

```bash
# Using Truvari
truvari bench \
    -b truth.vcf.gz \
    -c calls.vcf.gz \
    -o benchmark_results/ \
    --passonly \
    -r 500 \
    -p 0.7
```

## Common Issues

### High False Positives
- Filter by multiple caller support
- Require both PE and SR evidence
- Remove blacklist regions

### Missing Large SVs
- Check BAM insert size distribution
- Ensure sufficient coverage
- Consider long-read sequencing

### Breakpoint Imprecision
- Normal for short reads (~100bp uncertainty)
- Long reads improve resolution
- Use CIPOS/CIEND for confidence intervals

## Relationship with CNV Detection

| Tool | Detects |
|------|---------|
| Manta/Delly/LUMPY | Balanced SVs (INV, BND), small DEL/DUP |
| CNVkit/GATK CNV | Large copy number changes |
| Both together | Comprehensive SV detection |

For complete structural variation analysis, combine SV callers with CNV detection tools.

## Example Prompts

> "Call structural variants from my BAM file using multiple callers and merge the results"

> "Run somatic SV calling on my tumor-normal pair with Manta and Delly"

> "Filter my SV calls to keep only high-confidence variants with PE and SR support"

> "Annotate my structural variants with gene overlaps and pathogenicity predictions"
