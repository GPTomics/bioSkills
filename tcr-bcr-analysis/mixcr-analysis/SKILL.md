---
name: bio-tcr-bcr-analysis-mixcr-analysis
description: Perform V(D)J alignment and clonotype assembly from TCR-seq or BCR-seq data using MiXCR. Use when processing raw immune repertoire sequencing data to identify clonotypes and their frequencies.
tool_type: cli
primary_tool: MiXCR
---

# MiXCR Analysis

## Complete Workflow (Recommended)

```bash
# One-command analysis for common protocols
mixcr analyze amplicon \
    -s human \
    --starting-material rna \
    --5-end v-primers \
    --3-end c-primers \
    --adapters adapters-cdr3 \
    input_R1.fastq.gz input_R2.fastq.gz \
    output_prefix

# For 10x Genomics VDJ
mixcr analyze 10x-vdj-tcr \
    input_R1.fastq.gz input_R2.fastq.gz \
    output_prefix
```

## Step-by-Step Workflow

### Step 1: Align Reads

```bash
mixcr align \
    -s human \
    -p rna-seq \
    -OallowPartialAlignments=true \
    input_R1.fastq.gz input_R2.fastq.gz \
    alignments.vdjca
```

### Step 2: Assemble Clonotypes

```bash
# Assemble partial alignments (for fragmented data)
mixcr assemblePartial alignments.vdjca alignments_rescued.vdjca

# Extend alignments
mixcr extend alignments_rescued.vdjca alignments_extended.vdjca

# Assemble clonotypes
mixcr assemble alignments_extended.vdjca clones.clns
```

### Step 3: Export Results

```bash
# Export clonotypes to table
mixcr exportClones \
    -c TRB \
    -p fullImputed \
    clones.clns \
    clones.txt

# Export with specific columns
mixcr exportClones \
    -cloneId -count -fraction \
    -nFeature CDR3 -aaFeature CDR3 \
    -vGene -dGene -jGene \
    clones.clns \
    clones_custom.txt
```

## Preset Protocols

| Protocol | Use Case |
|----------|----------|
| `amplicon` | Targeted TCR/BCR amplicon sequencing |
| `shotgun` | RNA-seq with TCR/BCR extraction |
| `10x-vdj-tcr` | 10x Genomics TCR enrichment |
| `10x-vdj-bcr` | 10x Genomics BCR enrichment |

## Species Support

```bash
# Human
mixcr align -s human ...

# Mouse
mixcr align -s mmu ...

# Available: human, mmu, rat, rhesus, dog, pig, rabbit, chicken
```

## Output Format

Key columns in exported clones:

| Column | Description |
|--------|-------------|
| cloneId | Unique clone identifier |
| cloneCount | Number of reads |
| cloneFraction | Proportion of repertoire |
| nSeqCDR3 | Nucleotide CDR3 sequence |
| aaSeqCDR3 | Amino acid CDR3 sequence |
| allVHitsWithScore | V gene assignments |
| allDHitsWithScore | D gene assignments |
| allJHitsWithScore | J gene assignments |

## Quality Metrics

```bash
# Get alignment report
mixcr exportReports alignments.vdjca

# Key metrics:
# - Successfully aligned reads (>80% is good)
# - CDR3 found (>70% of aligned)
# - Clonotype count (varies by sample type)
```

## Parse MiXCR Output in Python

```python
import pandas as pd

def load_mixcr_clones(filepath):
    '''Load MiXCR clone table'''
    df = pd.read_csv(filepath, sep='\t')

    # Rename common columns
    df = df.rename(columns={
        'cloneCount': 'count',
        'cloneFraction': 'frequency',
        'aaSeqCDR3': 'cdr3_aa',
        'nSeqCDR3': 'cdr3_nt'
    })

    return df
```

## Related Skills

- **vdjtools-analysis** - Downstream diversity analysis
- **scirpy-analysis** - Single-cell VDJ integration
- **repertoire-visualization** - Visualize MiXCR output
