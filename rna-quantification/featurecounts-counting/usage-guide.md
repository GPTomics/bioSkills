# featureCounts - Usage Guide

## Overview
featureCounts is a read counting program from the Subread package that assigns aligned reads to genomic features such as genes, exons, or promoters.

## Prerequisites
```bash
# Conda (recommended)
conda install -c bioconda subread

# Ubuntu/Debian
sudo apt install subread

# From source
wget https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz
tar xzf subread-2.0.6-Linux-x86_64.tar.gz
export PATH=$PATH:$(pwd)/subread-2.0.6-Linux-x86_64/bin
```

## Quick Start
Tell your AI agent what you want to do:
- "Count reads in my BAM files using featureCounts"
- "Generate a gene count matrix from my aligned samples"
- "Run featureCounts on stranded paired-end data"

## Example Prompts
### Basic Counting
> "Run featureCounts on my aligned BAM files using the GRCh38 GTF annotation"

> "Count reads at the gene level for all samples in my data/ directory"

### Strand-Specific
> "Count reads with reverse-stranded library prep (TruSeq Stranded)"

> "Help me determine the strandedness of my RNA-seq data"

### Advanced Options
> "Count at the transcript level instead of gene level"

> "Run featureCounts allowing multi-mapping reads"

## What the Agent Will Do
1. Identify the input BAM files and annotation GTF
2. Determine library type (paired/single-end, stranded/unstranded)
3. Run featureCounts with appropriate parameters
4. Check the summary file for assignment rates
5. Extract the count matrix for downstream analysis

## Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-a` | Annotation file (GTF/GFF/SAF) | Required |
| `-o` | Output file | Required |
| `-t` | Feature type to count | exon |
| `-g` | Attribute for grouping | gene_id |
| `-s` | Strandedness (0/1/2) | 0 (unstranded) |
| `-p` | Paired-end mode | Off |
| `--countReadPairs` | Count fragments not reads | Off |
| `-T` | Number of threads | 1 |
| `-M` | Count multi-mappers | Off |
| `-O` | Count overlapping features | Off |

## Strandedness Guide

| Library Type | `-s` Value | Examples |
|--------------|------------|----------|
| Unstranded | 0 | Standard Illumina |
| Forward | 1 | Directional, some dUTP |
| Reverse | 2 | TruSeq Stranded, dUTP |

## Typical Workflows

### Standard RNA-seq (Illumina TruSeq Stranded)
```bash
featureCounts \
    -p --countReadPairs \
    -s 2 \
    -T 8 \
    -a Homo_sapiens.GRCh38.gtf \
    -o gene_counts.txt \
    sample1.bam sample2.bam sample3.bam
```

### Single-End Unstranded
```bash
featureCounts \
    -s 0 \
    -T 8 \
    -a annotation.gtf \
    -o gene_counts.txt \
    *.bam
```

### Transcript-Level Counting
```bash
featureCounts \
    -p --countReadPairs \
    -t exon \
    -g transcript_id \
    -O \
    -a annotation.gtf \
    -o transcript_counts.txt \
    *.bam
```

## Understanding Output

The main output file has these columns:
1. **Geneid** - Gene identifier from GTF
2. **Chr** - Chromosome(s) for the gene
3. **Start** - Start position(s)
4. **End** - End position(s)
5. **Strand** - Strand(s)
6. **Length** - Total exon length
7+ **Sample columns** - Raw counts per sample

## Quality Metrics

Check the `.summary` file for:
- **Assigned** - Reads successfully counted (target: >70%)
- **Unassigned_NoFeatures** - Reads not overlapping any feature
- **Unassigned_Ambiguity** - Reads overlapping multiple features
- **Unassigned_MultiMapping** - Multi-mapped reads (if not using `-M`)

## Tips
- Always use `-p --countReadPairs` for paired-end data - without this, each read in a pair is counted separately
- Match GTF to genome version - mismatched annotations cause low assignment rates
- Determine strandedness empirically using RSeQC's `infer_experiment.py`
- Process all samples together in a single run to produce an aligned count matrix
- Keep the summary file - it's essential for QC and troubleshooting
