# featureCounts Usage Guide

featureCounts is a read counting program from the Subread package that assigns aligned reads to genomic features such as genes, exons, or promoters.

## Installation

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

```bash
# Basic gene counting
featureCounts -a genes.gtf -o counts.txt aligned.bam

# Paired-end, reverse-stranded (Illumina TruSeq)
featureCounts -p --countReadPairs -s 2 -T 4 -a genes.gtf -o counts.txt *.bam
```

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

1. **Always use `-p --countReadPairs` for paired-end** - Without this, each read in a pair is counted separately
2. **Match GTF to genome version** - Mismatched annotations cause low assignment rates
3. **Determine strandedness empirically** - Use RSeQC's `infer_experiment.py`
4. **Process all samples together** - Single run produces aligned count matrix
5. **Keep the summary file** - Essential for QC and troubleshooting
