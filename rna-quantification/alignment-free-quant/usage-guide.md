# Alignment-Free Quantification Usage Guide

Salmon and kallisto perform transcript quantification directly from FASTQ files without traditional alignment, making them much faster than alignment-based methods.

## Installation

```bash
# Salmon
conda install -c bioconda salmon

# kallisto
conda install -c bioconda kallisto
```

## Quick Start

### Salmon

```bash
# Build index (one-time)
salmon index -t transcripts.fa -i salmon_index

# Quantify
salmon quant -i salmon_index -l A \
    -1 reads_R1.fq.gz -2 reads_R2.fq.gz \
    -o output_quant
```

### kallisto

```bash
# Build index (one-time)
kallisto index -i kallisto_index transcripts.fa

# Quantify
kallisto quant -i kallisto_index -o output_quant \
    reads_R1.fq.gz reads_R2.fq.gz
```

## Obtaining Transcriptomes

### Ensembl (Recommended)

```bash
# Human
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Mouse
wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
```

### GENCODE

```bash
# Human
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
```

## Building Decoy-Aware Salmon Index

For improved accuracy, include genomic decoy sequences:

```bash
# Get genome and transcriptome
wget genome.fa.gz
wget transcripts.fa.gz
gunzip *.gz

# Extract chromosome names as decoys
grep "^>" genome.fa | cut -d " " -f 1 | sed 's/>//g' > decoys.txt

# Concatenate (transcriptome first!)
cat transcripts.fa genome.fa > gentrome.fa

# Build index
salmon index -t gentrome.fa -d decoys.txt -i salmon_index -p 8
```

## Batch Processing Script

```bash
#!/bin/bash
INDEX="salmon_index"
THREADS=8

for r1 in *_R1.fastq.gz; do
    sample=$(basename $r1 _R1.fastq.gz)
    r2="${sample}_R2.fastq.gz"

    salmon quant -i $INDEX -l A \
        -1 $r1 -2 $r2 \
        -o ${sample}_quant \
        -p $THREADS \
        --gcBias --seqBias
done
```

## Understanding Output

### Salmon quant.sf

| Column | Description |
|--------|-------------|
| Name | Transcript ID |
| Length | Transcript length |
| EffectiveLength | Length adjusted for bias |
| TPM | Transcripts per million |
| NumReads | Estimated read count |

### kallisto abundance.tsv

| Column | Description |
|--------|-------------|
| target_id | Transcript ID |
| length | Transcript length |
| eff_length | Effective length |
| est_counts | Estimated counts |
| tpm | Transcripts per million |

## TPM vs Counts

- **TPM** - Normalized, comparable across samples, use for visualization
- **Counts** - Use with tximport for DESeq2/edgeR (they need raw counts)

## Tips

1. **Use decoy-aware Salmon index** for best accuracy
2. **Enable bias correction** with `--gcBias --seqBias` in Salmon
3. **Generate bootstraps** (`-b 100`) if using sleuth for DE
4. **Check mapping rates** - should be >70%
5. **Match transcriptome version** to your GTF annotation
