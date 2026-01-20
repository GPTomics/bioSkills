# Alignment-Free Quantification - Usage Guide

## Overview
Salmon and kallisto perform transcript quantification directly from FASTQ files without traditional alignment, making them much faster than alignment-based methods.

## Prerequisites
```bash
# Salmon
conda install -c bioconda salmon

# kallisto
conda install -c bioconda kallisto
```

## Quick Start
Tell your AI agent what you want to do:
- "Quantify my RNA-seq samples using Salmon"
- "Build a Salmon index with decoy sequences"
- "Run kallisto on my paired-end FASTQ files"

## Example Prompts
### Index Building
> "Build a decoy-aware Salmon index from the human transcriptome and genome"

> "Create a kallisto index from my transcripts.fa file"

### Quantification
> "Quantify paired-end reads for sample1 using Salmon with bias correction"

> "Run Salmon quant on all my FASTQ files in batch mode"

### Output Interpretation
> "Explain the columns in the Salmon quant.sf output file"

> "What's the difference between TPM and NumReads in Salmon output?"

## What the Agent Will Do
1. Download or locate the reference transcriptome (Ensembl/GENCODE)
2. Build the index (optionally with decoy sequences for Salmon)
3. Run quantification on each sample with appropriate parameters
4. Check mapping rates and quality metrics
5. Organize output files for downstream analysis with tximport

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
- Use decoy-aware Salmon index for best accuracy
- Enable bias correction with `--gcBias --seqBias` in Salmon
- Generate bootstraps (`-b 100`) if using sleuth for DE
- Check mapping rates - should be >70%
- Match transcriptome version to your GTF annotation
