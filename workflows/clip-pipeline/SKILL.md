---
name: bio-workflows-clip-pipeline
description: End-to-end CLIP-seq analysis from FASTQ to binding sites and motif enrichment. Use when analyzing protein-RNA interactions from CLIP-based methods.
tool_type: mixed
primary_tool: CLIPper
---

# CLIP-seq Pipeline

## Pipeline Overview

```
FASTQ → UMI extract → Trim → Align → Dedup → Peak call → Annotate → Motifs
```

## Step 1: Preprocessing

```bash
# Extract UMIs (if present)
umi_tools extract --stdin=reads.fastq.gz \
    --bc-pattern=NNNNNNNNNN \
    --stdout=extracted.fastq.gz

# Trim adapters
cutadapt -a AGATCGGAAGAG --minimum-length 20 \
    -o trimmed.fastq.gz extracted.fastq.gz
```

## Step 2: Alignment

```bash
STAR --genomeDir index \
    --readFilesIn trimmed.fastq.gz \
    --readFilesCommand zcat \
    --outFilterMismatchNmax 2 \
    --outSAMtype BAM SortedByCoordinate
```

## Step 3: Deduplication

```bash
umi_tools dedup -I aligned.bam -S dedup.bam \
    --output-stats=dedup_stats
```

## Step 4: Peak Calling

```bash
clipper -b dedup.bam -s hg38 -o peaks.bed
```

## Step 5: Annotation and Motifs

```bash
# Annotate peaks
bedtools intersect -a peaks.bed -b genes.gtf -wo > annotated.bed

# Motif analysis
bedtools getfasta -fi genome.fa -bed peaks.bed -fo peaks.fa
findMotifs.pl peaks.fa fasta motif_output -rna
```

## QC Checkpoints

1. **After dedup**: Check UMI duplication rate
2. **After peaks**: Verify peak count and width distribution
3. **After motifs**: Known RBP motif should be top hit

## Related Skills

- **clip-seq/** - Individual CLIP analysis skills
- **chip-seq/peak-calling** - Similar peak concepts
