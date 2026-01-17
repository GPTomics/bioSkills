# BWA-MEM2 Alignment Usage Guide

## Overview

bwa-mem2 is the successor to BWA-MEM, providing 2-3x faster alignment with nearly identical results. It's the standard choice for DNA short-read alignment including WGS, WES, and ChIP-seq.

## Installation

```bash
conda install -c bioconda bwa-mem2 samtools
```

## Quick Start

```bash
# 1. Index reference (once)
bwa-mem2 index reference.fa

# 2. Align paired-end reads
bwa-mem2 mem -t 8 reference.fa reads_R1.fq.gz reads_R2.fq.gz | \
    samtools sort -o aligned.bam -
samtools index aligned.bam
```

## Complete WGS Pipeline

```bash
# Variables
REF=reference.fa
R1=sample_R1.fq.gz
R2=sample_R2.fq.gz
SAMPLE=sample1
THREADS=16

# Index if needed
if [ ! -f ${REF}.bwt.2bit.64 ]; then
    bwa-mem2 index $REF
fi

# Align, sort, mark duplicates
bwa-mem2 mem -t $THREADS \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1" \
    $REF $R1 $R2 | \
    samtools fixmate -m -@ 4 - - | \
    samtools sort -@ 4 -T /tmp/${SAMPLE} - | \
    samtools markdup -@ 4 - ${SAMPLE}.markdup.bam

samtools index ${SAMPLE}.markdup.bam

# Alignment stats
samtools flagstat ${SAMPLE}.markdup.bam > ${SAMPLE}.flagstat.txt
```

## Read Groups

Read groups are essential for multi-sample analysis and GATK compatibility.

| Field | Description | Example |
|-------|-------------|---------|
| ID | Unique identifier | flowcell.lane |
| SM | Sample name | patient1 |
| PL | Platform | ILLUMINA |
| LB | Library | lib1 |
| PU | Platform unit | flowcell.lane.barcode |

```bash
# Full read group string
-R '@RG\tID:H0164.2\tSM:NA12878\tPL:ILLUMINA\tLB:Solexa-272222\tPU:H0164ALXX140820.2'
```

## Alignment Modes

### Default (Global)
Standard alignment allowing gaps:
```bash
bwa-mem2 mem reference.fa reads.fq
```

### Local-like Behavior
Soft clip low-quality ends:
```bash
bwa-mem2 mem -Y reference.fa reads.fq
```

### High Sensitivity
Lower seed length for divergent sequences:
```bash
bwa-mem2 mem -k 15 reference.fa reads.fq
```

## Memory Optimization

```bash
# Reduce memory with smaller batch size
bwa-mem2 mem -K 10000000 -t 4 reference.fa r1.fq r2.fq

# For limited memory systems, use fewer threads
bwa-mem2 mem -t 4 reference.fa r1.fq r2.fq
```

## Troubleshooting

### Index Issues
```bash
# If index incomplete, remove and rebuild
rm reference.fa.0123 reference.fa.amb reference.fa.ann reference.fa.bwt.2bit.64 reference.fa.pac
bwa-mem2 index reference.fa
```

### Out of Memory
- Reduce threads (-t)
- Reduce batch size (-K)
- Use original BWA-MEM as fallback

### Low Mapping Rate
- Check reference genome version matches reads
- Check read quality with FastQC
- Try more sensitive settings (-k 15)

## Performance Tips

1. Use SSD storage for reference and temporary files
2. Match thread count to available cores
3. Pipe directly to samtools to avoid intermediate SAM files
4. Use `-K` for reproducible results across different thread counts
