# BWA-MEM2 Alignment - Usage Guide

## Overview
bwa-mem2 is the successor to BWA-MEM, providing 2-3x faster alignment with nearly identical results. It's the standard choice for DNA short-read alignment including WGS, WES, and ChIP-seq.

## Prerequisites
```bash
conda install -c bioconda bwa-mem2 samtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Align my paired-end WGS reads to the human reference genome"
- "Index this reference genome for BWA alignment"
- "Run BWA-MEM2 with read groups for GATK compatibility"

## Example Prompts

### Basic Alignment
> "Align reads_R1.fq.gz and reads_R2.fq.gz to reference.fa using BWA-MEM2"

> "Index my reference genome for BWA-MEM2 alignment"

### WGS Pipeline
> "Run a complete WGS alignment pipeline with duplicate marking for sample NA12878"

> "Align my WGS reads and add proper read groups for downstream GATK analysis"

### Performance Optimization
> "Align reads with reduced memory usage for my 16GB system"

> "Run BWA-MEM2 with reproducible results across different thread counts"

## What the Agent Will Do
1. Index the reference genome if not already indexed
2. Align reads with appropriate settings (threads, read groups)
3. Sort and index the output BAM file
4. Optionally mark duplicates and generate alignment statistics

## Tips
- Always add read groups (-R) for multi-sample analysis and GATK compatibility
- Pipe directly to samtools to avoid large intermediate SAM files
- Use `-K` flag for reproducible results across different thread counts
- Check reference genome version matches your reads if mapping rate is low
- Use SSD storage for reference and temporary files for best performance

## Read Group Fields

| Field | Description | Example |
|-------|-------------|---------|
| ID | Unique identifier | flowcell.lane |
| SM | Sample name | patient1 |
| PL | Platform | ILLUMINA |
| LB | Library | lib1 |
| PU | Platform unit | flowcell.lane.barcode |

## Complete WGS Pipeline Example

```bash
REF=reference.fa
R1=sample_R1.fq.gz
R2=sample_R2.fq.gz
SAMPLE=sample1
THREADS=16

if [ ! -f ${REF}.bwt.2bit.64 ]; then
    bwa-mem2 index $REF
fi

bwa-mem2 mem -t $THREADS \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1" \
    $REF $R1 $R2 | \
    samtools fixmate -m -@ 4 - - | \
    samtools sort -@ 4 -T /tmp/${SAMPLE} - | \
    samtools markdup -@ 4 - ${SAMPLE}.markdup.bam

samtools index ${SAMPLE}.markdup.bam
samtools flagstat ${SAMPLE}.markdup.bam > ${SAMPLE}.flagstat.txt
```

## Troubleshooting

### Low Mapping Rate
- Check reference genome version matches reads
- Check read quality with FastQC
- Try more sensitive settings (-k 15)

### Out of Memory
- Reduce threads (-t)
- Reduce batch size (-K 10000000)
- Use original BWA-MEM as fallback
