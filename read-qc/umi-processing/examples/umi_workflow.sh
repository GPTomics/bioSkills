#!/bin/bash
# Reference: umi_tools 1.1+, STAR 2.7+, samtools 1.19+, Subread 2.0+ | Verify API if version differs
# UMI processing workflow for UMI-tagged RNA-seq. Order is mandatory: extract UMI (FASTQ) ->
# align -> dedup (dedup needs mapping coordinates). Only run this on libraries that HAVE UMIs;
# standard non-UMI bulk RNA-seq must NOT be deduplicated.
set -euo pipefail

SAMPLE=$1
R1=${SAMPLE}_R1.fastq.gz
R2=${SAMPLE}_R2.fastq.gz
STAR_INDEX=$2
GTF=$3

echo "=== UMI Processing Pipeline for ${SAMPLE} ==="

# Step 1: Extract UMIs (8bp at start of R1)
echo "Extracting UMIs..."
umi_tools extract \
    --stdin=${R1} \
    --read2-in=${R2} \
    --stdout=${SAMPLE}_R1_umi.fastq.gz \
    --read2-out=${SAMPLE}_R2_umi.fastq.gz \
    --bc-pattern=NNNNNNNN \
    --log=${SAMPLE}_extract.log

# Step 2: Align with STAR
echo "Aligning with STAR..."
STAR --runThreadN 8 \
    --genomeDir ${STAR_INDEX} \
    --readFilesIn ${SAMPLE}_R1_umi.fastq.gz ${SAMPLE}_R2_umi.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix ${SAMPLE}_ \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbGTFfile ${GTF} \
    --quantMode GeneCounts

# Step 3: Index BAM
echo "Indexing BAM..."
samtools index ${SAMPLE}_Aligned.sortedByCoord.out.bam

# Step 4: Deduplicate
echo "Deduplicating..."
umi_tools dedup \
    -I ${SAMPLE}_Aligned.sortedByCoord.out.bam \
    -S ${SAMPLE}_deduplicated.bam \
    --paired \
    --output-stats=${SAMPLE}_dedup \
    --log=${SAMPLE}_dedup.log

# Step 5: Index deduplicated BAM
samtools index ${SAMPLE}_deduplicated.bam

# Step 6: Count features. Since Subread 2.0.2, -p only declares paired-end input; counting
# FRAGMENTS (not reads) also requires --countReadPairs, else properly-paired fragments count ~2x.
echo "Counting features..."
featureCounts -T 8 -p --countReadPairs \
    -a ${GTF} \
    -o ${SAMPLE}_counts.txt \
    ${SAMPLE}_deduplicated.bam

# Report stats
echo "=== Statistics ==="
echo "Input reads: $(zcat ${R1} | wc -l | awk '{print $1/4}')"
echo "Aligned reads: $(samtools view -c ${SAMPLE}_Aligned.sortedByCoord.out.bam)"
echo "Deduplicated reads: $(samtools view -c ${SAMPLE}_deduplicated.bam)"

echo "Done: ${SAMPLE}_deduplicated.bam"
