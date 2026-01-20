#!/bin/bash
# Complete RNA-seq workflow: STAR + featureCounts + DESeq2
set -e

# Configuration
THREADS=8
GENOME_DIR="star_index"
GTF="genes.gtf"
SAMPLES="control_1 control_2 control_3 treated_1 treated_2 treated_3"
OUTDIR="results"

mkdir -p ${OUTDIR}/{trimmed,aligned,counts,qc}

echo "=== Step 1: Quality Control with fastp ==="
for sample in $SAMPLES; do
    echo "Processing ${sample}..."
    fastp \
        -i ${sample}_R1.fastq.gz \
        -I ${sample}_R2.fastq.gz \
        -o ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -O ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 35 \
        --html ${OUTDIR}/qc/${sample}_fastp.html \
        --json ${OUTDIR}/qc/${sample}_fastp.json \
        -w ${THREADS}
done

echo "=== Step 2: STAR Alignment ==="
for sample in $SAMPLES; do
    echo "Aligning ${sample}..."
    STAR \
        --genomeDir ${GENOME_DIR} \
        --readFilesIn ${OUTDIR}/trimmed/${sample}_R1.fq.gz ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix ${OUTDIR}/aligned/${sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --quantMode GeneCounts \
        --runThreadN ${THREADS}

    # Index BAM
    samtools index ${OUTDIR}/aligned/${sample}_Aligned.sortedByCoord.out.bam
done

echo "=== Step 3: featureCounts ==="
featureCounts \
    -T ${THREADS} \
    -p --countReadPairs \
    -s 2 \
    -a ${GTF} \
    -o ${OUTDIR}/counts/counts.txt \
    ${OUTDIR}/aligned/*_Aligned.sortedByCoord.out.bam

# Clean up column names in counts file
head -1 ${OUTDIR}/counts/counts.txt > ${OUTDIR}/counts/counts_clean.txt
tail -n +2 ${OUTDIR}/counts/counts.txt | \
    sed 's/_Aligned.sortedByCoord.out.bam//g' | \
    sed "s|${OUTDIR}/aligned/||g" >> ${OUTDIR}/counts/counts_clean.txt

echo "=== Step 4: Run DESeq2 (R script) ==="
cat << 'RSCRIPT' > ${OUTDIR}/run_deseq2.R
library(DESeq2)
library(apeglm)
library(ggplot2)

# Load counts
counts <- read.table('counts/counts_clean.txt', header = TRUE, row.names = 1, skip = 1)
counts <- counts[, 6:ncol(counts)]  # Remove annotation columns

# Sample metadata
samples <- colnames(counts)
conditions <- ifelse(grepl('control', samples), 'control', 'treated')
coldata <- data.frame(condition = factor(conditions), row.names = samples)

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10,]
dds$condition <- relevel(dds$condition, ref = 'control')
dds <- DESeq(dds)

res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
summary(res)

sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig), 'significant_genes.csv')
cat('Significant genes:', nrow(sig), '\n')
RSCRIPT

cd ${OUTDIR} && Rscript run_deseq2.R

echo "=== Pipeline Complete ==="
echo "Results in ${OUTDIR}/"
echo "  - QC reports: qc/"
echo "  - Alignments: aligned/"
echo "  - Counts: counts/"
echo "  - DE results: significant_genes.csv"
