#!/bin/bash
# Complete variant calling workflow: BWA-MEM2 + bcftools
set -e

# Configuration
THREADS=8
REF="reference.fa"
SAMPLES="sample1 sample2 sample3"
OUTDIR="results"

# Create directories
mkdir -p ${OUTDIR}/{trimmed,aligned,variants,qc}

echo "=== Variant Calling Pipeline ==="
echo "Reference: ${REF}"
echo "Samples: ${SAMPLES}"
echo "Threads: ${THREADS}"
echo ""

# Check reference is indexed
if [ ! -f "${REF}.bwt.2bit.64" ]; then
    echo "Indexing reference for bwa-mem2..."
    bwa-mem2 index ${REF}
fi

if [ ! -f "${REF}.fai" ]; then
    echo "Indexing reference for samtools..."
    samtools faidx ${REF}
fi

# Step 1: Quality Control
echo "=== Step 1: Quality Control ==="
for sample in $SAMPLES; do
    if [ ! -f "${OUTDIR}/trimmed/${sample}_R1.fq.gz" ]; then
        echo "QC: ${sample}"
        fastp \
            -i ${sample}_R1.fastq.gz \
            -I ${sample}_R2.fastq.gz \
            -o ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
            -O ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 50 \
            --html ${OUTDIR}/qc/${sample}_fastp.html \
            --json ${OUTDIR}/qc/${sample}_fastp.json \
            -w ${THREADS}
    else
        echo "Skipping QC for ${sample} (already done)"
    fi
done

# Step 2: Alignment
echo "=== Step 2: Alignment ==="
for sample in $SAMPLES; do
    if [ ! -f "${OUTDIR}/aligned/${sample}.markdup.bam" ]; then
        echo "Aligning: ${sample}"
        bwa-mem2 mem -t ${THREADS} \
            -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:lib1" \
            ${REF} \
            ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
            ${OUTDIR}/trimmed/${sample}_R2.fq.gz | \
        samtools view -@ ${THREADS} -bS - | \
        samtools fixmate -@ ${THREADS} -m - - | \
        samtools sort -@ ${THREADS} - | \
        samtools markdup -@ ${THREADS} - ${OUTDIR}/aligned/${sample}.markdup.bam

        samtools index ${OUTDIR}/aligned/${sample}.markdup.bam

        echo "Stats for ${sample}:"
        samtools flagstat ${OUTDIR}/aligned/${sample}.markdup.bam
        echo ""
    else
        echo "Skipping alignment for ${sample} (already done)"
    fi
done

# Step 3: Variant Calling
echo "=== Step 3: Variant Calling ==="
if [ ! -f "${OUTDIR}/variants/cohort.vcf.gz" ]; then
    echo "Calling variants jointly..."
    bcftools mpileup -Ou \
        -f ${REF} \
        --max-depth 250 \
        --min-MQ 20 \
        --min-BQ 20 \
        ${OUTDIR}/aligned/*.markdup.bam | \
    bcftools call -mv -Oz -o ${OUTDIR}/variants/cohort.vcf.gz

    bcftools index ${OUTDIR}/variants/cohort.vcf.gz
else
    echo "Skipping variant calling (already done)"
fi

# Step 4: Filtering
echo "=== Step 4: Filtering ==="
if [ ! -f "${OUTDIR}/variants/cohort.filtered.vcf.gz" ]; then
    echo "Applying quality filters..."
    bcftools filter -Oz \
        -e 'QUAL<20 || DP<10 || MQ<30' \
        -o ${OUTDIR}/variants/cohort.filtered.vcf.gz \
        ${OUTDIR}/variants/cohort.vcf.gz

    bcftools index ${OUTDIR}/variants/cohort.filtered.vcf.gz
else
    echo "Skipping filtering (already done)"
fi

# Step 5: Statistics
echo "=== Step 5: Statistics ==="
bcftools stats ${OUTDIR}/variants/cohort.filtered.vcf.gz > ${OUTDIR}/variants/stats.txt

# Summary
echo ""
echo "=== Pipeline Complete ==="
echo ""
echo "Results:"
echo "  QC reports: ${OUTDIR}/qc/"
echo "  Alignments: ${OUTDIR}/aligned/"
echo "  Raw VCF: ${OUTDIR}/variants/cohort.vcf.gz"
echo "  Filtered VCF: ${OUTDIR}/variants/cohort.filtered.vcf.gz"
echo "  VCF stats: ${OUTDIR}/variants/stats.txt"
echo ""

# Quick stats summary
echo "Variant Summary:"
bcftools stats ${OUTDIR}/variants/cohort.filtered.vcf.gz | grep -E "^SN"
echo ""
echo "Ti/Tv ratio:"
bcftools stats ${OUTDIR}/variants/cohort.filtered.vcf.gz | grep "TSTV"
