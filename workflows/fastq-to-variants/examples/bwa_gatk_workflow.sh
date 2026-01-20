#!/bin/bash
# Complete variant calling workflow: BWA-MEM2 + GATK HaplotypeCaller
set -e

# Configuration
THREADS=8
REF="reference.fa"
DBSNP="dbsnp.vcf.gz"
SAMPLES="sample1 sample2 sample3"
OUTDIR="results"

mkdir -p ${OUTDIR}/{trimmed,aligned,recal,gvcf,variants,qc}

echo "=== GATK Variant Calling Pipeline ==="
echo "Reference: ${REF}"
echo "dbSNP: ${DBSNP}"
echo "Samples: ${SAMPLES}"
echo ""

# Check reference files
if [ ! -f "${REF}.dict" ]; then
    echo "Creating sequence dictionary..."
    gatk CreateSequenceDictionary -R ${REF}
fi

if [ ! -f "${REF}.fai" ]; then
    echo "Indexing reference..."
    samtools faidx ${REF}
fi

if [ ! -f "${REF}.bwt.2bit.64" ]; then
    echo "Indexing reference for bwa-mem2..."
    bwa-mem2 index ${REF}
fi

# Step 1: QC
echo "=== Step 1: Quality Control ==="
for sample in $SAMPLES; do
    fastp \
        -i ${sample}_R1.fastq.gz \
        -I ${sample}_R2.fastq.gz \
        -o ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -O ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --detect_adapter_for_pe \
        --html ${OUTDIR}/qc/${sample}_fastp.html \
        -w ${THREADS}
done

# Step 2: Alignment
echo "=== Step 2: Alignment ==="
for sample in $SAMPLES; do
    bwa-mem2 mem -t ${THREADS} \
        -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
        ${REF} \
        ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        ${OUTDIR}/trimmed/${sample}_R2.fq.gz | \
    samtools view -@ ${THREADS} -bS - | \
    samtools fixmate -@ ${THREADS} -m - - | \
    samtools sort -@ ${THREADS} - | \
    samtools markdup -@ ${THREADS} - ${OUTDIR}/aligned/${sample}.markdup.bam

    samtools index ${OUTDIR}/aligned/${sample}.markdup.bam
done

# Step 3: Base Quality Score Recalibration
echo "=== Step 3: BQSR ==="
for sample in $SAMPLES; do
    echo "BQSR: ${sample}"

    gatk BaseRecalibrator \
        -R ${REF} \
        -I ${OUTDIR}/aligned/${sample}.markdup.bam \
        --known-sites ${DBSNP} \
        -O ${OUTDIR}/recal/${sample}_recal.table

    gatk ApplyBQSR \
        -R ${REF} \
        -I ${OUTDIR}/aligned/${sample}.markdup.bam \
        --bqsr-recal-file ${OUTDIR}/recal/${sample}_recal.table \
        -O ${OUTDIR}/recal/${sample}.recal.bam
done

# Step 4: HaplotypeCaller (GVCF mode)
echo "=== Step 4: HaplotypeCaller ==="
for sample in $SAMPLES; do
    echo "HaplotypeCaller: ${sample}"

    gatk HaplotypeCaller \
        -R ${REF} \
        -I ${OUTDIR}/recal/${sample}.recal.bam \
        -O ${OUTDIR}/gvcf/${sample}.g.vcf.gz \
        -ERC GVCF \
        --native-pair-hmm-threads ${THREADS}
done

# Step 5: Joint Genotyping
echo "=== Step 5: Joint Genotyping ==="

# Create sample map
> ${OUTDIR}/gvcf/sample_map.txt
for sample in $SAMPLES; do
    echo -e "${sample}\t${OUTDIR}/gvcf/${sample}.g.vcf.gz" >> ${OUTDIR}/gvcf/sample_map.txt
done

# GenomicsDBImport
gatk GenomicsDBImport \
    --sample-name-map ${OUTDIR}/gvcf/sample_map.txt \
    --genomicsdb-workspace-path ${OUTDIR}/genomicsdb \
    -L chr1 -L chr2 -L chr3  # Add all chromosomes or use interval list

# GenotypeGVCFs
gatk GenotypeGVCFs \
    -R ${REF} \
    -V gendb://${OUTDIR}/genomicsdb \
    -O ${OUTDIR}/variants/cohort.vcf.gz

# Step 6: Hard Filtering (for small cohorts; use VQSR for >30 samples)
echo "=== Step 6: Filtering ==="

# Filter SNPs
gatk SelectVariants \
    -V ${OUTDIR}/variants/cohort.vcf.gz \
    -select-type SNP \
    -O ${OUTDIR}/variants/cohort.snps.vcf.gz

gatk VariantFiltration \
    -R ${REF} \
    -V ${OUTDIR}/variants/cohort.snps.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "HighFS" \
    --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
    --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
    -O ${OUTDIR}/variants/cohort.snps.filtered.vcf.gz

# Filter Indels
gatk SelectVariants \
    -V ${OUTDIR}/variants/cohort.vcf.gz \
    -select-type INDEL \
    -O ${OUTDIR}/variants/cohort.indels.vcf.gz

gatk VariantFiltration \
    -R ${REF} \
    -V ${OUTDIR}/variants/cohort.indels.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "FS > 200.0" --filter-name "HighFS" \
    --filter-expression "SOR > 10.0" --filter-name "HighSOR" \
    -O ${OUTDIR}/variants/cohort.indels.filtered.vcf.gz

# Merge filtered variants
gatk MergeVcfs \
    -I ${OUTDIR}/variants/cohort.snps.filtered.vcf.gz \
    -I ${OUTDIR}/variants/cohort.indels.filtered.vcf.gz \
    -O ${OUTDIR}/variants/cohort.filtered.vcf.gz

echo "=== Pipeline Complete ==="
echo "Filtered VCF: ${OUTDIR}/variants/cohort.filtered.vcf.gz"
