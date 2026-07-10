#!/bin/bash
# Reference: BWA-MEM2 2.2.1+, GATK 4.5+, bcftools 1.19+, fastp 0.23+, samtools 1.19+ | Verify API if version differs
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
if [ ! -f "${REF%.*}.dict" ]; then
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
    # No `samtools collate` needed: fresh bwa-mem2 output is already name-grouped for fixmate -m.

    samtools index ${OUTDIR}/aligned/${sample}.markdup.bam
done

# Step 3: Base Quality Score Recalibration
# BQSR is the classic path shown here; it is honestly optional on modern binned-quality
# instruments (2-color NovaSeq). The alternative is HaplotypeCaller --dragen-mode (DRAGSTR),
# which models STR indel error internally and skips a separate BQSR step. See the SKILL body.
echo "=== Step 3: BQSR ==="
for sample in $SAMPLES; do
    echo "BQSR: ${sample}"

    # GATK best practice supplies dbSNP PLUS known indels (Mills + 1000G gold standard) as
    # --known-sites; dbSNP-only here is a simplification. Add: --known-sites ${MILLS_INDELS}
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

# Create sample map (`:` is a no-op command; `> file` alone truncates but reads as a stray redirect)
: > ${OUTDIR}/gvcf/sample_map.txt
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

# GATK canonical hard-filter starting points (lenient by design; SNP and indel thresholds differ
# because their error processes differ -- see variant-calling/filtering-best-practices). SNPs:
# QD<2 low confidence-per-depth, FS>60 strand bias, MQ<40 poor mapping, SOR>3 strand-odds bias.
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
    -O ${OUTDIR}/variants/cohort.merged.vcf.gz

# Normalize (left-align + split multiallelics) so the delivered VCF is the advertised NORMALIZED
# artifact (governing principle #2 / canonical step 6), before any annotation or benchmarking.
bcftools norm -m-any -f ${REF} -Oz \
    -o ${OUTDIR}/variants/cohort.norm.vcf.gz \
    ${OUTDIR}/variants/cohort.merged.vcf.gz

# Genotype-level pass after the site-level VariantFiltration (canonical step 7: site THEN genotype);
# null out low-confidence genotypes (-S .) without dropping the site.
bcftools filter -Oz -S . -e 'FMT/DP<8 | FMT/GQ<20' \
    -o ${OUTDIR}/variants/cohort.filtered.vcf.gz \
    ${OUTDIR}/variants/cohort.norm.vcf.gz
bcftools index ${OUTDIR}/variants/cohort.filtered.vcf.gz

echo "=== Pipeline Complete ==="
echo "Filtered VCF: ${OUTDIR}/variants/cohort.filtered.vcf.gz"
