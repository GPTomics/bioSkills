#!/bin/bash
# ONT structural variant calling with Sniffles2
set -e

THREADS=16
READS="nanopore_reads.fastq.gz"
REF="reference.fa"
SAMPLE="sample1"
OUTDIR="sv_results"

mkdir -p ${OUTDIR}/{qc,aligned,sv,annotation}

echo "=== ONT SV Calling Pipeline ==="

# Step 1: QC
echo "=== Step 1: Quality Control ==="
NanoPlot \
    --fastq ${READS} \
    --outdir ${OUTDIR}/qc \
    --threads ${THREADS} \
    --plots hex dot

echo "QC complete. Check ${OUTDIR}/qc/NanoStats.txt"

# Step 2: Alignment
echo "=== Step 2: Alignment ==="
minimap2 -ax map-ont \
    -t ${THREADS} \
    --MD \
    -Y \
    ${REF} \
    ${READS} | \
samtools sort -@ 4 -o ${OUTDIR}/aligned/${SAMPLE}.bam

samtools index ${OUTDIR}/aligned/${SAMPLE}.bam

echo "Alignment stats:"
samtools flagstat ${OUTDIR}/aligned/${SAMPLE}.bam | head -5

# Calculate coverage
avg_cov=$(samtools depth -a ${OUTDIR}/aligned/${SAMPLE}.bam | \
    awk '{sum+=$3; n++} END {printf "%.1f", sum/n}')
echo "Average coverage: ${avg_cov}x"

# Step 3: SV Calling
echo "=== Step 3: SV Calling ==="
sniffles \
    --input ${OUTDIR}/aligned/${SAMPLE}.bam \
    --vcf ${OUTDIR}/sv/${SAMPLE}.raw.vcf.gz \
    --reference ${REF} \
    --threads ${THREADS} \
    --minsvlen 50 \
    --output-rnames

# Step 4: Filtering
echo "=== Step 4: Filtering ==="

# Basic quality filter
bcftools view -i 'QUAL>=20' \
    ${OUTDIR}/sv/${SAMPLE}.raw.vcf.gz \
    -Oz -o ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz

bcftools index ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz

# Create separate files by SV type
for svtype in DEL INS DUP INV; do
    bcftools view -i "SVTYPE=\"${svtype}\"" \
        ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz \
        -Oz -o ${OUTDIR}/sv/${SAMPLE}.${svtype}.vcf.gz
done

# Step 5: Statistics
echo "=== Step 5: Statistics ==="
bcftools stats ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz > ${OUTDIR}/sv/stats.txt

echo ""
echo "=== Summary ==="
echo "Total SVs: $(bcftools view -H ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz | wc -l)"
echo "  DEL: $(bcftools view -H ${OUTDIR}/sv/${SAMPLE}.DEL.vcf.gz | wc -l)"
echo "  INS: $(bcftools view -H ${OUTDIR}/sv/${SAMPLE}.INS.vcf.gz | wc -l)"
echo "  DUP: $(bcftools view -H ${OUTDIR}/sv/${SAMPLE}.DUP.vcf.gz | wc -l)"
echo "  INV: $(bcftools view -H ${OUTDIR}/sv/${SAMPLE}.INV.vcf.gz | wc -l)"
echo ""
echo "Results: ${OUTDIR}/sv/${SAMPLE}.filtered.vcf.gz"
