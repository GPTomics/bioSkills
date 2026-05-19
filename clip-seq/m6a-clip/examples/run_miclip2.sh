#!/bin/bash
# Reference: miCLIP2 pipeline (Kortel 2021), m6Aboost 1.0+, PureCLIP 1.3.1+, samtools 1.19+ | Verify API if version differs
# miCLIP2 m6A detection pipeline: eCLIP-style processing + PureCLIP single-nt + m6Aboost ML classifier.
# Mettl3-KO calibration recommended; antibody false-positive rate is 30-50% without ML scoring.

R1=$1
R2=$2
STAR_INDEX=$3
GENOME=$4
SMINPUT_BAM=$5
OUT_PREFIX=${6:-"miclip2"}
THREADS=${7:-8}

# Step 1: eCLIP-style preprocessing (see clip-seq/clip-preprocessing for details)
umi_tools extract \
    --bc-pattern=NNNNNNNNNN \
    --stdin=$R1 --read2-in=$R2 \
    --stdout=${OUT_PREFIX}_R1.umi.fq.gz \
    --read2-out=${OUT_PREFIX}_R2.umi.fq.gz

cutadapt \
    -a AGATCGGAAGAGCACACGTCT \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 6 -m 18 -j $THREADS \
    -o ${OUT_PREFIX}_R1.trim.fq.gz \
    -p ${OUT_PREFIX}_R2.trim.fq.gz \
    ${OUT_PREFIX}_R1.umi.fq.gz \
    ${OUT_PREFIX}_R2.umi.fq.gz

# Step 2: STAR alignment (eCLIP parameter block)
STAR --runMode alignReads \
    --runThreadN $THREADS \
    --genomeDir $STAR_INDEX \
    --readFilesIn ${OUT_PREFIX}_R1.trim.fq.gz ${OUT_PREFIX}_R2.trim.fq.gz \
    --readFilesCommand zcat \
    --outFilterType BySJout \
    --outFilterMultimapNmax 1 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterScoreMinOverLread 0.66 \
    --outFilterMatchNminOverLread 0.66 \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUT_PREFIX}_

samtools index ${OUT_PREFIX}_Aligned.sortedByCoord.out.bam

# UMI dedup
umi_tools dedup \
    --stdin=${OUT_PREFIX}_Aligned.sortedByCoord.out.bam \
    --stdout=${OUT_PREFIX}_dedup.bam \
    --method=unique --paired
samtools index ${OUT_PREFIX}_dedup.bam

# Step 3: Single-nt CL site detection with PureCLIP
# m6A-specific HMM is not currently implemented; use general PureCLIP and post-filter
pureclip \
    -i ${OUT_PREFIX}_dedup.bam -bai ${OUT_PREFIX}_dedup.bam.bai \
    -g $GENOME \
    -ibam $SMINPUT_BAM -ibai ${SMINPUT_BAM}.bai \
    -o ${OUT_PREFIX}_sites.bed \
    -or ${OUT_PREFIX}_regions.bed \
    -nt $THREADS -dm 8

echo "PureCLIP sites: $(wc -l < ${OUT_PREFIX}_sites.bed)"

# Step 4: m6Aboost ML scoring
# Filters antibody false positives; trained on Mettl3-KO calibration data
# Without m6Aboost, miCLIP2 site list has 30-50% false positive rate
if command -v m6aboost > /dev/null; then
    m6aboost \
        --sites ${OUT_PREFIX}_sites.bed \
        --bam ${OUT_PREFIX}_dedup.bam \
        --genome $GENOME \
        --output ${OUT_PREFIX}_m6aboost.bed

    # Filter at m6Aboost score >= 0.5
    awk '$5 >= 0.5' ${OUT_PREFIX}_m6aboost.bed > ${OUT_PREFIX}_m6a_high.bed
    echo "High-confidence m6A sites (m6Aboost >= 0.5): $(wc -l < ${OUT_PREFIX}_m6a_high.bed)"
else
    echo "m6Aboost not installed; install from https://github.com/ZarnackGroup/m6Aboost"
    echo "Without ML filtering, expect 30-50% false-positive rate"
fi

# Step 5: DRACH context fraction (informational; 70-90% expected)
# Generate DRACH motif BED from genome first (use HOMER or custom regex), then:
# bedtools intersect -wa -u -s -a ${OUT_PREFIX}_m6a_high.bed -b drach_motifs.bed > ${OUT_PREFIX}_m6a_drach.bed
# DRACH_FRAC=$(echo "scale=4; $(wc -l < ${OUT_PREFIX}_m6a_drach.bed) / $(wc -l < ${OUT_PREFIX}_m6a_high.bed)" | bc)
# echo "DRACH-context fraction: $DRACH_FRAC (target 0.7-0.9)"

echo ""
echo "Cross-validate against orthogonal method:"
echo "  - GLORI (stoichiometric, antibody-free)"
echo "  - m6Anet (nanopore direct RNA)"
echo "  - DART-seq (in vivo enzyme-fusion)"
echo "Triangulation across 2+ methods = high-confidence m6A"
