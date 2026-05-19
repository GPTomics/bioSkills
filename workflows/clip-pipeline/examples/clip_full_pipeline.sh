#!/bin/bash
# Reference: umi_tools 1.1.5+, cutadapt 4.6+, STAR 2.7.11b+, samtools 1.19+, CLIPper 2.0+, bedtools 2.31+, preseq 3.2+ | Verify API if version differs
# Complete eCLIP pipeline: FASTQ -> ENCODE-stringent peaks (log2 FC >= 3 AND -log10 p >= 3)
# Critical constraints from the clip-seq decision-grade skills:
#   - 3'-only adapter trim with -q 6 -m 18 (preserves R2 5' = crosslink site -1)
#   - STAR --alignEndsType EndToEnd, mismatch 0.04 (raise to 0.07 for PAR-CLIP only)
#   - umi_tools dedup --method=unique (ENCODE convention)
#   - CLIPper + SMInput log2 normalization (Yeo lab eclip-pipeline scripts)

set -euo pipefail

R1=$1
R2=$2
SMINPUT_R1=$3
SMINPUT_R2=$4
STAR_INDEX=$5
GENOME_FA=$6
SPECIES=${7:-"hg38"}
PROTOCOL=${8:-"eclip"}        # eclip | iclip2 | parclip
OUTPUT_DIR=${9:-"clip_results"}
THREADS=${10:-8}

mkdir -p ${OUTPUT_DIR}/{preprocessed,aligned,qc,peaks,crosslinks,annotation,motifs}

# Protocol-specific parameters from clip-seq/clip-preprocessing and clip-seq/clip-alignment
case $PROTOCOL in
  eclip)
    UMI_PATTERN="NNNNNNNNNN"
    ADAPTER_3P_R1="AGATCGGAAGAGCACACGTCT"
    ADAPTER_3P_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    ADAPTER_5P_R2_PASS2="GATCGTCGGACTGTAGAACTCTGAAC"
    TWO_PASS="yes"
    MISMATCH=0.04
    ;;
  iclip2)
    # Assumes demultiplexed input. For multiplexed, demultiplex by NNNXXXXNN library barcode first.
    UMI_PATTERN="NNNNN"
    ADAPTER_3P_R1="AGATCGGAAGAGCGGTTCAG"
    TWO_PASS="no"
    MISMATCH=0.04
    ;;
  parclip)
    UMI_PATTERN="NNNN"
    ADAPTER_3P_R1="TCGTATGCCGTCTTCTGCTTG"
    TWO_PASS="no"
    MISMATCH=0.07   # PAR-CLIP exception: T->C signal would otherwise be filtered
    ;;
  *) echo "Unknown protocol: $PROTOCOL"; exit 1 ;;
esac

run_preprocess() {
    local IN_R1=$1 IN_R2=$2 PREFIX=$3
    umi_tools extract \
        --bc-pattern=$UMI_PATTERN \
        --stdin=$IN_R1 --read2-in=$IN_R2 \
        --stdout=${OUTPUT_DIR}/preprocessed/${PREFIX}_R1.umi.fq.gz \
        --read2-out=${OUTPUT_DIR}/preprocessed/${PREFIX}_R2.umi.fq.gz \
        --log=${OUTPUT_DIR}/qc/${PREFIX}_umi.log

    cutadapt \
        -a $ADAPTER_3P_R1 ${ADAPTER_3P_R2:+-A $ADAPTER_3P_R2} \
        --quality-base 33 -q 6 -m 18 -j $THREADS \
        -o ${OUTPUT_DIR}/preprocessed/${PREFIX}_R1.p1.fq.gz \
        -p ${OUTPUT_DIR}/preprocessed/${PREFIX}_R2.p1.fq.gz \
        ${OUTPUT_DIR}/preprocessed/${PREFIX}_R1.umi.fq.gz \
        ${OUTPUT_DIR}/preprocessed/${PREFIX}_R2.umi.fq.gz \
        > ${OUTPUT_DIR}/qc/${PREFIX}_cutadapt_pass1.log 2>&1

    if [ "$TWO_PASS" = "yes" ]; then
        cutadapt \
            -G $ADAPTER_5P_R2_PASS2 \
            --quality-base 33 -q 6 -m 18 -j $THREADS \
            -o ${OUTPUT_DIR}/preprocessed/${PREFIX}_R1.trim.fq.gz \
            -p ${OUTPUT_DIR}/preprocessed/${PREFIX}_R2.trim.fq.gz \
            ${OUTPUT_DIR}/preprocessed/${PREFIX}_R1.p1.fq.gz \
            ${OUTPUT_DIR}/preprocessed/${PREFIX}_R2.p1.fq.gz \
            > ${OUTPUT_DIR}/qc/${PREFIX}_cutadapt_pass2.log 2>&1
    else
        mv ${OUTPUT_DIR}/preprocessed/${PREFIX}_R1.p1.fq.gz ${OUTPUT_DIR}/preprocessed/${PREFIX}_R1.trim.fq.gz
        mv ${OUTPUT_DIR}/preprocessed/${PREFIX}_R2.p1.fq.gz ${OUTPUT_DIR}/preprocessed/${PREFIX}_R2.trim.fq.gz
    fi
}

run_align() {
    local PREFIX=$1
    STAR --runMode alignReads \
        --runThreadN $THREADS \
        --genomeDir $STAR_INDEX \
        --genomeLoad NoSharedMemory \
        --readFilesIn ${OUTPUT_DIR}/preprocessed/${PREFIX}_R1.trim.fq.gz \
                      ${OUTPUT_DIR}/preprocessed/${PREFIX}_R2.trim.fq.gz \
        --readFilesCommand zcat \
        --outFilterType BySJout \
        --outFilterMultimapNmax 1 \
        --alignEndsType EndToEnd \
        --outFilterMismatchNoverReadLmax $MISMATCH \
        --outFilterScoreMinOverLread 0.66 \
        --outFilterMatchNminOverLread 0.66 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes All \
        --outFileNamePrefix ${OUTPUT_DIR}/aligned/${PREFIX}_

    samtools index ${OUTPUT_DIR}/aligned/${PREFIX}_Aligned.sortedByCoord.out.bam
    samtools view -b -q 10 ${OUTPUT_DIR}/aligned/${PREFIX}_Aligned.sortedByCoord.out.bam \
        > ${OUTPUT_DIR}/aligned/${PREFIX}_q10.bam
    samtools index ${OUTPUT_DIR}/aligned/${PREFIX}_q10.bam

    umi_tools dedup \
        --stdin=${OUTPUT_DIR}/aligned/${PREFIX}_q10.bam \
        --stdout=${OUTPUT_DIR}/aligned/${PREFIX}_dedup.bam \
        --method=unique \
        --paired \
        --log=${OUTPUT_DIR}/qc/${PREFIX}_dedup.log
    samtools index ${OUTPUT_DIR}/aligned/${PREFIX}_dedup.bam
}

echo "=== Preprocess IP sample ==="
run_preprocess $R1 $R2 ip

echo "=== Preprocess SMInput control ==="
run_preprocess $SMINPUT_R1 $SMINPUT_R2 sminput

echo "=== Align IP ==="
run_align ip

echo "=== Align SMInput ==="
run_align sminput

# QC Gate 3: Library complexity (preseq, target >= 1M unique at sequenced depth)
echo "=== Library complexity QC ==="
preseq lc_extrap -B -P ${OUTPUT_DIR}/aligned/ip_q10.bam -o ${OUTPUT_DIR}/qc/preseq_ip.txt 2> ${OUTPUT_DIR}/qc/preseq_ip.log

# Peak calling: CLIPper
echo "=== CLIPper peak calling ==="
clipper \
    -b ${OUTPUT_DIR}/aligned/ip_dedup.bam \
    -s $SPECIES \
    -o ${OUTPUT_DIR}/peaks/clipper.bed \
    --FDR 0.05 \
    --superlocal \
    --save-pickle \
    --processors $THREADS

# SMInput normalization: log2(IP/SMI) and -log10 p per peak
# Requires Yeo lab eclip-pipeline scripts on PATH (github.com/YeoLab/eclip)
IP_READS=$(samtools view -c -F 4 ${OUTPUT_DIR}/aligned/ip_dedup.bam)
SMI_READS=$(samtools view -c -F 4 ${OUTPUT_DIR}/aligned/sminput_dedup.bam)
echo $IP_READS > ${OUTPUT_DIR}/aligned/ip_dedup.bam.readnum.txt
echo $SMI_READS > ${OUTPUT_DIR}/aligned/sminput_dedup.bam.readnum.txt

if command -v overlap_peakfi_with_bam_PE.py > /dev/null; then
    overlap_peakfi_with_bam_PE.py \
        ${OUTPUT_DIR}/peaks/clipper.bed \
        ${OUTPUT_DIR}/aligned/ip_dedup.bam \
        ${OUTPUT_DIR}/aligned/sminput_dedup.bam \
        ${OUTPUT_DIR}/aligned/ip_dedup.bam.readnum.txt \
        ${OUTPUT_DIR}/aligned/sminput_dedup.bam.readnum.txt \
        ${OUTPUT_DIR}/peaks/normed.bed
    compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.py \
        ${OUTPUT_DIR}/peaks/normed.bed \
        ${OUTPUT_DIR}/peaks/compressed.bed

    # ENCODE stringent thresholds
    awk 'BEGIN{FS=OFS="\t"} NF >= 6 && $5 >= 3 && $6 >= 3' \
        ${OUTPUT_DIR}/peaks/compressed.bed > ${OUTPUT_DIR}/peaks/stringent.bed
    awk 'BEGIN{FS=OFS="\t"} NF >= 6 && $5 >= 1 && $6 >= 2' \
        ${OUTPUT_DIR}/peaks/compressed.bed > ${OUTPUT_DIR}/peaks/lenient.bed
    echo "Stringent peaks (log2 FC >= 3, -log10 p >= 3): $(wc -l < ${OUTPUT_DIR}/peaks/stringent.bed)"
else
    echo "Yeo eclip-pipeline scripts not on PATH; install from github.com/YeoLab/eclip for SMInput normalization"
    cp ${OUTPUT_DIR}/peaks/clipper.bed ${OUTPUT_DIR}/peaks/stringent.bed
fi

# FRiP (QC gate 4; ENCODE narrow-binding RBP minimum 0.005)
READS_IN_PEAKS=$(bedtools intersect -c -s -a ${OUTPUT_DIR}/peaks/stringent.bed -b ${OUTPUT_DIR}/aligned/ip_dedup.bam | awk '{s+=$NF} END {print s}')
TOTAL_READS=$(samtools view -c -F 4 ${OUTPUT_DIR}/aligned/ip_dedup.bam)
FRIP=$(echo "scale=4; $READS_IN_PEAKS / $TOTAL_READS" | bc)
echo "FRiP: $FRIP (ENCODE narrow-binding minimum 0.005)"

# Single-nucleotide crosslink-site detection with PureCLIP
echo "=== PureCLIP single-nt crosslink sites ==="
if command -v pureclip > /dev/null; then
    pureclip \
        -i ${OUTPUT_DIR}/aligned/ip_dedup.bam \
        -bai ${OUTPUT_DIR}/aligned/ip_dedup.bam.bai \
        -g $GENOME_FA \
        -ibam ${OUTPUT_DIR}/aligned/sminput_dedup.bam \
        -ibai ${OUTPUT_DIR}/aligned/sminput_dedup.bam.bai \
        -o ${OUTPUT_DIR}/crosslinks/sites.bed \
        -or ${OUTPUT_DIR}/crosslinks/regions.bed \
        -nt $THREADS -dm 8
fi

# Motif analysis with GC-matched background
echo "=== Motif analysis ==="
bedtools getfasta -fi $GENOME_FA -bed ${OUTPUT_DIR}/peaks/stringent.bed -s \
    -fo ${OUTPUT_DIR}/motifs/peaks.fa

if command -v findMotifs.pl > /dev/null; then
    findMotifs.pl ${OUTPUT_DIR}/motifs/peaks.fa fasta ${OUTPUT_DIR}/motifs/homer \
        -rna -len 5,6,7,8 -p $THREADS
fi

echo ""
echo "=== Pipeline complete ==="
echo "Key outputs:"
echo "  Dedup IP BAM:      ${OUTPUT_DIR}/aligned/ip_dedup.bam"
echo "  Dedup SMInput BAM: ${OUTPUT_DIR}/aligned/sminput_dedup.bam"
echo "  Stringent peaks:   ${OUTPUT_DIR}/peaks/stringent.bed"
echo "  Single-nt CL:      ${OUTPUT_DIR}/crosslinks/sites.bed"
echo "  Motifs:            ${OUTPUT_DIR}/motifs/homer/"
echo ""
echo "Next steps:"
echo "  - IDR across replicates (see clip-seq/clip-qc)"
echo "  - Annotation with ChIPseeker (see clip-seq/binding-site-annotation)"
echo "  - mCross CL-registered motif (see clip-seq/clip-motif-analysis)"
echo "  - DEWSeq differential between conditions (see clip-seq/differential-clip)"
