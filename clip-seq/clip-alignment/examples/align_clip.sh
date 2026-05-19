#!/bin/bash
# Reference: STAR 2.7.11b+, samtools 1.19+, umi_tools 1.1.5+ | Verify API if version differs
# CLIP-seq alignment using ENCODE eCLIP parameter block.
# `--alignEndsType EndToEnd` is mandatory: soft-clip of 5' bases would destroy the truncation = crosslink -1 base.
# `--outFilterMismatchNoverReadLmax 0.04` is correct for eCLIP/iCLIP; PAR-CLIP needs 0.07 to retain T->C signal.

R1=$1
R2=$2                        # empty string for single-end (iCLIP/iCLIP2)
STAR_INDEX=$3
OUTPUT_PREFIX=${4:-"sample"}
PROTOCOL=${5:-"eclip"}       # eclip | iclip2 | parclip
ALLOW_MULTIMAP=${6:-"no"}    # yes for repeat-binding RBPs (then run CLAM downstream)
THREADS=${7:-8}

# Set mismatch ceiling per protocol
# 0.04 (4%) for iCLIP/eCLIP - tolerates RT errors but excludes degraded reads
# 0.07 (7%) for PAR-CLIP - must retain T->C reads with 20-50% per-T conversion
case $PROTOCOL in
  parclip)  MISMATCH=0.07 ;;
  *)        MISMATCH=0.04 ;;
esac

# Multi-mapper handling
# Default 1 (unique only) - ENCODE standard
# 100 with multNmax=-1 - emit all for CLAM EM rescue (MATR3, LINE-1 binders)
if [ "$ALLOW_MULTIMAP" = "yes" ]; then
    MULTI_FLAGS="--outFilterMultimapNmax 100 --outSAMmultNmax -1"
else
    MULTI_FLAGS="--outFilterMultimapNmax 1"
fi

# Build read input string
if [ -n "$R2" ]; then
    READS="--readFilesIn $R1 $R2"
else
    READS="--readFilesIn $R1"
fi

# Alignment
# --outFilterScoreMinOverLread 0.66 and --outFilterMatchNminOverLread 0.66:
#   ENCODE eCLIP stringency on top of mismatch ceiling; require >=66% matched bases
STAR --runMode alignReads \
    --runThreadN $THREADS \
    --genomeDir $STAR_INDEX \
    --genomeLoad NoSharedMemory \
    $READS \
    --readFilesCommand zcat \
    --outFilterType BySJout \
    $MULTI_FLAGS \
    --alignEndsType EndToEnd \
    --outFilterMismatchNoverReadLmax $MISMATCH \
    --outFilterScoreMinOverLread 0.66 \
    --outFilterMatchNminOverLread 0.66 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --outFileNamePrefix ${OUTPUT_PREFIX}_

samtools index ${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam

# MAPQ filter (255 = unique in STAR)
samtools view -b -q 10 ${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam > ${OUTPUT_PREFIX}_q10.bam
samtools index ${OUTPUT_PREFIX}_q10.bam

# UMI dedup
# --method=unique is the ENCODE convention; directional is more conservative but slower
DEDUP_FLAGS="--method=unique"
if [ -n "$R2" ]; then
    DEDUP_FLAGS="$DEDUP_FLAGS --paired"
fi

umi_tools dedup \
    --stdin=${OUTPUT_PREFIX}_q10.bam \
    --stdout=${OUTPUT_PREFIX}_dedup.bam \
    $DEDUP_FLAGS \
    --log=${OUTPUT_PREFIX}_dedup.log

samtools index ${OUTPUT_PREFIX}_dedup.bam

# Soft-clip sanity check - CLIP requires < 1% reads with S in CIGAR
TOTAL=$(samtools view -c ${OUTPUT_PREFIX}_dedup.bam)
SOFTCLIP=$(samtools view ${OUTPUT_PREFIX}_dedup.bam | awk '$6 ~ /S/' | wc -l)
PCT=$(echo "scale=4; 100 * $SOFTCLIP / $TOTAL" | bc)
echo ""
echo "Alignment complete: ${OUTPUT_PREFIX}_dedup.bam"
echo "Soft-clipped reads: $SOFTCLIP / $TOTAL ($PCT%)"
echo "  Expected < 1% for end-to-end CLIP alignment"
echo ""
echo "Library complexity: run preseq lc_extrap on ${OUTPUT_PREFIX}_q10.bam (pre-dedup)"
echo "Next: peak calling (see clip-seq/clip-peak-calling)"
