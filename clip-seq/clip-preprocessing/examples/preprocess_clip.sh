#!/bin/bash
# Reference: umi_tools 1.1.5+, cutadapt 4.6+, bowtie2 2.5.3+, samtools 1.19+ | Verify API if version differs
# eCLIP/iCLIP preprocessing aligned with ENCODE eCLIP standard.
# The 5' end of read 2 (paired-end) carries the crosslink -1 truncation; never trim it.

R1=$1
R2=$2
PROTOCOL=${3:-"eclip"}       # eclip | iclip2 | parclip
OUTPUT_PREFIX=${4:-"sample"}
RRNA_INDEX=${5:-""}          # bowtie2 index path for pre-mapping; "" to skip
THREADS=${6:-8}

case $PROTOCOL in
  eclip)
    # eCLIP: 10 nt UMI on R1 5' end
    UMI_PATTERN="NNNNNNNNNN"
    ADAPTER_3P_R1="AGATCGGAAGAGCACACGTCT"
    ADAPTER_3P_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    ADAPTER_5P_R2_PASS2="GATCGTCGGACTGTAGAACTCTGAAC"
    TWO_PASS="yes"
    ;;
  iclip2)
    # iCLIP/iCLIP2: NNNXXXXNN. Demultiplex first by the 4 X bases (NOT covered in this script).
    # After demultiplex, the surviving UMI is 5 nt (3 + 2 random) flanking the now-stripped barcode.
    UMI_PATTERN="NNNNN"
    ADAPTER_3P_R1="AGATCGGAAGAGCGGTTCAG"  # L3 adapter
    TWO_PASS="no"
    ;;
  parclip)
    # PAR-CLIP: 4 nt random barcode typical; T-to-C mutations are SIGNAL
    UMI_PATTERN="NNNN"
    ADAPTER_3P_R1="TCGTATGCCGTCTTCTGCTTG"
    TWO_PASS="no"
    ;;
  *)
    echo "Unknown protocol: $PROTOCOL"
    exit 1
    ;;
esac

echo "Protocol: $PROTOCOL"
echo "UMI pattern: $UMI_PATTERN"

# Step 1: UMI extraction
# Reason: random barcodes must be moved to read name BEFORE adapter trimming so that
# umi_tools dedup later can collapse PCR duplicates by (UMI, position, strand).
mkdir -p ${OUTPUT_PREFIX}_qc
if [ -n "$R2" ]; then
    umi_tools extract \
        --stdin=$R1 --read2-in=$R2 \
        --bc-pattern=$UMI_PATTERN \
        --stdout=${OUTPUT_PREFIX}_R1.umi.fq.gz \
        --read2-out=${OUTPUT_PREFIX}_R2.umi.fq.gz \
        --log=${OUTPUT_PREFIX}_qc/umi_extract.log
else
    umi_tools extract \
        --stdin=$R1 \
        --bc-pattern=$UMI_PATTERN \
        --stdout=${OUTPUT_PREFIX}_R1.umi.fq.gz \
        --log=${OUTPUT_PREFIX}_qc/umi_extract.log
fi

# Step 2: 3' adapter + light quality trim
# -q 6 is permissive on purpose; aggressive quality trimming chews back the 5' truncation base on R2
# -m 18 is mandatory; reads under 18 nt multi-map at > 50% and lose crosslink resolution
if [ -n "$R2" ]; then
    cutadapt \
        -a $ADAPTER_3P_R1 -A $ADAPTER_3P_R2 \
        --quality-base 33 -q 6 -m 18 \
        -j $THREADS \
        -o ${OUTPUT_PREFIX}_R1.trim1.fq.gz \
        -p ${OUTPUT_PREFIX}_R2.trim1.fq.gz \
        ${OUTPUT_PREFIX}_R1.umi.fq.gz ${OUTPUT_PREFIX}_R2.umi.fq.gz \
        > ${OUTPUT_PREFIX}_qc/cutadapt_pass1.log 2>&1
else
    cutadapt \
        -a $ADAPTER_3P_R1 \
        --quality-base 33 -q 6 -m 18 \
        -j $THREADS \
        -o ${OUTPUT_PREFIX}_R1.trim1.fq.gz \
        ${OUTPUT_PREFIX}_R1.umi.fq.gz \
        > ${OUTPUT_PREFIX}_qc/cutadapt_pass1.log 2>&1
fi

# Step 3 (eCLIP only): pass 2 strips read-through adapter from R2 5' end ONLY
# -G (uppercase) anchors on R2 5'; -g (lowercase) on R1 5' would destroy R1 truncation in single-end eCLIP
if [ "$TWO_PASS" = "yes" ] && [ -n "$R2" ]; then
    cutadapt \
        -G $ADAPTER_5P_R2_PASS2 \
        --quality-base 33 -q 6 -m 18 \
        -j $THREADS \
        -o ${OUTPUT_PREFIX}_R1.trim.fq.gz \
        -p ${OUTPUT_PREFIX}_R2.trim.fq.gz \
        ${OUTPUT_PREFIX}_R1.trim1.fq.gz ${OUTPUT_PREFIX}_R2.trim1.fq.gz \
        > ${OUTPUT_PREFIX}_qc/cutadapt_pass2.log 2>&1
else
    mv ${OUTPUT_PREFIX}_R1.trim1.fq.gz ${OUTPUT_PREFIX}_R1.trim.fq.gz
    [ -n "$R2" ] && mv ${OUTPUT_PREFIX}_R2.trim1.fq.gz ${OUTPUT_PREFIX}_R2.trim.fq.gz
fi

# Step 4 (optional): pre-map to rRNA. ENCODE eCLIP pipeline does this for speed.
# Reads aligning to rRNA are discarded; survivors feed the genome aligner.
if [ -n "$RRNA_INDEX" ]; then
    if [ -n "$R2" ]; then
        bowtie2 -x $RRNA_INDEX \
            -1 ${OUTPUT_PREFIX}_R1.trim.fq.gz -2 ${OUTPUT_PREFIX}_R2.trim.fq.gz \
            --un-conc-gz ${OUTPUT_PREFIX}_norrna_R%.fq.gz \
            -p $THREADS -S /dev/null \
            2> ${OUTPUT_PREFIX}_qc/bowtie2_rrna.log
        mv ${OUTPUT_PREFIX}_norrna_R1.fq.gz ${OUTPUT_PREFIX}_R1.final.fq.gz
        mv ${OUTPUT_PREFIX}_norrna_R2.fq.gz ${OUTPUT_PREFIX}_R2.final.fq.gz
    else
        bowtie2 -x $RRNA_INDEX \
            -U ${OUTPUT_PREFIX}_R1.trim.fq.gz \
            --un-gz ${OUTPUT_PREFIX}_R1.final.fq.gz \
            -p $THREADS -S /dev/null \
            2> ${OUTPUT_PREFIX}_qc/bowtie2_rrna.log
    fi
else
    mv ${OUTPUT_PREFIX}_R1.trim.fq.gz ${OUTPUT_PREFIX}_R1.final.fq.gz
    [ -n "$R2" ] && mv ${OUTPUT_PREFIX}_R2.trim.fq.gz ${OUTPUT_PREFIX}_R2.final.fq.gz
fi

# Step 5: report
echo ""
echo "Preprocessing complete"
ls -lh ${OUTPUT_PREFIX}_*final.fq.gz
echo ""
echo "Adapter trim retention (should be >= 70%):"
grep -E "passing filters|Pairs written" ${OUTPUT_PREFIX}_qc/cutadapt_pass1.log | tail -2
echo ""
echo "Next: align with STAR (see clip-seq/clip-alignment), then dedupe with:"
echo "  umi_tools dedup --stdin=aligned.bam --stdout=dedup.bam --method=unique --paired"
