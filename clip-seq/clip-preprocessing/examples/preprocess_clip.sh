#!/bin/bash
# CLIP-seq preprocessing with UMI handling

R1=$1
R2=$2
OUTPUT_PREFIX=${3:-"clip"}
UMI_LENGTH=${4:-10}

echo "CLIP-seq preprocessing"
echo "Input: $R1, $R2"
echo "UMI length: $UMI_LENGTH"

# Generate UMI pattern (e.g., NNNNNNNNNN for 10-nt UMI)
UMI_PATTERN=$(printf 'N%.0s' $(seq 1 $UMI_LENGTH))

# Step 1: Extract UMIs
# UMIs are typically in read 1 for eCLIP
# Pattern: N = UMI base, X = discard, C = cell barcode
echo "Extracting UMIs..."
umi_tools extract \
    --stdin=$R1 \
    --read2-in=$R2 \
    --bc-pattern=$UMI_PATTERN \
    --stdout=${OUTPUT_PREFIX}_R1_umi.fq.gz \
    --read2-out=${OUTPUT_PREFIX}_R2_umi.fq.gz \
    --log=${OUTPUT_PREFIX}_umi_extract.log

# Step 2: Adapter trimming
# eCLIP adapters (adjust for your protocol)
ADAPTER_R1="AGATCGGAAGAGCACACGTCT"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

echo "Trimming adapters..."
cutadapt \
    -a $ADAPTER_R1 \
    -A $ADAPTER_R2 \
    -m 18 \
    -j 4 \
    -o ${OUTPUT_PREFIX}_R1_trim.fq.gz \
    -p ${OUTPUT_PREFIX}_R2_trim.fq.gz \
    ${OUTPUT_PREFIX}_R1_umi.fq.gz \
    ${OUTPUT_PREFIX}_R2_umi.fq.gz \
    > ${OUTPUT_PREFIX}_cutadapt.log 2>&1

# Step 3: Second pass trimming (for read-through)
# Some CLIP protocols have adapter at both ends
echo "Second pass trimming..."
cutadapt \
    -g $ADAPTER_R1 \
    -m 18 \
    -j 4 \
    -o ${OUTPUT_PREFIX}_R1_final.fq.gz \
    ${OUTPUT_PREFIX}_R1_trim.fq.gz \
    >> ${OUTPUT_PREFIX}_cutadapt.log 2>&1

cutadapt \
    -g $ADAPTER_R2 \
    -m 18 \
    -j 4 \
    -o ${OUTPUT_PREFIX}_R2_final.fq.gz \
    ${OUTPUT_PREFIX}_R2_trim.fq.gz \
    >> ${OUTPUT_PREFIX}_cutadapt.log 2>&1

echo ""
echo "Preprocessing complete!"
echo "Output files:"
ls -lh ${OUTPUT_PREFIX}_*final.fq.gz

echo ""
echo "Trimming summary:"
grep -E "Total reads|Reads with adapters|Too short" ${OUTPUT_PREFIX}_cutadapt.log | head -6

echo ""
echo "Next step: align with STAR/bowtie2, then deduplicate with:"
echo "  umi_tools dedup --stdin=aligned.bam --stdout=deduped.bam"
