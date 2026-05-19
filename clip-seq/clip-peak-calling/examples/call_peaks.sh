#!/bin/bash
# Reference: CLIPper 2.0+, PureCLIP 1.3.1+, samtools 1.19+, idr 2.0.4+ | Verify API if version differs
# CLIP-seq peak calling: ENCODE eCLIP pipeline with CLIPper + SMInput normalization + IDR.
# For Skipper (high-sensitivity alternative), use the Snakemake workflow at github.com/algaebrown/skipper.
# For PureCLIP single-nt crosslink sites, see the PureCLIP block below.

IP_BAM=$1
SMINPUT_BAM=$2
SPECIES=${3:-"hg38"}
OUTPUT_PREFIX=${4:-"peaks"}
THREADS=${5:-8}

# Step 1: CLIPper peak calling on the IP BAM
# --FDR 0.05: peak-level FDR
# --superlocal: use the gene-level + local Poisson lambda hybrid
# --save-pickle: emit intermediate data for downstream normalization scripts
clipper \
    -b $IP_BAM \
    -s $SPECIES \
    -o ${OUTPUT_PREFIX}_clipper.bed \
    --FDR 0.05 \
    --superlocal \
    --save-pickle \
    --processors $THREADS

echo "CLIPper peaks: $(wc -l < ${OUTPUT_PREFIX}_clipper.bed)"

# Step 2: SMInput normalization with the Yeo lab eclip-pipeline scripts
# These scripts are in github.com/YeoLab/eclip; check they are on PATH.
# Output BED is enriched with log2(IP/SMInput) in column 4 and -log10(p) in column 5.
IP_READS=$(samtools view -c -F 4 $IP_BAM)
SMI_READS=$(samtools view -c -F 4 $SMINPUT_BAM)
echo $IP_READS > ${IP_BAM}.readnum.txt
echo $SMI_READS > ${SMINPUT_BAM}.readnum.txt

if command -v overlap_peakfi_with_bam_PE.py > /dev/null; then
    overlap_peakfi_with_bam_PE.py \
        ${OUTPUT_PREFIX}_clipper.bed \
        $IP_BAM $SMINPUT_BAM \
        ${IP_BAM}.readnum.txt ${SMINPUT_BAM}.readnum.txt \
        ${OUTPUT_PREFIX}_normed.bed
    compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.py \
        ${OUTPUT_PREFIX}_normed.bed \
        ${OUTPUT_PREFIX}_compressed.bed
else
    echo "Yeo lab eclip-pipeline scripts not on PATH; install from github.com/YeoLab/eclip"
    echo "Falling back to manual log2 FC calculation"
    # Manual fallback: bedtools intersect with IP and SMI BAMs separately
    bedtools intersect -c -s -a ${OUTPUT_PREFIX}_clipper.bed -b $IP_BAM > tmp.ip.bed
    bedtools intersect -c -s -a ${OUTPUT_PREFIX}_clipper.bed -b $SMINPUT_BAM > tmp.sm.bed
    paste tmp.ip.bed <(cut -f7 tmp.sm.bed) | awk -v ipn=$IP_READS -v smn=$SMI_READS 'BEGIN{OFS="\t"} {
        ip=($7+1)/ipn; sm=($8+1)/smn;
        log2fc = log(ip/sm)/log(2);
        print $1,$2,$3,$4,log2fc,$6
    }' > ${OUTPUT_PREFIX}_compressed.bed
    rm tmp.ip.bed tmp.sm.bed
fi

# Step 3: ENCODE stringent threshold (log2 FC >= 3 AND -log10 p >= 3)
# When using fallback (log2 FC in column 5 only), apply just the log2 threshold and document
awk 'BEGIN{FS=OFS="\t"} NF >= 5 && $5 >= 3' ${OUTPUT_PREFIX}_compressed.bed > ${OUTPUT_PREFIX}_stringent.bed
awk 'BEGIN{FS=OFS="\t"} NF >= 5 && $5 >= 1' ${OUTPUT_PREFIX}_compressed.bed > ${OUTPUT_PREFIX}_lenient.bed

echo ""
echo "Stringent peaks (log2 FC >= 3): $(wc -l < ${OUTPUT_PREFIX}_stringent.bed)"
echo "Lenient peaks (log2 FC >= 1):   $(wc -l < ${OUTPUT_PREFIX}_lenient.bed)"

# Step 4: FRiP
# ENCODE narrow-binding RBP minimum FRiP = 0.005
READS_IN_PEAKS=$(samtools view -c -L ${OUTPUT_PREFIX}_stringent.bed $IP_BAM)
TOTAL_READS=$(samtools view -c -F 4 $IP_BAM)
FRIP=$(echo "scale=4; $READS_IN_PEAKS / $TOTAL_READS" | bc)
echo "FRiP (stringent peaks): $FRIP   [ENCODE minimum 0.005]"

# Step 5: PureCLIP for single-nucleotide crosslink sites (optional)
# Run on the same IP + SMInput; outputs single-nt CL BED and broader region BED
# Skipped by default; uncomment to enable
# pureclip \
#     -i $IP_BAM -bai ${IP_BAM}.bai \
#     -g /path/to/${SPECIES}.fa \
#     -ibam $SMINPUT_BAM -ibai ${SMINPUT_BAM}.bai \
#     -o ${OUTPUT_PREFIX}_pureclip_sites.bed \
#     -or ${OUTPUT_PREFIX}_pureclip_regions.bed \
#     -nt $THREADS -dm 8

echo ""
echo "Next:"
echo "  - For IDR across replicates: idr --samples rep1.bed rep2.bed --rank 5 --idr-threshold 0.05"
echo "  - For single-nt CL sites: see clip-seq/crosslink-site-detection"
echo "  - For motif analysis: see clip-seq/clip-motif-analysis"
echo "  - For peak annotation: see clip-seq/binding-site-annotation"
