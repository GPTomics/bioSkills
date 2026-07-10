#!/bin/bash
# Reference: preseq 3.2+, picard 3.1+, samtools 1.19+, bedtools 2.31+, idr 2.0.4+, MultiQC 1.21+ | Verify API if version differs
# Comprehensive CLIP-seq library QC, ENCODE eCLIP convention.
# Five gates: preprocessing retention -> alignment rate -> library complexity -> FRiP -> IDR.

ALIGNED_BAM=$1                # pre-dedup BAM (preseq needs PCR duplicates)
DEDUP_BAM=$2                  # post-dedup BAM
PEAKS_BED=$3                  # stringent peaks (log2 FC >= 3)
SMINPUT_BAM=$4                # SMInput dedup BAM
REP2_PEAKS=${5:-""}           # second-replicate peaks for IDR (optional)
OUT_DIR=${6:-"qc"}

mkdir -p $OUT_DIR

echo "=== Gate 1: Preprocessing retention ==="
echo "(check cutadapt logs; target >= 70%)"
echo ""

echo "=== Gate 2: Alignment rate ==="
samtools flagstat $ALIGNED_BAM | tee ${OUT_DIR}/flagstat.txt
UNIQUE=$(samtools view -c -F 4 -q 255 $ALIGNED_BAM 2>/dev/null || samtools view -c -F 4 $ALIGNED_BAM)
TOTAL=$(samtools view -c $ALIGNED_BAM)
ALIGN_PCT=$(echo "scale=4; 100 * $UNIQUE / $TOTAL" | bc)
echo "Unique alignment: $ALIGN_PCT%"
echo "  Target: >= 60% (eCLIP), >= 70% (iCLIP/PAR-CLIP)"

echo ""
echo "=== Gate 3: Library complexity ==="
preseq lc_extrap -B -P $ALIGNED_BAM -o ${OUT_DIR}/preseq.txt 2> ${OUT_DIR}/preseq.log

# Extract predicted unique at 100M reads
EXP_100M=$(awk 'NR>1 && $1==1e+08' ${OUT_DIR}/preseq.txt | awk '{print $2}')
echo "Predicted unique fragments at 100M reads: $EXP_100M"
echo "  Target: >= 1e+07 (10M); < 3e+06 = library failed"

# Picard direct
picard EstimateLibraryComplexity \
    I=$ALIGNED_BAM \
    O=${OUT_DIR}/picard_complexity.txt 2> ${OUT_DIR}/picard.log
echo "Picard ESTIMATED_LIBRARY_SIZE (see ${OUT_DIR}/picard_complexity.txt)"

echo ""
echo "=== Gate 4: FRiP and IP enrichment ==="
DEDUP_TOTAL=$(samtools view -c -F 4 $DEDUP_BAM)
READS_IN_PEAKS=$(bedtools intersect -c -s -a $PEAKS_BED -b $DEDUP_BAM | awk '{s+=$NF} END {print s}')
FRIP=$(echo "scale=4; $READS_IN_PEAKS / $DEDUP_TOTAL" | bc)
echo "Dedup total reads: $DEDUP_TOTAL"
echo "Reads in stringent peaks: $READS_IN_PEAKS"
echo "FRiP: $FRIP"
echo "  Target: >= 0.005 (narrow-binding RBP)"
echo "  Atypical-binding RBPs (rare-transcript binders) exempt"

# IP vs SMInput enrichment. Comparing total library sizes measures sequencing depth, not enrichment;
# the informative quantity is the fraction of reads falling in peaks, depth-normalized in each library:
#   log2( (IP_in_peaks / IP_total) / (SMI_in_peaks / SMI_total) )
SMI_TOTAL=$(samtools view -c -F 4 $SMINPUT_BAM)
SMI_IN_PEAKS=$(bedtools intersect -c -s -a $PEAKS_BED -b $SMINPUT_BAM | awk '{s+=$NF} END {print s+0}')
echo "IP total reads:      $DEDUP_TOTAL  (in peaks: $READS_IN_PEAKS)"
echo "SMInput total reads: $SMI_TOTAL  (in peaks: $SMI_IN_PEAKS)"
if [ "$SMI_IN_PEAKS" -gt 0 ] && [ "$READS_IN_PEAKS" -gt 0 ]; then
    # bc -l: l() is natural log, so l(x)/l(2) is log2(x)
    ENRICH=$(echo "scale=4; l(($READS_IN_PEAKS / $DEDUP_TOTAL) / ($SMI_IN_PEAKS / $SMI_TOTAL)) / l(2)" | bc -l)
    echo "Global log2(IP/SMInput) in peaks: $ENRICH"
    awk -v e="$ENRICH" 'BEGIN { print (e+0 > 0.5) ? "  PASS: IP is enriched over SMInput" \
                                                   : "  FAIL: IP not enriched -- suspect a failed IP or a mismatched SMInput" }'
else
    ENRICH="NA"
    echo "Global log2(IP/SMInput) in peaks: NA (no reads in peaks in one library)"
fi

echo ""
echo "=== Gate 5: IDR replicate reproducibility ==="
if [ -n "$REP2_PEAKS" ]; then
    sort -k5,5gr $PEAKS_BED > ${OUT_DIR}/rep1.sorted.bed
    sort -k5,5gr $REP2_PEAKS > ${OUT_DIR}/rep2.sorted.bed

    idr --samples ${OUT_DIR}/rep1.sorted.bed ${OUT_DIR}/rep2.sorted.bed \
        --input-file-type bed --rank 5 \
        --output-file ${OUT_DIR}/idr.out \
        --idr-threshold 0.05 \
        --plot \
        --log-output-file ${OUT_DIR}/idr.log

    # ENCODE rule: Nt and Nself ratios both < 2
    NT=$(awk '$5 >= 540' ${OUT_DIR}/idr.out | wc -l)  # 540 = -log10(0.05) * 1000 from IDR encoding
    echo "IDR-passing peaks (true reps, threshold 0.05): $NT"
    echo "  Run pseudoreplicate IDR separately on each rep's split BAM"
    echo "  ENCODE rule: max(Nt, Nself) / min(Nt, Nself) <= 2"
else
    echo "Skipped (no second replicate provided)"
fi

echo ""
echo "=== Read distribution metagene ==="
if command -v read_distribution.py > /dev/null; then
    # Expects GENCODE BED12 reference - update path
    GENCODE_BED=${GENCODE_BED:-"gencode.v38.bed"}
    if [ -f "$GENCODE_BED" ]; then
        read_distribution.py -i $DEDUP_BAM -r $GENCODE_BED > ${OUT_DIR}/read_distribution.txt
        cat ${OUT_DIR}/read_distribution.txt
    else
        echo "Set GENCODE_BED to your annotation BED to run read_distribution.py"
    fi
else
    echo "RSeQC not installed; install with: conda install -c bioconda rseqc"
fi

echo ""
echo "=== rRNA contamination ==="
RRNA_READS=$(samtools idxstats $DEDUP_BAM 2>/dev/null | awk '$1 ~ /rRNA|45S|18S|28S|5_8S/ { sum+=$3 } END {print sum+0}')
RRNA_FRAC=$(echo "scale=4; $RRNA_READS / $DEDUP_TOTAL" | bc)
echo "rRNA reads: $RRNA_READS / $DEDUP_TOTAL ($RRNA_FRAC)"
echo "  Target: < 0.02 after rRNA pre-map; up to 0.30 without pre-map"

echo ""
echo "=== Summary ==="
echo "Aligned: $ALIGN_PCT% (target >= 60%)"
echo "Predicted unique at 100M: $EXP_100M (target >= 10M)"
echo "FRiP: $FRIP (target >= 0.005 narrow-binding)"
echo "IP/SMInput log2 enrichment in peaks: $ENRICH (target > 0.5)"
echo "rRNA fraction: $RRNA_FRAC (target < 0.02 post-prefilter)"
echo ""
echo "If multiple samples, run multiqc on the qc directory for a unified report:"
echo "  multiqc $OUT_DIR/"
