#!/bin/bash
# Reference: PureCLIP 1.3.1+, CTK 1.1.4+, samtools 1.19+, bedtools 2.31+ | Verify API if version differs
# Single-nucleotide crosslink-site detection for CLIP-seq.
# Method chosen per CLIP variant: PureCLIP (HMM, all variants) or CTK CITS (truncation, iCLIP/eCLIP) or CTK CIMS (mutation, HITS/PAR-CLIP).

DEDUP_BAM=$1
SMINPUT_BAM=$2
GENOME_FA=$3
PROTOCOL=${4:-"eclip"}             # eclip | iclip | hits | parclip
EXPRESSED_BED=${5:-""}             # restrict PureCLIP scope for convergence
OUT_PREFIX=${6:-"sample"}
THREADS=${7:-8}

case $PROTOCOL in
  eclip|iclip)
    # Truncation-based detection. PureCLIP is the HMM standard; CTK CITS is the empirical alternative.

    # PureCLIP: jointly models enrichment + truncation + CL motif sequence context
    # -dm 8: merge sites within 8 nt
    PUREFLAGS=""
    if [ -n "$EXPRESSED_BED" ]; then
        PUREFLAGS="-iv $EXPRESSED_BED"
    fi
    pureclip \
        -i $DEDUP_BAM -bai ${DEDUP_BAM}.bai \
        -g $GENOME_FA \
        -ibam $SMINPUT_BAM -ibai ${SMINPUT_BAM}.bai \
        -o ${OUT_PREFIX}_pureclip_sites.bed \
        -or ${OUT_PREFIX}_pureclip_regions.bed \
        -nt $THREADS -dm 8 \
        $PUREFLAGS

    echo "PureCLIP single-nt sites: $(wc -l < ${OUT_PREFIX}_pureclip_sites.bed)"
    ;;

  hits)
    # HITS-CLIP uses CIMS deletion mode. Requires deletion-tolerant alignment (BWA-aln upstream).
    if command -v parseAlignment.pl > /dev/null; then
        parseAlignment.pl --map-qual 1 --min-len 18 \
            --mutation-file ${OUT_PREFIX}_mut.txt \
            $DEDUP_BAM ${OUT_PREFIX}_tags.bed
        getMutationType.pl ${OUT_PREFIX}_tags.bed ${OUT_PREFIX}_mut.txt \
            -type del > ${OUT_PREFIX}_del.mut
        CIMS.pl ${OUT_PREFIX}_tags.bed ${OUT_PREFIX}_del.mut \
            -big -c ${OUT_PREFIX}_cims_cache \
            -p 0.01 ${OUT_PREFIX}_cims_del.bed
        echo "HITS-CLIP CIMS deletion sites: $(wc -l < ${OUT_PREFIX}_cims_del.bed)"
    else
        echo "CTK not on PATH; install from github.com/chaolinzhanglab/ctk"
        echo "Or use PureCLIP with --substitution-aware (alternative)"
    fi
    ;;

  parclip)
    # PAR-CLIP uses T->C substitution. Two complementary tools:
    #   PARalyzer: kernel-density clusters (Hafner 2010)
    #   CTK CIMS substitution: single-nt T->C positions
    # PureCLIP also handles PAR-CLIP via its general HMM.

    if command -v parseAlignment.pl > /dev/null; then
        parseAlignment.pl --map-qual 1 --min-len 18 \
            --mutation-file ${OUT_PREFIX}_mut.txt \
            $DEDUP_BAM ${OUT_PREFIX}_tags.bed

        getMutationType.pl ${OUT_PREFIX}_tags.bed ${OUT_PREFIX}_mut.txt \
            -type sub -nuc t -mut c > ${OUT_PREFIX}_t2c.mut

        CIMS.pl ${OUT_PREFIX}_tags.bed ${OUT_PREFIX}_t2c.mut \
            -big -c ${OUT_PREFIX}_cims_cache \
            -p 0.001 ${OUT_PREFIX}_cims_t2c.bed

        echo "PAR-CLIP CIMS T->C single-nt sites: $(wc -l < ${OUT_PREFIX}_cims_t2c.bed)"
    fi
    echo "For cluster-level PAR-CLIP, use PARalyzer with parameters from Hafner 2010 supplementary"
    ;;

  *)
    echo "Unknown protocol: $PROTOCOL"
    exit 1
    ;;
esac

echo ""
echo "Cross-validate: run a second method and check mCross motif is the same"
echo "Downstream:"
echo "  - mCross for motif registration: mCross -i sites.bed -g genome.fa -k 7"
echo "  - BEAPR allele-specific binding: bedtools intersect with het VCF"
echo "  - RBPNet deep learning: feed CL count distribution as training signal"
