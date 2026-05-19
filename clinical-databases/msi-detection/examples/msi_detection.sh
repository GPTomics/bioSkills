#!/bin/bash
# Reference: msisensor-pro 1.2+, msisensor 0.6+, samtools 1.19+ | Verify CLI flags if version differs.
#
# MSI detection workflow: tumor-only (MSIsensor-pro) and paired-normal (MSIsensor).
# FDA pembrolizumab MSI-H / dMMR pan-tumor approval (2017; Le 2015 NEJM seminal).
# Universal Lynch syndrome screening per NCCN / ACG (CRC <= 70 yr).

set -euo pipefail

MODE="${1:-tumor_only}"
TUMOR_BAM="${TUMOR_BAM:-tumor.bam}"
NORMAL_BAM="${NORMAL_BAM:-normal.bam}"
SAMPLE_ID="${SAMPLE_ID:-sample}"
REF_FA="${REF_FA:-GRCh38.fa}"
OUTPUT_DIR="${OUTPUT_DIR:-msi_output}"
THREADS="${THREADS:-16}"

mkdir -p "$OUTPUT_DIR"

# Step 1: One-time microsatellite list generation
if [ ! -f microsatellites.list ]; then
    echo "Generating microsatellite list from reference..."
    # -p min repeats; -m min length; -d output file
    msisensor-pro scan -d "$REF_FA" -o microsatellites.list -p 1 -m 5
fi

case "$MODE" in
    tumor_only)
        # MSIsensor-pro: no matched normal required
        # Requires pre-computed baseline from cohort of N normal samples
        # See msisensor-pro baseline for one-time baseline generation
        echo "=== MSI: Tumor-only (MSIsensor-pro) ==="
        if [ ! -f baseline.list ]; then
            echo "ERROR: baseline.list missing. Generate with:"
            echo "  msisensor-pro baseline -d microsatellites.list -i normal_samples.list -o baseline.list -b $THREADS"
            exit 1
        fi
        msisensor-pro pro \
            -d microsatellites.list \
            -t "$TUMOR_BAM" \
            -i "$SAMPLE_ID" \
            -o "$OUTPUT_DIR/${SAMPLE_ID}_msi" \
            -b "$THREADS" \
            --baseline baseline.list
        ;;
    paired)
        # MSIsensor: paired tumor-normal WES gold standard
        echo "=== MSI: Paired tumor-normal (MSIsensor) ==="
        msisensor msi \
            -d microsatellites.list \
            -n "$NORMAL_BAM" \
            -t "$TUMOR_BAM" \
            -o "$OUTPUT_DIR/${SAMPLE_ID}_msi" \
            -b "$THREADS"
        ;;
    ctdna)
        # MSIsensor-ct: cfDNA / liquid biopsy
        echo "=== MSI: ctDNA / liquid biopsy (MSIsensor-ct) ==="
        msisensor-ct \
            -d microsatellites.list \
            -t "$TUMOR_BAM" \
            -o "$OUTPUT_DIR/${SAMPLE_ID}_msi" \
            -b "$THREADS"
        ;;
    *)
        echo "Unknown mode: $MODE"
        echo "Usage: $0 [tumor_only|paired|ctdna]"
        exit 1
        ;;
esac

# Parse output
RESULT="$OUTPUT_DIR/${SAMPLE_ID}_msi"
if [ ! -f "$RESULT" ]; then
    echo "Error: result file missing"
    exit 1
fi

echo
echo "=== MSI Result ==="
column -t "$RESULT"

# Apply MSI-H classification
# FoCR-recommended threshold: >= 20% unstable loci
UNSTABLE_PCT=$(awk 'NR==2 {print $NF}' "$RESULT")
echo
echo "Unstable percentage: $UNSTABLE_PCT%"

if (( $(echo "$UNSTABLE_PCT >= 20" | bc -l) )); then
    echo "CLASSIFICATION: MSI-H"
    echo "  -> FDA pembrolizumab MSI-H / dMMR pan-tumor eligible (2017)"
    echo "  -> KEYNOTE-016/164/158 (Le 2015 NEJM seminal)"
    echo "  -> First-line CRC: KEYNOTE-177 (2020)"
    echo "  -> Confirm IHC (MLH1/MSH2/MSH6/PMS2 loss); apply Lynch workflow"
elif (( $(echo "$UNSTABLE_PCT >= 10" | bc -l) )); then
    echo "CLASSIFICATION: MSI-L (intermediate; clinically MSS per FDA)"
else
    echo "CLASSIFICATION: MSS"
fi

echo
echo "=== Lynch Syndrome Workflow ==="
echo "If MSI-H:"
echo "  1. Confirm with IHC (MLH1/MSH2/MSH6/PMS2)"
echo "  2. If MLH1 loss -> test MLH1 promoter methylation"
echo "     - Methylated -> sporadic (not Lynch)"
echo "     - Unmethylated -> Lynch suspect; germline testing"
echo "  3. If MSH2/6/PMS2 loss -> Lynch suspect; germline testing"
echo
echo "=== Reconciliation with other biomarkers ==="
echo "MSI-H + TMB-H: MSI-H is PRIMARY (Sha 2020 Cell Rep Med); TMB-H NOT additive."
echo "MSI-H + retained IHC: check POLE-exo via signatures (SBS10a/10b)"
echo "                      and consider rare MSH6-only subtype."
echo "MSS + TMB >= 100: investigate POLE-exo signature; ICI may still benefit."
