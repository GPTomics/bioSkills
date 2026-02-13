#!/bin/bash
# Reference: pandas 2.2+ | Verify API if version differs
# Non-coding RNA annotation with Infernal and tRNAscan-SE
set -euo pipefail

GENOME=$1
RFAM_CM=${2:-Rfam.cm}
RFAM_CLANIN=${3:-Rfam.clanin}
DOMAIN=${4:-B}   # B=bacteria, A=archaea, E=eukaryote
OUTDIR=${5:-ncrna_out}
THREADS=${6:-16}

mkdir -p $OUTDIR

echo "=== Non-Coding RNA Annotation ==="
echo "Genome: $GENOME"
echo "Domain: $DOMAIN"

# Verify Rfam database is pressed
if [ ! -f "${RFAM_CM}.i1m" ]; then
    echo "Pressing Rfam CM database..."
    cmpress $RFAM_CM
fi

# Run Infernal cmscan against Rfam
# --cut_ga: Rfam gathering threshold (family-specific, recommended)
# --rfam: Rfam-optimized speed settings
# --nohmmonly: CM-only scoring (more accurate for structured RNAs)
# --clanin: Resolve overlapping hits from same clan
echo ""
echo "Running Infernal cmscan..."
cmscan \
    --cut_ga \
    --rfam \
    --nohmmonly \
    --tblout ${OUTDIR}/infernal_hits.tbl \
    --fmt 2 \
    --cpu $THREADS \
    --clanin $RFAM_CLANIN \
    $RFAM_CM \
    $GENOME > ${OUTDIR}/infernal_hits.out

INFERNAL_HITS=$(grep -c -v '^#' ${OUTDIR}/infernal_hits.tbl || echo 0)
echo "Infernal hits: $INFERNAL_HITS"

# Run tRNAscan-SE
echo ""
echo "Running tRNAscan-SE (domain: $DOMAIN)..."
tRNAscan-SE \
    -${DOMAIN} \
    -o ${OUTDIR}/trnascan_results.txt \
    --gff ${OUTDIR}/trnascan.gff3 \
    --detail \
    --thread $THREADS \
    $GENOME

TRNA_COUNT=$(grep -c $'\ttRNA\t' ${OUTDIR}/trnascan.gff3 || echo 0)
echo "tRNAs found: $TRNA_COUNT"

# Run barrnap for rRNA (quick check)
echo ""
echo "Running barrnap for rRNA..."
BARRNAP_KINGDOM="bac"
if [ "$DOMAIN" = "E" ]; then BARRNAP_KINGDOM="euk"; fi
if [ "$DOMAIN" = "A" ]; then BARRNAP_KINGDOM="arc"; fi

barrnap --kingdom $BARRNAP_KINGDOM --threads $THREADS $GENOME > ${OUTDIR}/barrnap_rrna.gff3 2>/dev/null || echo "barrnap skipped"
RRNA_COUNT=$(grep -c 'rRNA' ${OUTDIR}/barrnap_rrna.gff3 || echo 0)
echo "rRNAs found: $RRNA_COUNT"

# Summary
echo ""
echo "=========================================="
echo "ncRNA Annotation Summary"
echo "=========================================="
echo "Infernal (Rfam) hits: $INFERNAL_HITS"
echo "tRNAs (tRNAscan-SE): $TRNA_COUNT"
echo "rRNAs (barrnap): $RRNA_COUNT"

# Count Infernal hits by type
echo ""
echo "Infernal hits by Rfam family (top 20):"
grep -v '^#' ${OUTDIR}/infernal_hits.tbl | awk '{print $3}' | sort | uniq -c | sort -rn | head -20

echo ""
echo "Results in: $OUTDIR"
echo "Run parse_ncrna.py to combine and categorize all results."
