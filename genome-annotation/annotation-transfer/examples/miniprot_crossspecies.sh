#!/bin/bash
# Reference: miniprot 0.13+ | Verify API if version differs
# Cross-species annotation transfer with MiniProt
set -euo pipefail

TARGET=$1
PROTEINS=$2
OUTDIR=${3:-miniprot_out}
THREADS=${4:-16}
# Max intron size: vertebrates 500000, insects 50000, fungi 5000, plants 100000
MAX_INTRON=${5:-200000}

mkdir -p $OUTDIR

echo "=== MiniProt Cross-Species Annotation Transfer ==="
echo "Target assembly: $TARGET"
echo "Protein sequences: $PROTEINS"
echo "Max intron size: $MAX_INTRON"

PROT_COUNT=$(grep -c '^>' $PROTEINS)
echo "Input proteins: $PROT_COUNT"

# Index target genome
echo ""
echo "Indexing target genome..."
miniprot -t $THREADS -d ${OUTDIR}/target.mpi $TARGET

# Align proteins to genome
# --gff: GFF3 output with gene models
# -G: Max intron size (species-dependent)
# Secondary alignments (paralogs/recent duplications) are emitted by default; tune with
# --outn (max alignments per protein) and --outs (min secondary-to-best score ratio).
# RAISING --outs is a STRICTER filter that suppresses secondaries -- do not raise it to "include" them.
echo ""
echo "Aligning proteins to genome..."
miniprot \
    -t $THREADS \
    --gff \
    -G $MAX_INTRON \
    --outn 5 \
    ${OUTDIR}/target.mpi \
    $PROTEINS > ${OUTDIR}/miniprot_alignments.gff

# Summary
ALIGNED=$(grep -c $'\tmRNA\t' ${OUTDIR}/miniprot_alignments.gff || echo 0)

echo ""
echo "=========================================="
echo "MiniProt Summary"
echo "=========================================="
echo "Input proteins: $PROT_COUNT"
echo "Aligned gene models: $ALIGNED"

if [ $PROT_COUNT -gt 0 ]; then
    RATE=$(echo "scale=1; $ALIGNED * 100 / $PROT_COUNT" | bc)
    echo "Alignment rate: ${RATE}%"
fi

echo ""
echo "Gene models: ${OUTDIR}/miniprot_alignments.gff"
echo ""
echo "These alignments can serve as evidence for BRAKER3 or as standalone annotation."
