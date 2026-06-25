#!/bin/bash
# Reference: ViennaRNA 2.6+ | Verify API if version differs
# Consensus RNA secondary structure prediction from a multiple sequence alignment.
# RNAalifold combines thermodynamics with covariation; a consensus structure is a
# HYPOTHESIS until covariation is statistically validated (see covariation-analysis / R-scape).

ALIGNMENT=$1
OUTPUT_PREFIX=${2:-"consensus"}

if [ -z "$ALIGNMENT" ]; then
    echo "Usage: $0 <alignment.sto|alignment.aln> [output_prefix]"
    echo ""
    echo "Format (Stockholm/Clustal/FASTA) is auto-detected; the alignment is a positional argument."
    exit 1
fi

echo "=== RNAalifold: consensus structure ==="
# --ribosum_scoring improves covariation detection; -d0 avoids dangle artifacts at gapped columns;
# --noPS suppresses the alirna.ps / alidot.ps PostScript files RNAalifold writes to the CWD by default.
RNAalifold \
    --ribosum_scoring \
    -d0 \
    -p \
    --noPS \
    "$ALIGNMENT" > "${OUTPUT_PREFIX}_rnaalifold.txt"

echo "Consensus structure and energy:"
head -3 "${OUTPUT_PREFIX}_rnaalifold.txt"

echo ""
echo "=== Structure Conservation Index (RNAz, if installed) ==="
# RNAz reads Clustal-W or MAF, NOT Stockholm -- convert first (esl-reformat clustal in.sto > in.aln).
# SCI ~ consensus MFE / mean single-sequence MFE; ~1.0 = a conserved structure.
# A negative z-score / high SCI means "more stable than random", NOT "the structure is correct" --
# statistical covariation support (R-scape) is the stronger evidence standard.
if command -v RNAz >/dev/null 2>&1; then
    ALN_CLUSTAL="$ALIGNMENT"
    if command -v esl-reformat >/dev/null 2>&1; then
        esl-reformat clustal "$ALIGNMENT" > "${OUTPUT_PREFIX}.aln" 2>/dev/null && ALN_CLUSTAL="${OUTPUT_PREFIX}.aln"
    fi
    RNAz "$ALN_CLUSTAL" 2>/dev/null | grep -E "SCI|z-score|SVM RNA-class probability" \
        || echo "  RNAz failed (it needs Clustal/MAF input, not Stockholm)."
else
    echo "  RNAz not installed; skipping SCI."
fi

echo ""
echo "Next: validate the consensus pairs with R-scape (covariation-analysis) before trusting them."
