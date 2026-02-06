#!/bin/bash
# Consensus RNA secondary structure prediction from a multiple sequence alignment.

ALIGNMENT=$1
OUTPUT_PREFIX=${2:-"consensus"}

if [ -z "$ALIGNMENT" ]; then
    echo "Usage: $0 <alignment.sto|alignment.aln> [output_prefix]"
    echo ""
    echo "Accepts Stockholm (.sto) or ClustalW (.aln) format alignments."
    exit 1
fi

echo "=== RNAalifold: Consensus structure ==="
RNAalifold \
    --ribosum_scoring \
    --cfactor 0.6 \
    --nfactor 0.5 \
    -p \
    --aln-stk="${OUTPUT_PREFIX}.sto" \
    "$ALIGNMENT" > "${OUTPUT_PREFIX}_rnaalifold.txt"

echo "MFE consensus structure:"
head -3 "${OUTPUT_PREFIX}_rnaalifold.txt"

echo ""
echo "=== SCI (Structure Conservation Index) ==="
# SCI > 0.5 suggests conserved structure; SCI > 1.0 indicates
# covariation supports the structure beyond thermodynamic stability
RNAz --no-shuffle "$ALIGNMENT" 2>/dev/null | grep -E "SCI|P-value|z-score"

echo ""
echo "Output files:"
echo "  ${OUTPUT_PREFIX}_rnaalifold.txt - Consensus structure and energy"
echo "  ${OUTPUT_PREFIX}.sto - Annotated Stockholm alignment"
echo "  alidot.ps - Base-pair probability dot plot (if generated)"
echo "  alirna.ps - Structure drawing (if generated)"
