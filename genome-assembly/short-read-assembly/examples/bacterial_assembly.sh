#!/bin/bash
# Bacterial genome assembly with SPAdes
set -euo pipefail

R1=$1
R2=$2
OUTDIR=$3
THREADS=${4:-16}

echo "=== Bacterial Genome Assembly ==="
echo "R1: $R1"
echo "R2: $R2"
echo "Output: $OUTDIR"

# Run SPAdes
spades.py \
    --isolate \
    --careful \
    -t $THREADS \
    -1 $R1 \
    -2 $R2 \
    -o $OUTDIR

# Check output
if [ -f "${OUTDIR}/scaffolds.fasta" ]; then
    echo ""
    echo "=== Assembly Statistics ==="
    echo "Scaffolds: $(grep -c '^>' ${OUTDIR}/scaffolds.fasta)"

    # Basic N50 calculation
    awk '/^>/{if(l) print l; l=0; next}{l+=length}END{print l}' ${OUTDIR}/scaffolds.fasta | \
        sort -rn | awk '{sum+=$1; lens[NR]=$1} END{
            for(i=1;i<=NR;i++){
                cumsum+=lens[i]
                if(cumsum>=sum/2){print "N50:", lens[i]; break}
            }
        }'

    echo ""
    echo "Assembly complete: ${OUTDIR}/scaffolds.fasta"
else
    echo "ERROR: Assembly failed"
    exit 1
fi
