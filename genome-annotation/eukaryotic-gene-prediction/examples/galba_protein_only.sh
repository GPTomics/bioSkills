#!/bin/bash
# Reference: galba 1.0.11+, busco 5.5+ | Verify API if version differs
# Eukaryotic gene prediction with GALBA using protein-only evidence
set -euo pipefail

GENOME=$1
PROTEINS=$2
OUTDIR=${3:-galba_out}
SPECIES=${4:-my_species}
THREADS=${5:-16}
LINEAGE=${6:-vertebrata_odb10}   # set to the DEEPEST applicable clade dataset, not the shallow eukaryota_odb10

echo "=== GALBA Protein-Only Gene Prediction ==="
echo "Genome: $GENOME"
echo "Proteins: $PROTEINS"
echo ""
echo "NOTE: GALBA works best with proteins from closely related species."
echo "If RNA-seq data is available, prefer BRAKER3 for better accuracy."

# Verify softmasking
LOWER=$(grep -v '^>' $GENOME | tr -cd 'a-z' | wc -c)
TOTAL=$(grep -v '^>' $GENOME | tr -cd 'a-zA-Z' | wc -c)
MASK_PCT=$(echo "scale=2; $LOWER * 100 / $TOTAL" | bc)
echo "Repeat masking: ${MASK_PCT}%"

if [ "$(echo "$MASK_PCT < 1" | bc)" -eq 1 ]; then
    echo "WARNING: Very low masking (<1%). Run RepeatMasker first."
    exit 1
fi

# Count input proteins
PROT_COUNT=$(grep -c '^>' $PROTEINS)
echo "Input proteins: $PROT_COUNT"

# Run GALBA
echo ""
echo "Running GALBA..."
# GALBA always treats the genome as soft-masked for repeats (no --softmasking flag, unlike BRAKER)
galba.pl \
    --genome=$GENOME \
    --prot_seq=$PROTEINS \
    --threads=$THREADS \
    --species=$SPECIES \
    --workingdir=$OUTDIR \
    --gff3

# Summary
echo ""
echo "=========================================="
echo "GALBA Prediction Summary"
echo "=========================================="
GENE_COUNT=$(grep -c $'\tgene\t' ${OUTDIR}/galba.gff3 || echo 0)
MRNA_COUNT=$(grep -c $'\tmRNA\t' ${OUTDIR}/galba.gff3 || echo 0)

echo "Genes: $GENE_COUNT"
echo "mRNAs: $MRNA_COUNT"

# BUSCO evaluation
echo ""
echo "Running BUSCO on predicted proteins..."
busco -i ${OUTDIR}/galba.aa -m proteins -l $LINEAGE -o ${OUTDIR}/busco_eval -c $THREADS --offline 2>/dev/null || echo "BUSCO skipped"

echo ""
echo "Results in: $OUTDIR"
