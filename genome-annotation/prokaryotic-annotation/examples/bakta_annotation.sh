#!/bin/bash
# Prokaryotic genome annotation with Bakta
set -euo pipefail

ASSEMBLY=$1
DB_PATH=${2:-/path/to/bakta_db}
OUTDIR=${3:-bakta_out}
PREFIX=${4:-my_genome}
LOCUS_TAG=${5:-MYORG}
THREADS=${6:-8}

mkdir -p $OUTDIR

echo "=== Bakta Prokaryotic Annotation ==="
echo "Assembly: $ASSEMBLY"
echo "Database: $DB_PATH"

# Download database if not present
if [ ! -d "$DB_PATH" ]; then
    echo "Downloading Bakta database (full)..."
    bakta_db download --output $DB_PATH --type full
fi

# Run Bakta annotation
echo ""
echo "Running Bakta..."
bakta \
    --db $DB_PATH \
    --output $OUTDIR \
    --prefix $PREFIX \
    --locus-tag $LOCUS_TAG \
    --threads $THREADS \
    $ASSEMBLY

# Verify output files
echo ""
echo "=========================================="
echo "Output Files"
echo "=========================================="
ls -lh ${OUTDIR}/${PREFIX}.*

# Annotation summary
echo ""
echo "=========================================="
echo "Annotation Summary"
echo "=========================================="
CDS_COUNT=$(grep -c $'\tCDS\t' ${OUTDIR}/${PREFIX}.gff3 || echo 0)
TRNA_COUNT=$(grep -c $'\ttRNA\t' ${OUTDIR}/${PREFIX}.gff3 || echo 0)
RRNA_COUNT=$(grep -c $'\trRNA\t' ${OUTDIR}/${PREFIX}.gff3 || echo 0)
HYPO_COUNT=$(grep -c 'hypothetical protein' ${OUTDIR}/${PREFIX}.tsv || echo 0)

echo "CDSs: $CDS_COUNT"
echo "tRNAs: $TRNA_COUNT"
echo "rRNAs: $RRNA_COUNT"
echo "Hypothetical proteins: $HYPO_COUNT"

if [ $CDS_COUNT -gt 0 ]; then
    HYPO_PCT=$(echo "scale=1; $HYPO_COUNT * 100 / $CDS_COUNT" | bc)
    echo "Hypothetical fraction: ${HYPO_PCT}%"
fi

# Run BUSCO on predicted proteins
echo ""
echo "Running BUSCO on predicted proteins..."
busco -i ${OUTDIR}/${PREFIX}.faa -m proteins -l bacteria_odb10 -o ${OUTDIR}/busco_proteins -c $THREADS --offline 2>/dev/null || echo "BUSCO skipped (database not available)"

echo ""
echo "Annotation complete. Results in: $OUTDIR"
