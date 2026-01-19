#!/bin/bash
# Assembly QC with QUAST and BUSCO
set -euo pipefail

ASSEMBLY=$1
OUTDIR=${2:-assembly_qc}
LINEAGE=${3:-bacteria_odb10}
THREADS=${4:-8}

mkdir -p $OUTDIR

echo "=== Assembly QC Pipeline ==="
echo "Assembly: $ASSEMBLY"
echo "Lineage: $LINEAGE"

# QUAST
echo ""
echo "Running QUAST..."
quast.py $ASSEMBLY -o ${OUTDIR}/quast -t $THREADS

# BUSCO
echo ""
echo "Running BUSCO..."
busco -i $ASSEMBLY -m genome -l $LINEAGE -o ${OUTDIR}/busco -c $THREADS --offline

# Print summaries
echo ""
echo "=========================================="
echo "QUAST Summary"
echo "=========================================="
grep -E "^(# contigs|Total length|Largest contig|N50|N90|L50|GC)" ${OUTDIR}/quast/report.txt

echo ""
echo "=========================================="
echo "BUSCO Summary"
echo "=========================================="
cat ${OUTDIR}/busco/short_summary*.txt | grep -E "^(\s+C:|Results)"

echo ""
echo "Full reports in: $OUTDIR"
