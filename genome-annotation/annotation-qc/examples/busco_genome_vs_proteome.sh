#!/bin/bash
# Reference: busco 5.5+, omark 0.3+ | Verify API if version differs
# The diagnostic fork: assembly-BUSCO vs proteome-BUSCO on the same assembly.
# A large gap (assembly high, proteome low) means the predictor missed present
# genes (fixable: evidence/masking/training). Both low means fix the assembly.
set -euo pipefail

GENOME=$1
PROTEINS=$2
LINEAGE=${3:-vertebrata_odb10}   # use the DEEPEST applicable clade dataset, not the shallow eukaryota_odb10
THREADS=${4:-16}

echo "=== Proteome BUSCO (the delivered annotation) ==="
busco -i "$PROTEINS" -m proteins -l "$LINEAGE" -o busco_proteins -c "$THREADS" --offline 2>/dev/null \
    || echo "BUSCO proteins skipped (lineage dataset not available offline)"

echo "=== Genome BUSCO (the assembly) ==="
busco -i "$GENOME" -m genome -l "$LINEAGE" -o busco_genome -c "$THREADS" --offline 2>/dev/null \
    || echo "BUSCO genome skipped (lineage dataset not available offline)"

echo ""
echo "Interpretation:"
echo "  proteome << genome   -> predictor missed present genes (fix evidence/masking/training)"
echo "  both low             -> fix the assembly (more data, purge_dups, decontaminate)"
echo "  Duplicated > 5-8% (no known WGD) -> purge_dups the assembly BEFORE annotating"

echo ""
echo "=== OMArk (proteome consistency + contamination; complements BUSCO) ==="
# Requires an OMAmer database (e.g. LUCA.h5). OMArk catches over-prediction/chimeras/contaminants BUSCO misses.
if [ -f "${OMAMER_DB:-LUCA.h5}" ]; then
    omamer search --db "${OMAMER_DB:-LUCA.h5}" --query "$PROTEINS" --out proteins.omamer
    omark -f proteins.omamer -d "${OMAMER_DB:-LUCA.h5}" -o omark_out
else
    echo "OMArk skipped (set OMAMER_DB to an OMAmer database to enable)"
fi
