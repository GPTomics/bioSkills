#!/bin/bash
# Reference: NCBI Datasets CLI 16.0+ | Verify API if version differs
# Bulk gene metadata + ortholog set via Datasets summary + dataformat TSV.

set -euo pipefail

SYMBOL="${1:-BRCA1}"
TAXON="${2:-Mammalia}"

echo "=== Gene records for ${SYMBOL} across ${TAXON} ==="
datasets summary gene symbol "${SYMBOL}" \
    --taxon "${TAXON}" \
    --as-json-lines \
  | dataformat tsv gene \
        --fields gene-id,symbol,taxname,description,chromosomes,nomenclature-authority-symbol \
  > "${SYMBOL}_${TAXON// /_}.tsv"

head "${SYMBOL}_${TAXON// /_}.tsv" | column -t -s $'\t'

echo
echo "=== NCBI ortholog set for ${SYMBOL} (one rep per species) ==="
datasets summary gene symbol "${SYMBOL}" --taxon human --ortholog --as-json-lines \
  | dataformat tsv gene \
        --fields gene-id,symbol,taxname,description \
  > "${SYMBOL}_orthologs.tsv"

head "${SYMBOL}_orthologs.tsv" | column -t -s $'\t'

echo
echo "Note: NCBI ortholog set is one rep per species. For tree-reconciled orthology"
echo "with co-orthologs and 1:many calls, use ortholog-inference (Compara, OMA)."

echo
echo "=== Inspect raw JSON-lines for fields available ==="
echo "datasets summary gene symbol ${SYMBOL} --taxon human --as-json-lines | jq -s 'first | keys' "
echo "Use field names from there in dataformat --fields"
