#!/bin/bash
# Reference: Foldseek 9+ | Verify API if version differs
# Structure-aware homology search via Foldseek; with ProstT5 path for sequence-only queries.

set -euo pipefail

QUERY="${1:-query.pdb}"          # or query.fa with ProstT5
DB_DIR="${2:-foldseek_dbs}"
TMP_DIR="${3:-foldseek_tmp}"

mkdir -p "${DB_DIR}" "${TMP_DIR}"

# Download AFDB Swiss-Prot subset (smaller, faster than full AlphaFoldDB).
if [ ! -d "${DB_DIR}/afdb_sp" ]; then
    echo "Downloading AlphaFoldDB Swiss-Prot subset..."
    foldseek databases Alphafold/Swiss-Prot "${DB_DIR}/afdb_sp" "${TMP_DIR}"
fi

if [[ "${QUERY}" == *.pdb ]] || [[ "${QUERY}" == *.cif ]]; then
    echo "=== Structure-to-structure search ==="
    foldseek easy-search "${QUERY}" "${DB_DIR}/afdb_sp" results.m8 "${TMP_DIR}" \
        --format-output query,target,fident,alnlen,evalue,bits,prob,qtmscore,ttmscore,lddt \
        --threads 8
else
    # Sequence-only path via ProstT5: predict 3Di alphabet directly from sequence, skip AF2.
    echo "=== Sequence-only search via ProstT5 ==="
    if [ ! -d "${DB_DIR}/prostt5" ]; then
        echo "Downloading ProstT5 weights..."
        foldseek databases ProstT5 "${DB_DIR}/prostt5" "${TMP_DIR}"
    fi
    foldseek easy-search "${QUERY}" "${DB_DIR}/afdb_sp" results.m8 "${TMP_DIR}" \
        --prostt5-model "${DB_DIR}/prostt5" \
        --format-output query,target,fident,alnlen,evalue,bits,prob,qtmscore,ttmscore,lddt \
        --threads 8
fi

echo
echo "=== Top hits (filtered prob > 0.9) ==="
echo "query target fident alnlen evalue bits prob qtmscore ttmscore lddt"
awk -F'\t' '$7 > 0.9' results.m8 | sort -k6,6 -gr | head -10

echo
echo "Note: fold-level similarity is necessary but not sufficient for homology."
echo "Cross-check with sequence evidence and conserved functional residues."
