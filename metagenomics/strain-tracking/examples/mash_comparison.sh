#!/bin/bash
# Reference: skani 0.2+, Mash 2.3+ | Verify API if version differs
# GENOME-vs-GENOME ANI for isolate/MAG comparison and dereplication - NOT in-situ strain resolution.
# ANI saturates above ~99.9% and cannot call same-vs-different strain; use inStrain/StrainPhlAn for that.
set -euo pipefail

INPUT_DIR=${1:-.}
OUTPUT_PREFIX=${2:-ani_results}

# skani is the modern fastANI replacement: faster and robust on FRAGMENTED MAGs (where fastANI degrades).
# Use -E for a sparse EDGE list (one row per pair: Ref Query ANI Af_ref Af_query Ref_name Query_name);
# the default `skani triangle` emits a lower-triangular MATRIX, which is not a pairwise table.
echo "=== skani all-vs-all ANI (edge list) ==="
skani triangle -E "${INPUT_DIR}"/*.fasta -o "${OUTPUT_PREFIX}_skani.tsv"

# ~95% ANI is the species boundary (Jain 2018). Above it, ANI cannot resolve strains - do not
# label high-ANI pairs "same strain"; that needs microdiversity-aware popANI/nGD.
# (-E emits a header row like skani dist; the numeric $3 >= 95 test below skips it - "ANI" coerces to 0.)
echo ""
echo "Pairs above the ~95% species boundary (same species, NOT necessarily same strain):"
awk -F'\t' '$3 >= 95 {print $1, $2, $3"%"}' "${OUTPUT_PREFIX}_skani.tsv" | head -20

# Mash remains useful for very fast distance triage of large genome sets:
# mash sketch -o "${OUTPUT_PREFIX}.msh" "${INPUT_DIR}"/*.fasta
# mash dist "${OUTPUT_PREFIX}.msh" "${OUTPUT_PREFIX}.msh" > "${OUTPUT_PREFIX}_mash.tsv"
echo "Results: ${OUTPUT_PREFIX}_skani.tsv (genome-vs-genome ANI; for strains use inStrain)."
