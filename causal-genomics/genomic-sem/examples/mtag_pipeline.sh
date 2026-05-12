#!/usr/bin/env bash
# Reference: MTAG 1.0.8+, Python 3.7+ | Verify API if version differs
#
# MTAG pipeline for multi-trait analysis of GWAS. Run on the same input sumstats
# used in the GenomicSEM common-factor GWAS for cross-method comparison.
#
# Expects each input sumstats file to have columns:
#   SNP, A1, A2, BETA (or Z), SE, P, N, EAF (or FRQ)
#
# MaxFDR > 5% per trait indicates MTAG assumptions are violated for that trait;
# prefer GenomicSEM with Q_SNP filtering in that regime (see SKILL.md reconciliation).

set -euo pipefail

MTAG_BIN="${MTAG_BIN:-./mtag/mtag.py}"
OUTPUT_DIR="${OUTPUT_DIR:-results/mtag}"

TRAIT1='data/trait1.mtag_in.txt'
TRAIT2='data/trait2.mtag_in.txt'
TRAIT3='data/trait3.mtag_in.txt'

mkdir -p "${OUTPUT_DIR}"

# Step 1: Run MTAG across all three traits
python "${MTAG_BIN}" \
    --sumstats "${TRAIT1},${TRAIT2},${TRAIT3}" \
    --out "${OUTPUT_DIR}/mtag" \
    --n_min 0 \
    --use_beta_se \
    --stream_stdout \
    --force

# Step 2: Check MaxFDR per trait (Turley 2018 threshold = 5%)
echo ''
echo '=== MaxFDR per trait ==='
grep -E '^.*maxFDR' "${OUTPUT_DIR}/mtag.log" || echo 'maxFDR not present; check MTAG version (>= 1.0.7 required)'

# Step 3: Extract genome-wide significant hits per trait
# MTAG writes per-trait sumstats as <out>_trait_<i>.txt; the mtag P column header is
# 'mtag_pval' (not the last field). Resolve its index from the header per file.
for i in 1 2 3; do
    OUT="${OUTPUT_DIR}/mtag_trait_${i}.txt"
    HITS="${OUTPUT_DIR}/mtag_trait_${i}.gws.tsv"
    awk 'NR==1 { for (k=1; k<=NF; k++) if ($k == "mtag_pval") pcol = k; print; next }
         pcol && $pcol != "" && $pcol+0 < 5e-08' "${OUT}" > "${HITS}"
    echo "Trait ${i} genome-wide significant SNPs (MTAG): $(($(wc -l < "${HITS}")-1))"
done

# Step 4: Compare with GenomicSEM common-factor SNPs
GSEM_FACTOR='results/common_factor_snps.tsv'
if [[ -f "${GSEM_FACTOR}" ]]; then
    echo ''
    echo '=== GenomicSEM vs MTAG overlap ==='
    for i in 1 2 3; do
        HITS="${OUTPUT_DIR}/mtag_trait_${i}.gws.tsv"
        OVERLAP=$(awk 'NR>1 {print $1}' "${HITS}" | \
            grep -F -f - <(awk 'NR>1 {print $1}' "${GSEM_FACTOR}") | wc -l)
        echo "Trait ${i} MTAG GWS that are also GenomicSEM common-factor SNPs: ${OVERLAP}"
    done
fi

# Step 5: Optional bivariate LDSC intercept check (sample overlap)
# (Run LDSC separately; intercept_2 > 0 indicates overlap)
echo ''
echo 'Done. Remember to inspect cross-trait LDSC intercept values for sample overlap.'
