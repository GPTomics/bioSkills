#!/bin/bash
# Reference: DIA-NN 1.9+ | Verify API if version differs
# DIA-NN predicted-library (directDIA) analysis: digest FASTA, predict a library,
# search, and filter the parquet report at the correct q-value level and context.
# Shows the commands with rationale; does not produce stray outputs on its own.

set -e

FASTA="uniprot_human_reviewed.fasta"
OUTPUT_DIR="diann_out"
THREADS=8

# QVALUE 0.01 = standard 1% precursor FDR (run context).
# MISSED_CLEAVAGES 1 = trypsin standard; higher expands the search/FDR burden.
# --mass-acc 0 = auto-optimize tolerances per file (do not hard-code ppm from another instrument).
# --reanalyse = two-pass global FDR / MBR that controls directDIA double-dipping.
# --fixed-mod UniMod:4,57.021464,C = fixed Cys carbamidomethylation (explicit form of the --unimod4 preset).
# --var-mod UniMod:35,15.994915,M = variable Met oxidation.
QVALUE=0.01
MISSED_CLEAVAGES=1

mkdir -p "$OUTPUT_DIR"

MZML_ARGS=""
for f in *.mzML; do
    MZML_ARGS="$MZML_ARGS --f $f"
done

echo "Running DIA-NN predicted-library analysis..."
diann \
    $MZML_ARGS \
    --lib "" --fasta "$FASTA" --fasta-search \
    --gen-spec-lib --predictor \
    --out "$OUTPUT_DIR/report.parquet" \
    --out-lib "$OUTPUT_DIR/report-lib.tsv" \
    --qvalue $QVALUE \
    --matrices \
    --mass-acc 0 \
    --reanalyse --smart-profiling \
    --cut "K*,R*" --missed-cleavages $MISSED_CLEAVAGES \
    --min-pep-len 7 --max-pep-len 30 \
    --fixed-mod UniMod:4,57.021464,C --var-mods 1 --var-mod UniMod:35,15.994915,M \
    --threads $THREADS

echo "Done. Main report (1.9+ default): $OUTPUT_DIR/report.parquet"
echo "Matrices: $OUTPUT_DIR/report.pg_matrix.tsv (verify *_matrix dotting vs installed version)"

# Filter the parquet report. Points to note for downstream code:
#   - read report.parquet, NOT report.tsv (1.9+ default).
#   - filter BOTH levels: Q.Value (precursor) AND PG.Q.Value (protein-group).
#   - for a cross-run matrix add Global.PG.Q.Value <= 0.01 (per-run FDR inflates across runs).
#   - the *_matrix.tsv files apply an extra 5% run-specific PG filter, so matrix count < report count is expected.
#   - convert DIA-NN's 0 (not-quantified) to NA before log2 / normalization.
