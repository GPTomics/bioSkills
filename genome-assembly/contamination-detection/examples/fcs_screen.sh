#!/bin/bash
# Reference: FCS-adaptor 0.5+, FCS-GX 0.5+ | Verify API if version differs
# Foreign-sequence screen for a SINGLE-ORGANISM assembly (eukaryote or isolate).
# This is PROBLEM 1. For PROBLEM 2 (MAG quality) use mag_qc_pipeline.sh instead.
# Crossing the two (FCS-GX flagged-length read as MIMAG contamination %) is a category error.
set -euo pipefail

ASSEMBLY=$1          # gzipped FASTA, e.g. assembly.fa.gz
TAX_ID=$2            # NCBI taxid of the SOURCE organism (NOT a contaminant)
GXDB="$3"            # path to the GX database, e.g. $GXDB_LOC/gxdb
OUTDIR=${4:-fcs_out}
LINEAGE=${5:-euk}    # --euk for eukaryotes, --prok for prokaryotes

mkdir -p "$OUTDIR"

# Step 1: FCS-adaptor FIRST - tiny DB, trivial RAM, unambiguous adaptor/vector hits.
echo "=== FCS-adaptor (run first; cheap) ==="
run_fcsadaptor.sh --fasta-input "$ASSEMBLY" --output-dir "$OUTDIR/adaptor" "--$LINEAGE"

# Step 2: FCS-GX - GenBank-mandatory cross-taxon screen.
# RAM WALL: the GX DB is ~470 GiB and wants a ~512 GiB-RAM host (or a tmpfs RAM disk / cloud VM).
# Underprovisioned it crawls (minutes become days), it does not error out.
echo "=== FCS-GX (needs ~512 GiB RAM; tax-id is the SOURCE organism) ==="
python3 ./fcs.py screen genome --fasta "$ASSEMBLY" --out-dir "$OUTDIR/gx" \
    --gx-db "$GXDB" --tax-id "$TAX_ID"

# Step 3: Apply ONLY auto-clean actions (EXCLUDE/TRIM/FIX). REVIEW items are left for the human -
# that tier is exactly where legitimate HGT/endosymbiont sequence hides (the tardigrade lesson).
REPORT="$OUTDIR/gx/$(basename "${ASSEMBLY%.gz}").${TAX_ID}.fcs_gx_report.txt"   # FCS-GX strips only .gz, keeps .fa
echo "=== Clean (EXCLUDE/TRIM/FIX only; REVIEW left for manual inspection) ==="
zcat "$ASSEMBLY" | python3 ./fcs.py clean genome \
    --action-report "$REPORT" \
    --output "$OUTDIR/clean.fasta" \
    --contam-fasta-out "$OUTDIR/contam.fasta"

echo "=== Done ==="
echo "Action report: $REPORT"
echo "Inspect REVIEW rows by hand. Before removing any foreign-flagged contig, check physical"
echo "integration (host flanks, host GC/coverage, host introns) - foreign-looking != contaminant."
echo "Next: build a BlobToolKit blob plot (GC x coverage x taxonomy) for the contigs in doubt."
