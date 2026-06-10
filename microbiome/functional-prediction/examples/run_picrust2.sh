#!/bin/bash
# Reference: PICRUSt2 2.5+ | Verify API if version differs
# Predict community functional POTENTIAL from 16S ASVs with PICRUSt2.
# This is PREDICTED gene-content potential interpolated from reference genomes,
# never measured gene content and never activity. NSTI gates and quantifies it.
set -euo pipefail

ASV_SEQS='asv_seqs.fna'        # representative ASV sequences (FASTA), from amplicon-processing
ASV_TABLE='asv_table.tsv'      # ASV abundance table (TSV or BIOM, samples as columns)
OUTPUT_DIR='picrust2_out'
THREADS=8
MAX_NSTI=2.0                   # PICRUSt2 default; ASVs ABOVE this are dropped before inference
                               # (their nearest sequenced genome is too distant to trust)
HSP_METHOD='mp'                # maximum parsimony, the recommended default; pic is faster but not recommended

picrust2_pipeline.py \
    -s "$ASV_SEQS" \
    -i "$ASV_TABLE" \
    -o "$OUTPUT_DIR" \
    -p "$THREADS" \
    --hsp_method "$HSP_METHOD" \
    --max_nsti "$MAX_NSTI" \
    --verbose

# Attach human-readable names to the MetaCyc pathway table
add_descriptions.py \
    -i "$OUTPUT_DIR/pathways_out/path_abun_unstrat.tsv.gz" \
    -m METACYC \
    -o "$OUTPUT_DIR/pathways_out/path_abun_described.tsv.gz"

# Mandatory: report the NSTI distribution and the fraction of reads dropped by the gate.
# The real quality file is marker_predicted_and_nsti.tsv.gz; the column is metadata_NSTI.
# A run that loses a large read fraction predicted function for a different community than was sampled.
python3 - "$OUTPUT_DIR/marker_predicted_and_nsti.tsv.gz" "$ASV_TABLE" "$MAX_NSTI" <<'PY'
import sys
import pandas as pd

nsti_path, table_path, max_nsti = sys.argv[1], sys.argv[2], float(sys.argv[3])
nsti = pd.read_csv(nsti_path, sep='\t').set_index('sequence')
counts = pd.read_csv(table_path, sep='\t', index_col=0)
reads_per_asv = counts.sum(axis=1)

dropped = nsti.index[nsti['metadata_NSTI'] > max_nsti]
reads_dropped_frac = reads_per_asv.reindex(dropped).sum() / reads_per_asv.sum()
print(f'mean NSTI {nsti.metadata_NSTI.mean():.3f}  median {nsti.metadata_NSTI.median():.3f}')
print(f'ASVs dropped at NSTI>{max_nsti}: {len(dropped)}/{len(nsti)}  reads dropped: {reads_dropped_frac:.1%}')
PY

echo 'Predicted POTENTIAL only - report NSTI, restrict claims to capacity, not activity.'
