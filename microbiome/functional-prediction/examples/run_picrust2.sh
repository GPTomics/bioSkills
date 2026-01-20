#!/bin/bash
# Run PICRUSt2 functional prediction pipeline

# Input files (exported from R/phyloseq)
ASV_SEQS="asv_seqs.fasta"
ASV_TABLE="asv_table.tsv"
OUTPUT_DIR="picrust2_output"
THREADS=4

# Check inputs exist
if [ ! -f "$ASV_SEQS" ]; then
    echo "Error: $ASV_SEQS not found"
    exit 1
fi

if [ ! -f "$ASV_TABLE" ]; then
    echo "Error: $ASV_TABLE not found"
    exit 1
fi

echo "Running PICRUSt2 pipeline..."
echo "  ASV sequences: $ASV_SEQS"
echo "  ASV table: $ASV_TABLE"
echo "  Output: $OUTPUT_DIR"
echo "  Threads: $THREADS"

# Run full pipeline
picrust2_pipeline.py \
    -s "$ASV_SEQS" \
    -i "$ASV_TABLE" \
    -o "$OUTPUT_DIR" \
    -p "$THREADS" \
    --stratified \
    --per_sequence_contrib \
    --verbose

# Check if successful
if [ $? -eq 0 ]; then
    echo ""
    echo "PICRUSt2 completed successfully!"
    echo ""
    echo "Output files:"
    echo "  KEGG orthologs: $OUTPUT_DIR/KO_metagenome_out/pred_metagenome_unstrat.tsv"
    echo "  EC numbers: $OUTPUT_DIR/EC_metagenome_out/pred_metagenome_unstrat.tsv"
    echo "  Pathways: $OUTPUT_DIR/pathways_out/path_abun_unstrat.tsv"
    echo ""

    # Add descriptions to pathway output
    echo "Adding pathway descriptions..."
    add_descriptions.py \
        -i "$OUTPUT_DIR/pathways_out/path_abun_unstrat.tsv" \
        -m METACYC \
        -o "$OUTPUT_DIR/pathways_out/path_abun_described.tsv"

    # Report NSTI statistics
    echo ""
    echo "NSTI statistics (quality metric - lower is better):"
    python3 -c "
import pandas as pd
nsti = pd.read_csv('$OUTPUT_DIR/marker_predicted_and_nsti.tsv.gz', sep='\t')
print(f'  Mean NSTI: {nsti[\"metadata_NSTI\"].mean():.3f}')
print(f'  Median NSTI: {nsti[\"metadata_NSTI\"].median():.3f}')
print(f'  ASVs with NSTI > 2 (unreliable): {(nsti[\"metadata_NSTI\"] > 2).sum()}/{len(nsti)}')
"
else
    echo "Error: PICRUSt2 pipeline failed"
    exit 1
fi
