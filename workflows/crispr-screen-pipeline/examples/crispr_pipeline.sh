#!/bin/bash
# Complete CRISPR screen analysis pipeline

# Configuration
LIBRARY="library.csv"
OUTPUT_PREFIX="screen_analysis"

# Step 1: Guide counting from FASTQ
echo "Step 1: Counting guides..."
mageck count \
    -l ${LIBRARY} \
    -n ${OUTPUT_PREFIX} \
    --sample-label Day0,Day14_Rep1,Day14_Rep2,Day14_Rep3 \
    --fastq Day0.fastq.gz Day14_Rep1.fastq.gz Day14_Rep2.fastq.gz Day14_Rep3.fastq.gz \
    --trim-5 0 \
    --pdf-report

# Step 2: MAGeCK RRA analysis
echo "Step 2: Running MAGeCK test..."
mageck test \
    -k ${OUTPUT_PREFIX}.count.txt \
    -t Day14_Rep1,Day14_Rep2,Day14_Rep3 \
    -c Day0 \
    -n ${OUTPUT_PREFIX}_rra \
    --pdf-report \
    --gene-lfc-method alphamedian

# Step 3: Extract hits
echo "Step 3: Calling hits..."
python3 << 'EOF'
import pandas as pd

gene_summary = pd.read_csv('screen_analysis_rra.gene_summary.txt', sep='\t')

# Negative selection hits (dropout)
neg_hits = gene_summary[(gene_summary['neg|fdr'] < 0.05) & (gene_summary['neg|lfc'] < -0.5)]
neg_hits = neg_hits.sort_values('neg|rank')

# Positive selection hits (enriched)
pos_hits = gene_summary[(gene_summary['pos|fdr'] < 0.05) & (gene_summary['pos|lfc'] > 0.5)]
pos_hits = pos_hits.sort_values('pos|rank')

print(f'Negative selection hits: {len(neg_hits)}')
print(f'Positive selection hits: {len(pos_hits)}')

neg_hits.to_csv('negative_hits.csv', index=False)
pos_hits.to_csv('positive_hits.csv', index=False)
EOF

echo "Pipeline complete!"
echo "Results: ${OUTPUT_PREFIX}_rra.gene_summary.txt"
echo "Hits: negative_hits.csv, positive_hits.csv"
