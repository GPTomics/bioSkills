#!/bin/bash
# Reference: ABC-Enhancer-Gene-Prediction 0.2.2+, samtools 1.19+, bedtools 2.31+, deepTools 3.5+ | Verify API if version differs
# ABC pipeline: ATAC + H3K27ac + Hi-C/Micro-C -> per-(enhancer, gene) regulatory scores.
# Reference: Fulco 2019 Nat Genet; Nasser 2021 Nature.

set -euo pipefail

ATAC_BAM=${1:-atac.dedup.bam}
H3K27AC_BAM=${2:-h3k27ac.dedup.bam}
HIC_DIR=${3:-hic_data/}                  # Cooler or hic format directory
GENOME_FA=${4:-hg38.fa}
SIZES=${5:-hg38.chrom.sizes}
GENE_BED=${6:-refseq_protein_coding.bed}
ABC_REPO=${7:-/path/to/ABC-Enhancer-Gene-Prediction}    # broadinstitute/ABC-Enhancer-Gene-Prediction
CELL_TYPE=${8:-K562}
OUTDIR=${9:-abc_out}
EFFECTIVE_GENOME=${10:-2913022398}        # hg38 100bp reads (deepTools)

mkdir -p $OUTDIR/{tracks,peaks,neighborhoods,predictions}

# 1. Generate normalized signal tracks
bamCoverage --bam $ATAC_BAM \
    --outFileName $OUTDIR/tracks/atac.bw \
    --binSize 50 --normalizeUsing RPGC \
    --effectiveGenomeSize $EFFECTIVE_GENOME \
    --numberOfProcessors 8

bamCoverage --bam $H3K27AC_BAM \
    --outFileName $OUTDIR/tracks/h3k27ac.bw \
    --binSize 50 --normalizeUsing RPGC \
    --effectiveGenomeSize $EFFECTIVE_GENOME \
    --numberOfProcessors 8

# 2. Define candidate enhancers (MACS narrowPeak, exclude promoters)
# (assume MACS3 already run upstream)
bedtools intersect -v -a atac_peaks.narrowPeak \
    -b promoter_regions.bed \
    > $OUTDIR/peaks/candidate_enhancers.bed

# 3. ABC neighborhoods: compute Activity per enhancer
# Path layout differs by ABC version: legacy = src/run.neighborhoods.py; Snakemake-based = workflow/scripts/run.neighborhoods.py
# Verify before running: `find $ABC_REPO -name run.neighborhoods.py`
python $ABC_REPO/workflow/scripts/run.neighborhoods.py \
    --candidate_enhancer_regions $OUTDIR/peaks/candidate_enhancers.bed \
    --genes $GENE_BED \
    --H3K27ac $OUTDIR/tracks/h3k27ac.bw \
    --DHS $OUTDIR/tracks/atac.bw \
    --chrom_sizes $SIZES \
    --ubiquitously_expressed_genes ubiquitously_expressed.txt \
    --cellType $CELL_TYPE \
    --outdir $OUTDIR/neighborhoods/

# 4. ABC predictions: Activity * Contact, normalize, threshold
python $ABC_REPO/workflow/scripts/predict.py \
    --enhancers $OUTDIR/neighborhoods/EnhancerList.txt \
    --genes $OUTDIR/neighborhoods/GeneList.txt \
    --HiCdir $HIC_DIR \
    --hic_resolution 5000 \
    --score_column ABC.Score \
    --threshold 0.02 \
    --cellType $CELL_TYPE \
    --outdir $OUTDIR/predictions/

# 5. Filter and report
echo "ABC predictions: $OUTDIR/predictions/EnhancerPredictionsAllPutative.txt"
echo "Above threshold (ABC.Score >= 0.02):"
awk -F'\t' 'NR > 1 && $11 >= 0.02' $OUTDIR/predictions/EnhancerPredictionsAllPutative.txt | wc -l

# Optional: cross-validation
# Compare against published CRISPRi-FlowFISH catalog (Fulco 2019 K562)
# wget https://www.engreitzlab.org/crispri-flowfish/K562_validated_pairs.txt
# Compute sensitivity / specificity at threshold 0.02
