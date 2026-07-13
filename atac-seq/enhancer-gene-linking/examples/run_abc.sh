#!/bin/bash
# Reference: ABC-Enhancer-Gene-Prediction 0.2.2+, samtools 1.19+, bedtools 2.31+, deepTools 3.5+ | Verify API if version differs
# ABC pipeline: ATAC + H3K27ac + Hi-C/Micro-C -> per-(enhancer, gene) regulatory scores.
# Reference: Fulco 2019 Nat Genet; Nasser 2021 Nature.

set -euo pipefail

ATAC_BAM=${1:-atac.dedup.bam}
H3K27AC_BAM=${2:-h3k27ac.dedup.bam}
HIC_DIR=${3:-hic_data/}                  # Cooler or hic format directory
# shellcheck disable=SC2034  # ABC consumes --chrom_sizes, not a genome FASTA; kept as a labelled slot so downstream arg positions match the usage string
GENOME_FA=${4:-hg38.fa}
SIZES=${5:-hg38.chrom.sizes}
GENE_BED=${6:-refseq_protein_coding.bed}
ABC_REPO=${7:-/path/to/ABC-Enhancer-Gene-Prediction}    # broadinstitute/ABC-Enhancer-Gene-Prediction
CELL_TYPE=${8:-K562}
OUTDIR=${9:-abc_out}
EFFECTIVE_GENOME=${10:-2913022398}        # hg38 total effective genome (deepTools, read-length-independent)

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
    --H3K27ac $H3K27AC_BAM \
    --DHS $ATAC_BAM \
    --chrom_sizes $SIZES \
    --chrom_sizes_bed ${SIZES}.bed \
    --ubiquitously_expressed_genes ubiquitously_expressed.txt \
    --cellType $CELL_TYPE \
    --outdir $OUTDIR/neighborhoods/

# 4. ABC predictions: Activity * Contact -- generates ALL unthresholded links
python $ABC_REPO/workflow/scripts/predict.py \
    --enhancers $OUTDIR/neighborhoods/EnhancerList.txt \
    --genes $OUTDIR/neighborhoods/GeneList.txt \
    --hic_file $HIC_DIR \
    --hic_type avg \
    --hic_resolution 5000 \
    --hic_pseudocount_distance 5000 \
    `# --hic_type choices: hic | juicebox | bedpe | avg -- must match the Hi-C input format` \
    `# --chrom_sizes and --hic_pseudocount_distance are both required=True in predict.py` \
    --chrom_sizes $SIZES \
    --score_column ABC.Score \
    --cellType $CELL_TYPE \
    --outdir $OUTDIR/predictions/

# 5. Threshold at ABC.Score >= 0.02 by column header (the ABC Snakemake pipeline runs
# filter_predictions.py with its full --output_* argument set; this is the standalone equivalent).
echo "ABC predictions (unthresholded): $OUTDIR/predictions/EnhancerPredictionsAllPutative.tsv.gz"
echo "Above threshold (ABC.Score >= 0.02):"
zcat $OUTDIR/predictions/EnhancerPredictionsAllPutative.tsv.gz | \
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="ABC.Score")c=i; next} $c>=0.02' | wc -l

# Optional: cross-validation
# Compare against published CRISPRi-FlowFISH catalog (Fulco 2019 K562)
# wget https://www.engreitzlab.org/crispri-flowfish/K562_validated_pairs.txt
# Compute sensitivity / specificity at threshold 0.02
