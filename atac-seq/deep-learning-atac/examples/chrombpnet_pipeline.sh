#!/bin/bash
# Reference: chrombpnet 0.1.7+, tangermeme 0.1+ | Verify API if version differs
# chromBPNet end-to-end: bias model -> accessibility model -> per-base profile + variant effect prediction.
# GPU required (~24 h on A100 per cell type).

set -euo pipefail

BAM=${1:-atac.dedup.bam}
PEAKS=${2:-peaks.narrowPeak}                       # Output of MACS3
NONPEAKS=${3:-nonpeaks.bed}                        # Background regions for bias model
GENOME=${4:-hg38.fa}
SIZES=${5:-hg38.chrom.sizes}
OUTDIR=${6:-chrombpnet_out}
DATA_TYPE=${7:-ATAC}                               # ATAC | DNASE
BIAS_SCALING=${8:-0.5}

mkdir -p $OUTDIR/{bias,model,splits,variants}

# 1. Generate train/val/test chromosome splits (output is JSON with chrom assignments)
# Train chroms = whatever is not in -tcr / -vcr (auto-inferred); `-tecr` flag does NOT exist.
chrombpnet prep splits \
    -c $SIZES \
    -tcr chr1 chr3 chr6 \
    -vcr chr8 chr20 \
    -op $OUTDIR/splits/fold_0

# 2. Train bias model from non-peak (background) regions
chrombpnet bias pipeline \
    -ibam $BAM \
    -d $DATA_TYPE \
    -g $GENOME -c $SIZES \
    -p $PEAKS -n $NONPEAKS \
    -fl $OUTDIR/splits/fold_0.json \
    -b $BIAS_SCALING \
    -o $OUTDIR/bias/

# 3. Train accessibility model with bias correction
chrombpnet pipeline \
    -ibam $BAM \
    -d $DATA_TYPE \
    -g $GENOME -c $SIZES \
    -p $PEAKS -n $NONPEAKS \
    -fl $OUTDIR/splits/fold_0.json \
    -b $OUTDIR/bias/bias.h5 \
    -o $OUTDIR/model/

echo "Trained model: $OUTDIR/model/chrombpnet_no_bias.h5"
echo "Bias model: $OUTDIR/bias/bias.h5"
echo "Bias-corrected bigWig: $OUTDIR/model/chrombpnet_no_bias_predictions/profile.bw"

# 4. Predict variant effects (optional; requires variants.tsv with chrom, pos, ref, alt)
# Uses kundajelab/variant-scorer (separate repo); `chrombpnet snp_score` is commented out in current chrombpnet
if [ -f "${9:-}" ]; then
    VARIANTS=$9
    VARIANT_SCORER=${VARIANT_SCORER:-/path/to/variant-scorer}    # Clone from kundajelab/variant-scorer
    python $VARIANT_SCORER/src/variant_scoring.py \
        --model $OUTDIR/model/chrombpnet_no_bias.h5 \
        --list $VARIANTS \
        --genome $GENOME \
        --chrom_sizes $SIZES \
        --out_prefix $OUTDIR/variants/predictions
    echo "Variant predictions: $OUTDIR/variants/predictions.variant_scores.tsv"
    echo "  Reports log2FC magnitudes; abs(log2FC) > 1 indicates strong-effect SNP"
fi

# 5. Validation: aggregate corrected profile at CTCF (sanity check; expect clean V)
echo "Next steps:"
echo "  - Run TOBIAS PlotAggregate on \$OUTDIR/model/.../profile.bw at CTCF motifs"
echo "  - For motif discovery: compute DeepLIFT contributions and run TF-MoDISco-lite"
echo "  - For variant scoring: tangermeme.variant_effect.substitution_effect on bulk SNP lists"
