#!/bin/bash
# Reference: Clair3 2.0+, bcftools 1.19+, samtools 1.19+ | Verify API if version differs
# Germline small-variant calling with Clair3. The MODEL must match the basecaller;
# there is no auto-detection and a mismatch silently degrades calls.
set -euo pipefail

BAM=${1:?Usage: $0 <aln.bam> <ref.fa> <model_path> [platform] [threads] [sample]}
REFERENCE=${2:?reference required}
MODEL_PATH=${3:?model_path required - a SPECIFIC model folder, e.g. .../r1041_e82_400bps_sup_v500}
PLATFORM=${4:-ont}        # ont | hifi | ilmn
THREADS=${5:-16}
SAMPLE=${6:-sample}

samtools quickcheck "$BAM" || { echo 'BAM invalid'; exit 1; }

# --include_all_ctgs is mandatory for non-human / draft references, or output is empty.
# --enable_phasing phases the final VCF (WhatsHap); --longphase_for_phasing swaps only the
# internal phaser to LongPhase.
run_clair3.sh \
    --bam_fn="$BAM" \
    --ref_fn="$REFERENCE" \
    --threads="$THREADS" \
    --platform="$PLATFORM" \
    --model_path="$MODEL_PATH" \
    --include_all_ctgs \
    --output="${SAMPLE}_clair3"

VCF="${SAMPLE}_clair3/merge_output.vcf.gz"
bcftools index -t "$VCF"

echo 'Variant counts:'
bcftools stats "$VCF" | grep '^SN'

# Benchmark with stratification (not just global F1) to expose ONT indel errors in
# homopolymer/STR strata - the residual error mode that whole-genome F1 hides:
#   hap.py giab_truth.vcf.gz "$VCF" -f giab_confident.bed -r "$REFERENCE" \
#       --engine=vcfeval --stratification giab_strat.tsv -o bench/${SAMPLE}
echo "Final VCF: $VCF"
