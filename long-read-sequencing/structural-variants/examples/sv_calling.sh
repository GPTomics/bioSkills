#!/bin/bash
# Reference: Sniffles 2.2+, cuteSV 2.1+, bcftools 1.19+ | Verify API if version differs
# Germline SV calling with Sniffles2 and cuteSV. The tandem-repeat BED is the single
# biggest FP-reduction lever in repeats; --reference makes insertions carry sequence.
set -euo pipefail

BAM=${1:?Usage: $0 <aln.bam> <ref.fa> <tandem_repeats.bed> [platform] [out_dir]}
REFERENCE=${2:?reference required}
TR_BED=${3:?tandem-repeat BED required (ships with Sniffles annotations/)}
PLATFORM=${4:-ont}        # ont | hifi | clr (sets the cuteSV parameter set)
OUTPUT_DIR=${5:-sv_calls}
mkdir -p "$OUTPUT_DIR"

# Sniffles2: support is coverage-derived (--minsupport auto), minsvlen default 35.
echo 'Calling SVs with Sniffles2...'
sniffles --input "$BAM" --vcf "${OUTPUT_DIR}/sniffles.vcf" \
    --reference "$REFERENCE" --tandem-repeats "$TR_BED"

# cuteSV parameters are NOT one-size-fits-all; pick the set by platform error rate.
case "$PLATFORM" in
    ont)  IB=100;  IR=0.3; DB=100;  DR=0.3 ;;
    hifi) IB=1000; IR=0.9; DB=1000; DR=0.5 ;;
    clr)  IB=100;  IR=0.3; DB=200;  DR=0.5 ;;
    *) echo "Unknown platform: $PLATFORM"; exit 1 ;;
esac

echo "Calling SVs with cuteSV ($PLATFORM parameters)..."
mkdir -p "${OUTPUT_DIR}/cutesv_work"
cuteSV "$BAM" "$REFERENCE" "${OUTPUT_DIR}/cutesv.vcf" "${OUTPUT_DIR}/cutesv_work" \
    --genotype \
    --max_cluster_bias_INS "$IB" --diff_ratio_merging_INS "$IR" \
    --max_cluster_bias_DEL "$DB" --diff_ratio_merging_DEL "$DR"

for vcf in sniffles cutesv; do
    echo "=== $vcf ==="
    bcftools stats "${OUTPUT_DIR}/${vcf}.vcf" | grep '^SN'
done
echo "Done. Benchmark with: truvari bench --base TRUTH --comp ${OUTPUT_DIR}/sniffles.vcf ... && truvari refine"
