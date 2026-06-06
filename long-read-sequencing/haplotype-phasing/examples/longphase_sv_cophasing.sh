#!/bin/bash
# Reference: longphase 1.7+, samtools 1.19+ | Verify API if version differs
# Whole-genome ONT phasing with LongPhase, co-phasing SNPs + indels + SVs (+ optional 5mC)
# into long blocks (~25 Mbp N50), then haplotagging. ~10x faster than WhatsHap on WGS.
set -euo pipefail

SNPS=${1:?Usage: $0 <snps.vcf> <svs.vcf> <aligned.bam> <ref.fa> [platform] [out_prefix]}
SVS=${2:?SV VCF required (e.g. from Sniffles2)}
BAM=${3:?aligned BAM required}
REFERENCE=${4:?reference required}
PLATFORM=${5:-ont}        # ont | pb (bare flags --ont / --pb, NOT --platform)
PREFIX=${6:-longphase}
THREADS=16

PLAT_FLAG="--ont"; [ "$PLATFORM" = pb ] && PLAT_FLAG="--pb"

# Co-phasing a phased SV bridges het-sparse gaps that SNPs alone cannot span -> longer blocks.
longphase phase -s "$SNPS" --indels --sv-file "$SVS" -b "$BAM" -r "$REFERENCE" \
    -o "$PREFIX" -t "$THREADS" $PLAT_FLAG
#   emits ${PREFIX}.vcf and ${PREFIX}_SV.vcf

# Haplotag with the phased small + SV VCFs; writes the same HP/PS read tags as WhatsHap.
longphase haplotag -s "${PREFIX}.vcf" --sv-file "${PREFIX}_SV.vcf" \
    -b "$BAM" -r "$REFERENCE" -o "${PREFIX}_haplotagged" -t "$THREADS"
samtools index "${PREFIX}_haplotagged.bam"

# Capture a count (grep -cm1); a `head | grep -q` pipeline would return 141 under
# `set -o pipefail` (samtools killed by SIGPIPE) and take the wrong branch.
if [ "$(samtools view "${PREFIX}_haplotagged.bam" | grep -cm1 'HP:i:' || true)" -ge 1 ]; then
    echo "HP tags present in ${PREFIX}_haplotagged.bam"
else
    echo 'WARNING: no HP tags - check phasing output.'
fi
