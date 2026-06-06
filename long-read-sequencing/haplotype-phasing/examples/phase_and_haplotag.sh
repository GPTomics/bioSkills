#!/bin/bash
# Reference: whatshap 2.3+, samtools 1.19+, htslib 1.19+ | Verify API if version differs
# Read-backed phasing with WhatsHap, then haplotagging. Phasing the VCF is NOT enough -
# downstream read-level tools need the BAM HP tag that only `haplotag` writes.
set -euo pipefail

VCF=${1:?Usage: $0 <het.vcf.gz> <aligned.bam> <ref.fa> [out_prefix]}
BAM=${2:?aligned BAM required}
REFERENCE=${3:?reference required}
PREFIX=${4:-phased}

# --reference enables realignment mode (rescues indel phasing on error-prone long reads);
# --indels also phases indels. Without --reference, phasing falls back to lower quality.
whatshap phase -o "${PREFIX}.vcf.gz" --reference "$REFERENCE" --indels "$VCF" "$BAM"
tabix -p vcf "${PREFIX}.vcf.gz"

# Separate, load-bearing step: write HP/PS read tags onto the BAM.
whatshap haplotag -o "${PREFIX}.haplotagged.bam" --reference "$REFERENCE" \
    --output-haplotag-list "${PREFIX}.htlist.tsv.gz" "${PREFIX}.vcf.gz" "$BAM"
samtools index "${PREFIX}.haplotagged.bam"

# Verify the tags survived - their absence is the silent failure mode.
# Capture a count (grep -cm1); `head | grep -q` would return 141 under `set -o pipefail`
# (samtools killed by SIGPIPE) and invert the check.
if [ "$(samtools view "${PREFIX}.haplotagged.bam" | grep -cm1 'HP:i:' || true)" -ge 1 ]; then
    echo 'HP tags present. Ready for modkit --partition-tag HP, Severus, IGV color-by-HP.'
else
    echo 'WARNING: no HP tags - check that phasing produced phased variants.'
fi

# Quality: ALWAYS report block N50 together with switch error (N50 alone is gameable).
whatshap stats --gtf "${PREFIX}.blocks.gtf" "${PREFIX}.vcf.gz"
echo "Compare to a trio/strand-seq truth with: whatshap compare --names truth,mine truth.vcf.gz ${PREFIX}.vcf.gz"
