#!/bin/bash
# Reference: modkit 0.3+, samtools 1.19+, htslib 1.19+ | Verify API if version differs
# Call 5mC from a Nanopore modBAM. The tags MUST be present and have survived alignment;
# modkit auto-thresholds at the 10th percentile of the ML distribution (NOT 0.5), and there
# is no --min-coverage flag (filter on Nvalid_cov, column 10, afterward).
set -euo pipefail

BAM=${1:?Usage: $0 <aligned_mod.bam> <reference.fa> [output_prefix] [min_cov]}
REFERENCE=${2:?reference required}
PREFIX=${3:-methylation}
MIN_COV=${4:-10}          # Nvalid_cov floor for confident single-site calls

# Methylation is a basecalling decision: no MM tags means re-basecall, not re-analyze.
# Capture a count (grep -cm1); `head | grep -q` would return 141 under `set -o pipefail`
# (samtools killed by SIGPIPE) and invert the check.
if [ "$(samtools view "$BAM" | grep -cm1 'MM:Z' || true)" -eq 0 ]; then
    echo 'ERROR: no MM/ML tags - this BAM has no methylation. Re-basecall from POD5 with a mods model.'
    exit 1
fi

# Pile up 5mC at CpG sites, strands combined. bedMethyl has 18 columns; col 10 = Nvalid_cov
# (the denominator), col 11 = percent_modified (0-100).
modkit pileup "$BAM" "${PREFIX}.bed" --ref "$REFERENCE" --cpg --combine-strands --threads 8
bgzip -f "${PREFIX}.bed" && tabix -p bed "${PREFIX}.bed.gz"

# Coverage filtering is post-hoc on Nvalid_cov, not a pileup flag.
zcat "${PREFIX}.bed.gz" | awk -v c="$MIN_COV" '$10 >= c' > "${PREFIX}.cov${MIN_COV}.bed"

echo "Sites with >= ${MIN_COV}x: $(wc -l < "${PREFIX}.cov${MIN_COV}.bed")"
echo "Mean percent modified (>= ${MIN_COV}x):"
awk '{s+=$11; n++} END {if(n) print s/n"%"}' "${PREFIX}.cov${MIN_COV}.bed"

# To compare against WGBS, combine 5mC+5hmC (bisulfite conflates them):
#   modkit pileup "$BAM" combined.bed --ref "$REFERENCE" --preset traditional
# For allele-specific methylation, phase+haplotag first, then:
#   modkit pileup haplotagged.bam asm/ --ref "$REFERENCE" --cpg --combine-strands --partition-tag HP
echo "bedMethyl: ${PREFIX}.bed.gz  (hand Nmod/Nvalid_cov to methylation-analysis/dmr-detection)"
