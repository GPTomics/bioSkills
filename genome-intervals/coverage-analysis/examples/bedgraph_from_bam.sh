#!/bin/bash
# Reference: mosdepth 0.3+, bedtools 2.31+, samtools 1.19+ | Verify API if version differs
# Coverage analysis: breadth-first (mosdepth), per-contig glance (samtools), and a bedGraph track (bedtools).
# Mark duplicates BEFORE running this -- depth on an un-deduped BAM is a vanity number.

# Usage: ./bedgraph_from_bam.sh alignments.bam [output_prefix]
BAM="${1:-alignments.bam}"
PREFIX="${2:-coverage}"

MIN_MAPQ=20        # drop multimappers/repeat reads; repeat coverage intentionally collapses (no neutral choice)
WINDOW=500         # mosdepth window size in bp; coarser = smaller output, finer = more resolution
CALLABLE_MIN=4     # min depth to call a base CALLABLE; tune to the variant caller
EXCESSIVE=150      # HIGH/excessive-depth ceiling; flags rDNA/artifact pileups

if [[ ! -f "$BAM" ]]; then
    echo "Usage: $0 <bam_file> [output_prefix]"
    echo "BAM file not found: $BAM"
    exit 1
fi

echo "=== mosdepth: median + breadth curve (the modern default) ==="
mosdepth --by "$WINDOW" -Q "$MIN_MAPQ" "$PREFIX" "$BAM"
echo "Wrote ${PREFIX}.mosdepth.summary.txt and ${PREFIX}.mosdepth.global.dist.txt"
echo "Per-chrom + total mean depth:"
cat "${PREFIX}.mosdepth.summary.txt"

echo -e "\n=== Breadth from the cumulative distribution (proportion >= depth) ==="
# global.dist.txt columns: chrom, depth, proportion_of_bases_at_least_this_depth ('total' = whole genome)
for D in 1 10 20 30; do
    awk -v d="$D" '$1=="total" && $2==d {printf "  breadth >= %2dx: %.1f%%\n", d, $3*100}' "${PREFIX}.mosdepth.global.dist.txt"
done
echo "  (median = the depth where proportion crosses 0.5)"
awk '$1=="total" && $3>=0.5 {m=$2} END {print "  median depth: " m "x"}' "${PREFIX}.mosdepth.global.dist.txt"

echo -e "\n=== mosdepth --quantize: callable-region BED ==="
# bins: [0,1)=NO_COVERAGE, [1,CALLABLE_MIN)=LOW, [CALLABLE_MIN,EXCESSIVE)=CALLABLE, [EXCESSIVE,inf)=HIGH
mosdepth --quantize "0:1:${CALLABLE_MIN}:${EXCESSIVE}:" "${PREFIX}.callable" "$BAM"
echo "Wrote ${PREFIX}.callable.quantized.bed.gz"

echo -e "\n=== samtools coverage: per-contig depth + breadth glance ==="
# 'coverage' column is BREADTH (% bases >=1x); 'meandepth' is depth
samtools coverage "$BAM"

echo -e "\n=== bedtools genomecov: bedGraph track (zeros included) ==="
bedtools genomecov -ibam "$BAM" -bga > "${PREFIX}.bedGraph"   # -bga marks zero-coverage gaps; -bg omits them
echo "Wrote ${PREFIX}.bedGraph ($(wc -l < "${PREFIX}.bedGraph") lines)"
echo "  (for short-insert/amplicon data use -pc for fragment coverage; for RNA-seq add -split)"

echo -e "\n=== Done ==="
ls -la "${PREFIX}"*
