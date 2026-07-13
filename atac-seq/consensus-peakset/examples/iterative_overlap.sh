#!/bin/bash
# Reference: bedtools 2.31+, samtools 1.19+, python 3.9+ | Verify API if version differs
# Corces 2018 iterative overlap removal: pool peaks -> re-center on summits -> 501bp fixed-width
# -> sort by signalValue -> greedy non-overlap -> blacklist filter -> SAF for featureCounts.

set -euo pipefail

PEAKS_GLOB=${1:-"peaks/per_sample/*.narrowPeak"}
GENOME_SIZES=${2:-hg38.chrom.sizes}
BLACKLIST=${3:-hg38-blacklist.v2.bed.gz}
HALF_WIDTH=${4:-250}                                  # Corces 2018 standard: 501 bp total
OUTDIR=${5:-consensus_out}

mkdir -p $OUTDIR

# 1. Pool all peaks; re-center on summit (col 2 + col 10); extend +/- HALF_WIDTH
echo "Pooling and re-centering..."
awk -v w=$HALF_WIDTH 'BEGIN{OFS="\t"}
    {summit = $2 + $10; print $1, summit - w, summit + w + 1, $4, $7, $6}' $PEAKS_GLOB | \
    awk '$2 >= 0' | \
    sort -k1,1 -k2,2n > $OUTDIR/pooled_recentered.bed

# Clamp to chromosome bounds
bedtools slop -i $OUTDIR/pooled_recentered.bed -g $GENOME_SIZES -b 0 > $OUTDIR/pooled_clamped.bed

# 2. Sort by signalValue (col 5; was narrowPeak col 7) descending
sort -k5,5gr $OUTDIR/pooled_clamped.bed > $OUTDIR/pooled_by_sig.bed

# 3. Iterative greedy non-overlap (Python because bedtools cluster doesn't preserve significance order)
OUTDIR="$OUTDIR" python3 - <<'PYEOF'
import os
indir = os.environ.get('OUTDIR', 'consensus_out')
kept = []
with open(f'{indir}/pooled_by_sig.bed') as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        chrom, start, end = cols[0], int(cols[1]), int(cols[2])
        if any(c == chrom and not (end <= s or start >= e) for c, s, e in kept):
            continue
        kept.append((chrom, start, end))
with open(f'{indir}/consensus_iterative.bed', 'w') as out:
    for c, s, e in sorted(kept):
        out.write(f'{c}\t{s}\t{e}\n')
print(f'Iterative overlap kept {len(kept)} peaks')
PYEOF

# 4. Blacklist filter
echo "Blacklist filtering..."
bedtools intersect -v -a $OUTDIR/consensus_iterative.bed -b $BLACKLIST | \
    sort -k1,1 -k2,2n > $OUTDIR/consensus_final.bed

echo "Final peakset: $(wc -l < $OUTDIR/consensus_final.bed) peaks at $((2 * HALF_WIDTH + 1)) bp fixed-width"

# 5. Convert to SAF format for featureCounts (BED start is 0-based half-open; SAF Start is 1-based inclusive, so +1)
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"}
     {print $1"_"$2"_"$3, $1, $2+1, $3, "+"}' $OUTDIR/consensus_final.bed \
    > $OUTDIR/consensus.saf

echo "Done. Use $OUTDIR/consensus_final.bed (BED) or $OUTDIR/consensus.saf (featureCounts SAF)."
echo "Next: featureCounts -F SAF -a $OUTDIR/consensus.saf -o counts.tsv -p --countReadPairs -T 8 *.bam"
