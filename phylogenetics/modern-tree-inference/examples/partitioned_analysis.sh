#!/bin/bash
# Reference: IQ-TREE 2.2+/2.3+ | Verify API if version differs
# Partitioned multi-gene analysis with per-partition models, merging, and
# concordance factors. The branch-length linkage flag is what people get wrong:
#   -p edge-linked PROPORTIONAL (default/recommended), -q edge-equal, -Q edge-unlinked.
# NOT spot-runnable offline: needs the iqtree2 binary and real alignments.
set -euo pipefail

CONCAT="${1:-concatenated.fasta}"
OUT="part_out"
mkdir -p "$OUT"

# Partition file: gene boundaries (models omitted so ModelFinder chooses + merges)
cat > "$OUT/partitions.nex" << 'EOF'
#nexus
begin sets;
    charset COI  = 1-657;
    charset CYTB = 658-1140;
    charset 16S  = 1141-1650;
    charset 28S  = 1651-2100;
end;
EOF

# -p           edge-linked proportional BLs (one rate multiplier per partition)
# -m MFP+MERGE per-partition models + greedy BIC merge of partitions that match
# -rcluster 10 relaxed clustering: only test top 10% most-similar pairs (tractable)
iqtree2 -s "$CONCAT" -p "$OUT/partitions.nex" -m MFP+MERGE -rcluster 10 \
        -B 1000 -bnni -alrt 1000 -T AUTO --seed 12345 --prefix "$OUT/part"

# Concordance factors expose nodes the concatenated bootstrap over-reads.
# -S builds one separate gene tree per partition (no concatenation); --gcf/--scfl
# then score the fixed concat tree against those per-locus trees and the sites.
iqtree2 -s "$CONCAT" -S "$OUT/partitions.nex" -m MFP -T AUTO --prefix "$OUT/loci"
iqtree2 -te "$OUT/part.treefile" -s "$CONCAT" \
        --gcf "$OUT/loci.treefile" --scfl 100 -T 4 --prefix "$OUT/concord"

echo "Partitioned tree: $OUT/part.treefile"
echo "Concordance:      $OUT/concord.cf.tree (UFBoot 100 + gCF ~33 => unresolved/ILS)"
grep -A 20 "Best-fit model" "$OUT/part.iqtree" || true
