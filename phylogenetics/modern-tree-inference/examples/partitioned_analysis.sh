#!/bin/bash
# Partitioned analysis for multi-gene datasets with IQ-TREE2

# Input files
CONCAT_ALIGNMENT="concatenated.fasta"
PARTITION_FILE="partitions.nex"

# Create example partition file
# Defines gene boundaries and optionally models per partition
cat > "$PARTITION_FILE" << 'EOF'
#nexus
begin sets;
    charset COI = 1-657;
    charset CYTB = 658-1140;
    charset 16S = 1141-1650;
    charset 28S = 1651-2100;
end;
EOF

# Option 1: Merged partition analysis
# -p: Find best partitioning scheme (may merge similar partitions)
# Recommended for most analyses - balances fit and complexity
iqtree2 -s "$CONCAT_ALIGNMENT" -p "$PARTITION_FILE" -m MFP -B 1000 -T AUTO --prefix merged

# Option 2: Edge-linked proportional (separate models, proportional branch lengths)
# -q: Keep partitions separate but link branch lengths proportionally
# Good when genes evolve at different rates but tree topology is shared
iqtree2 -s "$CONCAT_ALIGNMENT" -q "$PARTITION_FILE" -m MFP -B 1000 -T AUTO --prefix linked

# Option 3: Edge-unlinked (fully independent branch lengths per partition)
# -Q: Fully independent branch lengths per partition
# Use when genes may have very different evolutionary histories
iqtree2 -s "$CONCAT_ALIGNMENT" -Q "$PARTITION_FILE" -m MFP -B 1000 -T AUTO --prefix unlinked

# View partition schemes and models selected
grep -A 20 "Best-fit model" merged.iqtree
