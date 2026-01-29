#!/bin/bash
# IGV batch script for generating genome browser snapshots
#
# Usage: igv -b igv_batch.txt
# Requires: IGV installed and in PATH (or use full path to igv command)

# --- Create batch script ---
# Customize genome, files, and regions for your data

cat > igv_batch.txt << 'EOF'
# IGV Batch Commands
# Start new session
new

# Load reference genome
# Options: hg38, hg19, mm10, or path to .genome file
genome hg38

# Set snapshot output directory
snapshotDirectory ./igv_snapshots

# Load tracks
# Supported: BAM, bigWig, BED, VCF, GTF, etc.
load sample1.bam
load sample2.bam
load peaks.bed

# Configure display
# collapse/expand/squish for read display
collapse
maxPanelHeight 500

# Navigate to region and take snapshot
goto chr1:1000000-2000000
snapshot chr1_region1.png

# Another region
goto chr17:7565097-7590856
# This is the TP53 locus
snapshot tp53_locus.png

# Gene name navigation
goto MYC
snapshot myc_locus.png

# Sort reads by different criteria
# Options: base, strand, quality, sample, readGroup
sort base
snapshot myc_sorted.png

# Exit IGV
exit
EOF

echo "IGV batch script created: igv_batch.txt"
echo ""
echo "To run:"
echo "  1. Ensure your files (BAM, BED, etc.) are in the current directory"
echo "  2. Run: igv -b igv_batch.txt"
echo ""
echo "Common modifications:"
echo "  - Change 'genome hg38' to your reference"
echo "  - Add 'load' commands for your tracks"
echo "  - Add 'goto' and 'snapshot' for your regions"
echo ""
echo "To batch process regions from BED file:"

cat > process_regions.sh << 'SCRIPT'
#!/bin/bash
# Generate IGV batch from BED file of regions

BED_FILE="$1"
OUTPUT_DIR="${2:-igv_snapshots}"
GENOME="${3:-hg38}"

if [ -z "$BED_FILE" ]; then
    echo "Usage: $0 regions.bed [output_dir] [genome]"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Generate batch file
cat > igv_regions_batch.txt << EOF
new
genome $GENOME
snapshotDirectory $OUTPUT_DIR
load sample.bam
collapse
maxPanelHeight 500
EOF

while read -r line; do
    # Skip header/comment lines
    [[ "$line" =~ ^# ]] && continue
    [[ -z "$line" ]] && continue

    # Parse BED fields
    chr=$(echo "$line" | cut -f1)
    start=$(echo "$line" | cut -f2)
    end=$(echo "$line" | cut -f3)
    name=$(echo "$line" | cut -f4)

    # Use region as name if no name column
    [ -z "$name" ] && name="${chr}_${start}_${end}"

    echo "goto ${chr}:${start}-${end}" >> igv_regions_batch.txt
    echo "snapshot ${name}.png" >> igv_regions_batch.txt

done < "$BED_FILE"

echo "exit" >> igv_regions_batch.txt

echo "Generated: igv_regions_batch.txt"
echo "Run: igv -b igv_regions_batch.txt"
SCRIPT

chmod +x process_regions.sh
echo "Created: process_regions.sh"
echo "Usage: ./process_regions.sh regions.bed output_dir hg38"
