#!/bin/bash
# Generate coverage bedGraph from BAM files

# Usage: ./bedgraph_from_bam.sh alignments.bam [output_prefix]
BAM="${1:-alignments.bam}"
PREFIX="${2:-coverage}"

if [[ ! -f "$BAM" ]]; then
    echo "Usage: $0 <bam_file> [output_prefix]"
    echo "BAM file not found: $BAM"
    exit 1
fi

echo "Processing: $BAM"
echo "Output prefix: $PREFIX"

# Basic bedGraph (per-base coverage, non-zero only)
echo -e "\n=== Generating basic bedGraph ==="
bedtools genomecov -ibam "$BAM" -bg > "${PREFIX}.bedGraph"
echo "Created: ${PREFIX}.bedGraph"
echo "Lines: $(wc -l < "${PREFIX}.bedGraph")"

# bedGraph with zeros (all positions)
echo -e "\n=== Generating bedGraph with zeros ==="
bedtools genomecov -ibam "$BAM" -bga > "${PREFIX}.with_zeros.bedGraph"
echo "Created: ${PREFIX}.with_zeros.bedGraph"
echo "Lines: $(wc -l < "${PREFIX}.with_zeros.bedGraph")"

# Normalized bedGraph (CPM - counts per million)
echo -e "\n=== Generating CPM-normalized bedGraph ==="
TOTAL=$(samtools view -c "$BAM")
echo "Total reads: $TOTAL"
SCALE=$(echo "scale=10; 1000000/$TOTAL" | bc)
echo "Scale factor: $SCALE"
bedtools genomecov -ibam "$BAM" -bg -scale "$SCALE" > "${PREFIX}.cpm.bedGraph"
echo "Created: ${PREFIX}.cpm.bedGraph"

# Coverage histogram (genome-wide)
echo -e "\n=== Generating coverage histogram ==="
bedtools genomecov -ibam "$BAM" > "${PREFIX}.histogram.txt"
echo "Created: ${PREFIX}.histogram.txt"

# Summary statistics
echo -e "\n=== Coverage Statistics ==="
awk '$1 == "genome" {
    depth[NR] = $2;
    frac[NR] = $5;
    sum += $2 * $5;
}
END {
    print "Mean depth: " sum "x";
    # Find median (cumulative fraction >= 0.5)
    cumsum = 0;
    for (i in depth) {
        cumsum += frac[i];
        if (cumsum >= 0.5) {
            print "Median depth: " depth[i] "x";
            break;
        }
    }
}' "${PREFIX}.histogram.txt"

# Fraction of genome covered
echo -e "\n=== Genome Coverage ==="
awk '$1 == "genome" && $2 > 0 {sum += $5}
     END {print "Fraction covered (>0x): " sum * 100 "%"}' "${PREFIX}.histogram.txt"
awk '$1 == "genome" && $2 >= 10 {sum += $5}
     END {print "Fraction covered (>=10x): " sum * 100 "%"}' "${PREFIX}.histogram.txt"
awk '$1 == "genome" && $2 >= 30 {sum += $5}
     END {print "Fraction covered (>=30x): " sum * 100 "%"}' "${PREFIX}.histogram.txt"

echo -e "\n=== Done ==="
echo "Files created:"
ls -la ${PREFIX}*
