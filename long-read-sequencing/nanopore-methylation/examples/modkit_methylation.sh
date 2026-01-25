#!/bin/bash
# Nanopore Methylation Calling with modkit
# Demonstrates methylation analysis from ONT data

# =============================================================================
# Prerequisites
# =============================================================================

# Ensure BAM has methylation tags from dorado basecalling
# BAM should contain MM and ML tags

# Check if BAM has modification tags
echo "=== Checking BAM for modification tags ==="
samtools view -h input.bam | head -100 | grep -E "^@|MM:Z:|ML:B:"

# =============================================================================
# Basic Methylation Pileup
# =============================================================================

# Extract 5mC methylation at CpG sites
modkit pileup input.bam methylation.bed \
    --ref reference.fa \
    --cpg \
    --combine-strands \
    --threads 8

# Output: bedMethyl format
# chr  start  end  name  score  strand  thickStart  thickEnd  color  coverage  percent_modified

# =============================================================================
# Filtering Options
# =============================================================================

# Minimum coverage filter
modkit pileup input.bam methylation_cov10.bed \
    --ref reference.fa \
    --cpg \
    --combine-strands \
    --filter-threshold 0.5 \
    --min-coverage 10

# Specific regions only (faster)
modkit pileup input.bam promoters_meth.bed \
    --ref reference.fa \
    --cpg \
    --combine-strands \
    --include-bed promoter_regions.bed

# =============================================================================
# Sample Summary Statistics
# =============================================================================

echo "=== Methylation Summary ==="
modkit summary input.bam

# Output includes:
# - Total reads with modifications
# - Modification types (5mC, 6mA, etc.)
# - Fraction of bases modified

# =============================================================================
# Multiple Modification Types
# =============================================================================

# 5mC and 5hmC separately (if model supports)
modkit pileup input.bam meth_5mC.bed \
    --ref reference.fa \
    --mod-thresholds m:0.5 \
    --cpg

modkit pileup input.bam meth_5hmC.bed \
    --ref reference.fa \
    --mod-thresholds h:0.5 \
    --cpg

# =============================================================================
# Extract per-read methylation
# =============================================================================

# Get per-read modification calls (for visualization)
modkit extract input.bam per_read_meth.tsv \
    --ref reference.fa

# =============================================================================
# Compare Two Samples
# =============================================================================

# Generate bedMethyl for each sample
for sample in tumor normal; do
    modkit pileup "${sample}.bam" "${sample}_meth.bed" \
        --ref reference.fa \
        --cpg \
        --combine-strands \
        --min-coverage 10
done

# Find differentially methylated regions (requires R)
# See methylation-analysis/dmr-detection for downstream analysis

# =============================================================================
# Quality Metrics
# =============================================================================

echo "=== Quality Checks ==="

# Count CpG sites with sufficient coverage
echo "CpG sites with >=10x coverage:"
awk '$10 >= 10' methylation.bed | wc -l

# Distribution of methylation levels
echo "Methylation distribution:"
awk '$10 >= 10 {print int($11/10)*10}' methylation.bed | sort -n | uniq -c

# Average methylation
echo "Mean methylation (sites with >=10x):"
awk '$10 >= 10 {sum += $11; n++} END {print sum/n "%"}' methylation.bed

# Recommended thresholds:
# - Minimum coverage: 10x for reliable single-site calls
# - Probability threshold: 0.5 (default), increase for higher confidence
# - For population-level analysis: can use lower coverage
