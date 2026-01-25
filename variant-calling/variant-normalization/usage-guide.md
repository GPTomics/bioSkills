# Variant Normalization Usage Guide

## Overview

This guide covers normalizing VCF files for consistent variant representation.

## Prerequisites

- bcftools installed (`conda install -c bioconda bcftools`)
- Reference FASTA file (same one used for variant calling)
- Reference must be indexed (`samtools faidx reference.fa`)

## Quick Start

Tell your AI agent what you want to do:
- "Normalize my VCF by left-aligning indels and splitting multiallelic sites"
- "Prepare my VCF for annotation by normalizing variants against the reference"
- "Compare variants from two different callers after normalizing both VCFs"
- "Split multiallelic sites but keep SNPs and indels separate"

## Why Normalization Matters

### The Indel Representation Problem

The same indel can be written multiple ways:

```
# All represent the same 3bp deletion
Position 100:  ATCGA → A--GA  (ATG deleted)
Position 101:  TCGA  → T--A   (CGA deleted, shifted right)
Position 102:  CGA   → -GA    (C deleted, different start)
```

Different callers may choose different representations, making comparison impossible without normalization.

### Multiallelic Sites

A position with multiple alternate alleles:

```
chr1  100  .  A  G,T,C  30  PASS
```

Many tools require biallelic sites (one ALT per line). Splitting is part of normalization.

## Basic Normalization

### Left-Align Indels

```bash
bcftools norm -f reference.fa input.vcf.gz -Oz -o normalized.vcf.gz
bcftools index normalized.vcf.gz
```

This shifts indels to their leftmost position.

### Split Multiallelic Sites

```bash
bcftools norm -m-any input.vcf.gz -Oz -o split.vcf.gz
bcftools index split.vcf.gz
```

### Combined (Recommended)

```bash
bcftools norm -f reference.fa -m-any input.vcf.gz -Oz -o normalized.vcf.gz
bcftools index normalized.vcf.gz
```

## Understanding Multiallelic Splitting

### Before Splitting

```
chr1  100  .  A  G,T  30  PASS  GT  1/2
```

This represents a heterozygous site with two different alternate alleles.

### After Splitting (`-m-any`)

```
chr1  100  .  A  G  30  PASS  GT  1/0
chr1  100  .  A  T  30  PASS  GT  0/1
```

Two separate records, each with one alternate allele.

### Selective Splitting

```bash
# Split only SNPs (keep multiallelic indels)
bcftools norm -m-snps input.vcf.gz -Oz -o split_snps.vcf.gz

# Split only indels (keep multiallelic SNPs)
bcftools norm -m-indels input.vcf.gz -Oz -o split_indels.vcf.gz

# Split SNPs and indels separately (don't mix)
bcftools norm -m-both input.vcf.gz -Oz -o split_both.vcf.gz
```

## Joining Biallelic Sites

The reverse operation - combine biallelic sites at the same position:

```bash
bcftools norm -m+any input.vcf.gz -Oz -o joined.vcf.gz
```

### Before Joining

```
chr1  100  .  A  G  30  PASS
chr1  100  .  A  T  30  PASS
```

### After Joining

```
chr1  100  .  A  G,T  30  PASS
```

## Reference Allele Checking

### Check for Mismatches

```bash
bcftools norm -f reference.fa -c w input.vcf.gz > /dev/null
```

Options for `-c`:
- `w` - Warn on mismatch (default)
- `e` - Error on mismatch
- `x` - Exclude variants with mismatch
- `s` - Set REF to match reference

### Fix Mismatches

```bash
bcftools norm -f reference.fa -c s input.vcf.gz -Oz -o fixed.vcf.gz
```

### Remove Mismatches

```bash
bcftools norm -f reference.fa -c x input.vcf.gz -Oz -o clean.vcf.gz
```

## Removing Duplicates

After splitting or merging, duplicates may appear:

```bash
# Remove exact duplicates
bcftools norm -d exact input.vcf.gz -Oz -o deduped.vcf.gz

# Remove all duplicates at same position
bcftools norm -d all input.vcf.gz -Oz -o deduped.vcf.gz
```

Duplicate options:
- `exact` - Same CHROM, POS, REF, ALT
- `snps` - Duplicate SNPs only
- `indels` - Duplicate indels only
- `both` - Duplicate SNPs and indels
- `all` - Any duplicate at position
- `none` - Keep all (default)

## Atomizing Complex Variants

MNPs (multi-nucleotide polymorphisms) can be split into individual SNPs:

### Before Atomizing

```
chr1  100  .  ATG  GCA  30  PASS
```

### After Atomizing

```bash
bcftools norm --atomize input.vcf.gz -Oz -o atomized.vcf.gz
```

```
chr1  100  .  A  G  30  PASS
chr1  101  .  T  C  30  PASS
chr1  102  .  G  A  30  PASS
```

### Combined with Left-Alignment

```bash
bcftools norm -f reference.fa --atomize input.vcf.gz -Oz -o atomized.vcf.gz
```

## Complete Normalization Pipeline

### Standard Workflow

```bash
# 1. Normalize (left-align + split)
bcftools norm -f reference.fa -m-any input.vcf.gz -Oz -o step1.vcf.gz

# 2. Remove duplicates
bcftools norm -d exact step1.vcf.gz -Oz -o normalized.vcf.gz

# 3. Index
bcftools index normalized.vcf.gz

# Cleanup
rm step1.vcf.gz
```

### One-Step Pipeline

```bash
bcftools norm -f reference.fa -m-any -d exact input.vcf.gz -Oz -o normalized.vcf.gz
bcftools index normalized.vcf.gz
```

## Workflow: Compare Variant Callers

Before comparing VCFs from different callers, normalize both:

```bash
# Normalize GATK output
bcftools norm -f reference.fa -m-any gatk.vcf.gz -Oz -o gatk.norm.vcf.gz
bcftools index gatk.norm.vcf.gz

# Normalize bcftools output
bcftools norm -f reference.fa -m-any bcftools.vcf.gz -Oz -o bcftools.norm.vcf.gz
bcftools index bcftools.norm.vcf.gz

# Now compare
bcftools isec -p comparison gatk.norm.vcf.gz bcftools.norm.vcf.gz

# Count shared variants
echo "Shared variants: $(bcftools view -H comparison/0002.vcf | wc -l)"
```

## Workflow: Prepare for Annotation

Many annotation databases expect normalized variants:

```bash
# Full normalization
bcftools norm -f reference.fa -m-any -d exact variants.vcf.gz -Oz -o for_annotation.vcf.gz
bcftools index for_annotation.vcf.gz

# Now annotate with SnpEff, VEP, or database lookup
```

## Statistics

### Count Multiallelic Sites Before Normalization

```bash
# Count sites with multiple ALT alleles
bcftools view -H input.vcf.gz | awk -F'\t' '$5 ~ /,/' | wc -l
```

### Compare Before/After

```bash
echo "Before: $(bcftools view -H input.vcf.gz | wc -l) variants"

bcftools norm -f reference.fa -m-any input.vcf.gz -Oz -o normalized.vcf.gz

echo "After: $(bcftools view -H normalized.vcf.gz | wc -l) variants"
```

The count increases because multiallelic sites become multiple records.

## Troubleshooting

### "reference not found"

Need to provide reference FASTA:

```bash
bcftools norm -f reference.fa input.vcf.gz
```

### "REF does not match"

Wrong reference genome. Verify chromosomes match:

```bash
# Check VCF chromosomes
bcftools view -H input.vcf.gz | cut -f1 | sort -u

# Check reference chromosomes
grep "^>" reference.fa | head
```

### "not sorted"

Input must be sorted:

```bash
bcftools sort input.vcf.gz -Oz -o sorted.vcf.gz
bcftools norm -f reference.fa sorted.vcf.gz -Oz -o normalized.vcf.gz
```

### No Change After Normalization

Variants may already be normalized. Check:

```bash
bcftools norm -f reference.fa input.vcf.gz 2>&1 | grep "records modified"
```

## Example Prompts

> "Normalize my VCF by left-aligning indels and splitting multiallelic sites"

> "Prepare my VCF for annotation by normalizing variants against the reference"

> "Compare variants from two different callers after normalizing both VCFs"

> "Split multiallelic sites but keep SNPs and indels separate"

## See Also

- [bcftools norm documentation](http://www.htslib.org/doc/bcftools.html#norm)
- [Variant normalization paper](https://doi.org/10.1093/bioinformatics/btv112)
- **vcf-manipulation** - Compare normalized VCFs
- **variant-annotation** - Annotate after normalization
