# Consensus Sequences Usage Guide

## Overview

This guide covers generating sample-specific sequences by applying variants to a reference.

## Prerequisites

- bcftools installed (`conda install -c bioconda bcftools`)
- Reference FASTA (same one used for variant calling)
- VCF file must be indexed (`bcftools index input.vcf.gz`)

## Quick Start

Tell your AI agent what you want to do:
- "Generate a consensus FASTA by applying my VCF variants to the reference"
- "Extract both haplotypes for a specific gene region from my phased VCF"
- "Create consensus sequences for all samples in my multi-sample VCF"
- "Mask low-coverage regions with N characters in my consensus sequence"

## Understanding Consensus Generation

`bcftools consensus` modifies a reference sequence by applying variants from a VCF file:

```
Reference:  ATCGATCGATCG
VCF:        pos 5: T→C
            pos 8: C→G
Consensus:  ATCGACCGAGCG
```

This creates a sequence that represents what the sample's genome looks like.

## Basic Usage

### Generate Consensus

```bash
bcftools consensus -f reference.fa input.vcf.gz > consensus.fa
bcftools consensus -f reference.fa -o consensus.fa input.vcf.gz
```

### For Specific Sample

Multi-sample VCFs require specifying which sample:

```bash
bcftools consensus -f reference.fa -s sample_name input.vcf.gz > sample.fa
```

### List Available Samples

```bash
bcftools query -l input.vcf.gz
```

## Handling Heterozygous Sites

At heterozygous sites, you need to choose which allele to apply:

### Haplotype Options

| Option | Description |
|--------|-------------|
| `-H 1` | First haplotype (allele 1) |
| `-H 2` | Second haplotype (allele 2) |
| `-H A` | Apply all ALT alleles |
| `-H R` | Keep REF at heterozygous sites |
| `-I` | Use IUPAC ambiguity codes |

### Extract First Haplotype

```bash
bcftools consensus -f reference.fa -H 1 input.vcf.gz > haplotype1.fa
```

### Extract Second Haplotype

```bash
bcftools consensus -f reference.fa -H 2 input.vcf.gz > haplotype2.fa
```

### Use IUPAC Ambiguity Codes

For heterozygous sites, use ambiguity codes:

```bash
bcftools consensus -f reference.fa -I input.vcf.gz > consensus_iupac.fa
```

IUPAC codes:
- R = A or G
- Y = C or T
- M = A or C
- K = G or T
- W = A or T
- S = C or G

## Masking Regions

### Mask Low-Coverage Regions

Create a BED file of regions to mask:

```bash
# Find regions with depth < 10
samtools depth -a input.bam | \
    awk '$3<10 {print $1"\t"$2-1"\t"$2}' | \
    bedtools merge > low_coverage.bed

# Apply mask
bcftools consensus -f reference.fa -m low_coverage.bed input.vcf.gz > consensus.fa
```

### Custom Mask Character

Default mask character is N. Change with:

```bash
bcftools consensus -f reference.fa -m mask.bed -M X input.vcf.gz > consensus.fa
```

### Mark Missing Genotypes

```bash
bcftools consensus -f reference.fa -M N input.vcf.gz > consensus.fa
```

## Region Selection

### Specific Chromosome

```bash
bcftools consensus -f reference.fa -r chr1 input.vcf.gz > chr1_consensus.fa
```

### Specific Region

```bash
bcftools consensus -f reference.fa -r chr1:1000000-2000000 input.vcf.gz > region.fa
```

### Gene Region

```bash
# Extract gene coordinates from annotation
bcftools consensus -f reference.fa -r chr1:1234567-1245678 input.vcf.gz > gene.fa
```

## Filtering Variants Before Consensus

### Apply Only PASS Variants

```bash
bcftools view -f PASS input.vcf.gz | \
    bcftools consensus -f reference.fa > consensus.fa
```

### Apply Only High-Quality Variants

```bash
bcftools filter -i 'QUAL>=30 && INFO/DP>=10' input.vcf.gz | \
    bcftools consensus -f reference.fa > consensus.fa
```

### Apply Only SNPs

```bash
bcftools view -v snps input.vcf.gz | \
    bcftools consensus -f reference.fa > consensus_snps.fa
```

### Exclude Indels

```bash
bcftools view -V indels input.vcf.gz | \
    bcftools consensus -f reference.fa > consensus.fa
```

## Chain Files for Coordinate Mapping

When indels change sequence length, coordinate mapping is needed:

### Generate Chain File

```bash
bcftools consensus -f reference.fa -c chain.txt input.vcf.gz > consensus.fa
```

### Use Chain File

Chain files can be used with liftOver tools to convert coordinates:
- BED coordinates
- GTF/GFF annotations
- Variant positions

## Working with Multiple Samples

### Generate Consensus for Each Sample

```bash
for sample in $(bcftools query -l cohort.vcf.gz); do
    echo "Processing $sample..."
    bcftools consensus -f reference.fa -s "$sample" cohort.vcf.gz > "${sample}.fa"
done
```

### Generate Both Haplotypes

```bash
sample="sample1"
bcftools consensus -f reference.fa -s "$sample" -H 1 input.vcf.gz > "${sample}_hap1.fa"
bcftools consensus -f reference.fa -s "$sample" -H 2 input.vcf.gz > "${sample}_hap2.fa"
```

## Sequence Naming

### Default Names

Output uses reference chromosome names.

### Add Prefix

```bash
bcftools consensus -f reference.fa -p "sample1_" input.vcf.gz > sample1.fa
```

Results in: `>sample1_chr1`, `>sample1_chr2`, etc.

## Complete Workflows

### Phylogenetic Analysis Preparation

Generate consensus for all samples for multiple sequence alignment:

```bash
mkdir -p consensus

# Reference as first sequence
cp reference.fa consensus/reference.fa

# Each sample
for sample in $(bcftools query -l variants.vcf.gz); do
    bcftools view -s "$sample" variants.vcf.gz | \
        bcftools filter -i 'QUAL>=30' | \
        bcftools consensus -f reference.fa > "consensus/${sample}.fa"
done

# Combine all sequences
cat consensus/*.fa > all_samples.fa
```

### Viral Consensus with Depth Masking

```bash
# Create depth mask
samtools depth -a aligned.bam | \
    awk '$3<20 {print $1"\t"$2-1"\t"$2}' | \
    bedtools merge > low_depth.bed

# Generate consensus
bcftools filter -i 'QUAL>=30 && INFO/DP>=20' variants.vcf.gz | \
    bcftools view -f PASS | \
    bcftools consensus -f reference.fa -m low_depth.bed > consensus.fa

# Report coverage
masked=$(awk '{sum += $3-$2} END {print sum}' low_depth.bed)
total=$(samtools faidx reference.fa -n 1 | awk 'NR==2 {print length}')
echo "Masked $masked of $total bases ($(echo "scale=1; $masked*100/$total" | bc)%)"
```

### Gene-Specific Haplotypes

```bash
gene_region="chr1:1000000-1010000"
sample="patient_001"

# Haplotype 1
bcftools consensus -f reference.fa -r "$gene_region" -s "$sample" -H 1 \
    variants.vcf.gz > gene_hap1.fa

# Haplotype 2
bcftools consensus -f reference.fa -r "$gene_region" -s "$sample" -H 2 \
    variants.vcf.gz > gene_hap2.fa
```

## Verification

### Count Applied Variants

```bash
# Variants in VCF
bcftools view -H input.vcf.gz | wc -l

# Check for warnings during consensus
bcftools consensus -f reference.fa input.vcf.gz 2>&1 | grep -c "warning"
```

### Compare Consensus to Reference

```bash
# Simple character comparison
diff <(grep -v "^>" reference.fa | tr -d '\n') \
     <(grep -v "^>" consensus.fa | tr -d '\n') | head
```

## Troubleshooting

### "not indexed"

VCF file needs an index:

```bash
bcftools index input.vcf.gz
```

### "sequence not found"

Chromosome names don't match between reference and VCF:

```bash
# Check reference chromosomes
grep "^>" reference.fa | head

# Check VCF chromosomes
bcftools view -H input.vcf.gz | cut -f1 | sort -u | head
```

Fix by renaming:
```bash
bcftools annotate --rename-chrs rename.txt input.vcf.gz -Oz -o renamed.vcf.gz
```

### "REF does not match"

VCF was called with a different reference. The reference FASTA must exactly match what was used for variant calling.

### "overlapping records"

Variants at the same position or overlapping indels. This is usually handled automatically but check the output:

```bash
bcftools consensus -f reference.fa input.vcf.gz 2>&1 | grep -i "overlap"
```

## Example Prompts

> "Generate a consensus FASTA sequence by applying my VCF variants to the reference"

> "Extract both haplotypes for a specific gene region from my phased VCF"

> "Create consensus sequences for all samples in my multi-sample VCF"

> "Mask low-coverage regions in my consensus sequence with N characters"

## See Also

- [bcftools consensus documentation](http://www.htslib.org/doc/bcftools.html#consensus)
- **variant-calling** - Generate VCF files
- **filtering-best-practices** - Filter variants before consensus
- **variant-normalization** - Normalize variants first
