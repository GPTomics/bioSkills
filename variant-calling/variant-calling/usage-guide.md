# Variant Calling Usage Guide

## Overview

This guide covers calling germline SNPs and indels from aligned reads using bcftools mpileup/call, and choosing the right calling engine. bcftools is a position-based genotype-likelihood caller: fast, transparent, needs no training data, and is best suited for simple germline SNPs, non-model/organelle/microbial genomes, and rapid exploratory analysis. It is materially weaker on indels and in low-complexity, segmental-duplication, and MHC regions than local-reassembly callers, which discard the local alignment and reconstruct haplotypes. For highest indel/difficult-region accuracy on human data, hand off to DeepVariant or GATK HaplotypeCaller (DRAGEN-mode); for scalable cohort genotyping, use the GVCF joint-calling workflow.

## Prerequisites

- bcftools installed (`conda install -c bioconda bcftools`)
- Sorted, indexed BAM file
- Reference FASTA (same one used for alignment)
- Reference must be indexed (`samtools faidx reference.fa`)

## Quick Start

Tell your AI agent what you want to do:
- "Call SNPs and indels from my aligned BAM file"
- "Which caller should I use for my Illumina human WGS versus my bacterial genome?"
- "Call variants only in my target regions defined by a BED file"
- "Generate a VCF with depth and allelic depth annotations"
- "Call my mitochondrial or chrX/chrY variants with the right ploidy"

## Understanding the Workflow

Variant calling with bcftools is a two-step process:

1. **mpileup** - Generates genotype likelihoods at each position
2. **call** - Uses likelihoods to make variant calls

```
BAM + Reference -> mpileup -> genotype likelihoods -> call -> VCF
```

## Basic Variant Calling

### Simple Pipeline

```bash
bcftools mpileup -f reference.fa input.bam | bcftools call -mv -o variants.vcf
```

- `-f reference.fa` - Reference FASTA (required)
- `-m` - Use multiallelic caller (recommended)
- `-v` - Output variant sites only (skip reference calls)
- `-o variants.vcf` - Output file

### Compressed Output

For larger files, use compressed VCF:

```bash
bcftools mpileup -f reference.fa input.bam | bcftools call -mv -Oz -o variants.vcf.gz
bcftools index variants.vcf.gz
```

## Quality Filtering in mpileup

### Mapping Quality Filter

Exclude reads with low mapping quality:

```bash
bcftools mpileup -f reference.fa -q 20 input.bam | bcftools call -mv -o variants.vcf
```

### Base Quality Filter

Exclude bases with low quality scores:

```bash
bcftools mpileup -f reference.fa -Q 20 input.bam | bcftools call -mv -o variants.vcf
```

### Combined Quality Filters

```bash
bcftools mpileup -f reference.fa -q 20 -Q 20 input.bam | bcftools call -mv -o variants.vcf
```

## Adding Annotations

Annotations provide additional information about each variant.

### Basic Depth Annotation

```bash
bcftools mpileup -f reference.fa -a DP,AD input.bam | bcftools call -mv -o variants.vcf
```

- `DP` - Total read depth
- `AD` - Allelic depths (ref and alt counts)

### Full Annotation Set

```bash
bcftools mpileup -f reference.fa \
    -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,INFO/AD \
    input.bam | bcftools call -mv -o variants.vcf
```

- `FORMAT/DP` - Per-sample read depth
- `FORMAT/AD` - Per-sample allelic depths
- `FORMAT/ADF` - Forward strand allelic depths
- `FORMAT/ADR` - Reverse strand allelic depths
- `INFO/AD` - Total allelic depths across samples

## Region-Based Calling

### Single Region

```bash
bcftools mpileup -f reference.fa -r chr1:1000000-2000000 input.bam | \
    bcftools call -mv -o region.vcf
```

### Multiple Regions

```bash
bcftools mpileup -f reference.fa -r chr1:1000-2000,chr2:3000-4000 input.bam | \
    bcftools call -mv -o regions.vcf
```

### From BED File

```bash
bcftools mpileup -f reference.fa -R targets.bed input.bam | \
    bcftools call -mv -o targets.vcf
```

## Multi-Sample Calling

### Multiple BAM Files

```bash
bcftools mpileup -f reference.fa sample1.bam sample2.bam sample3.bam | \
    bcftools call -mv -o cohort.vcf
```

### BAM List File

Create a file with one BAM path per line:

```
# bams.txt
/path/to/sample1.bam
/path/to/sample2.bam
/path/to/sample3.bam
```

Then call:

```bash
bcftools mpileup -f reference.fa -b bams.txt | bcftools call -mv -o cohort.vcf
```

## Calling Models

### Multiallelic Caller (Recommended)

```bash
bcftools mpileup -f reference.fa input.bam | bcftools call -m -o variants.vcf
```

The multiallelic caller (`-m`) handles sites with multiple alternate alleles correctly.

### Consensus Caller (Legacy)

```bash
bcftools mpileup -f reference.fa input.bam | bcftools call -c -o variants.vcf
```

The consensus caller (`-c`) is legacy; use it only to reproduce older pipelines. The multiallelic caller is preferred for all new work.

## Ploidy Settings

### Haploid Organisms

```bash
bcftools mpileup -f reference.fa input.bam | \
    bcftools call -m --ploidy 1 -o variants.vcf
```

### Ploidy File for Mixed Samples

Create a ploidy file for complex cases (e.g., sex chromosomes):

```
# ploidy.txt
# CHROM  START  END  SEX  PLOIDY
chrX     1      -1   M    1
chrX     1      -1   F    2
chrY     1      -1   M    1
chrY     1      -1   F    0
*        1      -1   *    2
```

```bash
bcftools mpileup -f reference.fa input.bam | \
    bcftools call -m --ploidy-file ploidy.txt -o variants.vcf
```

## Output Options

### All Sites (Including Reference)

Remove `-v` to output all sites:

```bash
bcftools mpileup -f reference.fa input.bam | bcftools call -m -o all_sites.vcf
```

### Output Formats

| Flag | Format | Description |
|------|--------|-------------|
| `-Ov` | VCF | Uncompressed VCF (default) |
| `-Oz` | VCF.gz | Compressed VCF |
| `-Ou` | BCF | Uncompressed BCF |
| `-Ob` | BCF | Compressed BCF |

### Uncompressed BCF for Pipelines

For multi-step pipelines, use uncompressed BCF between steps:

```bash
bcftools mpileup -Ou -f reference.fa input.bam | \
    bcftools call -mv -Ou | \
    bcftools filter -Oz -o filtered.vcf.gz
```

## Performance Optimization

### Multi-Threading

```bash
bcftools mpileup -f reference.fa --threads 4 input.bam | \
    bcftools call -mv --threads 4 -o variants.vcf
```

### Parallel by Chromosome

```bash
for chr in chr1 chr2 chr3 chr4 chr5; do
    bcftools mpileup -Ou -f reference.fa -r "$chr" input.bam | \
        bcftools call -mv -Oz -o "${chr}.vcf.gz" &
done
wait

# Concatenate results
bcftools concat -Oz -o all.vcf.gz chr*.vcf.gz
bcftools index all.vcf.gz
```

### Max Depth Setting

For high-coverage data, set maximum depth to prevent memory issues:

```bash
bcftools mpileup -f reference.fa -d 1000 input.bam | bcftools call -mv -o variants.vcf
```

Default is 250. Increase for high-coverage samples.

## Complete Workflow Examples

### Standard SNP/Indel Calling

```bash
# Call variants with quality filters and annotations
bcftools mpileup -Ou -f reference.fa \
    -q 20 -Q 20 \
    -a FORMAT/DP,FORMAT/AD \
    input.bam | \
bcftools call -mv -Oz -o variants.vcf.gz

# Index for downstream use
bcftools index variants.vcf.gz

# Basic stats
bcftools stats variants.vcf.gz | head -50
```

### Multi-Sample Cohort Calling

```bash
# Call variants across samples
bcftools mpileup -Ou -f reference.fa \
    -q 20 -Q 20 \
    -a FORMAT/DP,FORMAT/AD \
    sample1.bam sample2.bam sample3.bam | \
bcftools call -mv -Oz -o cohort.vcf.gz

bcftools index cohort.vcf.gz

# Check sample names
bcftools query -l cohort.vcf.gz
```

### Targeted Region Calling

```bash
# Call variants in target regions only
bcftools mpileup -Ou -f reference.fa \
    -R targets.bed \
    -q 20 -Q 20 \
    -a FORMAT/DP,FORMAT/AD \
    input.bam | \
bcftools call -mv -Oz -o targets.vcf.gz

bcftools index targets.vcf.gz
```

## Troubleshooting

### "no FASTA reference"

Reference FASTA must be provided:

```bash
bcftools mpileup -f reference.fa input.bam | bcftools call -mv
```

### "reference mismatch"

The reference must match what was used for alignment. Check contig names:

```bash
# Check BAM header contigs
samtools view -H input.bam | grep @SQ

# Check reference contigs
grep "^>" reference.fa
```

### "no variants called"

Possible causes:
- Low coverage - check depth with `samtools depth`
- Quality thresholds too strict - try lower `-q` and `-Q` values
- Wrong reference genome

### Low Variant Quality

Add annotations to diagnose:

```bash
bcftools mpileup -f reference.fa -a DP,AD,ADF,ADR input.bam | \
    bcftools call -mv -o variants.vcf

# Check depths
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' variants.vcf
```

## Example Prompts

> "Call SNPs and indels from my aligned BAM file using bcftools"

> "Help me choose between bcftools, GATK HaplotypeCaller, and DeepVariant for my dataset"

> "Call variants only in my target regions defined by a BED file"

> "Generate a VCF with depth and allelic depth annotations"

> "Call my haploid bacterial genome, or set sex-chromosome ploidy for a human sample"

> "Compare my bcftools and GATK callsets fairly against a GIAB truth set"

## What the Agent Will Do

1. Verify the BAM file is sorted, indexed, and matches the reference
2. Select appropriate mpileup options (quality thresholds, annotations, regions)
3. Run bcftools mpileup piped to bcftools call with the multiallelic caller
4. Compress and index the output VCF
5. Report basic variant statistics (SNP/indel counts, Ti/Tv ratio)

## Tips

- Use `-Ou` (uncompressed BCF) between piped steps for speed
- Set `-d` (max depth) to 3-4x expected mean coverage to avoid memory issues
- Always request FORMAT/DP and FORMAT/AD annotations for downstream filtering
- For parallel processing, split by chromosome and concatenate results
- bcftools is fast but less accurate than GATK/DeepVariant for indels in homopolymers, segdups, and MHC; switch engines when those regions matter
- bcftools mpileup applies BAQ by default (downweights base qualities near likely misalignments to suppress false SNPs around indels); `-B` disables it, `-E` recomputes it
- Normalize with `bcftools norm -f ref.fa -m -any` before comparing, merging, or annotating any callset; un-normalized records mismatch spuriously
- Benchmark callers with hap.py + vcfeval inside the GIAB confident-region BED, reporting SNVs and indels separately - never `bcftools isec` on raw records

## Related Skills

- variant-calling/vcf-basics - View and query resulting VCF files
- variant-calling/filtering-best-practices - Filter variants by quality
- variant-calling/variant-normalization - Normalize indels after calling
- variant-calling/vcf-statistics - Ti/Tv, counts, and callset QC
- variant-calling/gatk-variant-calling - Local-reassembly calling with GATK HaplotypeCaller
- variant-calling/deepvariant - Deep-learning caller for best indel/difficult-region accuracy
- variant-calling/joint-calling - Scalable cohort genotyping via the GVCF workflow
