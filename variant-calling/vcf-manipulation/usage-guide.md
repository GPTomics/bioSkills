# VCF Manipulation Usage Guide

## Overview

This guide covers merging, concatenating, sorting, and comparing VCF files.

## Prerequisites

- bcftools installed (`conda install -c bioconda bcftools`)
- Input VCFs should be compressed (`.vcf.gz`) and indexed for most operations
- Index files with `bcftools index input.vcf.gz`

## Quick Start

Tell your AI agent what you want to do:
- "Merge VCF files from different samples into a single cohort VCF"
- "Concatenate per-chromosome VCFs into a genome-wide file"
- "Compare two variant callsets and find shared and unique variants"
- "Extract specific samples from a multi-sample VCF"

## Understanding Merge vs Concat

These two operations are commonly confused:

| Operation | Use When | Example |
|-----------|----------|---------|
| **merge** | Same variants, different samples | Combining per-sample VCFs from same caller |
| **concat** | Same samples, different regions | Combining per-chromosome VCFs |

### Merge Scenario

You have variants called separately for each sample:
```
sample1.vcf.gz  →  contains only sample1
sample2.vcf.gz  →  contains only sample2
sample3.vcf.gz  →  contains only sample3
```

Use `bcftools merge` to combine into one multi-sample VCF.

### Concat Scenario

You have variants called in parallel by chromosome:
```
chr1.vcf.gz  →  all samples, chromosome 1
chr2.vcf.gz  →  all samples, chromosome 2
chr3.vcf.gz  →  all samples, chromosome 3
```

Use `bcftools concat` to combine into one genome-wide VCF.

## Merging VCF Files

### Basic Merge

```bash
bcftools merge sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz -Oz -o merged.vcf.gz
bcftools index merged.vcf.gz
```

### Merge from File List

For many samples, create a file list:

```bash
# Create list of VCF files
ls *.vcf.gz > vcf_list.txt

# Merge all
bcftools merge -l vcf_list.txt -Oz -o merged.vcf.gz
```

### Handling Missing Data

When samples have variants at different positions:

```bash
# Default: missing genotypes shown as ./.
bcftools merge sample1.vcf.gz sample2.vcf.gz -Oz -o merged.vcf.gz

# Treat missing as homozygous reference
bcftools merge --missing-to-ref sample1.vcf.gz sample2.vcf.gz -Oz -o merged.vcf.gz
```

### Sample Name Conflicts

If VCFs have the same sample name:

```bash
# Force merge (keeps both with modified names)
bcftools merge --force-samples file1.vcf.gz file2.vcf.gz -Oz -o merged.vcf.gz
```

### Merge with Multi-threading

```bash
bcftools merge --threads 4 -l vcf_list.txt -Oz -o merged.vcf.gz
```

## Concatenating VCF Files

### Concatenate by Chromosome

```bash
# Files must be in order
bcftools concat chr1.vcf.gz chr2.vcf.gz chr3.vcf.gz -Oz -o genome.vcf.gz
bcftools index genome.vcf.gz
```

### Concatenate All Chromosomes

```bash
bcftools concat chr{1..22}.vcf.gz chrX.vcf.gz chrY.vcf.gz chrM.vcf.gz -Oz -o genome.vcf.gz
```

### From File List (Ordered)

```bash
# files.txt must list VCFs in genomic order
bcftools concat -f files.txt -Oz -o concatenated.vcf.gz
```

### Allow Overlapping Regions

When files may contain overlapping regions:

```bash
bcftools concat -a file1.vcf.gz file2.vcf.gz -Oz -o combined.vcf.gz
```

### Remove Duplicates While Concatenating

```bash
# Remove exact duplicate records
bcftools concat -a -d exact file1.vcf.gz file2.vcf.gz -Oz -o combined.vcf.gz

# Remove duplicate variants (same position and alleles)
bcftools concat -a -d all file1.vcf.gz file2.vcf.gz -Oz -o combined.vcf.gz
```

## Sorting VCF Files

### Basic Sort

```bash
bcftools sort input.vcf -Oz -o sorted.vcf.gz
bcftools index sorted.vcf.gz
```

### Large File Sort

Use temporary directory and memory limit:

```bash
bcftools sort -T /scratch/tmp -m 4G input.vcf.gz -Oz -o sorted.vcf.gz
```

## Intersecting VCF Files

### Compare Two Call Sets

```bash
bcftools isec -p comparison_dir caller1.vcf.gz caller2.vcf.gz
```

This creates:
- `comparison_dir/0000.vcf` - Variants only in caller1
- `comparison_dir/0001.vcf` - Variants only in caller2
- `comparison_dir/0002.vcf` - Shared variants (from caller1)
- `comparison_dir/0003.vcf` - Shared variants (from caller2)

### Output Compressed Files

```bash
bcftools isec -p comparison_dir -Oz caller1.vcf.gz caller2.vcf.gz
```

### Get Only Shared Variants

```bash
# Variants present in both files
bcftools isec -n=2 -w1 caller1.vcf.gz caller2.vcf.gz -Oz -o shared.vcf.gz
```

### Get Private Variants

```bash
# Variants only in file1 (not in file2)
bcftools isec -C caller1.vcf.gz caller2.vcf.gz -Oz -o only_in_caller1.vcf.gz
```

### Multi-File Comparison

```bash
# Find variants present in all 3 files
bcftools isec -n=3 -w1 file1.vcf.gz file2.vcf.gz file3.vcf.gz -Oz -o in_all_three.vcf.gz

# Find variants present in at least 2 of 3 files
bcftools isec -n+2 -w1 file1.vcf.gz file2.vcf.gz file3.vcf.gz -Oz -o in_at_least_two.vcf.gz
```

### Understanding -n Options

| Option | Meaning |
|--------|---------|
| `-n=2` | Present in exactly 2 files |
| `-n+2` | Present in 2 or more files |
| `-n-2` | Present in fewer than 2 files |
| `-n~11` | Present in file1 AND file2 |
| `-n~10` | Present in file1 but NOT file2 |
| `-n~01` | Present in file2 but NOT file1 |

## Subsetting VCF Files

### Extract Specific Samples

```bash
bcftools view -s sample1,sample2,sample3 input.vcf.gz -Oz -o subset.vcf.gz
```

### Exclude Specific Samples

```bash
bcftools view -s ^sample4,sample5 input.vcf.gz -Oz -o without.vcf.gz
```

### From Sample List File

```bash
# samples.txt: one sample name per line
bcftools view -S samples.txt input.vcf.gz -Oz -o subset.vcf.gz

# Exclude samples in list
bcftools view -S ^samples.txt input.vcf.gz -Oz -o without.vcf.gz
```

### Extract Region

```bash
bcftools view -r chr1:1000000-2000000 input.vcf.gz -Oz -o region.vcf.gz
```

### Extract Multiple Regions

```bash
# From BED file
bcftools view -R targets.bed input.vcf.gz -Oz -o targets.vcf.gz

# Multiple regions
bcftools view -r chr1:1000-2000,chr2:3000-4000 input.vcf.gz -Oz -o regions.vcf.gz
```

## Renaming Samples

### Create Rename File

```bash
# Format: old_name<tab>new_name
cat > rename.txt << EOF
SRR123456	patient_001
SRR123457	patient_002
SRR123458	patient_003
EOF
```

### Apply Rename

```bash
bcftools reheader -s rename.txt input.vcf.gz -o renamed.vcf.gz
```

### Verify Rename

```bash
# Before
bcftools query -l input.vcf.gz

# After
bcftools query -l renamed.vcf.gz
```

## Splitting VCF Files

### Split by Sample

```bash
for sample in $(bcftools query -l input.vcf.gz); do
    bcftools view -s "$sample" input.vcf.gz -Oz -o "${sample}.vcf.gz"
    bcftools index "${sample}.vcf.gz"
done
```

### Split by Chromosome

```bash
for chr in chr1 chr2 chr3 chr4 chr5; do
    bcftools view -r "$chr" input.vcf.gz -Oz -o "${chr}.vcf.gz"
    bcftools index "${chr}.vcf.gz"
done
```

## Complete Workflow Examples

### Merge Per-Sample VCFs into Cohort

```bash
# List all sample VCFs
ls samples/*.vcf.gz > sample_vcfs.txt

# Merge with missing handled as reference
bcftools merge -l sample_vcfs.txt --missing-to-ref -Oz -o cohort.vcf.gz
bcftools index cohort.vcf.gz

# Verify samples
bcftools query -l cohort.vcf.gz
```

### Combine Parallel Chromosome Calls

```bash
# After calling variants in parallel by chromosome
for chr in chr{1..22} chrX chrY chrM; do
    bcftools index "output/${chr}.vcf.gz"
done

# Concatenate all
bcftools concat output/chr*.vcf.gz -Oz -o genome.vcf.gz
bcftools index genome.vcf.gz
```

### Compare Caller Results

```bash
# Compare GATK vs bcftools
bcftools isec -p caller_comparison gatk.vcf.gz bcftools.vcf.gz

# Count variants in each category
echo "GATK only: $(bcftools view -H caller_comparison/0000.vcf | wc -l)"
echo "bcftools only: $(bcftools view -H caller_comparison/0001.vcf | wc -l)"
echo "Shared: $(bcftools view -H caller_comparison/0002.vcf | wc -l)"
```

### Extract Case vs Control Samples

```bash
# samples_case.txt and samples_control.txt contain sample names

# Extract cases
bcftools view -S samples_case.txt cohort.vcf.gz -Oz -o cases.vcf.gz
bcftools index cases.vcf.gz

# Extract controls
bcftools view -S samples_control.txt cohort.vcf.gz -Oz -o controls.vcf.gz
bcftools index controls.vcf.gz
```

## Troubleshooting

### "different samples" Error

You're using `concat` when you should use `merge`, or vice versa:
- `merge` = same positions, different samples
- `concat` = same samples, different positions

### "not sorted" Error

Input files must be sorted for `concat`:

```bash
bcftools sort input.vcf.gz -Oz -o sorted.vcf.gz
```

Or allow unsorted with `-a`:

```bash
bcftools concat -a file1.vcf.gz file2.vcf.gz -Oz -o combined.vcf.gz
```

### "index required" Error

Many operations need indexed input:

```bash
bcftools index input.vcf.gz
```

### "sample name conflict" Warning

Samples have the same name in different files:

```bash
# Either rename samples first
bcftools reheader -s rename.txt input.vcf.gz -o renamed.vcf.gz

# Or force merge
bcftools merge --force-samples file1.vcf.gz file2.vcf.gz -Oz -o merged.vcf.gz
```

## Example Prompts

> "Merge VCF files from different samples into a single cohort VCF"

> "Concatenate per-chromosome VCFs into a genome-wide file"

> "Compare two variant callsets and find shared and unique variants"

> "Extract specific samples from a multi-sample VCF"

## See Also

- [bcftools merge documentation](http://www.htslib.org/doc/bcftools.html#merge)
- [bcftools concat documentation](http://www.htslib.org/doc/bcftools.html#concat)
- [bcftools isec documentation](http://www.htslib.org/doc/bcftools.html#isec)
- **vcf-basics** - View and query VCF files
- **filtering-best-practices** - Filter variants before manipulation
