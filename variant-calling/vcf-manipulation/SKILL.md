---
name: bio-vcf-manipulation
description: Merge, concatenate, sort, intersect, and subset VCF files using bcftools. Use when combining variant files, comparing call sets, or restructuring VCF data.
tool_type: cli
primary_tool: bcftools
---

# VCF Manipulation

Merge, concat, sort, and compare VCF files using bcftools.

## Operations Overview

| Operation | Command | Use Case |
|-----------|---------|----------|
| Merge | `bcftools merge` | Combine samples from multiple VCFs |
| Concat | `bcftools concat` | Combine regions from multiple VCFs |
| Sort | `bcftools sort` | Sort unsorted VCF |
| Intersect | `bcftools isec` | Compare/intersect call sets |
| Subset | `bcftools view` | Extract samples or regions |

## bcftools merge

Combine multiple VCF files with **different samples** at the same positions.

### Basic Merge

```bash
bcftools merge sample1.vcf.gz sample2.vcf.gz -Oz -o merged.vcf.gz
```

### Merge Multiple Files

```bash
bcftools merge *.vcf.gz -Oz -o all_samples.vcf.gz
```

### Merge from File List

```bash
# files.txt: one VCF path per line
bcftools merge -l files.txt -Oz -o merged.vcf.gz
```

### Handle Missing Genotypes

```bash
# Output missing genotypes as ./. (default)
bcftools merge sample1.vcf.gz sample2.vcf.gz -Oz -o merged.vcf.gz

# Output missing as reference (0/0)
bcftools merge --missing-to-ref sample1.vcf.gz sample2.vcf.gz -Oz -o merged.vcf.gz
```

### Force Sample Names

When sample names conflict:

```bash
bcftools merge --force-samples sample1.vcf.gz sample2.vcf.gz -Oz -o merged.vcf.gz
```

### Merge Specific Regions

```bash
bcftools merge -r chr1:1000000-2000000 sample1.vcf.gz sample2.vcf.gz -Oz -o merged.vcf.gz
```

## bcftools concat

Combine VCF files with **same samples** from different regions.

### Concatenate Chromosomes

```bash
bcftools concat chr1.vcf.gz chr2.vcf.gz chr3.vcf.gz -Oz -o genome.vcf.gz
```

### Concatenate All Chromosomes

```bash
bcftools concat chr*.vcf.gz -Oz -o genome.vcf.gz
```

### From File List

```bash
# files.txt: one VCF path per line (in order)
bcftools concat -f files.txt -Oz -o concatenated.vcf.gz
```

### Allow Overlapping Regions

```bash
bcftools concat -a chr1_part1.vcf.gz chr1_part2.vcf.gz -Oz -o chr1.vcf.gz
```

### Remove Duplicates

```bash
bcftools concat -a -d all file1.vcf.gz file2.vcf.gz -Oz -o merged.vcf.gz
```

Options for `-d`:
- `snps` - Remove duplicate SNPs
- `indels` - Remove duplicate indels
- `both` - Remove duplicate SNPs and indels
- `all` - Remove all duplicates
- `exact` - Remove exact duplicates only

## bcftools sort

Sort VCF by chromosome and position.

### Basic Sort

```bash
bcftools sort input.vcf -Oz -o sorted.vcf.gz
```

### With Temporary Directory

For large files:

```bash
bcftools sort -T /tmp input.vcf.gz -Oz -o sorted.vcf.gz
```

### Memory Limit

```bash
bcftools sort -m 4G input.vcf.gz -Oz -o sorted.vcf.gz
```

## bcftools isec

Intersect and compare VCF files.

### Find Shared Variants

```bash
bcftools isec -p output_dir sample1.vcf.gz sample2.vcf.gz
```

Creates:
- `0000.vcf` - Private to sample1
- `0001.vcf` - Private to sample2
- `0002.vcf` - Shared (sample1 records)
- `0003.vcf` - Shared (sample2 records)

### Output Compressed

```bash
bcftools isec -p output_dir -Oz sample1.vcf.gz sample2.vcf.gz
```

### Intersection Only

```bash
bcftools isec -p output_dir -n=2 sample1.vcf.gz sample2.vcf.gz
# Only outputs variants present in exactly 2 files
```

### Comparison Options

| Flag | Description |
|------|-------------|
| `-n=2` | Present in exactly 2 files |
| `-n+2` | Present in 2 or more files |
| `-n-2` | Present in fewer than 2 files |
| `-n~11` | Boolean: file1 AND file2 |
| `-n~10` | Boolean: file1 AND NOT file2 |

### Two-File Intersection

```bash
# Variants in both files
bcftools isec -n=2 -w1 sample1.vcf.gz sample2.vcf.gz -Oz -o shared.vcf.gz

# Variants only in sample1
bcftools isec -n~10 -w1 sample1.vcf.gz sample2.vcf.gz -Oz -o only_sample1.vcf.gz
```

### Complement Mode

```bash
# Variants in file1 not in file2
bcftools isec -C sample1.vcf.gz sample2.vcf.gz -Oz -o unique.vcf.gz
```

## Subsetting VCF Files

### Extract Samples

```bash
bcftools view -s sample1,sample2 input.vcf.gz -Oz -o subset.vcf.gz
```

### Exclude Samples

```bash
bcftools view -s ^sample3 input.vcf.gz -Oz -o without_sample3.vcf.gz
```

### From Sample List File

```bash
# samples.txt: one sample name per line
bcftools view -S samples.txt input.vcf.gz -Oz -o subset.vcf.gz
```

### Extract Region

```bash
bcftools view -r chr1:1000000-2000000 input.vcf.gz -Oz -o region.vcf.gz
```

### Extract Multiple Regions

```bash
bcftools view -R regions.bed input.vcf.gz -Oz -o targets.vcf.gz
```

## Renaming Samples

### Single Sample

```bash
echo "old_name new_name" > rename.txt
bcftools reheader -s rename.txt input.vcf.gz -o renamed.vcf.gz
```

### Multiple Samples

```bash
# rename.txt format: old_name new_name
cat > rename.txt << EOF
sample1 patient_001
sample2 patient_002
sample3 patient_003
EOF

bcftools reheader -s rename.txt input.vcf.gz -o renamed.vcf.gz
```

## Splitting VCF Files

### Split by Sample

```bash
for sample in $(bcftools query -l input.vcf.gz); do
    bcftools view -s "$sample" input.vcf.gz -Oz -o "${sample}.vcf.gz"
done
```

### Split by Chromosome

```bash
for chr in $(bcftools view -h input.vcf.gz | grep "^##contig" | sed 's/.*ID=\([^,]*\).*/\1/'); do
    bcftools view -r "$chr" input.vcf.gz -Oz -o "${chr}.vcf.gz"
done
```

### Split Multiallelic Sites

```bash
bcftools norm -m-any input.vcf.gz -Oz -o split.vcf.gz
```

## Common Workflows

### Merge Cohort VCFs

```bash
# Create file list
ls *.vcf.gz > files.txt

# Merge all samples
bcftools merge -l files.txt -Oz -o cohort.vcf.gz
bcftools index cohort.vcf.gz
```

### Combine Chromosome VCFs

```bash
# After parallel variant calling by chromosome
bcftools concat chr{1..22}.vcf.gz chrX.vcf.gz chrY.vcf.gz -Oz -o genome.vcf.gz
bcftools index genome.vcf.gz
```

### Compare Two Callers

```bash
# Find variants called by both GATK and bcftools
bcftools isec -p comparison gatk.vcf.gz bcftools.vcf.gz

# Count results
wc -l comparison/*.vcf
```

### Extract Passing Variants

```bash
bcftools view -f PASS input.vcf.gz -Oz -o pass_only.vcf.gz
bcftools index pass_only.vcf.gz
```

## cyvcf2 Python Operations

**Note:** True VCF merging (combining samples at matching positions) is complex.
Use `bcftools merge` for production work. cyvcf2 is better for filtering/querying.

### Concatenate Records (Not True Merge)

```python
from cyvcf2 import VCF, Writer

# WARNING: This concatenates records, not a true merge
# For actual merging of samples, use bcftools merge
vcf1 = VCF('file1.vcf.gz')
writer = Writer('combined.vcf', vcf1)

for variant in vcf1:
    writer.write_record(variant)

writer.close()
vcf1.close()
```

### Find Shared Positions

```python
from cyvcf2 import VCF

# Load positions from first VCF
vcf1_positions = set()
for variant in VCF('sample1.vcf.gz'):
    vcf1_positions.add((variant.CHROM, variant.POS))

# Check second VCF
shared = 0
unique = 0
for variant in VCF('sample2.vcf.gz'):
    if (variant.CHROM, variant.POS) in vcf1_positions:
        shared += 1
    else:
        unique += 1

print(f'Shared: {shared}')
print(f'Unique to sample2: {unique}')
```

## Quick Reference

| Task | Command |
|------|---------|
| Merge samples | `bcftools merge s1.vcf.gz s2.vcf.gz -Oz -o merged.vcf.gz` |
| Concat regions | `bcftools concat chr1.vcf.gz chr2.vcf.gz -Oz -o all.vcf.gz` |
| Sort VCF | `bcftools sort input.vcf -Oz -o sorted.vcf.gz` |
| Intersect | `bcftools isec -p dir a.vcf.gz b.vcf.gz` |
| Extract samples | `bcftools view -s sample1 input.vcf.gz` |
| Rename samples | `bcftools reheader -s names.txt input.vcf.gz` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `different samples` | merge vs concat confusion | Use merge for samples, concat for regions |
| `not sorted` | Unsorted input to concat | Sort first or use `-a` flag |
| `sample name conflict` | Duplicate sample names | Use `--force-samples` |
| `index required` | Missing index for merge/isec | Run `bcftools index` first |

## Related Skills

- vcf-basics - View and query VCF files
- filtering-best-practices - Filter variants before manipulation
- variant-normalization - Normalize before comparing
- vcf-statistics - Compare statistics after manipulation
