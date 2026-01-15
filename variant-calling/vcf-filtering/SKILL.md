---
name: bio-vcf-filtering
description: Filter VCF variants by quality, depth, allele frequency, and other criteria using bcftools filter and view expressions. Use when removing low-quality variants or selecting specific variant types.
tool_type: cli
primary_tool: bcftools
---

# VCF Filtering

Filter variants using bcftools filter and view.

## Two Approaches

| Command | Purpose |
|---------|---------|
| `bcftools filter` | Apply soft/hard filters using expressions |
| `bcftools view` | Select variants by type, region, sample |

## bcftools filter

### Basic Quality Filter

```bash
bcftools filter -e 'QUAL<30' input.vcf.gz -o filtered.vcf
```

### Soft Filter (Mark, Don't Remove)

```bash
bcftools filter -s 'LowQual' -e 'QUAL<30' input.vcf.gz -o marked.vcf
```

- `-s` sets FILTER column instead of removing
- Variants failing filter get "LowQual" in FILTER column
- Passing variants remain "PASS" or "."

### Hard Filter (Remove Variants)

```bash
bcftools filter -e 'QUAL<30' input.vcf.gz -o filtered.vcf
# Variants with QUAL<30 are removed
```

### Include Instead of Exclude

```bash
bcftools filter -i 'QUAL>=30' input.vcf.gz -o filtered.vcf
# Only keeps variants with QUAL>=30
```

## Common Filter Expressions

### Quality Score

```bash
bcftools filter -e 'QUAL<30' input.vcf.gz
```

### Read Depth

```bash
# INFO depth
bcftools filter -e 'INFO/DP<10' input.vcf.gz

# FORMAT depth (any sample)
bcftools filter -e 'FORMAT/DP<10' input.vcf.gz

# FORMAT depth (all samples)
bcftools filter -e 'MIN(FORMAT/DP)<10' input.vcf.gz
```

### Allele Frequency

```bash
# Minimum allele frequency
bcftools filter -e 'INFO/AF<0.05' input.vcf.gz

# Maximum allele frequency (rare variants)
bcftools filter -i 'INFO/AF<0.01' input.vcf.gz
```

### Strand Bias

```bash
bcftools filter -e 'INFO/FS>60' input.vcf.gz
```

### Mapping Quality

```bash
bcftools filter -e 'INFO/MQ<40' input.vcf.gz
```

### Combined Filters

```bash
bcftools filter -e 'QUAL<30 || INFO/DP<10 || INFO/MQ<40' input.vcf.gz
```

## Expression Syntax

### Operators

| Operator | Meaning |
|----------|---------|
| `<`, `<=`, `>`, `>=` | Comparison |
| `=`, `==` | Equals |
| `!=` | Not equals |
| `&&`, `\|\|` | AND, OR |
| `!` | NOT |

### Aggregate Functions

| Function | Description |
|----------|-------------|
| `MIN(x)` | Minimum across samples |
| `MAX(x)` | Maximum across samples |
| `AVG(x)` | Average across samples |
| `SUM(x)` | Sum across samples |

### Special Values

| Value | Description |
|-------|-------------|
| `.` | Missing value |
| `PASS` | Filter status |

### Check for Missing Values

```bash
# Exclude missing QUAL
bcftools filter -e 'QUAL="."' input.vcf.gz

# Include only if DP exists
bcftools filter -i 'INFO/DP!="."' input.vcf.gz
```

## bcftools view Filtering

### Filter by Variant Type

```bash
# SNPs only
bcftools view -v snps input.vcf.gz -o snps.vcf.gz

# Indels only
bcftools view -v indels input.vcf.gz -o indels.vcf.gz

# Both SNPs and indels
bcftools view -v snps,indels input.vcf.gz

# Exclude SNPs
bcftools view -V snps input.vcf.gz -o no_snps.vcf.gz
```

### Variant Types

| Type | Description |
|------|-------------|
| `snps` | Single nucleotide polymorphisms |
| `indels` | Insertions and deletions |
| `mnps` | Multi-nucleotide polymorphisms |
| `other` | Other variant types |

### Filter by Region

```bash
bcftools view -r chr1:1000000-2000000 input.vcf.gz -o region.vcf.gz

# Multiple regions
bcftools view -r chr1:1000-2000,chr2:3000-4000 input.vcf.gz
```

### Filter by Samples

```bash
# Include samples
bcftools view -s sample1,sample2 input.vcf.gz -o subset.vcf.gz

# Exclude samples
bcftools view -s ^sample3,sample4 input.vcf.gz -o subset.vcf.gz
```

### Filter by Genotype

```bash
# Sites with at least one non-reference genotype
bcftools view -c 1 input.vcf.gz -o polymorphic.vcf.gz

# Sites where all samples are non-reference
bcftools view -C 0 input.vcf.gz -o fixed_alt.vcf.gz
```

### Apply FILTER Column

```bash
# Keep only PASS variants
bcftools view -f PASS input.vcf.gz -o pass_only.vcf.gz

# Keep PASS and missing filter
bcftools view -f .,PASS input.vcf.gz
```

## Multi-Step Filtering

### Pipeline Approach

```bash
bcftools filter -e 'QUAL<30' input.vcf.gz | \
    bcftools filter -e 'INFO/DP<10' | \
    bcftools view -v snps -Oz -o filtered_snps.vcf.gz
```

### Named Soft Filters

```bash
bcftools filter -s 'LowQual' -e 'QUAL<30' input.vcf.gz | \
    bcftools filter -s 'LowDepth' -e 'INFO/DP<10' -Oz -o marked.vcf.gz

# Later, remove all filtered variants
bcftools view -f PASS marked.vcf.gz -Oz -o pass_only.vcf.gz
```

## cyvcf2 Python Filtering

### Basic Filter

```python
from cyvcf2 import VCF, Writer

vcf = VCF('input.vcf.gz')
writer = Writer('filtered.vcf', vcf)

for variant in vcf:
    if variant.QUAL and variant.QUAL >= 30:
        writer.write_record(variant)

writer.close()
vcf.close()
```

### Multi-Criteria Filter

```python
from cyvcf2 import VCF, Writer

vcf = VCF('input.vcf.gz')
writer = Writer('filtered.vcf', vcf)

for variant in vcf:
    qual = variant.QUAL or 0
    dp = variant.INFO.get('DP', 0)

    if qual >= 30 and dp >= 10:
        writer.write_record(variant)

writer.close()
vcf.close()
```

### Filter by Genotype

```python
from cyvcf2 import VCF, Writer

vcf = VCF('input.vcf.gz')
writer = Writer('filtered.vcf', vcf)

for variant in vcf:
    # Keep if at least one sample is heterozygous
    # gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
    if 1 in variant.gt_types:
        writer.write_record(variant)

writer.close()
vcf.close()
```

### Filter by Sample Depth

```python
from cyvcf2 import VCF, Writer
import numpy as np

vcf = VCF('input.vcf.gz')
writer = Writer('filtered.vcf', vcf)

for variant in vcf:
    depths = variant.format('DP')
    if depths is not None and np.min(depths) >= 10:
        writer.write_record(variant)

writer.close()
vcf.close()
```

## Common Filter Sets

### Stringent Quality Filter

```bash
bcftools filter -e 'QUAL<50 || INFO/DP<20 || INFO/MQ<50' input.vcf.gz
```

### GATK-Style Hard Filters (SNPs)

```bash
bcftools filter -s 'GATK_SNP' \
    -e 'TYPE="snp" && (QUAL<30 || INFO/DP<10 || INFO/MQ<40 || INFO/FS>60)' \
    input.vcf.gz
```

### GATK-Style Hard Filters (Indels)

```bash
bcftools filter -s 'GATK_INDEL' \
    -e 'TYPE="indel" && (QUAL<30 || INFO/DP<10 || INFO/FS>200)' \
    input.vcf.gz
```

### Population Genetics Filter

```bash
bcftools filter -i 'INFO/AF>=0.05 && INFO/AF<=0.95' input.vcf.gz
```

## Quick Reference

| Task | Command |
|------|---------|
| Quality filter | `bcftools filter -e 'QUAL<30' in.vcf.gz` |
| Depth filter | `bcftools filter -e 'INFO/DP<10' in.vcf.gz` |
| SNPs only | `bcftools view -v snps in.vcf.gz` |
| Indels only | `bcftools view -v indels in.vcf.gz` |
| PASS only | `bcftools view -f PASS in.vcf.gz` |
| Soft filter | `bcftools filter -s 'LowQ' -e 'QUAL<30' in.vcf.gz` |
| Region | `bcftools view -r chr1:1-1000 in.vcf.gz` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `no such INFO tag` | Tag not in VCF | Check header with `bcftools view -h` |
| `syntax error` | Invalid expression | Check operator syntax (`\|\|` not `or`) |
| `empty output` | Filter too strict | Relax thresholds |

## Related Skills

- **variant-calling** - Generate VCF files
- **vcf-basics** - View and query VCF files
- **vcf-statistics** - Evaluate filter effects
- **variant-normalization** - Normalize before filtering
