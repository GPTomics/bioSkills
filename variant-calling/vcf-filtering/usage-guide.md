# VCF Filtering Usage Guide

This guide covers filtering variants by quality, type, and other criteria.

## Prerequisites

- bcftools installed (`conda install -c bioconda bcftools`)
- cyvcf2 installed for Python filtering (`pip install cyvcf2`)
- Input VCF/BCF file (compressed and indexed recommended)

## Understanding Filter Approaches

bcftools offers two main filtering commands:

1. **bcftools filter** - Expression-based filtering on annotations
2. **bcftools view** - Selection by type, region, sample

Both can exclude variants (hard filter) or mark them (soft filter).

## Soft vs Hard Filtering

### Hard Filtering

Removes variants that fail criteria:

```bash
bcftools filter -e 'QUAL<30' input.vcf.gz -o filtered.vcf
# Variants with QUAL<30 are gone
```

### Soft Filtering

Marks variants in FILTER column but keeps them:

```bash
bcftools filter -s 'LowQual' -e 'QUAL<30' input.vcf.gz -o marked.vcf
# Variants with QUAL<30 have FILTER="LowQual"
# Passing variants have FILTER="PASS" or "."
```

Soft filtering is useful when you want to:
- Keep all variants for review
- Apply multiple filter names
- Delay the decision about which filters to enforce

Later, convert soft to hard filter:

```bash
bcftools view -f PASS marked.vcf.gz -o pass_only.vcf
```

## Filter Expressions

### Basic Syntax

Include variants matching expression:
```bash
bcftools filter -i 'QUAL>=30' input.vcf.gz
```

Exclude variants matching expression:
```bash
bcftools filter -e 'QUAL<30' input.vcf.gz
```

Both commands produce the same result for this example.

### Accessing Fields

```bash
# QUAL column directly
QUAL<30

# INFO fields
INFO/DP<10
INFO/AF<0.05
INFO/MQ<40

# FORMAT fields (sample-level)
FORMAT/DP<10
FORMAT/GQ<20
```

### Operators

| Operator | Meaning | Example |
|----------|---------|---------|
| `<` | Less than | `QUAL<30` |
| `<=` | Less or equal | `QUAL<=30` |
| `>` | Greater than | `INFO/DP>100` |
| `>=` | Greater or equal | `INFO/DP>=10` |
| `=` or `==` | Equals | `TYPE="snp"` |
| `!=` | Not equals | `FILTER!="PASS"` |
| `&&` | AND | `QUAL>=30 && INFO/DP>=10` |
| `\|\|` | OR | `QUAL<30 \|\| INFO/DP<10` |
| `!` | NOT | `!INFO/DB` |

### Aggregate Functions

For multi-sample VCFs, aggregate across samples:

| Function | Description | Example |
|----------|-------------|---------|
| `MIN(x)` | Minimum value | `MIN(FORMAT/DP)<10` |
| `MAX(x)` | Maximum value | `MAX(FORMAT/GQ)<20` |
| `AVG(x)` | Average value | `AVG(FORMAT/DP)<15` |
| `SUM(x)` | Sum | `SUM(FORMAT/DP)<50` |

Example:
```bash
# Keep variants where ALL samples have DP>=10
bcftools filter -e 'MIN(FORMAT/DP)<10' input.vcf.gz
```

### Missing Values

Check for missing annotations:

```bash
# Exclude variants with missing QUAL
bcftools filter -e 'QUAL="."' input.vcf.gz

# Include only variants with DP annotation
bcftools filter -i 'INFO/DP!="."' input.vcf.gz
```

## Filtering by Variant Type

### Using bcftools view

```bash
# SNPs only
bcftools view -v snps input.vcf.gz -o snps.vcf.gz

# Indels only
bcftools view -v indels input.vcf.gz -o indels.vcf.gz

# SNPs and indels (exclude other types)
bcftools view -v snps,indels input.vcf.gz

# Exclude indels
bcftools view -V indels input.vcf.gz -o no_indels.vcf.gz
```

### Using Expressions

```bash
# SNPs only
bcftools filter -i 'TYPE="snp"' input.vcf.gz

# Indels only
bcftools filter -i 'TYPE="indel"' input.vcf.gz
```

## Filtering by Region

### Single Region

```bash
bcftools view -r chr1:1000000-2000000 input.vcf.gz
```

### Multiple Regions

```bash
bcftools view -r chr1:1000-2000,chr2:3000-4000 input.vcf.gz
```

### From BED File

```bash
bcftools view -R targets.bed input.vcf.gz
```

### Exclude Regions

```bash
bcftools view -T ^exclude.bed input.vcf.gz
```

## Filtering by Sample

### Include Specific Samples

```bash
bcftools view -s sample1,sample2 input.vcf.gz
```

### From Sample List File

```bash
# samples.txt: one sample name per line
bcftools view -S samples.txt input.vcf.gz
```

### Exclude Samples

```bash
bcftools view -s ^sample3,sample4 input.vcf.gz
```

## Filtering by Genotype

### Sites with Alternate Alleles

```bash
# At least one sample has non-ref genotype
bcftools view -c 1 input.vcf.gz

# At least 5 alternate alleles across samples
bcftools view -c 5 input.vcf.gz
```

### Expression-Based Genotype Filtering

```bash
# Keep heterozygous sites
bcftools filter -i 'COUNT(GT="het")>0' input.vcf.gz

# Keep sites where all samples are called
bcftools filter -e 'COUNT(GT="mis")>0' input.vcf.gz
```

## Common Filter Recipes

### Quality Control Filter

Standard quality filters for most analyses:

```bash
bcftools filter -e 'QUAL<30 || INFO/DP<10' input.vcf.gz -Oz -o qc.vcf.gz
```

### Stringent Filter

For high-confidence variant sets:

```bash
bcftools filter -e 'QUAL<50 || INFO/DP<20 || INFO/MQ<50' input.vcf.gz -Oz -o stringent.vcf.gz
```

### SNP Hard Filters (GATK-Style)

```bash
bcftools view -v snps input.vcf.gz | \
    bcftools filter -e 'QUAL<30 || INFO/DP<10 || INFO/MQ<40 || INFO/FS>60' \
    -Oz -o filtered_snps.vcf.gz
```

### Indel Hard Filters (GATK-Style)

```bash
bcftools view -v indels input.vcf.gz | \
    bcftools filter -e 'QUAL<30 || INFO/DP<10 || INFO/FS>200' \
    -Oz -o filtered_indels.vcf.gz
```

### Common Variants Only

For population genetics (exclude rare variants):

```bash
bcftools filter -i 'INFO/AF>=0.05 && INFO/AF<=0.95' input.vcf.gz
```

### Rare Variants Only

For disease studies:

```bash
bcftools filter -i 'INFO/AF<0.01' input.vcf.gz
```

## Multi-Step Filtering Pipeline

### Sequential Filters

```bash
bcftools filter -e 'QUAL<30' input.vcf.gz | \
    bcftools filter -e 'INFO/DP<10' | \
    bcftools view -v snps | \
    bcftools view -f PASS \
    -Oz -o final.vcf.gz

bcftools index final.vcf.gz
```

### Named Soft Filters

Apply multiple named filters, then decide which to enforce:

```bash
# Apply soft filters
bcftools filter -s 'LowQual' -e 'QUAL<30' input.vcf.gz | \
    bcftools filter -s 'LowDepth' -e 'INFO/DP<10' | \
    bcftools filter -s 'LowMQ' -e 'INFO/MQ<40' \
    -Oz -o annotated.vcf.gz

# Check filter distribution
bcftools query -f '%FILTER\n' annotated.vcf.gz | sort | uniq -c

# Apply specific filters
bcftools view -f PASS,LowMQ annotated.vcf.gz  # Keep PASS and LowMQ
bcftools view -f PASS annotated.vcf.gz         # PASS only
```

## Python Filtering with cyvcf2

### Basic Quality Filter

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

def passes_filters(variant, min_qual=30, min_dp=10, min_mq=40):
    qual = variant.QUAL or 0
    dp = variant.INFO.get('DP', 0)
    mq = variant.INFO.get('MQ', 0)
    return qual >= min_qual and dp >= min_dp and mq >= min_mq

vcf = VCF('input.vcf.gz')
writer = Writer('filtered.vcf', vcf)

for variant in vcf:
    if passes_filters(variant):
        writer.write_record(variant)

writer.close()
vcf.close()
```

### Filter SNPs Only

```python
from cyvcf2 import VCF, Writer

vcf = VCF('input.vcf.gz')
writer = Writer('snps.vcf', vcf)

for variant in vcf:
    if variant.is_snp and (variant.QUAL or 0) >= 30:
        writer.write_record(variant)

writer.close()
vcf.close()
```

### Count Filtered Variants

```python
from cyvcf2 import VCF

total = 0
passed = 0

for variant in VCF('input.vcf.gz'):
    total += 1
    if (variant.QUAL or 0) >= 30:
        passed += 1

print(f'Total: {total}')
print(f'Passed: {passed}')
print(f'Filtered: {total - passed}')
print(f'Pass rate: {passed/total*100:.1f}%')
```

## Evaluating Filter Impact

### Before and After Counts

```bash
# Before filtering
bcftools view -H input.vcf.gz | wc -l

# After filtering
bcftools filter -e 'QUAL<30' input.vcf.gz | bcftools view -H | wc -l
```

### Compare Statistics

```bash
# Original stats
bcftools stats input.vcf.gz > original.stats

# Filtered stats
bcftools filter -e 'QUAL<30' input.vcf.gz | bcftools stats > filtered.stats

# Compare
diff original.stats filtered.stats
```

## Troubleshooting

### "no such INFO tag"

The tag doesn't exist in the VCF. Check available tags:

```bash
bcftools view -h input.vcf.gz | grep "^##INFO"
```

### "syntax error in expression"

Common issues:
- Use `||` not `or`
- Use `&&` not `and`
- Quote strings: `TYPE="snp"` not `TYPE=snp`

### Empty Output

Filter too strict. Check how many variants pass:

```bash
bcftools filter -i 'QUAL>=30' input.vcf.gz | bcftools view -H | wc -l
```

### Filter Not Working

Make sure the field has values:

```bash
bcftools query -f '%INFO/DP\n' input.vcf.gz | head
```

## See Also

- [bcftools filter documentation](http://www.htslib.org/doc/bcftools.html#filter)
- [bcftools view documentation](http://www.htslib.org/doc/bcftools.html#view)
- **vcf-basics** - Understanding VCF format
- **vcf-statistics** - Evaluating filter effects
