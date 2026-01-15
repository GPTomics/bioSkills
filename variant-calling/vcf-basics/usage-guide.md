# VCF/BCF Basics Usage Guide

This guide covers viewing, querying, and understanding variant files.

## Prerequisites

- bcftools installed (`conda install -c bioconda bcftools`)
- cyvcf2 installed (`pip install cyvcf2`)
- htslib for bgzip/tabix (`conda install -c bioconda htslib`)

## Understanding VCF Format

VCF (Variant Call Format) is the standard format for storing genetic variants.

### File Structure

A VCF file has three sections:

1. **Meta-information lines** (##) - Describe the file and define fields
2. **Header line** (#CHROM) - Column names
3. **Data lines** - One variant per line

### Example VCF

```
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO            FORMAT      SAMPLE1     SAMPLE2
chr1    1000    rs123   A       G       50      PASS    DP=100;AF=0.5   GT:DP:GQ    0/1:50:99   0/0:50:99
chr1    2000    .       CT      C       30      PASS    DP=80;AF=0.25   GT:DP:GQ    0/1:40:80   0/1:40:85
```

### Fixed Columns

| Column | Description | Example |
|--------|-------------|---------|
| CHROM | Chromosome | chr1 |
| POS | 1-based position | 1000 |
| ID | Variant ID or `.` | rs123 |
| REF | Reference allele | A |
| ALT | Alternate allele(s) | G |
| QUAL | Quality score | 50 |
| FILTER | PASS or filter name | PASS |
| INFO | Variant annotations | DP=100;AF=0.5 |
| FORMAT | Sample field format | GT:DP:GQ |

### INFO Field

Semicolon-separated key=value pairs describing the variant:

```
DP=100;AF=0.5;MQ=60
```

Common INFO fields:
- `DP` - Total read depth
- `AF` - Allele frequency
- `MQ` - Mapping quality
- `AN` - Total alleles
- `AC` - Allele count

### FORMAT and Sample Fields

FORMAT defines the order of sample values. Sample columns contain colon-separated values:

```
FORMAT      SAMPLE1
GT:DP:GQ    0/1:50:99
```

Common FORMAT fields:
- `GT` - Genotype
- `DP` - Read depth
- `GQ` - Genotype quality
- `AD` - Allelic depths
- `PL` - Phred-scaled likelihoods

## Viewing VCF Files

### Basic Viewing

```bash
# View entire file
bcftools view input.vcf.gz | less

# View first 20 variants
bcftools view -H input.vcf.gz | head -20

# View header only
bcftools view -h input.vcf.gz
```

### View Specific Region

```bash
# Region query (requires index)
bcftools view input.vcf.gz chr1:1000000-2000000

# Multiple regions
bcftools view -r chr1:1000-2000,chr2:3000-4000 input.vcf.gz
```

### View Specific Samples

```bash
# Include samples
bcftools view -s sample1,sample2 input.vcf.gz

# Exclude samples
bcftools view -s ^sample3,sample4 input.vcf.gz

# List samples
bcftools query -l input.vcf.gz
```

## Querying VCF Files

### Basic Query

```bash
# Extract CHROM, POS, REF, ALT
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' input.vcf.gz
```

### Query with INFO Fields

```bash
# Include depth and allele frequency
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\n' input.vcf.gz
```

### Query Sample Data

```bash
# Get genotypes for all samples
bcftools query -f '%CHROM\t%POS[\t%GT]\n' input.vcf.gz

# Get sample name with genotype
bcftools query -f '%CHROM\t%POS[\t%SAMPLE=%GT]\n' input.vcf.gz

# Get depth per sample
bcftools query -f '%CHROM\t%POS[\t%DP]\n' input.vcf.gz
```

### Query with Header

```bash
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\n' input.vcf.gz
```

### Query to TSV

```bash
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP[\t%GT]\n' \
    input.vcf.gz > variants.tsv
```

## Format Conversion

### Compression

VCF files should be compressed with bgzip (not gzip):

```bash
# Compress VCF
bgzip input.vcf

# Decompress
bgzip -d input.vcf.gz

# Keep original
bgzip -k input.vcf
```

### VCF to BCF

BCF (binary VCF) is faster to process:

```bash
bcftools view -Ob -o output.bcf input.vcf.gz
bcftools index output.bcf
```

### BCF to VCF

```bash
bcftools view -Oz -o output.vcf.gz input.bcf
```

## Indexing

Index files enable random access:

```bash
# CSI index (default, supports large chromosomes)
bcftools index input.vcf.gz

# TBI index (tabix format)
bcftools index -t input.vcf.gz

# Force reindex
bcftools index -f input.vcf.gz
```

## Working with cyvcf2

### Basic Iteration

```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')

for variant in vcf:
    print(f'{variant.CHROM}:{variant.POS} {variant.REF}>{",".join(variant.ALT)}')

vcf.close()
```

### Using Context Manager Pattern

```python
from cyvcf2 import VCF

def process_vcf(path):
    vcf = VCF(path)
    try:
        for variant in vcf:
            yield variant
    finally:
        vcf.close()

for v in process_vcf('input.vcf.gz'):
    print(v.CHROM, v.POS)
```

### Accessing Variant Properties

```python
from cyvcf2 import VCF

for variant in VCF('input.vcf.gz'):
    # Basic properties
    chrom = variant.CHROM
    pos = variant.POS
    ref = variant.REF
    alts = variant.ALT  # List of alternate alleles
    qual = variant.QUAL
    filter_status = variant.FILTER  # None if PASS

    # Variant type
    var_type = variant.var_type  # 'snp', 'indel', 'mnp', etc.
    is_snp = variant.is_snp
    is_indel = variant.is_indel

    # Allele info
    is_transition = variant.is_transition
    num_called = variant.num_called
    num_het = variant.num_het
    num_hom_ref = variant.num_hom_ref
    num_hom_alt = variant.num_hom_alt

    break
```

### Accessing INFO Fields

```python
from cyvcf2 import VCF

for variant in VCF('input.vcf.gz'):
    # Get INFO field (returns None if not present)
    dp = variant.INFO.get('DP')
    af = variant.INFO.get('AF')
    mq = variant.INFO.get('MQ')

    # Check if field exists
    if 'DP' in variant.INFO:
        print(f'Depth: {variant.INFO["DP"]}')
```

### Accessing Genotypes

```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
samples = vcf.samples

for variant in vcf:
    # gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
    gt_types = variant.gt_types

    # gt_bases: actual genotype strings like 'A/G'
    gt_bases = variant.gt_bases

    # gt_phases: phasing information
    gt_phases = variant.gt_phases

    for sample, gt_type, gt_base in zip(samples, gt_types, gt_bases):
        print(f'{sample}: {gt_base} (type={gt_type})')

    break
```

### Accessing FORMAT Fields

```python
from cyvcf2 import VCF

for variant in VCF('input.vcf.gz'):
    # Returns numpy array (samples x values)
    depths = variant.format('DP')
    gqs = variant.format('GQ')
    ads = variant.format('AD')  # Allelic depths

    if depths is not None:
        print(f'Depths: {depths.flatten()}')
```

### Region Queries

```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')

# Query specific region
for variant in vcf('chr1:1000000-2000000'):
    print(f'{variant.CHROM}:{variant.POS}')

# Multiple regions
for region in ['chr1:1000-2000', 'chr2:3000-4000']:
    for variant in vcf(region):
        print(variant.POS)
```

### Writing VCF

```python
from cyvcf2 import VCF, Writer

vcf = VCF('input.vcf.gz')
writer = Writer('output.vcf', vcf)  # Inherits header

for variant in vcf:
    # Filter and write
    if variant.QUAL and variant.QUAL > 30:
        writer.write_record(variant)

writer.close()
vcf.close()
```

## Common Tasks

### Count Variants

```bash
# Total variants
bcftools view -H input.vcf.gz | wc -l

# SNPs only
bcftools view -v snps input.vcf.gz | bcftools view -H | wc -l

# Indels only
bcftools view -v indels input.vcf.gz | bcftools view -H | wc -l
```

### List Samples

```bash
bcftools query -l input.vcf.gz
```

### Extract Variant IDs

```bash
bcftools query -f '%ID\n' input.vcf.gz | grep -v '^\.$'
```

### Get Allele Frequencies

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' input.vcf.gz
```

## Troubleshooting

### "no BGZF EOF marker"

File compressed with gzip instead of bgzip:

```bash
# Fix: recompress
gunzip input.vcf.gz
bgzip input.vcf
```

### "index required"

Missing index file:

```bash
bcftools index input.vcf.gz
```

### "sample not found"

Check available samples:

```bash
bcftools query -l input.vcf.gz
```

## See Also

- [VCF specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
- [bcftools documentation](http://www.htslib.org/doc/bcftools.html)
- [cyvcf2 documentation](https://brentp.github.io/cyvcf2/)
