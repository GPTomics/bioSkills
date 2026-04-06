---
name: bio-vcf-basics
description: View, query, and understand VCF/BCF variant files using bcftools and cyvcf2. Use when inspecting variants, extracting specific fields, or understanding VCF format structure.
tool_type: cli
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# VCF/BCF Basics

View and query variant files using bcftools and cyvcf2.

## Format Overview

| Format | Description | Use Case |
|--------|-------------|----------|
| VCF | Text format, human-readable | Debugging, small files |
| VCF.gz | Compressed VCF (bgzip) | Standard distribution |
| BCF | Binary VCF | Fast processing, large files |

## VCF Format Structure

```
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1
chr1    1000    rs123   A       G       30      PASS    DP=50   GT:DP   0/1:25
```

### Header Lines (##)
- `##fileformat` - VCF version
- `##INFO` - INFO field definitions
- `##FORMAT` - FORMAT field definitions
- `##FILTER` - Filter definitions
- `##contig` - Reference contigs
- `##reference` - Reference genome

### Column Header (#CHROM)
Fixed columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
Followed by sample columns

### Data Columns

| Column | Description |
|--------|-------------|
| CHROM | Chromosome |
| POS | 1-based position of the first base in REF |
| ID | Variant identifier (e.g., rs number) or `.` if novel |
| REF | Reference allele |
| ALT | Alternate allele(s), comma-separated. `*` indicates an overlapping deletion at this site |
| QUAL | Phred-scaled probability that a variant exists at this site (site-level, NOT per-sample) |
| FILTER | PASS or semicolon-separated filter names. `.` means filters were not applied |
| INFO | Semicolon-separated key=value pairs (site-level annotations) |
| FORMAT | Colon-separated format keys defining per-sample field order |
| SAMPLE | Colon-separated values matching FORMAT order |

## Critical Field Interpretation

Understanding what each field actually measures -- and what it does not -- is essential for filtering and interpretation decisions.

### QUAL vs GQ: Different Questions

| Metric | Scope | Question Answered | When Most Useful |
|--------|-------|-------------------|------------------|
| QUAL | Site-level | "Is there any variant here at all?" | Multi-sample calling where site confidence matters |
| GQ | Sample-level | "Is this specific genotype assignment correct?" | Per-sample genotype confidence |
| QD | Site-level | QUAL normalized by allele depth | Preferred over raw QUAL for filtering (depth-independent) |

QUAL can be high even when an individual sample's genotype is uncertain. Conversely, a sample may have a confident genotype (high GQ) at a site with moderate QUAL.

### AD vs DP Discrepancy

The sum of AD (allelic depth) values is often less than DP (total depth). This is expected behavior, not an error:
- **DP** counts all reads spanning the position, including uninformative reads (low base quality, ambiguous alignment)
- **AD** counts only reads that confidently support a specific allele
- INFO/DP (site-level) differs from FORMAT/DP (per-sample) -- site DP is the sum across all samples

### Key INFO Annotations for Filtering

| Annotation | Meaning | What It Detects |
|-----------|---------|-----------------|
| QD | QUAL / allele depth | Low values suggest variant quality not supported by reads |
| FS | Fisher strand bias (phred-scaled) | Variant reads predominantly on one strand (artifact) |
| SOR | Strand odds ratio | Same as FS but handles high-depth sites better |
| MQ | Root mean square mapping quality | Low values indicate reads map ambiguously (paralogous regions) |
| MQRankSum | MQ difference: ref vs alt reads | Very negative = alt reads map much worse than ref (suspicious) |
| ReadPosRankSum | Read position: ref vs alt reads | Very negative = variant only at read ends (misalignment artifact) |

### PL (Phred-Scaled Likelihoods)

PL encodes the relative likelihood of each possible genotype. For a biallelic site: `PL = [P(0/0), P(0/1), P(1/1)]`. The most likely genotype always has PL=0; GQ equals the difference between the lowest and second-lowest PL values.

For multiallelic sites with n alleles, PL contains n*(n+1)/2 values covering all possible diploid genotypes.

## bcftools view

**Goal:** View, subset, and convert VCF/BCF files from the command line.

**Approach:** Use bcftools view with flags for header control, region selection, sample extraction, and format conversion.

**"Show me what's in this VCF file"** → Display VCF contents with optional filtering by header, region, or sample.

### View VCF
```bash
bcftools view input.vcf.gz | head
```

### View Header Only
```bash
bcftools view -h input.vcf.gz
```

### View Without Header
```bash
bcftools view -H input.vcf.gz | head
```

### View Specific Region
```bash
bcftools view input.vcf.gz chr1:1000000-2000000
```

### View Specific Samples
```bash
bcftools view -s sample1,sample2 input.vcf.gz
```

### Exclude Samples
```bash
bcftools view -s ^sample3 input.vcf.gz
```

## bcftools query

**Goal:** Extract specific fields from a VCF in a custom tabular format.

**Approach:** Use bcftools query with format specifiers for CHROM, POS, INFO, and FORMAT fields.

**"Extract positions and genotypes from my VCF"** → Pull specific columns from variant records into a flat text format.

Extract specific fields in custom format.

### Basic Query
```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' input.vcf.gz
```

### Query with INFO Fields
```bash
bcftools query -f '%CHROM\t%POS\t%INFO/DP\t%INFO/AF\n' input.vcf.gz
```

### Query with Sample Fields
```bash
bcftools query -f '%CHROM\t%POS[\t%GT]\n' input.vcf.gz
```

### Query Specific Samples
```bash
bcftools query -f '%CHROM\t%POS[\t%SAMPLE=%GT]\n' -s sample1,sample2 input.vcf.gz
```

### Include Header
```bash
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\n' input.vcf.gz
```

### Common Format Specifiers

| Specifier | Description |
|-----------|-------------|
| `%CHROM` | Chromosome |
| `%POS` | Position |
| `%ID` | Variant ID |
| `%REF` | Reference allele |
| `%ALT` | Alternate allele |
| `%QUAL` | Quality score |
| `%FILTER` | Filter status |
| `%INFO/TAG` | INFO field value |
| `%TYPE` | Variant type (snp, indel, etc.) |
| `[%GT]` | Genotype (per sample) |
| `[%DP]` | Depth (per sample) |
| `[%SAMPLE]` | Sample name |
| `\n` | Newline |
| `\t` | Tab |

## Format Conversion

**Goal:** Convert between VCF, compressed VCF, and BCF formats.

**Approach:** Use bcftools view with output format flags (-Ov, -Oz, -Ob) and bgzip/index for compression and indexing.

### VCF to BCF
```bash
bcftools view -Ob -o output.bcf input.vcf.gz
```

### BCF to VCF
```bash
bcftools view -Ov -o output.vcf input.bcf
```

### Compress VCF (bgzip)
```bash
bgzip input.vcf
# Creates input.vcf.gz
```

### Index VCF/BCF
```bash
bcftools index input.vcf.gz
# Creates input.vcf.gz.csi

bcftools index -t input.vcf.gz
# Creates input.vcf.gz.tbi (tabix index)
```

## Output Format Options

| Flag | Format |
|------|--------|
| `-Ov` | Uncompressed VCF |
| `-Oz` | Compressed VCF (bgzip) |
| `-Ou` | Uncompressed BCF |
| `-Ob` | Compressed BCF |

## Genotype Encoding

| Genotype | Meaning |
|----------|---------|
| `0/0` | Homozygous reference |
| `0/1` | Heterozygous |
| `1/1` | Homozygous alternate |
| `1/2` | Heterozygous for two different ALT alleles (compound het at multiallelic site) |
| `./.` | Missing genotype (no confident call) |
| `0\|1` | Phased heterozygous (allele before `\|` is on haplotype 1) |

### Phased vs Unphased

- `/` separates **unphased** alleles -- the two chromosomal copies are known, but which came from which parent is not
- `|` separates **phased** alleles -- haplotype assignment is known (e.g., from read-backed phasing, trio analysis, or long-read sequencing)
- Phasing matters for compound heterozygosity: two variants on the same gene are pathogenic only if they are on different haplotypes (in *trans*), not the same haplotype (in *cis*)

### Multiallelic Genotypes

At multiallelic sites (e.g., ALT = G,T), allele indices reference the comma-separated ALT list: 0=REF, 1=first ALT, 2=second ALT. Genotype `1/2` means one copy of each ALT allele. Splitting multiallelic sites into biallelic records with `bcftools norm -m-` converts `1/2` into two `0/1` records, losing compound heterozygosity information -- see variant-normalization for caveats.

## cyvcf2 Python Alternative

**Goal:** Read, query, and write VCF files programmatically in Python.

**Approach:** Use cyvcf2's VCF reader to iterate variants, access properties/INFO/FORMAT fields, and write filtered output with Writer.

**"Parse this VCF in Python"** → Open VCF with cyvcf2 and iterate variant records with attribute-style access to fields.

### Open and Iterate
```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
for variant in vcf:
    print(f'{variant.CHROM}:{variant.POS} {variant.REF}>{variant.ALT[0]}')
```

### Access Variant Properties
```python
from cyvcf2 import VCF

for variant in VCF('input.vcf.gz'):
    print(f'Chrom: {variant.CHROM}')
    print(f'Pos: {variant.POS}')
    print(f'ID: {variant.ID}')
    print(f'Ref: {variant.REF}')
    print(f'Alt: {variant.ALT}')  # List
    print(f'Qual: {variant.QUAL}')
    print(f'Filter: {variant.FILTER}')
    print(f'Type: {variant.var_type}')  # snp, indel, etc.
    break
```

### Access INFO Fields
```python
from cyvcf2 import VCF

for variant in VCF('input.vcf.gz'):
    dp = variant.INFO.get('DP')
    af = variant.INFO.get('AF')
    print(f'{variant.CHROM}:{variant.POS} DP={dp} AF={af}')
```

### Access Genotypes
```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
samples = vcf.samples  # List of sample names

for variant in vcf:
    gts = variant.gt_types  # 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
    for sample, gt in zip(samples, gts):
        gt_str = ['HOM_REF', 'HET', 'UNKNOWN', 'HOM_ALT'][gt]
        print(f'{sample}: {gt_str}')
    break
```

### Access Sample Fields
```python
from cyvcf2 import VCF

for variant in VCF('input.vcf.gz'):
    depths = variant.format('DP')  # numpy array
    gqs = variant.format('GQ')     # Genotype quality
    print(f'Depths: {depths}')
```

### Fetch Region
```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
for variant in vcf('chr1:1000000-2000000'):
    print(f'{variant.CHROM}:{variant.POS}')
```

### Get Header Info
```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
print(f'Samples: {vcf.samples}')
print(f'Contigs: {vcf.seqnames}')

# INFO field definitions
for info in vcf.header_iter():
    if info['HeaderType'] == 'INFO':
        print(f'{info["ID"]}: {info["Description"]}')
```

### Write VCF
```python
from cyvcf2 import VCF, Writer

vcf = VCF('input.vcf.gz')
writer = Writer('output.vcf', vcf)

for variant in vcf:
    if variant.QUAL > 30:
        writer.write_record(variant)

writer.close()
vcf.close()
```

## Quick Reference

| Task | bcftools | cyvcf2 |
|------|----------|--------|
| View VCF | `bcftools view file.vcf.gz` | `VCF('file.vcf.gz')` |
| View header | `bcftools view -h file.vcf.gz` | `vcf.header_iter()` |
| Get region | `bcftools view file.vcf.gz chr1:1-1000` | `vcf('chr1:1-1000')` |
| Query fields | `bcftools query -f '%CHROM\t%POS\n'` | Loop with properties |
| Count variants | `bcftools view -H file.vcf.gz \| wc -l` | `sum(1 for _ in vcf)` |
| VCF to BCF | `bcftools view -Ob -o out.bcf in.vcf.gz` | Use Writer |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `no BGZF EOF marker` | Not bgzipped | Use `bgzip` not `gzip` |
| `index required` | Missing index for region query | Run `bcftools index` |
| `sample not found` | Wrong sample name | Check with `bcftools query -l` |

## Related Skills

- variant-calling - Generate VCF from alignments
- filtering-best-practices - Filter variants by quality/criteria
- vcf-manipulation - Merge, concat, intersect VCFs
- alignment-files/pileup-generation - Generate pileup for calling
