---
name: bio-variant-annotation
description: Add annotations to VCF files and predict variant consequences using bcftools annotate and csq. Use when adding database information, custom annotations, or determining functional impact.
tool_type: cli
primary_tool: bcftools
---

# Variant Annotation

Add annotations and predict consequences using bcftools.

## Annotation Types

| Tool | Purpose |
|------|---------|
| `bcftools annotate` | Add/remove INFO and FORMAT fields |
| `bcftools csq` | Predict functional consequences |

## bcftools annotate

### Add Annotations from Database

```bash
bcftools annotate -a dbsnp.vcf.gz -c ID input.vcf.gz -Oz -o annotated.vcf.gz
```

### Annotation Columns (`-c`)

| Option | Description |
|--------|-------------|
| `ID` | Copy ID column |
| `QUAL` | Copy QUAL column |
| `FILTER` | Copy FILTER column |
| `INFO` | Copy all INFO fields |
| `INFO/TAG` | Copy specific INFO field |
| `+INFO/TAG` | Add to existing values |

### Add rsIDs from dbSNP

```bash
bcftools annotate -a dbsnp.vcf.gz -c ID input.vcf.gz -Oz -o with_rsids.vcf.gz
```

### Add Multiple Annotations

```bash
bcftools annotate -a database.vcf.gz -c ID,INFO/AF,INFO/CAF input.vcf.gz -Oz -o annotated.vcf.gz
```

### Add INFO Field from BED

```bash
# BED with 4th column as annotation
bcftools annotate -a regions.bed.gz -c CHROM,FROM,TO,INFO/REGION \
    -h <(echo '##INFO=<ID=REGION,Number=1,Type=String,Description="Region name">') \
    input.vcf.gz -Oz -o annotated.vcf.gz
```

### Add INFO Field from TAB File

```bash
# Tab file: CHROM POS VALUE
bcftools annotate -a annotations.tab.gz -c CHROM,POS,INFO/SCORE \
    -h <(echo '##INFO=<ID=SCORE,Number=1,Type=Float,Description="Custom score">') \
    input.vcf.gz -Oz -o annotated.vcf.gz
```

## Removing Annotations

### Remove INFO Fields

```bash
bcftools annotate -x INFO/DP,INFO/MQ input.vcf.gz -Oz -o clean.vcf.gz
```

### Remove All INFO Fields

```bash
bcftools annotate -x INFO input.vcf.gz -Oz -o minimal.vcf.gz
```

### Remove FORMAT Fields

```bash
bcftools annotate -x FORMAT/AD,FORMAT/DP input.vcf.gz -Oz -o clean.vcf.gz
```

### Remove ID Column

```bash
bcftools annotate -x ID input.vcf.gz -Oz -o no_ids.vcf.gz
```

### Keep Only Specific Fields

```bash
# Remove all INFO except DP
bcftools annotate -x ^INFO/DP input.vcf.gz -Oz -o minimal.vcf.gz
```

## Setting Annotations

### Set ID from Fields

```bash
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' input.vcf.gz -Oz -o with_ids.vcf.gz
```

### Set ID with rsID Fallback

```bash
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' input.vcf.gz -Oz -o with_ids.vcf.gz
# '+' means only set if ID is missing
```

## Header Modifications

### Add Header Line

```bash
bcftools annotate -h new_headers.txt input.vcf.gz -Oz -o annotated.vcf.gz
```

Header file format:
```
##INFO=<ID=CUSTOM,Number=1,Type=Integer,Description="Custom field">
##FILTER=<ID=MyFilter,Description="Custom filter">
```

### Rename Chromosomes

```bash
# rename.txt: old_name new_name
bcftools annotate --rename-chrs rename.txt input.vcf.gz -Oz -o renamed.vcf.gz
```

Example rename.txt:
```
1	chr1
2	chr2
X	chrX
```

## bcftools csq

Predict functional consequences using GFF annotation.

### Basic Consequence Prediction

```bash
bcftools csq -f reference.fa -g genes.gff3.gz input.vcf.gz -Oz -o consequences.vcf.gz
```

### Output Format

Adds `BCSQ` INFO field with format:
```
consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change
```

### Consequence Types

| Consequence | Description |
|-------------|-------------|
| `synonymous` | No amino acid change |
| `missense` | Amino acid change |
| `stop_gained` | Introduces stop codon |
| `stop_lost` | Removes stop codon |
| `frameshift` | Changes reading frame |
| `splice_donor` | Affects splice donor |
| `splice_acceptor` | Affects splice acceptor |
| `intron` | Intronic variant |
| `intergenic` | Between genes |
| `5_prime_utr` | 5' UTR variant |
| `3_prime_utr` | 3' UTR variant |

### Phase-Aware Predictions

```bash
bcftools csq -f reference.fa -g genes.gff3.gz -p s input.vcf.gz -Oz -o phased.vcf.gz
```

Phase options (`-p`):
- `a` - Take first allele, ignore phase
- `m` - Merge haplotypes
- `r` - Require all sites phased
- `R` - Require at least one site phased
- `s` - Skip unphased sites

### Local Phasing Window

```bash
bcftools csq -f reference.fa -g genes.gff3.gz -l input.vcf.gz -Oz -o local_phased.vcf.gz
```

## Extracting Consequences

### Query BCSQ Field

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' consequences.vcf.gz
```

### Filter by Consequence

```bash
# Keep only missense variants
bcftools view -i 'INFO/BCSQ~"missense"' consequences.vcf.gz -Oz -o missense.vcf.gz
```

### Split BCSQ for Analysis

```bash
bcftools query -f '%CHROM\t%POS\t%INFO/BCSQ\n' consequences.vcf.gz | \
    awk -F'|' '{print $1"\t"$2"\t"$3}'
```

## Annotation Databases

### Common Databases

| Database | Content | Format |
|----------|---------|--------|
| dbSNP | rsIDs, frequencies | VCF |
| gnomAD | Population frequencies | VCF |
| ClinVar | Clinical significance | VCF |
| COSMIC | Cancer mutations | VCF |

### Download and Prepare

```bash
# dbSNP (example)
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all.vcf.gz
tabix -p vcf common_all.vcf.gz
```

### Annotate with gnomAD

```bash
bcftools annotate -a gnomad.vcf.gz \
    -c INFO/AF,INFO/AF_popmax \
    input.vcf.gz -Oz -o annotated.vcf.gz
```

### Annotate with ClinVar

```bash
bcftools annotate -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNDN \
    input.vcf.gz -Oz -o annotated.vcf.gz
```

## Common Workflows

### Add rsIDs and Filter

```bash
bcftools annotate -a dbsnp.vcf.gz -c ID input.vcf.gz | \
    bcftools view -i 'ID!="."' -Oz -o known_variants.vcf.gz
```

### Annotate with Population Frequency

```bash
bcftools annotate -a gnomad.vcf.gz -c INFO/AF input.vcf.gz | \
    bcftools filter -i 'INFO/AF<0.01' -Oz -o rare.vcf.gz
```

### Full Annotation Pipeline

```bash
# Normalize first
bcftools norm -f reference.fa -m-any input.vcf.gz | \
# Add rsIDs
bcftools annotate -a dbsnp.vcf.gz -c ID | \
# Add population frequencies
bcftools annotate -a gnomad.vcf.gz -c INFO/AF | \
# Predict consequences
bcftools csq -f reference.fa -g genes.gff3.gz -Oz -o annotated.vcf.gz
```

## cyvcf2 Annotation Access

### Read Annotations

```python
from cyvcf2 import VCF

for variant in VCF('annotated.vcf.gz'):
    rsid = variant.ID
    af = variant.INFO.get('AF')
    bcsq = variant.INFO.get('BCSQ')

    if bcsq:
        parts = bcsq.split('|')
        consequence = parts[0]
        gene = parts[1] if len(parts) > 1 else None
        print(f'{variant.CHROM}:{variant.POS} {consequence} {gene}')
```

### Filter by Annotation

```python
from cyvcf2 import VCF, Writer

vcf = VCF('annotated.vcf.gz')
writer = Writer('rare.vcf', vcf)

for variant in vcf:
    af = variant.INFO.get('AF', 1.0)
    if af < 0.01:
        writer.write_record(variant)

writer.close()
vcf.close()
```

## Quick Reference

| Task | Command |
|------|---------|
| Add rsIDs | `bcftools annotate -a dbsnp.vcf.gz -c ID in.vcf.gz` |
| Add INFO field | `bcftools annotate -a db.vcf.gz -c INFO/AF in.vcf.gz` |
| Remove field | `bcftools annotate -x INFO/DP in.vcf.gz` |
| Set ID | `bcftools annotate --set-id '%CHROM\_%POS' in.vcf.gz` |
| Rename chroms | `bcftools annotate --rename-chrs map.txt in.vcf.gz` |
| Consequences | `bcftools csq -f ref.fa -g genes.gff in.vcf.gz` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `no such tag` | INFO field doesn't exist | Check with `bcftools view -h` |
| `chromosome not found` | Annotation file has different chroms | Use `--rename-chrs` |
| `index required` | Annotation file not indexed | Run `tabix -p vcf file.vcf.gz` |

## Related Skills

- **variant-normalization** - Normalize before annotating
- **vcf-filtering** - Filter by annotations
- **vcf-basics** - Query annotated fields
- **database-access/entrez-fetch** - Download annotation databases
