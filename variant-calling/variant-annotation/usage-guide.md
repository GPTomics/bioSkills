# Variant Annotation Usage Guide

## Overview

This guide covers adding annotations to VCF files and predicting functional consequences.

## Prerequisites

- bcftools installed (`conda install -c bioconda bcftools`)
- Annotation databases (dbSNP, gnomAD, ClinVar, etc.)
- For consequence prediction: reference FASTA and GFF3 gene annotation

## Quick Start

Tell your AI agent what you want to do:
- "Add rsIDs from dbSNP to my VCF file"
- "Annotate my variants with allele frequencies from gnomAD"
- "Predict functional consequences of my variants using bcftools csq"
- "Add ClinVar clinical significance annotations to my VCF"

## bcftools annotate Overview

`bcftools annotate` adds or removes annotations from VCF files. Common uses:

- Add rsIDs from dbSNP
- Add population frequencies from gnomAD
- Add clinical significance from ClinVar
- Remove unwanted fields to reduce file size
- Rename chromosomes

## Adding Annotations from VCF

### Add rsIDs from dbSNP

```bash
bcftools annotate -a dbsnp.vcf.gz -c ID input.vcf.gz -Oz -o with_rsids.vcf.gz
bcftools index with_rsids.vcf.gz
```

The `-c ID` tells bcftools to copy the ID column from the annotation file.

### Add INFO Fields

```bash
# Add allele frequency
bcftools annotate -a gnomad.vcf.gz -c INFO/AF input.vcf.gz -Oz -o with_af.vcf.gz

# Add multiple fields
bcftools annotate -a gnomad.vcf.gz -c INFO/AF,INFO/AC,INFO/AN input.vcf.gz -Oz -o annotated.vcf.gz
```

### Add All INFO Fields

```bash
bcftools annotate -a database.vcf.gz -c INFO input.vcf.gz -Oz -o annotated.vcf.gz
```

### Column Specification Options

| Option | Description |
|--------|-------------|
| `ID` | Variant identifier |
| `QUAL` | Quality score |
| `FILTER` | Filter status |
| `INFO` | All INFO fields |
| `INFO/TAG` | Specific INFO field |
| `+INFO/TAG` | Append to existing field |
| `-INFO/TAG` | Remove field |

## Adding Annotations from BED/TAB Files

### BED File Annotation

```bash
# regions.bed.gz has 4 columns: CHROM, START, END, REGION_NAME

# First add header definition
echo '##INFO=<ID=REGION,Number=1,Type=String,Description="Genomic region">' > header.txt

bcftools annotate \
    -a regions.bed.gz \
    -c CHROM,FROM,TO,INFO/REGION \
    -h header.txt \
    input.vcf.gz -Oz -o annotated.vcf.gz
```

### TAB File Annotation

```bash
# scores.tab.gz has columns: CHROM, POS, SCORE

echo '##INFO=<ID=SCORE,Number=1,Type=Float,Description="Conservation score">' > header.txt

bcftools annotate \
    -a scores.tab.gz \
    -c CHROM,POS,INFO/SCORE \
    -h header.txt \
    input.vcf.gz -Oz -o annotated.vcf.gz
```

## Removing Annotations

### Remove Specific Fields

```bash
# Remove INFO fields
bcftools annotate -x INFO/DP,INFO/MQ input.vcf.gz -Oz -o clean.vcf.gz

# Remove FORMAT fields
bcftools annotate -x FORMAT/AD,FORMAT/PL input.vcf.gz -Oz -o clean.vcf.gz

# Remove ID column
bcftools annotate -x ID input.vcf.gz -Oz -o no_ids.vcf.gz
```

### Remove All of a Type

```bash
# Remove all INFO fields
bcftools annotate -x INFO input.vcf.gz -Oz -o minimal.vcf.gz

# Remove all FORMAT fields except GT
bcftools annotate -x FORMAT input.vcf.gz -Oz -o gt_only.vcf.gz
```

### Keep Only Specific Fields

```bash
# Remove all INFO except DP and AF
bcftools annotate -x ^INFO/DP,INFO/AF input.vcf.gz -Oz -o minimal.vcf.gz
```

## Setting Variant IDs

### Create IDs from Fields

```bash
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' input.vcf.gz -Oz -o with_ids.vcf.gz
```

Example output: `chr1_12345_A_G`

### Only Set Missing IDs

```bash
# '+' prefix means only set if ID is currently '.'
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' input.vcf.gz -Oz -o with_ids.vcf.gz
```

### ID Format Specifiers

| Specifier | Value |
|-----------|-------|
| `%CHROM` | Chromosome |
| `%POS` | Position |
| `%REF` | Reference allele |
| `%ALT` | Alternate allele |
| `%FIRST_ALT` | First alternate allele |

## Chromosome Renaming

### UCSC to Ensembl Style

```bash
# rename.txt: old_name<tab>new_name
cat > rename.txt << EOF
chr1	1
chr2	2
chr3	3
chrX	X
chrY	Y
chrM	MT
EOF

bcftools annotate --rename-chrs rename.txt input.vcf.gz -Oz -o renamed.vcf.gz
```

### Ensembl to UCSC Style

```bash
cat > rename.txt << EOF
1	chr1
2	chr2
3	chr3
X	chrX
Y	chrY
MT	chrM
EOF

bcftools annotate --rename-chrs rename.txt input.vcf.gz -Oz -o renamed.vcf.gz
```

## Consequence Prediction with bcftools csq

### Requirements

1. Reference FASTA
2. GFF3 gene annotation (from Ensembl, NCBI, etc.)

### Basic Usage

```bash
bcftools csq -f reference.fa -g genes.gff3.gz input.vcf.gz -Oz -o consequences.vcf.gz
```

### Output Format

Adds `BCSQ` INFO field:
```
BCSQ=missense|GENE1|ENST00000123456|protein_coding|+|123A>V|c.368C>T
```

Fields are pipe-separated:
1. Consequence type
2. Gene name
3. Transcript ID
4. Biotype
5. Strand
6. Amino acid change
7. DNA change

### Common Consequences

| Consequence | Impact |
|-------------|--------|
| `synonymous` | Low - same amino acid |
| `missense` | Moderate - different amino acid |
| `stop_gained` | High - premature stop |
| `frameshift` | High - reading frame shift |
| `splice_donor` | High - affects splicing |
| `splice_acceptor` | High - affects splicing |
| `intron` | Low - intronic |
| `intergenic` | Low - between genes |

### Query Consequences

```bash
# List all consequences
bcftools query -f '%CHROM\t%POS\t%INFO/BCSQ\n' consequences.vcf.gz | head

# Extract just consequence and gene
bcftools query -f '%INFO/BCSQ\n' consequences.vcf.gz | cut -d'|' -f1,2 | sort | uniq -c
```

### Filter by Consequence

```bash
# Keep high-impact variants
bcftools view -i 'INFO/BCSQ~"stop_gained" || INFO/BCSQ~"frameshift"' \
    consequences.vcf.gz -Oz -o high_impact.vcf.gz
```

## Annotation Database Setup

### dbSNP

```bash
# Download
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

# Annotate
bcftools annotate -a GCF_000001405.40.gz -c ID input.vcf.gz -Oz -o with_rsids.vcf.gz
```

### gnomAD

gnomAD files are large. Download only needed chromosomes:

```bash
# Download chromosome 1
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr1.vcf.bgz
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr1.vcf.bgz.tbi

# Annotate
bcftools annotate -a gnomad.genomes.v4.0.sites.chr1.vcf.bgz \
    -c INFO/AF,INFO/AF_popmax \
    -r chr1 \
    input.vcf.gz -Oz -o with_gnomad.vcf.gz
```

### ClinVar

```bash
# Download
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# Annotate
bcftools annotate -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNDN \
    input.vcf.gz -Oz -o with_clinvar.vcf.gz
```

## Complete Annotation Workflows

### Basic Pipeline

```bash
# 1. Add rsIDs
bcftools annotate -a dbsnp.vcf.gz -c ID input.vcf.gz | \
# 2. Add population frequencies
bcftools annotate -a gnomad.vcf.gz -c INFO/AF | \
# 3. Output
bgzip -c > annotated.vcf.gz

bcftools index annotated.vcf.gz
```

### Clinical Variant Analysis

```bash
# Normalize
bcftools norm -f reference.fa -m-any input.vcf.gz | \
# Add clinical significance
bcftools annotate -a clinvar.vcf.gz -c INFO/CLNSIG,INFO/CLNDN | \
# Predict consequences
bcftools csq -f reference.fa -g genes.gff3.gz | \
# Filter pathogenic
bcftools view -i 'INFO/CLNSIG~"Pathogenic"' -Oz -o pathogenic.vcf.gz
```

### Rare Variant Analysis

```bash
# Annotate with population frequency
bcftools annotate -a gnomad.vcf.gz -c INFO/AF input.vcf.gz | \
# Filter rare variants (AF < 1%)
bcftools filter -i 'INFO/AF<0.01 || INFO/AF="."' | \
# Predict consequences
bcftools csq -f reference.fa -g genes.gff3.gz -Oz -o rare_consequences.vcf.gz
```

## Troubleshooting

### "no such tag"

The annotation file doesn't have the requested field:

```bash
# Check available fields
bcftools view -h annotation.vcf.gz | grep "^##INFO"
```

### "chromosome not found"

Chromosome names don't match between files:

```bash
# Check chromosomes in each file
bcftools view -H input.vcf.gz | cut -f1 | sort -u | head
bcftools view -H annotation.vcf.gz | cut -f1 | sort -u | head

# Rename if needed
bcftools annotate --rename-chrs rename.txt input.vcf.gz
```

### "index required"

Annotation file needs an index:

```bash
tabix -p vcf annotation.vcf.gz
```

### Slow Annotation

For large files, restrict to regions:

```bash
bcftools annotate -a database.vcf.gz -c INFO/AF \
    -r chr1,chr2,chr3 \
    input.vcf.gz -Oz -o annotated.vcf.gz
```

## Example Prompts

> "Add rsIDs from dbSNP to my VCF file"

> "Annotate my variants with allele frequencies from gnomAD"

> "Predict functional consequences of my variants using bcftools csq"

> "Add ClinVar clinical significance annotations to my VCF"

## See Also

- [bcftools annotate documentation](http://www.htslib.org/doc/bcftools.html#annotate)
- [bcftools csq documentation](http://www.htslib.org/doc/bcftools.html#csq)
- **variant-normalization** - Normalize before annotating
- **filtering-best-practices** - Filter by annotations
