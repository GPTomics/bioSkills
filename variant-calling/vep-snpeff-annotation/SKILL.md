---
name: bio-variant-calling-vep-snpeff-annotation
description: Comprehensive variant effect prediction with VEP (Ensembl), SnpEff, and ANNOVAR. Provides detailed functional annotations including gene impact, protein changes, population frequencies, and clinical significance beyond basic bcftools csq.
tool_type: cli
primary_tool: VEP
---

# VEP/SnpEff Annotation

Comprehensive variant annotation tools that provide detailed functional predictions, population frequencies, and clinical annotations.

## Tool Comparison

| Tool | Strengths | Speed | Database |
|------|-----------|-------|----------|
| VEP | Ensembl integration, plugins | Moderate | Ensembl |
| SnpEff | Fast, simple | Fast | Custom/Ensembl |
| ANNOVAR | Flexible databases | Moderate | Multiple |

## Ensembl VEP

### Installation

```bash
# Conda (recommended)
conda install -c bioconda ensembl-vep

# Download cache for offline use
vep_install -a cf -s homo_sapiens -y GRCh38 --CONVERT
```

### Basic Annotation

```bash
vep -i input.vcf -o output.vep.txt --cache --offline
```

### VCF Output

```bash
vep -i input.vcf -o output.vcf --vcf --cache --offline
```

### Common Options

```bash
vep -i input.vcf -o output.vcf \
    --vcf \
    --cache \
    --offline \
    --species homo_sapiens \
    --assembly GRCh38 \
    --everything \
    --fork 4
```

### --everything Flag

Enables commonly used annotations:
- `--sift b` - SIFT predictions
- `--polyphen b` - PolyPhen predictions
- `--ccds` - CCDS identifiers
- `--hgvs` - HGVS nomenclature
- `--symbol` - Gene symbols
- `--numbers` - Exon/intron numbers
- `--domains` - Protein domains
- `--regulatory` - Regulatory features
- `--canonical` - Canonical transcript
- `--protein` - Protein sequence
- `--biotype` - Transcript biotype
- `--af` - 1000 Genomes frequencies
- `--af_gnomade` - gnomAD exome frequencies
- `--af_gnomadg` - gnomAD genome frequencies
- `--max_af` - Maximum allele frequency
- `--pubmed` - PubMed IDs
- `--uniprot` - UniProt identifiers

### Filtering by Impact

```bash
# Only output high/moderate impact
vep -i input.vcf -o output.vcf --vcf \
    --cache --offline \
    --pick \
    --filter_common \
    --filter "IMPACT in HIGH,MODERATE"
```

### Pick One Consequence

```bash
# Pick most severe consequence per variant
vep -i input.vcf -o output.vcf --vcf \
    --cache --offline \
    --pick

# Pick per gene (one consequence per gene)
vep -i input.vcf -o output.vcf --vcf \
    --cache --offline \
    --per_gene --pick_order canonical,biotype,rank
```

### Plugins

```bash
# CADD scores
vep -i input.vcf -o output.vcf --vcf \
    --cache --offline \
    --plugin CADD,whole_genome_SNVs.tsv.gz

# dbNSFP (multiple predictors)
vep -i input.vcf -o output.vcf --vcf \
    --cache --offline \
    --plugin dbNSFP,dbNSFP4.3a.gz,ALL

# Multiple plugins
vep -i input.vcf -o output.vcf --vcf \
    --cache --offline \
    --plugin CADD,cadd.tsv.gz \
    --plugin dbNSFP,dbnsfp.gz,SIFT_score,Polyphen2_HDIV_score \
    --plugin SpliceAI,spliceai.vcf.gz
```

### Custom Annotations

```bash
# Add custom VCF annotations
vep -i input.vcf -o output.vcf --vcf \
    --cache --offline \
    --custom gnomad.vcf.gz,gnomAD,vcf,exact,0,AF

# Add BED annotations
vep -i input.vcf -o output.vcf --vcf \
    --cache --offline \
    --custom regions.bed.gz,Regions,bed,overlap
```

### VEP Output Fields

| Field | Description |
|-------|-------------|
| Consequence | SO term (e.g., missense_variant) |
| IMPACT | HIGH, MODERATE, LOW, MODIFIER |
| SYMBOL | Gene symbol |
| Gene | Ensembl gene ID |
| Feature | Transcript ID |
| BIOTYPE | Transcript biotype |
| EXON | Exon number |
| HGVSc | HGVS coding change |
| HGVSp | HGVS protein change |
| cDNA_position | Position in cDNA |
| CDS_position | Position in CDS |
| Protein_position | Position in protein |
| Amino_acids | Ref/Alt amino acids |
| Codons | Ref/Alt codons |
| SIFT | SIFT prediction |
| PolyPhen | PolyPhen prediction |

## SnpEff

### Installation

```bash
conda install -c bioconda snpeff

# List available databases
snpEff databases | grep -i "homo_sapiens"

# Download database
snpEff download GRCh38.105
```

### Basic Annotation

```bash
snpEff ann GRCh38.105 input.vcf > output.vcf
```

### Common Options

```bash
snpEff ann \
    -v \
    -stats stats.html \
    -csvStats stats.csv \
    GRCh38.105 \
    input.vcf > output.vcf
```

### Filter by Impact

```bash
# Annotate then filter
snpEff ann GRCh38.105 input.vcf | \
    SnpSift filter "(ANN[*].IMPACT = 'HIGH')" > high_impact.vcf

# Or filter during annotation
snpEff ann -no-downstream -no-intergenic -no-intron \
    GRCh38.105 input.vcf > coding.vcf
```

### SnpEff ANN Field

Format: `Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA_pos|CDS_pos|Protein_pos|Distance|Notes`

```bash
# Extract ANN field
bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' output.vcf | \
    awk -F'|' '{print $1"\t"$2"\t"$3"\t"$4}'
```

### SnpEff Impact Categories

| Impact | Examples |
|--------|----------|
| HIGH | Stop gained, frameshift, splice donor/acceptor |
| MODERATE | Missense, inframe indel |
| LOW | Synonymous, splice region |
| MODIFIER | Intron, intergenic, UTR |

### Add Database Annotations (SnpSift)

```bash
# dbSNP
SnpSift annotate dbsnp.vcf.gz input.vcf > annotated.vcf

# ClinVar
SnpSift annotate clinvar.vcf.gz input.vcf > annotated.vcf

# dbNSFP
SnpSift dbnsfp -db dbNSFP4.3a.txt.gz input.vcf > annotated.vcf

# Chain multiple
snpEff ann GRCh38.105 input.vcf | \
    SnpSift annotate dbsnp.vcf.gz | \
    SnpSift annotate clinvar.vcf.gz > fully_annotated.vcf
```

### SnpSift Filtering

```bash
# Filter by expression
SnpSift filter "(QUAL >= 30) & (DP >= 10)" input.vcf > filtered.vcf

# Filter by ClinVar significance
SnpSift filter "(exists CLNSIG) & (CLNSIG has 'Pathogenic')" input.vcf > pathogenic.vcf

# Filter by allele frequency
SnpSift filter "(AF < 0.01)" input.vcf > rare.vcf
```

## ANNOVAR

### Installation

```bash
# Download from ANNOVAR website (requires registration)
# https://annovar.openbioinformatics.org/

# Download databases
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad30_genome humandb/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20230416 humandb/
```

### Convert VCF to ANNOVAR Format

```bash
convert2annovar.pl -format vcf4 input.vcf > input.avinput
```

### Gene-Based Annotation

```bash
annotate_variation.pl -geneanno -dbtype refGene -buildver hg38 \
    input.avinput humandb/
```

### Table Annotation

```bash
table_annovar.pl input.vcf humandb/ \
    -buildver hg38 \
    -out annotated \
    -remove \
    -protocol refGene,gnomad30_genome,clinvar_20230416,dbnsfp42a \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
```

### Common ANNOVAR Databases

| Database | Type | Content |
|----------|------|---------|
| refGene | Gene | RefSeq gene annotations |
| ensGene | Gene | Ensembl gene annotations |
| gnomad30_genome | Filter | gnomAD v3 frequencies |
| clinvar | Filter | ClinVar pathogenicity |
| dbnsfp42a | Filter | dbNSFP predictions |
| cosmic70 | Filter | COSMIC mutations |

## Python: Parse Annotated VCF

### Parse VEP Annotations

```python
from cyvcf2 import VCF

def parse_vep_csq(csq_string, csq_header):
    '''Parse VEP CSQ field into dict.'''
    fields = csq_header.split('|')
    values = csq_string.split('|')
    return dict(zip(fields, values))

vcf = VCF('vep_output.vcf')
csq_header = None
for h in vcf.header_iter():
    if h['HeaderType'] == 'INFO' and h['ID'] == 'CSQ':
        csq_header = h['Description'].split('Format: ')[1].rstrip('"')
        break

for variant in vcf:
    csq = variant.INFO.get('CSQ')
    if csq:
        for transcript in csq.split(','):
            parsed = parse_vep_csq(transcript, csq_header)
            if parsed.get('IMPACT') in ('HIGH', 'MODERATE'):
                print(f"{variant.CHROM}:{variant.POS} {parsed['SYMBOL']} {parsed['Consequence']}")
```

### Parse SnpEff Annotations

```python
from cyvcf2 import VCF

def parse_snpeff_ann(ann_string):
    '''Parse SnpEff ANN field.'''
    fields = ['Allele', 'Annotation', 'Impact', 'Gene_Name', 'Gene_ID',
              'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank',
              'HGVS_c', 'HGVS_p', 'cDNA_pos', 'CDS_pos', 'Protein_pos', 'Distance']
    values = ann_string.split('|')
    return dict(zip(fields, values[:len(fields)]))

for variant in VCF('snpeff_output.vcf'):
    ann = variant.INFO.get('ANN')
    if ann:
        for transcript in ann.split(','):
            parsed = parse_snpeff_ann(transcript)
            if parsed['Impact'] == 'HIGH':
                print(f"{variant.CHROM}:{variant.POS} {parsed['Gene_Name']} {parsed['Annotation']}")
```

## Complete Annotation Pipeline

```bash
#!/bin/bash
set -euo pipefail

INPUT=$1
REFERENCE=$2
VEP_CACHE=$3
OUTPUT_PREFIX=$4

# Normalize variants
bcftools norm -f $REFERENCE -m-any $INPUT -Oz -o ${OUTPUT_PREFIX}_norm.vcf.gz
bcftools index ${OUTPUT_PREFIX}_norm.vcf.gz

# VEP annotation
vep -i ${OUTPUT_PREFIX}_norm.vcf.gz \
    -o ${OUTPUT_PREFIX}_vep.vcf \
    --vcf \
    --cache --offline --dir_cache $VEP_CACHE \
    --assembly GRCh38 \
    --everything \
    --pick \
    --fork 4

bgzip ${OUTPUT_PREFIX}_vep.vcf
bcftools index ${OUTPUT_PREFIX}_vep.vcf.gz

# Filter high/moderate impact
bcftools view -i 'INFO/CSQ~"HIGH" || INFO/CSQ~"MODERATE"' \
    ${OUTPUT_PREFIX}_vep.vcf.gz \
    -Oz -o ${OUTPUT_PREFIX}_filtered.vcf.gz

echo "Done: ${OUTPUT_PREFIX}_filtered.vcf.gz"
```

## Annotation Field Reference

### Pathogenicity Predictors

| Predictor | Deleterious | Benign |
|-----------|-------------|--------|
| SIFT | < 0.05 | >= 0.05 |
| PolyPhen-2 (HDIV) | > 0.957 (probably), > 0.453 (possibly) | <= 0.453 |
| CADD | > 20 (top 1%), > 30 (top 0.1%) | < 10 |
| REVEL | > 0.5 | < 0.5 |

### Clinical Significance (ClinVar)

| Code | Meaning |
|------|---------|
| Pathogenic | Disease-causing |
| Likely_pathogenic | Probably disease-causing |
| Uncertain_significance | VUS |
| Likely_benign | Probably not disease-causing |
| Benign | Not disease-causing |

## Related Skills

- variant-annotation - Basic bcftools annotate/csq
- variant-normalization - Normalize before annotating
- vcf-filtering - Filter by annotations
- database-access - Download annotation databases
