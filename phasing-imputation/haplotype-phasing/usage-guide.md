# Haplotype Phasing Usage Guide

## Overview

Haplotype phasing determines which alleles are inherited together on each chromosome. This is essential for:
- Genotype imputation
- Haplotype-based association tests
- Population genetics (haplotype frequency, identity by descent)
- Compound heterozygosity detection

## Choosing a Phasing Tool

| Tool | Best For | Memory | Speed |
|------|----------|--------|-------|
| Beagle 5.4 | General purpose, <50k samples | Medium | Fast |
| SHAPEIT5 | Large biobanks (>50k samples) | Low | Very fast |
| Eagle2 | Phasing with ref panel | Medium | Fast |

## Beagle Workflow

```bash
# Download Beagle
wget https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar

# Download genetic maps
wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip plink.GRCh38.map.zip -d genetic_maps/

# Prepare input
bcftools view -m2 -M2 -v snps input.vcf.gz -Oz -o biallelic.vcf.gz
bcftools index biallelic.vcf.gz

# Phase per chromosome
for chr in {1..22}; do
    java -Xmx16g -jar beagle.jar \
        gt=biallelic.vcf.gz \
        chrom=chr${chr} \
        map=genetic_maps/plink.chr${chr}.GRCh38.map \
        out=phased_chr${chr} \
        nthreads=8
done

# Merge chromosomes
ls phased_chr*.vcf.gz > filelist.txt
bcftools concat -f filelist.txt -Oz -o phased.vcf.gz
bcftools index phased.vcf.gz
```

## SHAPEIT5 Workflow (Large Datasets)

```bash
# Install
conda install -c bioconda shapeit5

# For datasets >10,000 samples, use two-stage approach
# Stage 1: Phase common variants (MAF > 0.1%)
shapeit5_phase_common \
    --input input.bcf \
    --map genetic_map.txt \
    --output phased_common.bcf \
    --thread 16 \
    --log phase_common.log \
    --filter-maf 0.001

# Stage 2: Phase rare variants using common as scaffold
shapeit5_phase_rare \
    --input input.bcf \
    --scaffold phased_common.bcf \
    --map genetic_map.txt \
    --output phased_all.bcf \
    --thread 16 \
    --log phase_rare.log
```

## Using Reference Panel for Better Phasing

```bash
# Beagle with 1000 Genomes reference
java -Xmx16g -jar beagle.jar \
    gt=study.vcf.gz \
    ref=1000GP.chr22.vcf.gz \
    map=plink.chr22.GRCh38.map \
    out=phased_with_ref \
    nthreads=8

# Benefits:
# - Better phasing accuracy for rare variants
# - Required if imputing from same reference
```

## Input QC Before Phasing

```bash
# 1. Check for strand issues
bcftools +fixref input.vcf.gz -Oz -o fixed.vcf.gz -- -f reference.fa

# 2. Filter to biallelic SNPs
bcftools view -m2 -M2 -v snps fixed.vcf.gz -Oz -o biallelic.vcf.gz

# 3. Remove problematic regions
bcftools view -T ^exclude_regions.bed biallelic.vcf.gz -Oz -o filtered.vcf.gz

# 4. Check sample missingness
bcftools stats filtered.vcf.gz | grep "PSC"
```

## Evaluating Phasing Quality

### With Trio Data
```bash
# If you have parent-offspring trios
# Compare phased haplotypes to Mendelian inheritance

# Extract trio
bcftools view -s child,mother,father phased.vcf.gz -Oz -o trio.vcf.gz

# Count Mendelian errors
bcftools +mendelian trio.vcf.gz -p child,mother,father
```

### Switch Error Rate
```python
def calculate_switch_error(truth_vcf, test_vcf):
    '''Calculate switch error rate between phased VCFs.'''
    # Implementation depends on having truth haplotypes
    # Switch error = rate of switches needed to match truth
    pass
```

## Common Issues

### High Missingness
```bash
# Beagle can handle missing data, but better to filter extreme cases
bcftools view -i 'F_MISSING < 0.1' input.vcf.gz -Oz -o low_missing.vcf.gz
```

### Memory Errors
```bash
# Increase heap size
java -Xmx32g -jar beagle.jar ...

# Or phase in smaller chunks
java -jar beagle.jar gt=input.vcf.gz chrom=chr22:1-25000000 ...
```

### Multiallelic Sites
```bash
# Split multiallelic (Beagle handles biallelic best)
bcftools norm -m- input.vcf.gz -Oz -o split.vcf.gz
```

## Performance Tips

1. **Use genetic maps** - significantly improves accuracy
2. **Process chromosomes in parallel** - each chr independent
3. **Match reference genome build** - GRCh37 vs GRCh38
4. **Filter before phasing** - remove low-quality variants
5. **Use reference panel** - especially for rare variants

## Output Format

Phased genotypes use `|` separator:
```
# Unphased
GT: 0/1  (unknown which allele on which haplotype)

# Phased
GT: 0|1  (REF on haplotype 1, ALT on haplotype 2)
GT: 1|0  (ALT on haplotype 1, REF on haplotype 2)
```
