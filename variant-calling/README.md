# variant-calling

## Overview

Variant calling and VCF/BCF file manipulation using bcftools and cyvcf2. Covers calling SNPs/indels from alignments, structural variant detection, filtering, normalization, annotation, and downstream analysis.

**Tool type:** cli | **Primary tools:** bcftools, cyvcf2, Manta, Delly

## Skills

| Skill | Description |
|-------|-------------|
| vcf-basics | View, query, understand VCF/BCF format structure |
| variant-calling | Call SNPs/indels from BAM files using mpileup/call |
| gatk-variant-calling | GATK HaplotypeCaller, GVCF workflow, VQSR filtering |
| structural-variant-calling | Call SVs (DEL, DUP, INV, INS, BND) with Manta/Delly |
| vcf-filtering | Filter variants by quality, depth, type, expressions |
| vcf-manipulation | Merge, concat, sort, intersect VCF files |
| variant-normalization | Left-align indels, split multiallelic sites |
| variant-annotation | Add annotations, predict functional consequences |
| vep-snpeff-annotation | Comprehensive annotation with VEP, SnpEff, ANNOVAR |
| vcf-statistics | Generate quality metrics, Ti/Tv ratio, concordance |
| consensus-sequences | Apply variants to reference FASTA |

## Example Prompts

- "Call variants from my aligned BAM file"
- "Run GATK HaplotypeCaller on my sample"
- "Joint genotype my cohort with GATK"
- "Call structural variants with Manta"
- "Detect deletions and inversions with Delly"
- "Merge SV calls from multiple callers"
- "View the first 20 variants in my VCF"
- "Extract chromosome 1 variants to a new file"
- "List all sample names in this VCF"
- "Filter variants with QUAL < 30"
- "Keep only SNPs with depth >= 10"
- "Extract PASS variants only"
- "Get rare variants with AF < 0.01"
- "Merge VCF files from different samples"
- "Compare variants between two callers"
- "Find variants shared between two files"
- "Concatenate per-chromosome VCFs"
- "Normalize indels to left-aligned representation"
- "Split multiallelic sites to biallelic"
- "Add rsIDs from dbSNP"
- "Predict functional consequences"
- "Calculate transition/transversion ratio"
- "Generate consensus sequence from variants"
- "Annotate variants with VEP"
- "Run SnpEff on my VCF"
- "Add CADD scores to my variants"

## Requirements

```bash
# bcftools
conda install -c bioconda bcftools

# cyvcf2
pip install cyvcf2

# SV callers
conda install -c bioconda manta delly smoove survivor

# Annotation tools
conda install -c bioconda ensembl-vep snpeff

# GATK
conda install -c bioconda gatk4
```

## Related Skills

- **alignment-files** - Prepare BAM files for variant calling
- **copy-number** - CNV detection (complementary to SV calling)
- **long-read-sequencing** - Long-read SV detection
- **population-genetics** - Population-level analysis of variants
- **database-access** - Download reference databases (dbSNP, gnomAD)
- **sequence-io** - Work with consensus FASTA output
