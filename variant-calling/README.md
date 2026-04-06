# variant-calling

## Overview

Variant calling and VCF/BCF file manipulation. Covers germline SNP/indel calling (bcftools, GATK HaplotypeCaller, DeepVariant), structural variant detection (Manta, Delly, GRIDSS), filtering (VQSR, hard filters, allele-specific), normalization, annotation (VEP, SnpEff, ANNOVAR), clinical interpretation (ACMG/ClinVar), and VCF utilities.

**Tool type:** mixed | **Primary tools:** bcftools, GATK, DeepVariant, VEP, Manta, Delly

## Skills

| Skill | Description |
|-------|-------------|
| vcf-basics | View, query, understand VCF/BCF format structure and field interpretation (QUAL vs GQ, AD vs DP) |
| variant-calling | Call SNPs/indels from BAM files using bcftools mpileup/call with caller selection guidance |
| gatk-variant-calling | GATK HaplotypeCaller with DRAGEN-GATK mode, GVCF workflow, VQSR and hard filtering |
| deepvariant | Deep learning variant calling with CNN-based model selection (WGS/WES/PacBio/ONT) |
| joint-calling | Multi-sample joint genotyping with GenomicsDB/CombineGVCFs and GLnexus, cohort size guidance |
| structural-variant-calling | Call SVs with Manta/Delly/GRIDSS consensus approach, short-read vs long-read decision framework |
| filtering-best-practices | VQSR decision tree, hard filter rationale, allele-specific filtering, common pitfalls |
| vcf-manipulation | Merge, concat, sort, intersect VCF files with normalization-before-comparison guidance |
| variant-normalization | Left-align indels, split multiallelic sites with pipeline order and splitting caveats |
| variant-annotation | VEP/SnpEff/ANNOVAR with MANE Select transcripts, tool concordance data, AlphaMissense |
| clinical-interpretation | ClinVar/ACMG classification rules, pathogenicity score interpretation, re-analysis guidance |
| vcf-statistics | QC metrics with expected ranges, stratified evaluation, population-scale QC |
| consensus-sequences | Apply variants to reference FASTA with phasing requirements and diploid considerations |

## Example Prompts

- "Call variants from my aligned BAM file"
- "Run GATK HaplotypeCaller on my sample"
- "Joint genotype my cohort with GATK"
- "Call structural variants with Manta"
- "Detect deletions and inversions with Delly"
- "Merge SV calls from multiple callers"
- "View the first 20 variants in my VCF"
- "Filter variants with QUAL < 30"
- "Keep only SNPs with depth >= 10"
- "Extract PASS variants only"
- "Get rare variants with AF < 0.01"
- "Merge VCF files from different samples"
- "Normalize indels to left-aligned representation"
- "Add rsIDs from dbSNP"
- "Annotate variants with VEP"
- "Run SnpEff on my VCF"
- "Add CADD scores to my variants"
- "Generate consensus sequence from variants"

## Requirements

```bash
# Core tools
conda install -c bioconda bcftools htslib samtools
pip install cyvcf2

# Variant callers
conda install -c bioconda gatk4
# DeepVariant: use Docker (google/deepvariant:1.6.1)

# SV callers
conda install -c bioconda manta delly smoove survivor
# GRIDSS: requires Java 11+ and R

# Annotation tools
conda install -c bioconda ensembl-vep snpeff

# Joint calling
# GLnexus: use Docker (quay.io/mlin/glnexus:v1.4.1)
```

## Related Skills

- **alignment-files** - Prepare BAM files for variant calling
- **copy-number** - CNV detection (complementary to SV calling)
- **long-read-sequencing** - Long-read SV detection
- **population-genetics** - Population-level analysis of variants
- **database-access** - Download reference databases (dbSNP, gnomAD)
