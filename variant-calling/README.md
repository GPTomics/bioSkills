# variant-calling

## Overview

Variant calling and VCF/BCF file manipulation. Covers germline SNP/indel calling (bcftools, GATK HaplotypeCaller, DeepVariant), structural variant detection (Manta, Delly, GRIDSS), filtering (VQSR, hard filters, allele-specific), normalization, annotation (VEP, SnpEff, ANNOVAR), clinical interpretation (ACMG/ClinVar), and VCF utilities.

**Tool type:** mixed | **Primary tools:** bcftools, GATK, DeepVariant, VEP, Manta, Delly

## Skills

| Skill | Description |
|-------|-------------|
| vcf-basics | View, query, and interpret VCF/BCF: QUAL/GQ/PL-GL levels, AD/DP and allele balance, GT phasing/ploidy/PS, missing-vs-hom-ref, INFO/FORMAT Number A/R/G, symbolic alleles/END, and the gVCF <NON_REF> model |
| variant-calling | Call germline SNPs/indels from a BAM with bcftools mpileup/call; engine-selection hub (bcftools vs GATK vs DeepVariant vs DRAGEN) with ploidy and normalization guidance |
| gatk-variant-calling | GATK HaplotypeCaller: local-reassembly + PairHMM mechanism, the -ERC GVCF model, BQSR-vs-DRAGSTR/--dragen-mode, and ploidy/mitochondria/contamination edge cases |
| deepvariant | CNN image-classification germline calling: platform model selection (WGS/WES/PACBIO/ONT_R104/HYBRID), three-stage pipeline, DeepTrio de-novo, GLnexus joint calling; no post-hoc GATK filters/BQSR |
| joint-calling | Cohort joint genotyping (GATK GenomicsDBImport/CombineGVCFs/GenotypeGVCFs, GLnexus for DeepVariant); the merge-vs-genotype trap, N+1, spanning-deletion allele, biobank scaling |
| structural-variant-calling | Call SVs by reconstructing 4 signals (RP/SR/depth/assembly): caller blind-spots, SVLEN-sign/BND VCF traps, force-genotyping cohorts, Truvari benchmarking, long-read switch |
| filtering-best-practices | Site vs genotype filtering: VQSR/VETS/NVScoreVariants vs hard filters by cohort size, SNP/indel thresholds, the hom-alt missing-annotation trap, Ti/Tv validation |
| vcf-manipulation | Merge, concat, sort, intersect, and subset VCFs with bcftools; normalize-before-combine, the single-sample-merge 0/0 trap, and merge vs concat vs isec selection |
| variant-normalization | Left-align/trim to parsimonious canonical form, atomize MNPs, split multiallelics; pipeline order, vt-vs-bcftools discordance, the HGVS 3'-rule clash, and the silent homopolymer annotation miss |
| variant-annotation | VEP/SnpEff/ANNOVAR non-determinism, MANE Select vs --pick danger, HGVS 3'-rule clash, NMD/PVS1, one calibrated predictor (REVEL/AlphaMissense/SpliceAI), grpmax filtering AF |
| clinical-interpretation | Current ACMG/AMP (graded PVS1, PM2_Supporting, calibrated PP3/BP4, Bayesian points) + somatic tiers/oncogenicity, ClinVar-star and grpmax-FAF evidence rules, VCEP override, VUS reanalysis |
| vcf-statistics | Interpret QC metrics (Ti/Tv, het/hom by ancestry, HWE excess-het, contamination, relatedness) and screen cohorts for swaps/sex with bcftools stats, vcftools, somalier, peddy, KING |
| consensus-sequences | Apply variants to a reference FASTA (bcftools) or build viral consensus (iVar), with phasing/chimera checks, no-coverage masking, and diploid/symbolic-SV boundaries |

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
