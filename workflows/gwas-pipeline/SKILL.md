---
name: bio-workflows-gwas-pipeline
description: End-to-end GWAS workflow from VCF to association results. Covers PLINK QC, population structure correction, and association testing for case-control or quantitative traits. Use when running genome-wide association studies.
tool_type: mixed
primary_tool: PLINK2
workflow: true
depends_on:
  - population-genetics/plink-basics
  - phasing-imputation/reference-panels
  - phasing-imputation/haplotype-phasing
  - phasing-imputation/genotype-imputation
  - phasing-imputation/imputation-qc
  - population-genetics/population-structure
  - population-genetics/association-testing
  - population-genetics/rare-variant-association
  - population-genetics/linkage-disequilibrium
qc_checkpoints:
  - after_qc: "Sample/variant call rates >95%, HWE p>1e-6"
  - after_imputation: "INFO/R2 or DR2 filtered (MAF-stratified), cases+controls imputed together, dosages carried forward"
  - after_structure: "No population stratification bias"
  - after_association: "Lambda ~1.0, expected QQ plot"
---

## Version Compatibility

Reference examples tested with: ggplot2 3.5+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# GWAS Pipeline

**"Run a GWAS from my genotype data"** -> Orchestrate sample/variant QC (PLINK2), population stratification (PCA), association testing (linear/logistic regression), multiple testing correction, and Manhattan/QQ plot visualization.

Complete workflow for genome-wide association studies from genotype data to significant associations.

## Workflow Overview

```
VCF/PLINK files
    |
    v
[1. QC Filtering] ------> Sample and variant QC
    |
    v
[2. Phase + Impute] ----> Align to panel, phase, impute to dosages, filter by R2 (-> phasing-imputation)
    |
    v
[3. LD Pruning] --------> Independent variants for PCA
    |
    v
[4. Population Structure] --> PCA for covariates
    |
    v
[5. Association Testing] --> Logistic/linear regression on dosages
    |
    v
[6. Results] -----------> Manhattan plot, QQ plot
    |
    v
Significant associations
```

## Step 1: Data Import and QC

### Convert VCF to PLINK

```bash
# VCF to PLINK binary format
plink2 --vcf input.vcf.gz \
    --make-bed \
    --out study

# Or with phenotype/covariate files
plink2 --vcf input.vcf.gz \
    --pheno phenotypes.txt \
    --make-bed \
    --out study
```

QC order is load-bearing (-> population-genetics/plink-basics): variant missingness runs FIRST, in its own invocation, so a sample is not dropped for missingness driven by variants slated for removal.

### Variant call-rate (first)

```bash
plink2 --bfile study --missing --out study_stats
plink2 --bfile study --geno 0.05 --make-bed --out study_var_qc   # variant missingness BEFORE sample missingness
```

### Sample QC

```bash
# Sample missingness AFTER variant missingness.
plink2 --bfile study_var_qc --mind 0.05 --make-bed --out study_sample_qc

# Sex check: split the pseudoautosomal region first or male PAR het reads as a sex error (see plink-basics).
plink2 --bfile study_sample_qc --split-par hg38 --check-sex --out study_sex_check

# KING-robust relatedness - structure-robust and IBD-free (this is the point of KING; no --genome needed).
plink2 --bfile study_sample_qc --king-cutoff 0.0884 --make-bed --out study_unrelated
```

### Variant QC: MAF and controls-only HWE

```bash
# Controls-only HWE with mid-p. plink2 is NOT controls-only by default, so gate to controls explicitly:
# a true risk variant depletes heterozygotes in cases and would fail a case-inclusive HWE test (see plink-basics).
plink2 --bfile study_unrelated --keep-if "PHENO1 == control" --hwe 1e-6 midp --write-snplist --out hwe_pass
plink2 --bfile study_unrelated --maf 0.01 --extract hwe_pass.snplist --make-bed --out study_qc
plink2 --bfile study_qc --freq --out study_qc
```

**QC Checkpoint:**
- Variant call rate >95% (applied before sample missingness)
- Sample call rate >95%
- MAF >1%
- HWE applied to controls only with mid-p (p>1e-6)

## Step 2: Phasing and Imputation

Most GWAS impute the QC'd array genotypes up to a dense reference panel before association, to increase variant density and harmonize across platforms and studies. This stage is owned by the phasing-imputation skills; the workflow only orchestrates the handoff. The decisions that matter here, in order:

1. **Select and prepare the panel** (-> phasing-imputation/reference-panels). Match the panel ancestry to the cohort (TOPMed for diverse/admixed, HRC or 1000G for European, HGDP+1kGP for diverse-and-downloadable), reconcile genome build, and run the strand/allele harmonization check; a flipped palindromic SNP or build mismatch corrupts results without erroring.
2. **Phase, then impute** (-> phasing-imputation/haplotype-phasing, phasing-imputation/genotype-imputation). Pre-phase the QC'd VCF (Eagle2/SHAPEIT5) and impute against the panel per chromosome (Beagle/Minimac4/IMPUTE5), or upload to the Michigan/TOPMed server for controlled-access panels. Impute cases and controls TOGETHER; separate imputation manufactures false associations. The output carries dosages (DS), not hard calls.
3. **Filter by quality** (-> phasing-imputation/imputation-qc). Drop variants below an INFO/R2/DR2 cutoff paired with a MAF floor, MAF-stratified, because a flat cutoff is a hidden rare-variant filter. Carry dosages forward.

**Goal:** Increase variant density and harmonize across platforms by imputing the QC'd genotypes against a dense ancestry-matched panel, carrying dosages (not hard calls) into association.

**Approach:** Export the QC'd genotypes to VCF, hand off to the phasing-imputation skills for strand/build harmonization, phasing, and per-chromosome imputation (cases and controls together), then filter on the engine's quality field plus a MAF floor.

```bash
# Convert QC'd PLINK back to VCF, align to the panel, phase + impute (see phasing-imputation skills for the full commands)
plink2 --bfile study_qc --export vcf bgz --out study_qc
# ... reference-panels: strand/build harmonization; haplotype-phasing: phase; genotype-imputation: impute to dosages ...
# Post-imputation quality filter on the engine's field (DR2 Beagle / R2 Minimac), with a MAF floor
bcftools view -e 'INFO/DR2<0.3 || INFO/AF<0.01 || INFO/AF>0.99' imputed.vcf.gz -Oz -o imputed.qc.vcf.gz
```

**QC Checkpoint:** cases and controls imputed together; INFO/R2 filtered (MAF-stratified) with a MAF floor; dosages (DS), not hard calls, carried into association.

## Step 3: LD Pruning for PCA

```bash
# Exclude long-range-LD regions and inversions FIRST (MHC, 8p23.1, 17q21.31, LCT) - their internal r2 is
# high and real, so a window prune keeps them and a PC then tracks the inversion (-> population-structure).
plink2 --bfile study_qc --exclude range longrange_ld.txt --make-bed --out study_noLR

# Identify independent variants (r2 0.1 matches the linkage-disequilibrium / population-structure default).
plink2 --bfile study_noLR --indep-pairwise 50 5 0.1 --out pruned
plink2 --bfile study_noLR --extract pruned.prune.in --make-bed --out study_pruned
```

## Step 4: Population Structure (PCA)

```bash
# Calculate principal components
plink2 --bfile study_pruned \
    --pca 10 \
    --out study_pca

# The eigenvec file contains PCs for use as covariates
```

### Visualize PCA

```r
library(ggplot2)

# Load PCA results
pca <- read.table('study_pca.eigenvec', header = FALSE)
colnames(pca) <- c('FID', 'IID', paste0('PC', 1:10))

# Load phenotype for coloring
pheno <- read.table('phenotypes.txt', header = TRUE)
pca <- merge(pca, pheno, by = c('FID', 'IID'))

# Plot
ggplot(pca, aes(x = PC1, y = PC2, color = as.factor(PHENO))) +
    geom_point(alpha = 0.5) +
    labs(title = 'PCA of Study Samples', color = 'Phenotype') +
    theme_minimal()
ggsave('pca_plot.pdf', width = 8, height = 6)
```

## Step 5: Association Testing

Run the association on imputed DOSAGES, not hard calls, so the imputation uncertainty is propagated (PLINK2 reads dosages with `--vcf imputed.qc.vcf.gz dosage=DS`, or use a `.pgen` built from dosages). The examples below use the QC'd best-guess genotypes for brevity; substitute the dosage input for an imputed analysis.

The engine choice is set by sample structure, not convenience (-> population-genetics/association-testing). PC-covariate GLM (below) is valid only for unrelated samples whose confounding is continuous ancestry; any related, family-based, or fine-scale-structured cohort needs a linear mixed model (BOLT-LMM, SAIGE, regenie) with leave-one-chromosome-out, because PCs cannot remove a covariance structure. Use SPA (SAIGE) or Firth when the case:control ratio is more extreme than ~1:10 or minor-allele counts are low. For aggregating rare variants by gene rather than testing them one at a time, branch to population-genetics/rare-variant-association (burden/SKAT/SKAT-O via SAIGE-GENE+/regenie).

### Case-Control (Binary Trait)

```bash
# Logistic regression with PCA covariates
plink2 --bfile study_qc \
    --pheno phenotypes.txt \
    --covar study_pca.eigenvec \
    --covar-col-nums 3-12 \
    --glm firth-fallback hide-covar \
    --out gwas_results

# Binary --glm defaults to firth-fallback -> results in gwas_results.PHENO.glm.logistic.hybrid
```

### Quantitative Trait

```bash
# Linear regression
plink2 --bfile study_qc \
    --pheno phenotypes.txt \
    --pheno-name BMI \
    --covar study_pca.eigenvec \
    --covar-col-nums 3-12 \
    --glm hide-covar \
    --out gwas_bmi

# Results in gwas_bmi.BMI.glm.linear
```

### With Additional Covariates

```bash
# Include age, sex, and PCs
plink2 --bfile study_qc \
    --pheno phenotypes.txt \
    --covar covariates.txt \
    --covar-name AGE,SEX,PC1-PC10 \
    --glm hide-covar \
    --out gwas_adjusted
```

## Step 6: Results Visualization

### Manhattan Plot

```r
library(qqman)

# Load results
results <- read.table('gwas_results.PHENO.glm.logistic.hybrid', header = TRUE)
results <- results[!is.na(results$P),]

# Manhattan plot
png('manhattan.png', width = 1200, height = 600)
manhattan(results, chr = 'X.CHROM', bp = 'POS', snp = 'ID', p = 'P',
          suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8))
dev.off()

# QQ plot
png('qq_plot.png', width = 600, height = 600)
qq(results$P)
dev.off()
```

### Calculate Genomic Inflation

```r
# Lambda (genomic inflation factor)
chisq <- qchisq(1 - results$P, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
cat('Lambda:', round(lambda, 3), '\n')
# Lambda should be close to 1.0 (1.0-1.1 acceptable)
```

### Extract Significant Hits

```bash
# Select the P column by header, not a fixed index (firth-fallback adds columns and shifts positions).
awk 'NR==1{for(i=1;i<=NF;i++) if($i=="P") p=i; print; next} $p<5e-8' \
    gwas_results.PHENO.glm.logistic.hybrid > significant_hits.txt
awk 'NR==1{for(i=1;i<=NF;i++) if($i=="P") p=i; print; next} $p<1e-5' \
    gwas_results.PHENO.glm.logistic.hybrid > suggestive_hits.txt
```

## Parameter Recommendations

| Step | Parameter | Value |
|------|-----------|-------|
| Sample QC | --mind | 0.05 |
| Variant QC | --geno | 0.05 |
| Variant QC | --maf | 0.01 |
| Variant QC | --hwe | 1e-6 |
| LD pruning | --indep-pairwise | 50 5 0.1 (after long-range-LD exclusion) |
| PCA | --pca | 10 |
| Significance | p-value | 5e-8 |

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| Elevated lambda | Polygenicity OR confounding - lambda alone cannot tell them apart | Check the LDSC intercept (~1 = polygenic inflation, not confounding); do NOT reflexively genomic-control, which over-corrects true signal |
| Residual structure after PCs | Relatedness/fine-scale structure (a covariance PCs cannot remove) | Switch to a linear mixed model (BOLT-LMM/SAIGE/regenie) with LOCO |
| Inflation at low MAC or extreme case:control ratio | Score/Wald test anti-conservative | Use SPA (SAIGE) or Firth regression |
| No significant hits | Low power | Increase sample size, meta-analysis; for rare variants aggregate by gene (rare-variant-association) |
| QQ deviation at low end | Batch effects / differential missingness | Check technical artifacts and case/control missingness skew |

## Complete Pipeline Script

```bash
#!/bin/bash
set -e

INPUT_VCF="genotypes.vcf.gz"
PHENO="phenotypes.txt"
OUTDIR="gwas_results"

mkdir -p ${OUTDIR}

# Step 1: Convert, then QC in order (variant missingness, then sample, then MAF + controls-only HWE).
plink2 --vcf ${INPUT_VCF} --make-bed --out ${OUTDIR}/raw
plink2 --bfile ${OUTDIR}/raw --geno 0.05 --make-bed --out ${OUTDIR}/var_qc
plink2 --bfile ${OUTDIR}/var_qc --mind 0.05 --king-cutoff 0.0884 --make-bed --out ${OUTDIR}/samp_qc
plink2 --bfile ${OUTDIR}/samp_qc --keep-if "PHENO1 == control" --hwe 1e-6 midp \
    --write-snplist --out ${OUTDIR}/hwe_pass
plink2 --bfile ${OUTDIR}/samp_qc --maf 0.01 --extract ${OUTDIR}/hwe_pass.snplist \
    --make-bed --out ${OUTDIR}/qc

# Step 2: LD pruning (exclude long-range-LD regions first; r2 0.1)
plink2 --bfile ${OUTDIR}/qc --exclude range longrange_ld.txt --make-bed --out ${OUTDIR}/noLR
plink2 --bfile ${OUTDIR}/noLR --indep-pairwise 50 5 0.1 --out ${OUTDIR}/pruned
plink2 --bfile ${OUTDIR}/noLR --extract ${OUTDIR}/pruned.prune.in \
    --make-bed --out ${OUTDIR}/pruned_set

# Step 3: PCA
plink2 --bfile ${OUTDIR}/pruned_set --pca 10 --out ${OUTDIR}/pca

# Step 4: Association on the full QC'd set (binary --glm defaults to firth-fallback -> .glm.logistic.hybrid)
plink2 --bfile ${OUTDIR}/qc --pheno ${PHENO} \
    --covar ${OUTDIR}/pca.eigenvec --covar-col-nums 3-12 \
    --glm firth-fallback hide-covar --out ${OUTDIR}/gwas

echo "=== GWAS Complete ==="
echo "Results: ${OUTDIR}/gwas.*.glm.*"
```

## Related Skills

- database-access/ensembl-rest - VEP annotation for top GWAS variants (per-variant); local VEP for >1K
- database-access/biomart-queries - Bulk SNP-to-gene mapping via BioMart
- population-genetics/plink-basics - PLINK file formats and commands
- phasing-imputation/reference-panels - Select and prepare the reference panel; strand/build harmonization
- phasing-imputation/haplotype-phasing - Pre-phase the QC'd genotypes before imputation
- phasing-imputation/genotype-imputation - Impute to dosages against the panel
- phasing-imputation/imputation-qc - Filter imputed variants by R2 and MAF before association
- population-genetics/population-structure - PCA and admixture
- population-genetics/association-testing - Single-variant models (PC-GLM vs LMM, SPA/Firth) on dosages
- population-genetics/rare-variant-association - Gene-based burden/SKAT/SKAT-O for rare variants
- population-genetics/linkage-disequilibrium - LD concepts
