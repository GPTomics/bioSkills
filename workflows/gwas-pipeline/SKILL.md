---
name: bio-workflows-gwas-pipeline
description: Orchestrates the GWAS pipeline from genotypes to association results, chaining PLINK2 QC (variant-then-sample missingness, controls-only HWE, KING relatedness), panel harmonization + joint phasing/imputation to dosages, long-range-LD-excluded PCA, and an engine chosen by sample structure (PLINK2-GLM / regenie / SAIGE / BOLT-LMM), with LDSC-intercept diagnostics. Use when committing the genome build + ancestry-matched imputation panel once (ancestry match > panel size), running the strand/allele harmonization gate (drop intermediate-frequency palindromes), imputing cases+controls TOGETHER on dosages, excluding long-range-LD regions before PCA, choosing an LMM when relatedness/structure is present (PCs cannot remove a covariance), or separating polygenicity from confounding via the LDSC intercept. Hands mechanism to the population-genetics and phasing-imputation component skills; not a re-teach of any single step.
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

Reference examples tested with: PLINK 2.0 (alpha 5+), regenie 3.4+, SAIGE 1.3+, Eagle 2.4+ / SHAPEIT5, Minimac4 / Beagle 5.4, bcftools 1.19+, LDSC 1.0, qqman 0.1.9+, ggplot2 3.5+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: regenie/SAIGE run a two-step design (step 1 whole-genome ridge / null GLMM on LD-pruned common markers; step 2 tests the imputed set with LOCO) — step-1 and step-2 sample IDs + variance-ratio file MUST match or calibration silently fails. Binary PLINK2 `--glm` defaults to `firth-fallback` (output `.glm.logistic.hybrid`). Sequencing/WGS and non-EUR LD shift the genome-wide threshold off the 5e-8 EUR-array folklore. Confirm in-tool before quoting.

# GWAS Pipeline

**"Run a GWAS from my genotype data"** -> Orchestrate sample/variant QC (PLINK2), population stratification (PCA), association testing (linear/logistic regression), multiple testing correction, and Manhattan/QQ plot visualization.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

A GWAS result is decided at four seams, not inside the association test.

1. **The genome build + ancestry-matched imputation panel is a made-once commitment.** The panel fixes what can be imputed (HRC MAF floor ~5e-4; TOPMed imputes far rarer) AND the build (HRC/1000G-P3 are GRCh37; TOPMed/1000G-NYGC are GRCh38). Ancestry match > panel size: a bigger ancestry-mismatched panel imputes WORSE because there are no matching haplotypes to copy.
2. **Strand/allele harmonization is a made-once GATE, and the silent corruptor.** Align every study variant's alleles to the panel REF/ALT, fix strand, and DROP unresolvable palindromes (A/T, C/G at MAF>0.4). A flipped palindrome or a build mismatch flips BETA and cancels/manufactures signal WITHOUT erroring — caught only on the AF-concordance plot (points on the y=1-x anti-diagonal), not by a crash.
3. **QC runs before association and in a defensible internal order; the engine is chosen by structure, not convenience.** Variant missingness before sample missingness; HWE in CONTROLS only (a true non-additive risk variant deviates in cases); impute cases+controls TOGETHER; exclude long-range-LD regions before PCA. PC-covariate GLM is valid only for unrelated, continuous-ancestry cohorts — relatedness/structure needs an LMM (regenie/SAIGE/BOLT) with LOCO, because PCs remove a mean shift, not a covariance.
4. **Refuse the bare lambda.** lambda_GC alone cannot separate polygenicity from confounding — read it WITH the LDSC intercept (~1 = polygenic, fine; >1 = confounding). Do not reflexively genomic-control (it over-corrects true signal). Filter INFO/R2 MAF-stratified (a flat cutoff is a hidden rare-variant filter).

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Genome build + imputation panel (ancestry-matched) | What can be imputed + the build; ancestry mismatch imputes worse regardless of size |
| Strand/allele harmonization | A flipped palindrome / build mismatch flips BETA silently; drop MAF>0.4 palindromes |
| Ancestry/LD reference for structure | Flows into PCA and (if LMM) the GRM; must match the cohort |
| Association engine (by structure) | PC-GLM (unrelated) vs LMM+LOCO (related/structured); SPA/Firth for imbalance/low MAC |

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
    --pheno phenotypes.txt \
    --make-bed \
    --out study   # load --pheno at conversion so PHENO1 is in the .fam for the controls-only HWE gate below
```

QC order is critical (-> population-genetics/plink-basics): variant missingness runs FIRST, in its own invocation, so a sample is not dropped for missingness driven by variants slated for removal.

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

The engine choice is set by sample structure, not convenience (-> population-genetics/association-testing). PC-covariate GLM (below) is valid only for unrelated samples whose confounding is continuous ancestry; any related, family-based, or fine-scale-structured cohort needs a linear mixed model with leave-one-chromosome-out, because PCs cannot remove a covariance structure.

| Situation | Engine | Hand off to |
|-----------|--------|-------------|
| Unrelated, continuous-ancestry confounding only | PLINK2 `--glm` (PC covariates) | population-genetics/association-testing |
| Relatedness / family / fine-scale structure | LMM (regenie / SAIGE / BOLT-LMM) with LOCO | population-genetics/association-testing |
| Biobank binary, case:control worse than ~1:10 or low MAC | SAIGE (SPA) or regenie (`--firth`/`--spa`) | population-genetics/association-testing |
| Rare-variant, aggregate by gene | burden/SKAT/SKAT-O via SAIGE-GENE+/regenie | population-genetics/rare-variant-association |
| Threshold | 5e-8 is EUR-common-array folklore; ~1e-8-3e-8 (African LD), ~5e-9 (WGS/rare) | population-genetics/association-testing |

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
# comment.char='' is REQUIRED: PLINK2's header starts with '#CHROM', which the default
# comment.char='#' would eat (dropping the header and the first variant).
results <- read.table('gwas_results.PHENO.glm.logistic.hybrid', header = TRUE, comment.char = '')
results <- results[!is.na(results$P),]
# qqman's manhattan needs a NUMERIC chromosome: recode X/Y/MT (e.g. X->23) or filter to autosomes first.
results$X.CHROM <- suppressWarnings(as.integer(sub('^chr', '', results$X.CHROM)))
results <- results[!is.na(results$X.CHROM),]

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

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| "0 variants matched" the panel | chr-naming (`1` vs `chr1`) or build mismatch | `bcftools annotate --rename-chrs`; match study build to the panel before any check |
| Flipped BETA / cancelled or manufactured signal | Strand-flip / allele-swap at intermediate-frequency palindromes | Run the Rayner check + AF-concordance FreqPlot; drop MAF>0.4 palindromes |
| Case-control association that is imputation artifact | Cases and controls imputed SEPARATELY | Impute jointly; separate imputation makes error differ between arms |
| Disproportionate rare-variant loss | Flat INFO/R2 cutoff | MAF-stratified R2 thresholds + a MAF floor |
| A PC tracks an inversion, not ancestry | Long-range-LD region kept before pruning | Exclude MHC/8p23.1/17q21.31/LCT BEFORE `--indep-pairwise` |
| Elevated lambda | Polygenicity OR confounding (lambda alone cannot separate) | Read the LDSC intercept (~1 = polygenic, fine); do NOT reflexively genomic-control |
| Residual structure after PCs | Relatedness/fine-scale structure (a covariance PCs cannot remove) | LMM (regenie/SAIGE/BOLT) with LOCO |
| Inflation at low MAC / extreme case:control ratio | Score/Wald anti-conservative | SPA (SAIGE) or Firth |
| chrX mis-analyzed | chrX coded as an autosome | split-PAR, sex covariate, explicit X handling |

## References

- McCarthy S, Das S, Kretzschmar W, et al (2016) A reference panel of 64,976 haplotypes for genotype imputation. *Nature Genetics* 48:1279-1283. DOI 10.1038/ng.3643. (HRC.)
- Taliun D, Harris DN, Kessler MD, et al (2021) Sequencing of 53,831 diverse genomes from the NHLBI TOPMed Program. *Nature* 590:290-299. DOI 10.1038/s41586-021-03205-y. (TOPMed diversity.)
- Sheng X, Xia L, Cahoon JL, et al (2023) Inverted genomic regions between reference genome builds. *HGG Advances* 4:100159. DOI 10.1016/j.xhgg.2022.100159. (liftover/BBIS strand danger.)
- Mbatchou J, Barnard L, Backman J, et al (2021) Computationally efficient whole-genome regression for quantitative and binary traits. *Nature Genetics* 53:1097-1103. DOI 10.1038/s41588-021-00870-7. (regenie.)

## Complete Pipeline Script

```bash
#!/bin/bash
set -e

INPUT_VCF="genotypes.vcf.gz"
PHENO="phenotypes.txt"
OUTDIR="gwas_results"

mkdir -p ${OUTDIR}

# Step 1: Convert, then QC in order (variant missingness, then sample, then MAF + controls-only HWE).
plink2 --vcf ${INPUT_VCF} --pheno ${PHENO} --make-bed --out ${OUTDIR}/raw   # --pheno embeds PHENO1 in .fam for the controls-only HWE gate
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
- workflows/causal-genomics-pipeline - Downstream: the emitted sumstats (CHR/POS/EA/OA/EAF/BETA/SE/P/N, build documented) are its input contract
