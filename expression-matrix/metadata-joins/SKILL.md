---
name: bio-expression-matrix-metadata-joins
description: Aligns sample metadata with count matrices and constructs design matrices for downstream DE, handling the alphabetical-reference-level trap (relevel BEFORE DESeq), LRT reduced-model rules, the interaction-term resultsNames trap, continuous-covariate scaling and splines, repeated measures via duplicateCorrelation or dream, high-cardinality categorical pseudo-singular designs, sample swap detection via XIST/RPS4Y1 expression and somalier/NGSCheckMate genotypes, SABV (sex-as-biological-variable) mandate, Simpson's-paradox collapsing of technical replicates, and the `~ 0 + group` parameterization for clean contrasts. Use when building a design matrix, troubleshooting reversed fold-change direction, encoding paired or repeated-measures designs, detecting sample swaps, deciding sex-as-covariate, or aggregating technical replicates.
tool_type: mixed
primary_tool: pandas
---

## Version Compatibility

Reference examples tested with: pandas 2.2+, DESeq2 1.42+, edgeR 4.0+, limma 3.58+, variancePartition / dream 1.32+, somalier 0.2.18+ (CLI), pyensembl 2.3+, anndata 0.10+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metadata Joins

**"Align my sample metadata with my count matrix and build a design"** -> Reconcile sample identifiers, validate the join, set factor levels explicitly to control fold-change direction, encode the experimental structure (paired, repeated measures, interaction) in the design formula, and detect swaps before downstream DE.

## The Single Most Important Modern Insight -- Alphabetical reference levels invert fold changes silently

DESeq2 picks the reference level alphabetically if not told otherwise. With `condition = c('Treated', 'Untreated')`, `T < U`, so `Treated` becomes the reference. The reported `log2FoldChange` is then `Untreated vs Treated` -- **the opposite of what the methods section says**. No error; the volcano plot looks plausible; the gene list is correct but with reversed sign.

```r
dds$condition <- relevel(dds$condition, ref = 'Untreated')
coldata$condition <- factor(coldata$condition, levels = c('Untreated', 'Treated'))
```

Set BEFORE `DESeq()`; relevel-then-DESeq again to take effect.

A second insight: in interaction designs `~ genotype * treatment`, `results(dds, name='treatment_drug_vs_vehicle')` returns the drug effect IN THE WT REFERENCE only, not the marginal/average effect. The interaction coefficient `genotypeKO.treatmentdrug` is the DIFFERENCE in drug effect between KO and WT (the difference of differences), NOT the drug effect in KO. The cleanest fix is to use `~ 0 + group` with `group = paste(genotype, treatment, sep='_')` so every comparison of interest is a single contrast.

## Algorithmic Taxonomy

| Pattern | Encoding | Tests | Caveat |
|---------|----------|-------|--------|
| Simple two-group | `~ condition` | `results(name='condition_treated_vs_control')` | Set reference level explicitly |
| With known batch | `~ batch + condition` | `results(name='condition_...')` | Variable of interest LAST is the convention; not required |
| Paired (pre/post per subject) | `~ subject + condition` | `results(name='condition_...')` | Pairing FIRST; subject absorbs baseline variability |
| Repeated measures (>2 time per subject) | DREAM `~ condition + (1|subject)` | `topTable(coef='condition')` | Mixed model; uses voom + lmer |
| Interaction (2x2) | `~ A * B` (expands to A + B + A:B) | `results(name='A.B')` for diff-in-diff; combined factor for cleaner contrasts | resultsNames trap |
| Multi-group all pairwise | `~ 0 + group` + `makeContrasts` | Contrasts named directly | DESeq2 needs intercept; works for edgeR/limma |
| Multi-level "any change" | LRT `reduced = ~ 1` | `results(dds)` from LRT fit | padj is omnibus; LFC is one specific level |

## Decision Tree by Scenario

| Scenario | Recommended approach |
|----------|---------------------|
| Tumor / normal paired by patient | `~ patient + tissue`; pairing FIRST -- absorbs subject variability |
| Pre / post drug, same patient (2 time points) | `~ subject + condition` |
| Longitudinal, 4+ time points per subject | DREAM mixed model with `(1 \| subject)` random effect |
| Known batch (sequencing run, library prep date) | `~ batch + condition`; do NOT subtract |
| Many batches (>10 levels, n=20-40) | `~ condition + (1 \| batch)` via DREAM; or aggregate batches |
| Continuous covariate (age, RIN) | Center first; include linearly OR via `ns(x, df=3)` if non-linear |
| Mixed-sex cohort | Include sex unless sex-specific; report sex-stratified sensitivity |
| Multi-group, want pairwise contrasts | `~ 0 + group` + `makeContrasts` |
| Interaction question (does effect of A differ across B?) | Either `~ A * B` (handle resultsNames trap) or combined factor `~ 0 + group` |
| Many technical reps within bio reps | `duplicateCorrelation` (limma) OR aggregate via `collapseReplicates` (DESeq2) |
| Cohort >=20 samples | Run somalier or NGSCheckMate genotype check at the matrix-build step |

## Basic Join

**Goal:** Align count matrix columns with metadata rows; remove samples present in only one source; verify alignment before downstream use.

**Approach:** Intersect sample identifiers; reorder both data sources; assert match.

```python
import pandas as pd

counts = pd.read_csv('counts.tsv', sep='\t', index_col=0)
metadata = pd.read_csv('metadata.csv', index_col=0)

common = counts.columns.intersection(metadata.index)
only_counts = set(counts.columns) - set(metadata.index)
only_meta = set(metadata.index) - set(counts.columns)
if only_counts: print(f'In counts not metadata: {only_counts}')
if only_meta: print(f'In metadata not counts: {only_meta}')

counts = counts[common]
metadata = metadata.loc[common]
assert all(counts.columns == metadata.index)
```

```r
common <- intersect(colnames(counts), rownames(coldata))
counts <- counts[, common]
coldata <- coldata[common, , drop = FALSE]
stopifnot(all(colnames(counts) == rownames(coldata)))
```

When sample names differ by formatting (underscore vs dash, BAM suffix, case), try systematic transformations -- replace `_` <-> `-`, strip `.bam`, lower-case, take prefix before `_` -- before giving up.

## Reference Level (Critical -- Set Before DESeq)

**Goal:** Pick the baseline against which fold changes are reported, controlling sign direction.

**Approach:** `relevel()` or construct the factor with explicit `levels=`. Set BEFORE `DESeq()` runs.

```r
coldata$condition <- factor(coldata$condition, levels = c('control', 'treated'))

coldata$condition <- relevel(coldata$condition, ref = 'control')
```

```python
metadata['condition'] = pd.Categorical(metadata['condition'],
                                        categories=['control', 'treated'],
                                        ordered=True)
```

With factor levels explicit, the LFC reads `treated / control` -- treated up means LFC > 0.

A reminder for edgeR / limma users: `makeContrasts(Treated - Control)` is explicit and immune to the reference-level trap. Same caution applies less because the user names the contrast.

## Paired Designs

```r
design = ~ patient + tissue
```

Pairing variable FIRST (convention). Patient absorbs inter-subject baseline variability, dramatically increasing power for the tissue effect.

For DESeq2:
```r
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ patient + tissue)
dds <- DESeq(dds)
res <- results(dds, name = 'tissue_tumor_vs_normal')
```

For edgeR:
```r
design <- model.matrix(~ patient + tissue, coldata)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 'tissuetumor')
```

Common mistake: writing `~ tissue + patient`. Numerically the model is the same; the convention of pairing-first improves readability and matches the natural mental model.

## Interaction Terms -- the resultsNames Trap

```r
design = ~ genotype + treatment + genotype:treatment
dds <- DESeqDataSetFromMatrix(counts, coldata, design = design)
dds <- DESeq(dds)
resultsNames(dds)
```

Output names (after relevel):
```
"Intercept"
"genotype_KO_vs_WT"
"treatment_drug_vs_vehicle"
"genotypeKO.treatmentdrug"
```

| Question | Wrong answer | Right answer |
|----------|--------------|--------------|
| Drug effect averaged over genotypes | `results(name='treatment_drug_vs_vehicle')` -- this is drug effect in WT only | Combined factor `~ 0 + group`; or contrast that explicitly averages |
| Is drug effect different between genotypes? | n/a | `results(name='genotypeKO.treatmentdrug')` -- the interaction coefficient IS the difference of differences |
| Drug effect in KO | `results(name='treatment_drug_vs_vehicle')` -- WRONG, this is drug in WT | `results(contrast=list(c('treatment_drug_vs_vehicle', 'genotypeKO.treatmentdrug')))` (sum of main + interaction) |

The cleaner alternative for designs with many contrasts of interest:

```r
coldata$group <- factor(paste(coldata$genotype, coldata$treatment, sep = '_'))
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ 0 + group)
dds <- DESeq(dds)

res_drug_in_ko <- results(dds, contrast = c('group', 'KO_drug', 'KO_vehicle'))
res_drug_in_wt <- results(dds, contrast = c('group', 'WT_drug', 'WT_vehicle'))
res_diff       <- results(dds, contrast = list(c('groupKO_drug', 'groupWT_vehicle'),
                                                c('groupKO_vehicle', 'groupWT_drug')))
```

`~ 0 + group` parameterization is the long-standing edgeR / limma recommendation (Smyth and Robinson User's Guides) for any design with multiple pairwise contrasts of interest.

## LRT and the Reduced Model

```r
dds <- DESeq(dds, test = 'LRT', reduced = ~ batch)
```

Reduced model drops the term being tested. With `design = ~ batch + condition` and `reduced = ~ batch`, the LRT tests condition. With interaction designs:

```r
dds <- DESeq(dds, test = 'LRT', reduced = ~ genotype + treatment)
```

(Tests the interaction.)

The reduced model must be NESTED in the full model (every term in reduced must appear in full).

## Continuous Covariates

| Encoding | Assumption | When |
|----------|------------|------|
| Linear (`+ age`) | log-expression linear in age | Limited range, biologically linear |
| Centered linear (`+ I(age - mean(age))`) | As linear, interpretable intercept | Standard for age-RIN-day covariates |
| Natural spline (`+ ns(age, df=3)`) | Smooth nonlinear | Wide age range with non-monotonic effects |
| Polynomial (`+ poly(age, 2)`) | Quadratic; orthogonal polynomials | Limited use; splines usually better |

```r
coldata$age_c <- coldata$age - mean(coldata$age)
design = ~ age_c + RIN + condition

library(splines)
design = ~ ns(age, df = 3) + condition
```

DO NOT include library size as a covariate -- it is handled by size factors / normalization factors internally.

## Repeated Measures -- duplicateCorrelation vs dream

**Goal:** Correctly model within-subject correlation when the same subject contributes multiple samples.

**Approach:** For technical reps within bio reps OR paired pre/post: `duplicateCorrelation` (limma) is adequate. For >2 time points per subject or random slopes: DREAM (`variancePartition`).

```r
library(limma)
library(edgeR)

v <- voom(y, design)
corfit <- duplicateCorrelation(v, design, block = coldata$donor)
v <- voom(y, design, block = coldata$donor, correlation = corfit$consensus)
corfit <- duplicateCorrelation(v, design, block = coldata$donor)
fit <- lmFit(v, design, block = coldata$donor, correlation = corfit$consensus)
fit <- eBayes(fit, robust = TRUE)
```

The double pass is intentional: estimate correlation, re-voom with correlation, re-estimate, fit. Limma's `duplicateCorrelation` assumes ONE within-subject correlation across all genes -- approximation.

For proper per-gene mixed models:

```r
library(variancePartition)

form <- ~ condition + (1 | donor)
vobj <- voomWithDreamWeights(y, form, coldata)
fitmm <- dream(vobj, form, coldata)
fitmm <- eBayes(fitmm)
tt <- topTable(fitmm, coef = 'condition')
```

See `differential-expression/timeseries-de` for full longitudinal designs.

Before committing to a design, `variancePartition::fitExtractVarPartModel(vobj, form, coldata)` quantifies the fraction of expression variance explained by each covariate, gene-by-gene. If a candidate "nuisance" covariate explains <1% of variance across most genes, it can usually be dropped; if a known biological factor explains <5% and isn't of direct interest, model it as a random effect rather than fixed.

## Sample Swap Detection (Mandatory for Cohort >= 20)

### Sex check via XIST and chrY expression

**Goal:** Detect mislabeled samples by checking that gene-expression sex matches reported sex.

**Approach:** Compare XIST expression (high in XX, low/absent in XY) against chrY-gene expression (DDX3Y, RPS4Y1, UTY, KDM5D, EIF1AY -- high in XY, absent in XX).

```python
import pandas as pd

def sex_check(counts, metadata, sex_column='sex'):
    y_genes = ['DDX3Y', 'RPS4Y1', 'UTY', 'KDM5D', 'EIF1AY']
    y_avail = [g for g in y_genes if g in counts.index]
    if 'XIST' not in counts.index or not y_avail:
        return None
    predicted = pd.Series('unknown', index=counts.columns)
    predicted[counts.loc[y_avail].sum() > counts.loc['XIST']] = 'M'
    predicted[counts.loc['XIST'] > counts.loc[y_avail].sum()] = 'F'
    if sex_column in metadata.columns:
        mis = predicted != metadata[sex_column]
        if mis.any():
            print(f'SEX MISMATCHES: {list(metadata.index[mis])}')
    return predicted
```

```r
sex_check <- function(counts, coldata, sex_col = 'sex') {
    y_genes <- c('DDX3Y', 'RPS4Y1', 'UTY', 'KDM5D', 'EIF1AY')
    y_expr <- colSums(counts[intersect(y_genes, rownames(counts)), , drop = FALSE])
    predicted <- ifelse(y_expr > counts['XIST', ], 'M', 'F')
    mis <- predicted != coldata[[sex_col]]
    if (any(mis)) cat('Sex mismatches:', colnames(counts)[mis], '\n')
    predicted
}
```

CAVEAT: tumors with X loss, sex chromosome aneuploidies, HeLa (XXX with mixed inactivation) muddy this. Genotype-based methods are more robust.

### Genotype-based: somalier and NGSCheckMate

```bash
somalier extract -d extracted/ --sites sites.GRCh38.vcf.gz \
    -f reference.fa sample.bam
somalier relate --infer extracted/*.somalier
```

```bash
ncm_fastq.py -l fastq_list.txt -O outdir -bed common_sites.bed
```

Somalier (Pedersen et al. 2020 *Genome Med* 12:62) extracts a few thousand SNP sketches per sample (sub-second per sample) and computes pairwise relatedness from BAM/CRAM/VCF. NGSCheckMate (Lee et al. 2017 *NAR* 45:e103) computes VAF correlation across a common-SNP panel; works on FASTQ/BAM/VCF including RNA-seq.

For any cohort >=20 samples, run one of these at the matrix-build step. Catching a swap in raw data is cheap; finding it after DE is expensive.

## SABV -- Sex as Biological Variable

NIH 2016+ requires sex consideration in vertebrate animal and human studies. Mauvais-Jarvis F et al. 2020 *Lancet* 396:565 reviews effect-size differences across diseases.

Practical implication:
- Always include sex in the metadata, even when not in the model.
- For sex-balanced cohorts, `~ sex + condition` rarely hurts and captures real biology.
- For sex-confounded cohorts (all-male disease cohort vs mixed-sex control), can't rescue but documents the limitation.
- For chrX/chrY analyses, sex MUST be in the model OR the analysis is uninterpretable.
- Report DE counts overall AND sex-stratified.

## Simpson's Paradox -- Collapsing Technical Replicates

**Goal:** Aggregate technical replicates from the same subject correctly; never treat them as independent biological replicates.

**Approach:** Sum (not average) technical replicates of the same subject BEFORE downstream DE.

```r
library(DESeq2)
dds_collapsed <- collapseReplicates(dds, groupby = dds$subject)
```

```python
counts_per_subject = counts.T.groupby(metadata['subject']).sum().T
metadata_per_subject = metadata.drop_duplicates(subset='subject').set_index('subject')
```

Why sum and not average? Reads add. Two technical replicates yielding 1M reads each are equivalent to one library yielding 2M reads. Averaging would understate the effective library size.

Treating technical replicates as independent biological samples is the cardinal sin: it inflates the apparent sample size and deflates standard errors. With 4 patients x 3 tech reps = 12 samples, naive DE assumes 12 independent observations; the truth is closer to 4. p-values are compressed ~3x.

## High-Cardinality Categorical Covariates

**Goal:** Handle batch / lane / well covariates with many levels without making the design matrix singular.

**Approach:** Aggregate to fewer levels, model as random effect via DREAM, or drop if confounded.

```r
ct <- table(coldata$condition, coldata$batch)
ct

ad <- alias(model.matrix(~ batch + condition, coldata))
ad$Complete
```

For a batch with 30 levels in n=40 samples: 30 batch coefficients + condition + intercept = 32 parameters for 40 observations. Symptoms: `Matrix not positive definite` (DESeq2), degenerate p-values (limma).

Fixes:
- Aggregate batches (sequencing run, sequencing pool, library prep date as proxy).
- Random effect: `~ condition + (1 | batch)` via DREAM borrows information across batch levels via shrinkage.
- Drop the covariate if confounded with condition (`alias()` reveals collinearity).

## ~ 0 + group Parameterization

```r
design_default <- model.matrix(~ group, coldata)
# columns: (Intercept), groupB, groupC  -- A is reference

design_nointercept <- model.matrix(~ 0 + group, coldata)
# columns: groupA, groupB, groupC  -- each column is mean of that group
```

With `~ 0 + group`, every contrast reads as `B - A`:

```r
library(limma)

con <- makeContrasts(BvsA = groupB - groupA,
                     CvsA = groupC - groupA,
                     BvsC = groupB - groupC,
                     levels = design_nointercept)
fit <- glmQLFit(y, design_nointercept, robust = TRUE)
qlf <- glmQLFTest(fit, contrast = con[, 'BvsA'])
```

DESeq2 needs an intercept internally, so `~ 0 + group` works directly with edgeR/limma but DESeq2 uses `contrast=` to achieve the same effect.

## Create DESeq2 / edgeR / AnnData Containers

```r
dds <- DESeqDataSetFromMatrix(as.matrix(counts), coldata, design = ~ batch + condition)

y <- DGEList(counts = as.matrix(counts), group = coldata$condition)
y$samples <- cbind(y$samples, coldata)

adata <- ad.AnnData(X = t(as.matrix(counts)), obs = coldata, var = data.frame(row.names = rownames(counts)))
```

AnnData convention is cells (samples) in rows -- transpose from the typical R genes-in-rows convention.

## Per-Method Failure Modes

### Fold-change direction reversed

**Trigger:** Methods says "treated vs control"; published volcano shows expected up-genes on the LEFT.

**Mechanism:** Factor levels left at alphabetical default; `c('Treated','Untreated')` -> `T < U` -> Treated is reference -> LFC is Untreated/Treated.

**Symptom:** Known up-regulated genes appear down; reviewer questions direction.

**Fix:** `relevel(coldata$condition, ref = 'Untreated')` BEFORE `DESeq()`. Re-run.

### Interaction coefficient mistaken for main effect

**Trigger:** `~ A * B` design; `results(name='B_drug_vs_vehicle')` reported as "drug effect"; reviewer asks about genotype-specific effect.

**Mechanism:** With interaction, `B_drug_vs_vehicle` is drug effect IN THE A REFERENCE LEVEL only, not averaged across A.

**Symptom:** Drug effect doesn't match the marginal estimate from a separate `~ B`-only fit.

**Fix:** Use `~ 0 + group` with combined factor; OR extract per-stratum results explicitly using contrasts that sum main + interaction.

### Pseudoreplication -- 12 samples, only 3 subjects

**Trigger:** 3 subjects x 4 conditions = 12 samples; vanilla DESeq2 with `~ condition`; many DE genes.

**Mechanism:** Same subject contributes multiple observations; not independent. Effective sample size for testing condition is ~3, not 12.

**Symptom:** p-value histogram anti-conservative; replication low.

**Fix:** Include subject in design (`~ subject + condition`) OR use DREAM with random subject.

### Sample swap caught in DE results, not metadata

**Trigger:** PCA shows "control" sample clustering with treated; investigation reveals it was mislabeled at thaw.

**Mechanism:** Manual sample tracking is error-prone; cohort >=20 inevitably has swaps.

**Symptom:** One sample dramatically off its group cluster; DE gene list dominated by sample-specific effects.

**Fix:** Run somalier or NGSCheckMate at the matrix-build step, BEFORE DE. Catching a swap early is cheap; finding it post-DE is expensive.

### Sex confounded, chrY genes dominate top hits

**Trigger:** Mixed-sex cohort; sex not in design; PCA shows clear sex split on PC1.

**Mechanism:** Sex distribution differs between groups; "treatment effect" partially captures sex.

**Symptom:** Top DE genes are DDX3Y, RPS4Y1, UTY (chrY) and XIST -- not biology of interest.

**Fix:** Add sex to design (`~ sex + condition`). For sex-specific analyses, stratify and report each separately.

### duplicateCorrelation single-pass

**Trigger:** limma user wrote a one-pass `duplicateCorrelation` + lmFit; QC reviewer asks why no re-voom.

**Mechanism:** The proper pattern is: voom -> dupCor -> re-voom WITH correlation -> dupCor again -> lmFit. The first voom doesn't know about block structure; re-voom with correlation gets better weights.

**Symptom:** Slightly inflated DE counts vs the two-pass pattern.

**Fix:** Implement the two-pass pattern per the limma User's Guide (section 9.7).

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| Sample names don't match between counts and metadata | Underscore/dash inconsistency, BAM suffix, case | `fuzzy_match_samples()` or manual normalization; report what failed |
| DESeq2 design not full rank | Confounded covariates | `alias(design)$Complete` to identify; aggregate or drop |
| `Matrix not positive definite` | High-cardinality batch with few samples | Aggregate batches or use random effect via DREAM |
| LFC direction reversed | Alphabetical reference level | `relevel()` before `DESeq()` |
| `results(name='...')` returns drug effect in WT only | Interaction design and naming trap | Use combined factor `~ 0 + group`; or contrast summing main + interaction |
| Inflated DE list with 12 samples from 3 subjects | Pseudoreplication | Include subject; or use DREAM |
| Sex effect appears as treatment effect | Sex not in design | Add `~ sex + condition` |
| Sample distance heatmap shows mixing groups | Likely swap | Run somalier or NGSCheckMate |

## References

- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15(12):550. doi:10.1186/s13059-014-0550-8
- Robinson MD, McCarthy DJ, Smyth GK. 2010. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* 26(1):139-140. doi:10.1093/bioinformatics/btp616
- Ritchie ME et al. 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 43(7):e47. doi:10.1093/nar/gkv007
- Smyth GK, Michaud J, Scott HS. 2005. Use of within-array replicate spots for assessing differential expression in microarray experiments. *Bioinformatics* 21(9):2067-2075. doi:10.1093/bioinformatics/bti270
- Hoffman GE, Roussos P. 2021. dream: powerful differential expression analysis for repeated measures designs. *Bioinformatics* 37(2):192-201. doi:10.1093/bioinformatics/btaa687
- Hoffman GE, Schadt EE. 2016. variancePartition: interpreting drivers of variation in complex gene expression studies. *BMC Bioinformatics* 17:483. doi:10.1186/s12859-016-1323-z
- Pedersen BS, Bhetariya PJ, Brown J, Kravitz SN, Marth GT, Jensen RL, Bronner MP, Underhill HR, Quinlan AR. 2020. Somalier: rapid relatedness estimation for cancer and germline studies using efficient genome sketches. *Genome Medicine* 12:62. doi:10.1186/s13073-020-00761-2
- Lee S, Lee S, Ouellette S, Park W-Y, Lee EA, Park PJ. 2017. NGSCheckMate: software for validating sample identity in next-generation sequencing studies within and across data types. *Nucleic Acids Res* 45(11):e103. doi:10.1093/nar/gkx193
- Mauvais-Jarvis F, Bairey Merz N, Barnes PJ et al. 2020. Sex and gender: modifiers of health, disease, and medicine. *Lancet* 396(10250):565-582. doi:10.1016/S0140-6736(20)31561-0

## Related Skills

- counts-ingest - Building count matrices before metadata join
- gene-id-mapping - Annotating result tables with symbols
- normalization - Reference for downstream normalization choice
- sparse-handling - AnnData obs metadata convention
- differential-expression/deseq2-basics - Where the design formula matters; relevel; interactions
- differential-expression/edger-basics - edgeR design matrix conventions
- differential-expression/batch-correction - Design-inclusion vs subtraction; confounding
- differential-expression/timeseries-de - Repeated measures; DREAM mixed model details
- single-cell/preprocessing - scRNA-seq sample metadata handling
