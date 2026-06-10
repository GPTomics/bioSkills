---
name: bio-methylation-ewas-design
description: Designs and defends an epigenome-wide association study (EWAS) on 450K/EPIC array or bisulfite methylation - the layer deciding whether a hit is credible. Covers the confounding hierarchy (cell composition covariates as the dominant confounder, batch/Sentrix chip/array position, age/sex, smoking AHRR cg05575921, ancestry/mQTL, reverse causation), chip randomization (no-rescue theorem), surrogate variable analysis sva/SmartSVA, ComBat, RUVm, over-correction, genomic inflation lambda vs GWAS genomic control, BACON bias/inflation, genome-wide significance threshold 450K/EPIC, FWER vs FDR, pwrEWAS power, meta-analysis, EWAS Catalog/Atlas, methylation risk scores. Use when designing an EWAS, choosing a covariate set, randomizing a plate layout, interpreting lambda, applying BACON, setting a threshold, powering a study, or using an MRS. For the per-site test see differential-cpg-testing; for cell fractions see cell-type-deconvolution; for causal mQTL orientation see causal-genomics/mendelian-randomization.
tool_type: mixed
primary_tool: meffil
---

## Version Compatibility

Reference examples tested with: meffil 1.3+, sva 3.50+, bacon 1.30+, limma 3.58+, missMethyl 1.36+, pwrEWAS 1.16+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Two versions decide everything downstream. The ARRAY (450K vs EPIC v1 ~850K vs EPIC v2 ~935K) sets the CpG universe and therefore the effective number of tests and the genome-wide threshold - confirm against the current manifest before quoting a threshold. The CELL-TYPE REFERENCE panel must match the tissue and age (cord blood, adult blood, and tumor need different references); a mismatched reference silently biases the cell-fraction covariates. BACON slot names, pwrEWAS argument names, and sva method options have shifted across Bioconductor releases - check `?function` on the installed build.

# EWAS Design

**"Run an EWAS for my phenotype"** -> Build the covariate model and the chip randomization BEFORE the per-CpG test, because the confounders are larger than the signal - the p-value is the last and least interesting thing.
- R: `meffil.ewas(beta, variable, covariates=data.frame(age, sex, cellprops, chip, plate))` then BACON-correct the test statistics.

Scope: the design-and-inference layer of an EWAS - confounding, randomization, batch/SVA/RUV correction, genomic inflation/BACON, genome-wide thresholds, power, meta-analysis, replication, and methylation risk scores. The per-CpG test mechanics (beta vs M-value, limma) -> differential-cpg-testing. Cell-fraction estimation algorithms (Houseman/IDOL/EpiDISH) -> cell-type-deconvolution. Normalization internals (funnorm/noob) -> array-preprocessing. mQTL Mendelian randomization for causal orientation -> causal-genomics/mendelian-randomization. MRS weight-learning (penalized regression, nested CV) -> the machine-learning category.

## The Single Most Important Modern Insight -- An EWAS Hit Is a Cell-Composition Difference Until Proven Otherwise

The defining feature of the field is that the confounders are LARGER than the signal: a true effect is typically a sub-2% absolute methylation difference, while cell mix, batch, age, sex, smoking, and ancestry each move methylation 10-50%. An EWAS is therefore won or lost BEFORE the per-CpG test - in the chip randomization, the covariate set, and the replication cohort. Four corollaries, each of which a common misuse violates:

1. **Composition before regulation.** Whole blood is a leukocyte mixture; a CpG can be 90% methylated in granulocytes and 10% in lymphocytes. Any phenotype that shifts the mixture (age, inflammation, infection, stress, smoking, most diseases) produces methylation differences that are composition artifacts, not within-cell regulatory change. Cell-proportion covariates are mandatory, not optional; the default suspicion for any unadjusted hit is "a blood-count difference."
2. **Confounder > signal.** Design (randomization, matching) is the primary defense; statistical adjustment is the backstop. No analysis rescues a design that confounded chip/plate with phenotype at the bench.
3. **Blood is the lamppost, not the keys.** Blood is convenient; the disease tissue (brain, adipose, tumor) is usually inaccessible and weakly correlated (blood-brain mean r ~0.15, Hannon 2015). A blood EWAS for a brain trait tests a confounded surrogate, and a cross-sectional design cannot orient cause vs consequence (reverse causation).
4. **Discovery is a hypothesis; replication is the finding.** With tiny effects and pervasive confounding, a single-cohort genome-wide-significant CpG is a lead, not a result. The EWAS Catalog/Atlas exist so a "novel" hit can be checked against the generic age/smoking/cell-comp CpGs everyone finds.

Organize the analysis around defending these four, not around listing `sva`/`ComBat`/`bacon` functions.

## The Confounding Hierarchy (ordered by damage, not convenience)

| Rank | Confounder | Why it dominates | Defense | Owner |
|------|-----------|------------------|---------|-------|
| 1 | Cell composition | Each leukocyte has a radically different methylome; any mixture shift is a fake signal | Estimate cell fractions and ALWAYS adjust | -> cell-type-deconvolution (estimation); decide here |
| 2 | Batch / Sentrix chip / array position / plate / scan date | EPIC runs 8 samples per BeadChip (450K: 12); position and chip imprint methylation | Randomize at design time; then SVA/RUVm/ComBat + chip/position covariates | here + array-preprocessing |
| 3 | Age and sex | Age is the most reproducible methylome correlate (the clock); sex drives X/Y and X-inactivation | Always include; analyze sex chromosomes separately; sex is a mislabel QC check | here |
| 4 | Smoking | Largest reproducible blood exposure signature; confounds half of all health phenotypes | Adjust (methylation smoking score > self-report); AHRR cg05575921 is the positive control | here |
| 5 | Genetic / ancestry / mQTL | Many CpGs are under genetic control; population structure confounds like GWAS | Genetic PCs; analyze within ancestry; drop SNP-affected/cross-reactive probes | here + array-preprocessing |
| 6 | Reverse causation / tissue relevance | Cross-sectional methylation may be a consequence, not a cause; blood != disease tissue | Prospective/longitudinal or MZ-discordant design; cross-tissue concordance | here; causal -> causal-genomics/mendelian-randomization |

Default covariate set for a blood EWAS: age, sex, cell proportions, chip/position (or SVs), and where relevant genetic PCs and a smoking score. The omission of any one is the most common EWAS error.

## Design Beats Correction -- The No-Rescue Theorem

If all cases ran on chip A and all controls on chip B, batch and phenotype are inseparable - ComBat will remove the batch AND the signal, or preserve a spurious one. No statistical adjustment recovers a confounded design (general principle: experimental-design/batch-design).

**Goal:** Make technical batch orthogonal to phenotype before any sample touches the array.

**Approach:** Randomize (or block) sample-to-chip, sample-to-position, and sample-to-plate assignment so case/control, age, and sex are balanced across chips; reserve no chip for one group. This single uncorrectable decision matters more than every analysis choice that follows.

```r
# Stratified randomization of samples to 8-position EPIC BeadChips (450K: 12), balancing case/control per chip
meta$chip <- NA_integer_
for (grp in split(seq_len(nrow(meta)), meta$case)) {
  meta$chip[sample(grp)] <- rep(seq_len(ceiling(length(grp) / 4)), each = 4, length.out = length(grp))
}
# Then interleave the two case groups across chips so no chip is single-group; verify balance:
with(meta, table(chip, case))   # every chip should hold both groups
```

## Batch / SVA / RUV Correction (and the over-correction trap)

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| sva | Leek & Storey 2007 *PLoS Genet* 3:e161 | latent surrogate variables built orthogonal to the variable of interest | reference-free soak-up of cell mix + unknown batch |
| SmartSVA | Chen 2017 *BMC Genomics* 18:413 | order-of-magnitude-faster SVA with explicit convergence | the de-facto EWAS SVA for large cohorts |
| ComBat | Johnson 2007 *Biostatistics* 8:118 | empirical-Bayes removal of a SPECIFIED batch (chip/plate) | known, labelled batch; protect biology via `mod=` |
| RUVm | Maksimovic 2015 *Nucleic Acids Res* 43:e106 | two-stage RUV-inverse using Illumina's ~600 negative control probes | array EWAS; estimates unwanted variation from control probes |
| funnorm | Fortin 2014 *Genome Biol* 15:503 | control-probe PCA at the normalization stage | fix unwanted variation BEFORE testing (-> array-preprocessing) |

The central tension: cell-composition and technical variation MUST be removed, but aggressive correction removes REAL biology when the unwanted variation overlaps the phenotype. Every tool above can erase a true effect.

**The positive-control check (the pragmatic referee).** Monitor a known signal through correction. In a blood EWAS containing smokers, smoking->cg05575921 (AHRR, ~18% hypomethylation, Joehanes 2016) MUST survive. If adding surrogate variables makes the QQ plot look clean BUT kills AHRR, the pipeline is over-corrected. A clean QQ with a dead positive control is a broken pipeline, not a good one. Do not keep adding SVs until lambda hits 1.0.

ComBat and limma/EWAS modeling run on M-values (logit of beta), not betas (betas are bounded [0,1] and heteroscedastic); effect sizes are reported back on the beta / delta-beta scale for interpretability. Pass the biological covariate to ComBat's `mod=` so it is protected. sva: pass the FULL model (`mod`, including the variable of interest) AND the null model (`mod0`) so SVs are orthogonal to the phenotype; choose the number with `num.sv(method='be')`.

## Genomic Inflation and BACON -- Why GWAS Intuition Fails

Lambda (the genomic inflation factor) = median observed chi-square / expected median. In GWAS, lambda > 1 signals stratification and is corrected by genomic control (divide all statistics by lambda). This reasoning is WRONG for EWAS:

1. EWAS statistics are routinely BOTH inflated AND deflated for non-GWAS reasons - residual cell-composition variation, the strong correlation among CpGs, un-modeled technical variation, and the fact that a strong exposure (smoking) genuinely associates with a large fraction of the genome (real signal that legitimately inflates lambda). A lambda of 1.2 may be real biology in a smoking EWAS and cell-composition leakage in an under-powered case-control study - the SAME number means different things.
2. Genomic control assumes a single multiplicative inflation on a mostly-null genome and ignores BIAS (a systematic mean shift off zero). EWAS violates all three assumptions.

**BACON (van Iterson 2017 *Genome Biol* 18:19)** fits a Bayesian Gaussian mixture to the observed test statistics, estimates the empirical-null distribution, and reports BOTH a bias (mean shift) AND an inflation (scale) - then standardizes statistics against that empirical null without assuming the genome is mostly null. Report lambda but correct with BACON; show QQ plots before and after. Caveat: BACON can over-deflate a highly polygenic exposure if the alternative component is large - apply it with the QQ plot, not as a reflex.

## Genome-Wide Thresholds and the FWER-vs-FDR Decision

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| 450K: P < 2.4e-7 | Saffari 2018 *Genet Epidemiol* 42:20 | empirical effective-test FWER 5%; ~210,000 effective tests (< probe count due to correlation) |
| EPIC v1: P < 9e-8 | Mansell 2019 *BMC Genomics* 20:366 | null-simulation FWER 5% for the ~850K EPIC array |
| Pragmatic 1e-7 | community | round number between the two array-specific values |
| FDR (BH) q < 0.05 | Benjamini-Hochberg | more powerful; discovery / exposure scans; report as a sensitivity layer |

Naive Bonferroni on the probe count is over-conservative because CpGs are correlated (co-methylation, shared mQTLs); raw p is anti-conservative because there are ~850K tests. Use the array-specific FWER threshold for the cross-study-comparable headline claim and the EWAS Catalog; use BH-FDR for hypothesis-generating discovery. BH's independence assumption is imperfect under probe correlation but robust to positive dependence. Region-level (DMR) multiple-testing is owned by dmr-detection; mechanical p.adjust is owned by differential-cpg-testing. EPIC v2 (~935K probes, renamed/dropped vs v1) may shift the threshold - verify against the manifest.

## Power and Meta-Analysis

**Goal:** Size an EWAS for a realistic (tiny) effect using site-specific methylation variance, not a single assumed sd.

**Approach:** pwrEWAS (Graw 2019 *BMC Bioinformatics* 20:218) simulates DNAm using empirical per-CpG variance from a reference dataset, so power reflects the real site-specific variance structure (bimodal sites near 0/1 behave differently from intermediate sites). Specify target delta-beta, sample size, and the genome-wide threshold; expect to need hundreds-to-thousands of samples for a 1-2% effect - which is why meta-analysis dominates.

EWAS meta-analysis is the field standard for power: each cohort runs an identical pre-specified pipeline (meffil enables this), BACON-corrects per cohort, uploads summary statistics, and a central site does fixed-effect inverse-variance pooling with heterogeneity testing (Cochran's Q, I^2). meffil's distinguishing feature is distributed normalization - cohorts normalize locally without sharing individual data, reducing meta-analysis heterogeneity - plus automated selection of the number of normalization PCs.

## Study Designs

| Design | Controls for | Cost |
|--------|--------------|------|
| Case-control (cross-sectional) | nothing inherently; adjust on age/sex/ancestry/cell-mix, randomize chips | reverse causation, composition confounding |
| Longitudinal / prospective | reverse causation (methylation measured before outcome); within-person change via mixed model | needs follow-up; repeated samples |
| MZ-discordant twin | genetic + shared-environment confounding by design (paired within-pair analysis) | discordant pairs are rare -> power-limited |
| Exposure EWAS | smoking is the template/positive-control | exposure misclassification (self-report) |
| Meta-analysis | low power of single cohorts | requires harmonized pipelines + BACON per cohort |

## Replication and Look-Up

Replication is the real significance bar. After discovery, triage every hit against both databases - a "novel" CpG that is in fact a top smoking/age/blood-cell CpG is almost certainly residual confounding.

- **EWAS Catalog** (Battram 2022 *Wellcome Open Res* 7:41; ewascatalog.org) - published associations at P < 1e-4 (so a Catalog "hit" is a lookup, not a genome-wide claim) plus de-novo EWAS; check whether a CpG was reported for any trait.
- **EWAS Atlas / Open Platform** (Li 2019 *Nucleic Acids Res* 47:D983) - curated associations with a trait-ENRICHMENT tool for interpreting a CpG set.

## Methylation Risk Scores (MRS)

An EWAS produces per-CpG associations; an MRS turns many of them into one number - a weighted CpG sum, MRS = sum_j w_j * beta_ij, with weights from a training EWAS or a penalized fit. It is the methylation analogue of a polygenic risk score (PRS), and an epigenetic clock is the age/health special case (-> epigenetic-clocks).

**The load-bearing distinction from a PRS.** A PRS sums germline variants - fixed at conception, antecedent, plausibly causal. An MRS sums methylation - modifiable, tissue/time-specific, and frequently a CONSEQUENCE of the exposure/trait rather than a cause. The methylation smoking score (Elliott 2014 *Clin Epigenetics* 6:4; generalized by Sugden 2019) does not predict a propensity to smoke; it MEASURES the footprint smoking left. So an MRS is reverse-causal-by-default: report it as a predictive BIOMARKER / objective exposure proxy, and reserve "risk" and "cause" for prospectively- or MR-supported claims. PRS = germline cause; MRS = state consequence.

The most useful design role is as a better-measured confounder: self-reported smoking is biased and coarse, so adjusting an EWAS for the methylation smoking score controls residual smoking confounding far better whenever the phenotype is smoking-correlated (caveat: if smoking is on the causal path, over-adjusting via the score removes real signal - the same confounder-vs-mediator tension as cell composition). Other DNAm scores exist for BMI/alcohol/education (McCartney 2018) and circulating proteins (EpiScores, Gadd 2022 *eLife* 11:e71802). Portability fails across array (450K vs EPIC drop CpGs), tissue, and ancestry - report how many score CpGs are present on the array. Defer MRS weight-learning (penalized regression, nested CV, calibration, leakage) to the machine-learning category; an MRS validated in its own training cohort is not validated.

## Per-Method Failure Modes

### EWAS without cell-composition covariates
**Trigger:** running the regression on whole-blood betas with no cell-fraction adjustment. **Mechanism:** any phenotype that shifts the leukocyte mixture moves methylation 10-50%. **Symptom:** the top hits are known granulocyte/lymphocyte-proportion CpGs. **Fix:** estimate cell proportions (-> cell-type-deconvolution) and always include them; treat any unadjusted hit as composition until shown otherwise.

### Chip/plate confounded with phenotype, then "fixed" by ComBat
**Trigger:** all cases on chip A, controls on chip B. **Mechanism:** batch and phenotype are inseparable. **Symptom:** ComBat removes the signal or fabricates one; no error. **Fix:** randomize sample-to-chip/position at design time - this is uncorrectable after the bench.

### GWAS genomic control applied to EWAS
**Trigger:** dividing all test statistics by lambda. **Mechanism:** assumes single multiplicative inflation on a mostly-null genome and ignores bias. **Symptom:** real polygenic signal (smoking) deflated or residual confounding under-corrected. **Fix:** use BACON to estimate empirical-null bias AND inflation.

### Over-correction with too many surrogate variables
**Trigger:** adding SVs/PCs until lambda hits 1.0. **Mechanism:** latent factors correlated with the phenotype deflate true signal. **Symptom:** clean QQ but the positive control (AHRR) is gone. **Fix:** monitor cg05575921 through correction; stop when adding SVs starts killing it.

### Bonferroni-on-probe-count or raw-p reporting
**Trigger:** 0.05/850000, or reporting nominal p. **Mechanism:** probe correlation makes Bonferroni too strict; ignoring 850K tests is too loose. **Symptom:** missed real hits or a flood of false ones. **Fix:** array-specific FWER threshold (2.4e-7 / 9e-8) for claims, BH-FDR for discovery.

### Blood EWAS interpreted as disease-tissue mechanism
**Trigger:** a blood EWAS for a brain/adipose/tumor phenotype read as tissue biology. **Mechanism:** blood-brain methylation r ~0.15. **Symptom:** hits with no plausible blood mechanism. **Fix:** justify tissue relevance (BECon/blood-brain tools); frame blood as a biomarker, not mechanism, unless the CpG is cross-tissue concordant.

### MRS treated as a germline-like risk score
**Trigger:** describing a DNAm score for a disease as "risk" like a PRS. **Mechanism:** methylation is modifiable and frequently downstream. **Symptom:** a "DNAm risk score" that is actually the footprint of disease/treatment/behavior. **Fix:** call it a predictive biomarker; orient cause vs consequence only with prospective or MR evidence.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| 450K genome-wide P < 2.4e-7 | Saffari 2018 *Genet Epidemiol* 42:20 | empirical effective-test FWER 5% (~210,000 effective tests) |
| EPIC genome-wide P < 9e-8 | Mansell 2019 *BMC Genomics* 20:366 | null-simulation FWER 5% for ~850K EPIC |
| BH-FDR q < 0.05 | Benjamini-Hochberg | discovery / exposure scans; report alongside FWER |
| AHRR cg05575921 ~18% hypomethylation in smokers | Joehanes 2016 *Circ Cardiovasc Genet* 9:436 | the canonical positive control; must survive correction |
| EWAS Catalog inclusion P < 1e-4 | Battram 2022 *Wellcome Open Res* 7:41 | a Catalog entry is a lookup, not a genome-wide claim |
| Detect 1-2% delta-beta -> hundreds-to-thousands of samples | pwrEWAS, Graw 2019 | tiny site-specific effects drive the field to meta-analysis |
| Run ComBat/limma on M-values, report delta-beta | Du 2010 *BMC Bioinformatics* 11:587 | betas are bounded/heteroscedastic; M-values stabilize variance |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Top hits are blood-cell CpGs | no cell-composition covariates | estimate and adjust cell fractions |
| ComBat removed the signal | batch confounded with phenotype | randomize at design time; cannot fix after |
| Lambda misread as pass/fail | GWAS dogma applied to EWAS | pair lambda with QQ shape + positive control; correct with BACON |
| Clean QQ but AHRR gone | over-correction with too many SVs | stop adding SVs when the positive control dies |
| Threshold flood or famine | Bonferroni-on-probe-count or raw p | use array-specific FWER / BH-FDR |
| "Novel" CpG is a known generic hit | no Catalog/Atlas triage | look up the CpG before claiming novelty |
| MRS called a risk score | PRS connotations on modifiable methylation | report as biomarker; orient cause only with prospective/MR design |

## References

- van Iterson M, van Zwet EW, Heijmans BT; BIOS Consortium. 2017. Controlling bias and inflation in epigenome- and transcriptome-wide association studies using the empirical null distribution. *Genome Biol* 18:19.
- Saffari A, Silver MJ, Zavattari P, et al. 2018. Estimation of a significance threshold for epigenome-wide association studies. *Genet Epidemiol* 42:20-33.
- Mansell G, Gorrie-Stone TJ, Bao Y, et al. 2019. Guidance for DNA methylation studies: statistical insights from the Illumina EPIC array. *BMC Genomics* 20:366.
- Leek JT, Storey JD. 2007. Capturing heterogeneity in gene expression studies by surrogate variable analysis. *PLoS Genet* 3:e161.
- Chen J, Behnam E, Huang J, et al. 2017. Fast and robust adjustment of cell mixtures in epigenome-wide association studies with SmartSVA. *BMC Genomics* 18:413.
- Johnson WE, Li C, Rabinovic A. 2007. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics* 8:118-127.
- Maksimovic J, Gagnon-Bartsch JA, Speed TP, Oshlack A. 2015. Removing unwanted variation in a differential methylation analysis of Illumina HumanMethylation450 array data. *Nucleic Acids Res* 43:e106.
- Min JL, Hemani G, Davey Smith G, Relton C, Suderman M. 2018. Meffil: efficient normalization and analysis of very large DNA methylation datasets. *Bioinformatics* 34:3983-3989.
- Graw S, Henn R, Thompson JA, Koestler DC. 2019. pwrEWAS: a user-friendly tool for comprehensive power estimation for epigenome wide association studies (EWAS). *BMC Bioinformatics* 20:218.
- Joehanes R, Just AC, Marioni RE, et al. 2016. Epigenetic Signatures of Cigarette Smoking. *Circ Cardiovasc Genet* 9:436-447.
- Elliott HR, Tillin T, McArdle WL, et al. 2014. Differences in smoking associated DNA methylation patterns in South Asians and Europeans. *Clin Epigenetics* 6:4.
- Gadd DA, Hillary RF, McCartney DL, et al. 2022. Epigenetic scores for the circulating proteome as tools for disease prediction. *eLife* 11:e71802.
- McCartney DL, Hillary RF, Stevenson AJ, et al. 2018. Epigenetic prediction of complex traits and death. *Genome Biol* 19:136.
- Battram T, Yousefi P, Crawford G, et al. 2022. The EWAS Catalog: a database of epigenome-wide association studies. *Wellcome Open Res* 7:41.
- Li M, Zou D, Li Z, et al. 2019. EWAS Atlas: a curated knowledgebase of epigenome-wide association studies. *Nucleic Acids Res* 47:D983-D988.
- Hannon E, Lunnon K, Schalkwyk L, Mill J. 2015. Interindividual methylomic variation across blood, cortex, and cerebellum. *Epigenetics* 10:1024-1032.
- Michels KB, Binder AM, Dedeurwaerder S, et al. 2013. Recommendations for the design and analysis of epigenome-wide association studies. *Nat Methods* 10:949-955.

## Related Skills

- differential-cpg-testing - The per-site test this design layer feeds
- cell-type-deconvolution - Cell-fraction covariates (the dominant confounder)
- array-preprocessing - Normalization choice (funnorm) vs model-level batch correction
- array-qc-filtering - Probe filtering and chip/position batch diagnosis
- causal-genomics/mendelian-randomization - mQTL-based causal orientation (reverse causation)
- experimental-design/batch-design - General randomization and batch-design principles
- clinical-biostatistics/multiplicity-graphical - FWER for confirmatory trials (contrast with discovery FDR)
- workflows/methylation-pipeline - End-to-end pipeline
