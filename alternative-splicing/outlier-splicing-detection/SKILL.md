---
name: bio-outlier-splicing-detection
description: Detects aberrant splicing in single rare-disease patients vs a control panel using FRASER 2.0 (Bioconductor; Beta-binomial autoencoder on Intron Jaccard Index, default delta cutoff 0.1, q hyperparameter), OUTRIDER (gene-level outlier expression via autoencoder denoising), LeafcutterMD (Dirichlet-multinomial outlier mode of LeafCutter for annotation-free junctions), and DROP (Snakemake pipeline integrating FRASER2 + OUTRIDER + monoallelic expression for clinical diagnostics). The statistical model is fundamentally different from differential splicing — single-sample-vs-cohort outlier detection rather than two-group comparison. Standard tool in EU rare-disease (Solve-RD) and NIH UDN programs. Use when applying RNA-seq to undiagnosed Mendelian disease, validating predicted splice variants in clinical samples, or detecting cryptic splicing in disease tissue.
tool_type: r
primary_tool: FRASER
---

## Version Compatibility

Reference examples tested with: FRASER 2.0 (>=1.99.0), OUTRIDER 1.20+, LeafcutterMD via leafcutter 0.2.9+, DROP 1.4+, R 4.4+, BiocManager 1.30+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Outlier Splicing Detection

For clinical RNA-seq diagnostics in rare disease, the question is not "what differs between groups?" but "what is aberrant in this single patient relative to a panel of unaffected samples?". The statistical framework is **single-sample-vs-cohort outlier detection**, fundamentally different from two-group differential splicing. Tools in this space are designed for clinical Mendelian diagnostic settings.

## Tool Taxonomy

| Tool | Statistic | Test target | Fails when |
|------|-----------|-------------|------------|
| FRASER 2.0 | Beta-binomial autoencoder on Intron Jaccard Index | Splicing outliers (per-sample, per-junction) | Cohort <20 samples; tissue mismatch |
| OUTRIDER | Autoencoder-denoised expression Z-score | Gene-level expression outliers (LoF, monoallelic) | Cohort <20 samples |
| LeafcutterMD | Dirichlet-multinomial outlier mode | Annotation-free intron usage | Beta-binomial fits poorly OR few controls |
| DROP | Snakemake pipeline | All of above + monoallelic expression | Pipeline complexity for small projects |

Core reference: **FRASER 2.0** for splicing outliers, **OUTRIDER** for expression outliers, **DROP** to combine. Standard tool in EU rare-disease programs (Solve-RD) and NIH UDN.

## Decision Tree by Diagnostic Scenario

| Scenario | Recommended approach |
|----------|----------------------|
| Single rare-disease patient + panel of n>=50 controls | FRASER 2.0 (Intron Jaccard Index) |
| Single patient + small panel (n=20-50) | FRASER 2.0 with auxiliary GTEx controls; tune q carefully |
| Patient + cohort <20 | Insufficient for outlier detection; consider differential or recruit more samples |
| Outlier expression suspected (loss of function, monoallelic) | OUTRIDER on same cohort |
| Annotation-free outlier (cryptic exon, novel junction) | LeafcutterMD |
| Integrated diagnostic pipeline (splicing + expression + MAE) | DROP |
| TDP-43 ALS post-mortem brain (cryptic exons) | FRASER 2.0; expect UNC13A, STMN2, ATG4B |
| SF3B1-mutant cancer sample | FRASER 2.0 with cohort-matched RNA-seq; expect cryptic 3'ss |
| Familial dysautonomia (ELP1) | FRASER 2.0 in fibroblast/iPSC; CNS tissue gives strongest signal |
| Stargardt deep-intronic ABCA4 | FRASER 2.0 in retina-relevant tissue |
| Solid tumor splicing biomarker | Differential splicing (n>=10 vs cohort) — see differential-splicing skill |
| RNA validation of SpliceAI hit | FRASER 2.0 + cross-reference with predicted variant location |

## When to Use Outlier vs Differential

**Outlier regime** (this skill):
- Single patient or small case series vs control panel
- Question: "What is aberrant in this patient?"
- Statistical model: single-sample p-value vs cohort distribution

**Differential regime** (differential-splicing skill):
- Two well-defined groups, n>=3 each
- Question: "What differs between groups?"
- Statistical model: two-group LRT or related

If you have n>=10 patients with a shared phenotype, prefer **differential** (more power); if single patient or heterogeneous case series, use **outlier**.

## FRASER 2.0 Workflow

**Goal:** Detect aberrant splicing in patient samples vs cohort using the Intron Jaccard Index.

**Approach:** Count split reads per junction, compute Intron Jaccard Index per intron, fit a Beta-binomial autoencoder to estimate expected values, then flag outliers by p-value and delta.

```r
library(FRASER); library(BiocParallel)

bam_files <- list.files('bams/', pattern='.bam$', full.names=TRUE)
sample_table <- data.frame(
    sampleID = gsub('.bam', '', basename(bam_files)),
    bamFile = bam_files,
    pairedEnd = TRUE
)

settings <- FraserDataSet(
    colData = sample_table,
    workingDir = 'fraser_workdir',
    name = 'rare_disease_cohort'
)

settings <- countRNAData(settings, BPPARAM = MulticoreParam(8))
fds <- calculatePSIValues(settings)

fds <- filterExpressionAndVariability(
    fds,
    minDeltaPsi = 0.0,
    minExpressionInOneSample = 20,
    quantile = 0.05,
    quantileMinExpression = 1
)

fitMetrics(fds) <- 'jaccard'
currentType(fds) <- 'jaccard'  # canonical setter for active metric in FRASER 2.0
fds <- FRASER(
    fds,
    q = c(jaccard = 10),
    BPPARAM = MulticoreParam(8)
)

results <- results(
    fds,
    psiType = 'jaccard',
    padjCutoff = 0.05,
    deltaPsiCutoff = 0.1
)

patient_results <- results[results$sampleID == 'PATIENT_001', ]
patient_results <- patient_results[order(patient_results$padjust), ]
```

**FRASER 2.0 changes vs FRASER 1.x:**
- Default `psiType` changed from three metrics (psi5, psi3, theta) to single **Intron Jaccard Index**
- Default `deltaPsiCutoff` dropped from 0.3 to **0.1**
- Pseudocount and filtering parameter optimization
- Bioconductor package version >=1.99.0 == FRASER 2.0

`q = 10` is the autoencoder dimension hyperparameter. **Tune via `optimHyperParams(fds, type='jaccard')` for cohort-specific optimum** — too low: confounders not removed; too high: real signal absorbed.

## OUTRIDER for Gene-Level Outlier Expression

**Goal:** Detect genes with aberrantly high or low expression in patient samples.

**Approach:** Autoencoder denoising of expression matrix; outliers identified by Z-score and adjusted p-value.

```r
library(OUTRIDER); library(BiocParallel)

countTable <- read.table('counts.tsv', header=TRUE, row.names=1)
ods <- OutriderDataSet(countData = countTable)

ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
ods <- estimateBestQ(ods, BPPARAM = MulticoreParam(8))
ods <- OUTRIDER(ods, BPPARAM = MulticoreParam(8))

res <- results(ods, padjCutoff = 0.05, zScoreCutoff = 0)
patient_outliers <- res[res$sampleID == 'PATIENT_001', ]
```

OUTRIDER (Brechtmann 2018 *Am J Hum Genet*) catches loss-of-function alleles producing transcript collapse, monoallelic effects, and tissue-inappropriate expression — complements splice outlier detection.

## LeafcutterMD for Annotation-Free Outlier Intron Usage

**Goal:** Detect outlier intron usage relative to a control panel without annotation dependence.

**Approach:** Run LeafCutter in MD (Mahalanobis Distance) mode against the control panel.

```bash
for bam in *.bam; do
    regtools junctions extract -a 8 -m 50 -s XS "$bam" -o "${bam%.bam}.junc"
done

ls *.junc > juncfiles.txt
python leafcutter_cluster_regtools.py -j juncfiles.txt -o leafcutter -m 50 -l 500000

leafcutterMD.R \
    --num_threads 4 \
    --output_prefix patient_outlier \
    leafcutter_perind_numers.counts.gz
```

LeafcutterMD (Jenkinson 2020 *Bioinformatics*) reports per-sample p-values per intron-cluster; useful when FRASER's Beta-binomial model fits poorly or when novel-junction sensitivity matters.

## DROP Pipeline (Integrated Workflow)

**Goal:** Run FRASER2 + OUTRIDER + monoallelic expression in a unified Snakemake pipeline for clinical diagnostics.

**Approach:** DROP is distributed via **bioconda** (not PyPI). Install in a dedicated environment, then configure with patient + control sample sheets; pipeline handles QC, alignment, counting, autoencoding, and reporting.

```bash
# Install via bioconda (DROP is not on PyPI)
mamba create -n drop_env -c conda-forge -c bioconda drop --override-channels
conda activate drop_env

drop init my_diagnostic_run
cd my_diagnostic_run

# Edit config.yaml:
#  - sample_table: samples.tsv (patient + controls)
#  - aberrantSplicing: enabled
#  - aberrantExpression: enabled
#  - mae: enabled (monoallelic expression)

snakemake --cores 16 --use-conda
```

DROP (Yepez 2021 *Nat Protocols*) is the standard tool in EU rare-disease genome+RNA-seq programs (Solve-RD) and the NIH UDN. v1.4+ uses FRASER 2.0. The **MAE module** uses a custom z-score test on heterozygous SNPs from RNA-seq (allele-specific expression) — useful for catching dominant-negative or monoallelic LoF that splicing/expression outliers miss. Cohort >=30 samples recommended for confident outlier detection.

## Variant + Outlier Integration

**Goal:** Connect a candidate splice-altering DNA variant to RNA-level confirmation.

**Approach:** Cross-reference SpliceAI hits with FRASER2 outliers in the same sample.

```r
library(dplyr)

variants <- read.table('spliceai_hits.tsv', header=TRUE, sep='\t')
fraser_hits <- read.table('fraser_results.tsv', header=TRUE, sep='\t')

confirmed <- variants %>%
    filter(delta_max >= 0.2) %>%
    inner_join(
        fraser_hits %>% filter(sampleID == 'PATIENT_001', padjust < 0.05),
        by = c('chrom' = 'seqnames'),
        relationship = 'many-to-many'
    ) %>%
    filter(abs(pos - start) < 1000 | abs(pos - end) < 1000)
```

A SpliceAI hit + concordant FRASER2 outlier in the patient = strong PS3 functional evidence in the ACMG framework. This integration is the highest-value clinical pipeline step — converts a computational PP3 to functional PS3.

## Cohort Size and Power

| Cohort size | Power | Comment |
|-------------|-------|---------|
| n < 20 | Marginal | High FDR; consider GTEx tissue-matched controls as auxiliary |
| n = 20-50 | Acceptable | FRASER autoencoder can fit; tune q carefully |
| n >= 50 | Recommended | Standard clinical diagnostic cohort size |
| n >= 100 | Optimal | Tissue-matched and batch-matched gives best calibration |

GTEx-derived tissue-matched controls can supplement small in-house cohorts but introduce batch effects; use only when in-house n < 30 and document the pooling strategy.

## Tissue Choice for Mendelian RNA-seq

| Tissue | Pros | Cons | Genes captured |
|--------|------|------|----------------|
| Whole blood (PAXgene) | Easy, standard | Globin contamination; many disease genes silent | ~70-80% of clinical genes |
| Fibroblast (skin biopsy) | Reasonable expression | Requires culture; senescence variability | ~75-85% |
| Muscle biopsy | Best for muscular dystrophy | Invasive | ~85-90% for muscle disorders |
| iPSC-derived neuron / cardiomyocyte | Disease-relevant tissue | Cost, variability | ~95% if differentiation works |
| Urine sediment | Non-invasive | Low yield | ~50-60% |

For UDN-style cases: blood first, then fibroblast if blood lacks expression of candidate gene. Critical: **a negative blood RNA-seq doesn't rule out a candidate gene that's silent in blood** — verify gene expression with GTEx tissue panel before committing to the tissue.

## Hyperparameter Tuning

```r
fds <- estimateBestQ(fds, type='jaccard', useOHT=TRUE)
plotEncDimSearch(fds, type='jaccard')
```

The encoding dimension `q` should be where the loss curve plateaus. Too low: confounders not removed; too high: real signal absorbed.

For typical 50-100 sample cohorts, q=8-15 is the usual range. For very small cohorts (n=20-30), q=5-8.

## Per-Tool Failure Modes

### FRASER 2.0: Q Hyperparameter Mistuning

**Trigger:** Default q=10 used without tuning; or wrong q for cohort size.

**Mechanism:** Q is the autoencoder bottleneck dimension; too small → confounders leak into outlier signal; too large → real biological signal absorbed by autoencoder.

**Symptom:** Either no significant outliers (q too high) or many spurious calls clustering by batch (q too low).

**Fix:** Run `estimateBestQ(fds, type='jaccard', useOHT=TRUE)` (Optimal Hard Thresholding default; very fast); use the returned `bestQ(fds)` value. For exhaustive search, pass `useOHT=FALSE, q_param=c(2, 5, 10, 15)` and inspect `plotEncDimSearch` for the plateau.

### FRASER 2.0: Tissue Mismatch

**Trigger:** Patient sample from different tissue than majority of controls.

**Mechanism:** FRASER autoencoder learns tissue-specific expression patterns; tissue-mismatched patient appears as global outlier.

**Symptom:** Hundreds of "significant" outliers in patient; not biologically interpretable.

**Fix:** Strict tissue matching; if controls are mixed-tissue, use only controls from patient's tissue.

### OUTRIDER: Few Controls

**Trigger:** Cohort <20 samples.

**Mechanism:** Autoencoder needs sufficient samples to learn expression covariance; fails to fit at very small cohort sizes.

**Symptom:** Convergence warnings; uncalibrated p-values.

**Fix:** Pool with GTEx auxiliary controls; or use simpler outlier methods (z-score on log-CPM).

### LeafcutterMD: Cluster Count Limits

**Trigger:** Very few clusters in patient sample (low coverage or filtered out).

**Mechanism:** LeafcutterMD computes Mahalanobis distance over clusters; few observations → unstable distance.

**Symptom:** Inflated or deflated p-values; few significant calls.

**Fix:** Increase coverage; relax filtering (`-m 10` instead of 50); or switch to FRASER2.

### DROP: Snakemake Pipeline Failures

**Trigger:** Missing dependencies or incompatible R/Bioconductor versions.

**Mechanism:** DROP orchestrates many tools; version mismatches cascade through pipeline.

**Symptom:** Snakemake step fails partway through; cryptic R errors.

**Fix:** Use `--use-conda` flag for environment isolation; pin versions in environment yamls.

## Reconciliation: When Outlier Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| FRASER2 sig, OUTRIDER not | Splicing change without expression collapse | Standard splicing outlier; report |
| OUTRIDER sig, FRASER2 not | Expression LoF without splicing change | Likely promoter / regulatory; not splicing |
| Both sig at same gene | LoF allele triggering NMD on splicing-disrupted transcript | Strong combined evidence; expect downstream |
| LeafcutterMD sig, FRASER2 not | Novel cryptic event not in annotation | High-priority novel finding; investigate |
| All tools null but biology suggests change | Underpowered cohort or wrong tissue | Verify gene expression in tissue; recruit larger cohort |

## Disease-Specific Expectations

| Condition | Expected outlier signature | Tissue |
|-----------|----------------------------|--------|
| ALS / FTD (TDP-43 loss) | Cryptic exons in UNC13A, STMN2, ATG4B | Post-mortem brain ONLY |
| SF3B1-mutant MDS / CLL / uveal melanoma | Aberrant 3'ss ~10-30nt upstream of canonical | Bone marrow / tumor tissue |
| Spinal muscular atrophy (untreated SMN2) | SMN exon 7 skipping | Fibroblast / iPSC-MN |
| Familial dysautonomia (ELP1 c.2204+6T>C) | ELP1 exon 20 skipping (>=99% in CNS, partial elsewhere) | iPSC-neuron > fibroblast > blood |
| Deep-intronic CFTR / USH2A / CEP290 | Pseudoexon inclusion | Cognate disease tissue (lung / retina) |
| Duchenne muscular dystrophy (DMD) | Out-of-frame exon skipping pattern | Muscle biopsy |
| Stargardt (ABCA4) deep-intronic | Pseudoexon in retina | Retinal organoid / iPSC-RPE |

For each, the gene must be expressed in the queried tissue. Verify with GTEx before assuming negative result rules out the gene.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `FRASER: cohort too small` | n<10 | Pool with auxiliary controls; or recruit more patients |
| `FRASER: countRNAData failed on chromosome X` | BAM index missing or corrupted | Re-index BAMs; check `samtools idxstats` |
| `estimateBestQ: convergence not reached` | Default q range insufficient | Expand `q_param=c(2,5,10,15,20)`; or use `useOHT=TRUE` for the fast deterministic alternative |
| `OUTRIDER: estimateBestQ slow` | Default range tries q=2-30 | Restrict to expected range (`q=4:15`) |
| `DROP: snakemake job failed at FRASER` | DROP-FRASER version mismatch | Update DROP to latest; verify FRASER 2.0 compatibility |
| `LeafcutterMD: insufficient clusters` | Cluster filter too strict | Lower `-m` minimum cluster reads |
| `Variant integration: chrom format mismatch` | VCF uses 1, FRASER uses chr1 (or vice versa) | Standardize with `bcftools annotate --rename-chrs` |

## Common Pitfalls

- **Using bulk differential-splicing tools for n=1 vs cohort** — rMATS, leafcutter (regular), SUPPA2 are not designed for this. Use FRASER2 / LeafcutterMD.
- **Ignoring tissue choice** — clinical gene expression varies dramatically across blood / fibroblast / muscle. A negative blood RNA-seq doesn't rule out a candidate gene that's silent in blood.
- **Forgetting batch effects** — combining in-house and external (GTEx) controls introduces sequencing batch confounding; use ComBat or include batch as covariate.
- **Skipping the variant + outlier integration** — RNA-only outlier without DNA variant suggests cellular state or technical artifact; DNA-only prediction without RNA confirmation is supporting only (PP3, not PS3).
- **Treating all FRASER2 outliers as pathogenic** — many are benign tissue-specific variation. Filter against gnomAD splice constraint and disease gene panels.
- **Q hyperparameter not tuned** — default `q=10` works for ~50-100 sample cohorts; tune for outliers.
- **Wrong default delta cutoff for FRASER 1.x vs 2.0** — 1.x default 0.3, 2.0 default 0.1; document which version.

## Quality Thresholds

| Metric | Recommendation | Source |
|--------|----------------|--------|
| Cohort size | n>=50 (ideal); n>=20 (minimum) | Solve-RD / UDN convention |
| FRASER 2.0 padj | < 0.05 | Standard |
| FRASER 2.0 delta Jaccard | >= 0.1 (default in v2.0) | Scheller 2023 *AJHG* |
| OUTRIDER padj | < 0.05 | Brechtmann 2018 *AJHG* |
| OUTRIDER zScore | abs >= 2 | Conservative |
| Sequencing depth | >=50M PE reads/sample | Standard for AS analysis |
| Tissue match between patient and controls | Required | Critical for FRASER2 calibration |
| Batch match | Strongly recommended | Reduces autoencoder confounding |

## Related Skills

- splice-variant-prediction - SpliceAI / Pangolin for in-silico prediction; integration target
- differential-splicing - When testing multiple patients vs controls (>=10 vs cohort)
- splicing-qc - Library / depth / tissue prerequisites
- variant-calling/clinical-interpretation - ACMG/AMP framework integration
- workflows/clinical-trial-pipeline - Trial-grade RNA-seq diagnostics

## References

- Mertes et al 2021 *Nat Commun* - FRASER 1.x
- Scheller et al 2023 *Am J Hum Genet* - FRASER 2.0 (Intron Jaccard Index)
- Brechtmann et al 2018 *Am J Hum Genet* - OUTRIDER
- Jenkinson et al 2020 *Bioinformatics* - LeafcutterMD
- Yepez et al 2021 *Nat Protocols* - DROP pipeline
- Cummings et al 2017 *Sci Transl Med* - RNA-seq for muscular dystrophy diagnostics
- Kremer et al 2017 *Nat Commun* - RNA-seq for mitochondrial disease
- Brown et al 2022 *Nature* - UNC13A cryptic exon (TDP-43 / ALS)
- Klim et al 2019 *Nat Neurosci* - STMN2 cryptic splicing (ALS)
- Darman et al 2015 *Cell Rep* - SF3B1 cryptic 3'ss
- Walker et al 2023 *Am J Hum Genet* - ClinGen SVI splicing recommendations
