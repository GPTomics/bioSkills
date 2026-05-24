---
name: bio-crispr-screens-mageck-analysis
description: Analyzes pooled CRISPR screens with MAGeCK (Li et al 2014), covering count generation (mageck count), the RRA two-condition workflow (mageck test using alpha-RRA over per-sgRNA negative-binomial p-values), the MLE multi-condition workflow (mageck mle with explicit design matrix and beta-score output), normalization choice (median vs total vs control-sgRNA vs spike-in), sgRNA efficiency injection, paired-sample testing, time-course design, drug-screen versus dropout-screen design matrices, MAGeCKFlute and MAGeCK-VISPR downstream visualization, and decision logic for when to use MAGeCK vs JACKS / BAGEL2 / drugZ / Chronos. Use when running a fresh CRISPR screen analysis, picking RRA vs MLE for the experimental design, choosing a normalization method from QC signatures, debugging MLE convergence failure or NaN beta scores, comparing MAGeCK output across tools, or building a batch-aware multi-cell-line / multi-condition MLE design matrix.
tool_type: cli
primary_tool: MAGeCK
---

## Version Compatibility

Reference examples tested with: MAGeCK 0.5.9+, MAGeCKFlute 2.0+ (R/Bioconductor), MAGeCK-VISPR 0.5.6+, pandas 2.2+, numpy 1.26+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mageck --version`, `mageck count --help`, `mageck test --help`, `mageck mle --help`
- R: `packageVersion('MAGeCKFlute')`, `?FluteRRA`, `?FluteMLE`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## MAGeCK CRISPR Screen Analysis

**"Run MAGeCK on my pooled CRISPR screen"** -> Count sgRNAs from FASTQ, normalize across samples, and rank genes by enrichment or depletion using either the robust rank aggregation (RRA) test for two-condition designs or the maximum-likelihood (MLE) model with explicit design matrix for multi-condition / time-course / drug screens.

- CLI: `mageck count` -> `mageck test` for two-condition RRA
- CLI: `mageck mle` for multi-condition / time-course / multi-cell-line MLE
- R: `MAGeCKFlute::FluteRRA()` / `FluteMLE()` for downstream visualization and pathway analysis
- Python: `mageck-vispr` for interactive QC + result dashboard

## RRA vs MLE Decision Tree

| Experimental design | Recommended | Why |
|---------------------|-------------|-----|
| Two conditions (e.g. drug vs vehicle, treated vs untreated), single cell line, no covariates | `mageck test` (RRA) | RRA is more robust to outlier sgRNAs; faster; default for most published screens |
| Time series (Day 0 -> Day 7 -> Day 14 -> Day 21) | `mageck mle` | RRA cannot model multiple timepoints jointly; MLE estimates per-condition beta scores |
| Multi-cell-line panel (e.g. DepMap-style 5-50 lines) | `mageck mle` with cell-line covariate, or Chronos | MLE handles >2 conditions; Chronos preferred at DepMap scale |
| Paired samples (each replicate matched donor/cell prep) | `mageck mle` with paired design | RRA does not support pairing |
| Combinatorial (treatment x cell line x time) | `mageck mle` with full factorial design | RRA only handles 1 factor |
| Drug screen (vehicle vs drug, multiple doses) | `mageck mle` with dose covariate OR drugZ (preferred for chemogenomic) | drugZ optimized for chemogenomic; see [[drugz-chemogenomic]] |
| Essentiality (Day 0 -> endpoint, single condition) | `mageck test` for simple dropout; `mageck mle` if multi-cell-line | RRA suffices; or BAGEL2 for Bayesian essentiality; see [[bagel-essentiality]] |
| Multi-batch / multi-screen joint analysis | `mageck mle` with batch covariate, JACKS for guide efficacy, or Chronos | See [[batch-correction]] |

**Fails when:**
- Using RRA for time-series and treating each timepoint as a separate test: false-discovery inflation from un-modeled temporal correlation. Use MLE.
- Using MLE without explicit design matrix in a screen with severe batch effects: beta scores include batch variance. Add batch covariate.
- Using MAGeCK for drug screens without vehicle control: comparing drug vs Day 0 conflates drug effect with proliferation. Compare drug vs vehicle, not Day 0. See [[drugz-chemogenomic]].

## The RRA Algorithm (under the hood)

**Why this matters for postdoc-level use:** RRA is not "non-parametric"; it tests whether per-gene sgRNA p-value ranks are more clustered toward the extremes than expected under the uniform null. The chain:

1. `mageck count` normalizes raw counts (median by default) and outputs `*.count_normalized.txt`.
2. `mageck test` computes a per-sgRNA p-value under a NB model fitted to per-gene variance (mean-variance trend stabilized via empirical Bayes shrinkage; same family as edgeR/DESeq2 but simpler dispersion estimation).
3. Per-sgRNA p-values are ranked; ranks are normalized to percentile (rho).
4. For each gene with k sgRNAs, the alpha-RRA score is the minimum over `Pr(beta(i, k-i+1) <= rho_i)` for i in 1..k -- the probability of observing the i-th smallest rank in a uniform sample.
5. Gene-level p-value is computed by permuting sgRNA-to-gene assignment to derive an empirical null over alpha-RRA scores.
6. Multiple-testing correction is Benjamini-Hochberg.

**Critical assumption:** RRA assumes most sgRNAs are non-changing (used to estimate the NB dispersion). If >40% of sgRNAs change, median normalization fails and the dispersion estimate is biased. Symptom: every gene appears significant. Fix: use control-sgRNA normalization (`--norm-method control`) with non-targeting controls as the reference.

## The MLE Model (under the hood)

**Why this matters:** MLE assumes the per-sgRNA log-count is Negative-Binomial with mean determined by a linear combination of condition betas plus an sgRNA-efficiency term (if supplied). The chain:

1. Define a design matrix where rows are samples and columns are conditions; entries are 1 if the sample belongs to the condition, 0 otherwise. A baseline column (all 1s) is required.
2. Per-gene, the model is `log(count_ij) = baseline_i + sum_c (beta_gc * design_jc) + log(sgrna_efficiency_i) + log(size_factor_j)` where i is sgRNA, j is sample, c is condition.
3. The optimizer finds the per-gene beta_gc maximizing the NB likelihood; per-gene Wald-statistic gives a p-value per condition.
4. Each beta represents the log-fold-change in that condition relative to baseline, accounting for sgRNA efficiency.

**Critical assumption:** Beta scores are interpretable only when the design matrix is correctly specified. A common error is omitting batch as a covariate -- the resulting betas absorb batch variance. The fix is to add batch columns; the betas then estimate biology after batch adjustment.

**Convergence failure:** MLE NaN beta scores indicate optimizer divergence -- usually because a gene has too few sgRNAs with non-zero counts in the relevant conditions. The output includes such genes with NaN; do not interpret them as zero effect.

## Count sgRNAs from FASTQ

**Goal:** Quantify sgRNA representation from raw sequencing data.

**Approach:** Run `mageck count` to map FASTQ reads to the sgRNA library reference, producing a normalized count matrix and QC summary. Set `--norm-method median` (default; robust to outliers) or `--norm-method control` (when many sgRNAs change).

```bash
mageck count \
    --list-seq library.csv \                       # sgRNA library: header sgRNA,Gene,Sequence
    --sample-label Plasmid,Day0,Veh_r1,Veh_r2,Drug_r1,Drug_r2 \
    --fastq Plasmid.fq.gz Day0.fq.gz Veh_r1.fq.gz Veh_r2.fq.gz Drug_r1.fq.gz Drug_r2.fq.gz \
    --norm-method median \                         # see normalization decision below
    --output-prefix screen \
    --trim-5 CACCG                                 # 5'-adapter (lentiGuide-Puro convention)

# Outputs:
#   screen.count.txt           raw counts
#   screen.count_normalized.txt normalized counts (median-scaled)
#   screen.countsummary.txt    per-sample QC: Gini, reads, mapping rate, % zero-count
#   screen.log                  per-FASTQ mapping stats
```

**Library file format (tab- or comma-separated; first row is header):**

```
sgRNA,Gene,Sequence
BRCA1_1,BRCA1,ATGGATTTATCTGCTCTTCG
BRCA1_2,BRCA1,CAGCAGATACTTGATGCATC
NTC_0001,NonTargeting_0001,GACGCATCGAATCAATAGCC
```

## Normalization Decision

| `--norm-method` | When to use | Mechanism | Fails when |
|------------------|-------------|-----------|------------|
| `median` (default) | Standard screen, <40% guides change | Scale each sample to a common median of all guides | Heavy selection (>40% guides change) inflates median, biasing scaling |
| `total` | Equal sequencing depth assumed, no outliers | Scale to total reads | Sensitive to PCR jackpots / outlier high-count guides |
| `none` | Already-normalized inputs (DEPRECATED for raw FASTQ counting) | Skips scaling | Almost never appropriate; only for already-normalized inputs |
| `control` | Heavy selection screens; library has ≥500 NTCs | Scale each sample so non-targeting controls have constant median | Requires `--control-sgrna ntcs.txt` listing NTC sgRNAs |

**Diagnostic:** Run with `median` first. If essentialome PR-AUC against CEGv2 is high and Gini is in range, accept. If PR-AUC <0.5 despite reasonable Gini, retry with `control` -- this isolates whether median normalization was masking the signal.

## MAGeCK Test (RRA for Two-Condition)

**Goal:** Identify genes significantly enriched or depleted between treatment and control.

**Approach:** Compute per-sgRNA NB-model p-values, rank, apply alpha-RRA to get per-gene scores, permute for gene-level FDR. Outputs separate columns for positive and negative selection.

```bash
mageck test \
    --count-table screen.count.txt \
    --treatment-id Drug_r1,Drug_r2 \
    --control-id Veh_r1,Veh_r2 \
    --norm-method median \
    --output-prefix drug_vs_veh \
    --gene-lfc-method median \                     # alternative: mean (less robust)
    --sort-criteria pos                            # rank for positive selection
# Outputs: drug_vs_veh.gene_summary.txt (gene-level), drug_vs_veh.sgrna_summary.txt
```

**Gene-summary columns:**

| Column | Meaning |
|--------|---------|
| `id` | Gene symbol |
| `num` | Number of sgRNAs in library for this gene |
| `neg|score`, `pos|score` | alpha-RRA score, negative- or positive-selection direction |
| `neg|p-value`, `pos|p-value` | Permutation p-value |
| `neg|fdr`, `pos|fdr` | BH-corrected FDR |
| `neg|rank`, `pos|rank` | Rank by score (1 = top hit in that direction) |
| `neg|lfc`, `pos|lfc` | Median log-fold-change of sgRNAs in that direction |

**Interpretation rule:** A gene is a strong essential if `neg|fdr < 0.05` AND `neg|lfc < -1`. A gene is a strong resistance hit if `pos|fdr < 0.05` AND `pos|lfc > 1`. Genes with `fdr < 0.05` but `|lfc| < 0.5` are statistically significant but biologically weak -- worth flagging for orthogonal validation.

## MAGeCK MLE (Multi-Condition)

**Goal:** Estimate gene effects across complex experimental designs with multiple conditions, time points, batches, or covariates.

**Approach:** Specify a design matrix mapping samples to conditions, then run `mageck mle` which fits a NB GLM with sgRNA-efficiency term and outputs per-gene per-condition beta scores.

```bash
# Design matrix: design.txt (tab-separated)
# Samples must match sample-label in mageck count output
# baseline column must be present and all 1
cat > design.txt <<EOF
Samples	baseline	day7	day14	day21
Day0	1	0	0	0
Day7_r1	1	1	0	0
Day7_r2	1	1	0	0
Day14_r1	1	0	1	0
Day14_r2	1	0	1	0
Day21_r1	1	0	0	1
Day21_r2	1	0	0	1
EOF

mageck mle \
    --count-table screen.count.txt \
    --design-matrix design.txt \
    --output-prefix timecourse_mle \
    --norm-method median \
    --sgrna-efficiency efficiency.txt \            # OPTIONAL: from JACKS or library design
    --sgrna-eff-name-column 1 \
    --sgrna-eff-score-column 2 \
    --max-sgrnapergene-permutation 10
# Outputs: timecourse_mle.gene_summary.txt with beta scores per condition + Wald p-values
```

**Gene-summary columns (MLE):**

| Column | Meaning |
|--------|---------|
| `Gene` | Gene symbol |
| `sgRNA` | Number of sgRNAs |
| `<condition>|beta` | Effect-size estimate (log-fold-change relative to baseline) |
| `<condition>|z` | Wald z-statistic |
| `<condition>|p-value` | Two-sided p-value |
| `<condition>|fdr` | BH-corrected FDR |
| `<condition>|wald-fdr` | Wald-statistic-based FDR (alternative) |

**Interpretation rule:** Beta scores are log2-fold-changes; a beta of -1 in condition day21 means sgRNAs are 2-fold depleted in day-21 relative to baseline. NaN betas indicate convergence failure (typically a gene with too few non-zero counts in that condition); exclude from interpretation.

## Sample MAGeCK Test for Drug Screen with sgRNA Efficiency

```bash
mageck test \
    --count-table screen.count.txt \
    --treatment-id Drug_r1,Drug_r2,Drug_r3 \
    --control-id Veh_r1,Veh_r2,Veh_r3 \
    --norm-method control \                         # NTCs as normalization reference
    --control-sgrna ntcs.txt \                      # one NTC sgRNA per line
    --gene-lfc-method median \
    --variance-estimation-samples Plasmid,Day0 \    # estimate NB dispersion from early samples
    --sgrna-efficiency efficiency.txt \
    --sgrna-eff-name-column 1 \
    --sgrna-eff-score-column 2 \
    --output-prefix drug_screen_normalized
```

**Reading order:** Run `mageck count` -> screen-qc skill for QC -> `mageck test` or `mageck mle` -> `MAGeCKFlute` for visualization -> downstream pathway analysis.

## Time-Course Analysis

**Goal:** Identify genes with consistent depletion or enrichment across multiple time points.

**Approach:** Run `mageck mle` with a design matrix where each timepoint is a separate column, then test for monotonic trends in per-condition beta scores.

```python
import pandas as pd
import numpy as np

def time_course_consistency(mle_results, conditions=['day7', 'day14', 'day21']):
    '''Identify genes with monotonic beta trends across timepoints.
    Returns genes where all betas same sign and trend is monotone.'''
    beta_cols = [f'{c}|beta' for c in conditions]
    fdr_cols = [f'{c}|fdr' for c in conditions]
    df = mle_results[['Gene'] + beta_cols + fdr_cols].copy()
    df['all_negative'] = (df[beta_cols] < 0).all(axis=1)
    df['all_positive'] = (df[beta_cols] > 0).all(axis=1)
    df['monotone'] = df[beta_cols].apply(lambda x: (np.diff(x) <= 0).all() or (np.diff(x) >= 0).all(), axis=1)
    df['any_sig'] = (df[fdr_cols] < 0.05).any(axis=1)
    return df[df['monotone'] & df['any_sig']].sort_values(beta_cols[-1])
```

## Visualizing Results

**Goal:** Generate publication-grade volcano plot and rank plot of MAGeCK output.

**Approach:** Load gene_summary.txt, plot `-log10(fdr)` vs LFC, color by significance, annotate top hits.

```python
import matplotlib.pyplot as plt
import numpy as np

def volcano(gene_summary_path, direction='neg', fdr_threshold=0.05, lfc_threshold=1.0):
    '''direction: "neg" for dropout, "pos" for enrichment.'''
    df = pd.read_csv(gene_summary_path, sep='\t')
    lfc_col = f'{direction}|lfc'
    fdr_col = f'{direction}|fdr'
    fig, ax = plt.subplots(figsize=(9, 7))
    sig = (df[fdr_col] < fdr_threshold) & (np.abs(df[lfc_col]) > lfc_threshold)
    ax.scatter(df.loc[~sig, lfc_col], -np.log10(df.loc[~sig, fdr_col].clip(lower=1e-10)),
                c='lightgray', alpha=0.4, s=10)
    ax.scatter(df.loc[sig, lfc_col], -np.log10(df.loc[sig, fdr_col].clip(lower=1e-10)),
                c='red' if direction == 'neg' else 'blue', alpha=0.7, s=18)
    top = df.loc[sig].nsmallest(15, fdr_col)
    for _, r in top.iterrows():
        ax.annotate(r['id'], (r[lfc_col], -np.log10(max(r[fdr_col], 1e-10))), fontsize=7)
    ax.axhline(-np.log10(fdr_threshold), ls='--', c='black', lw=0.5)
    ax.axvline(0, c='gray', lw=0.5)
    ax.set_xlabel(f'{direction} LFC')
    ax.set_ylabel(f'-log10({direction} FDR)')
    return fig
```

## MAGeCKFlute Integration (R)

**Goal:** Use the canonical pathway-analysis dashboard for downstream visualization.

**Approach:** Load MAGeCK output in R, run FluteRRA or FluteMLE, which produces volcano, rank, square-plot, and KEGG/Reactome enrichment.

```r
library(MAGeCKFlute)
# After mageck test
FluteRRA(gene_summary = "drug_vs_veh.gene_summary.txt",
         sgrna_summary = "drug_vs_veh.sgrna_summary.txt",
         organism = "hsa",
         outdir = "flute_output/")

# After mageck mle
FluteMLE(gene_summary = "timecourse_mle.gene_summary.txt",
         treatname = "day21", ctrlname = "baseline",
         organism = "hsa",
         outdir = "flute_mle_output/")
```

## MAGeCK-VISPR Interactive Dashboard

```bash
mageck-vispr init my_screen
# Edit config.yaml to point to counts + library
mageck-vispr run --cores 8
# Generates interactive HTML at output/index.html
```

## Failure Modes

### NaN beta scores in MLE output

**Trigger:** A gene has fewer than 2 non-zero-count sgRNAs in the relevant condition.
**Mechanism:** MLE optimizer cannot fit beta when likelihood is flat (insufficient evidence).
**Symptom:** `mle.gene_summary.txt` shows NaN in specific gene-condition cells.
**Fix:** Filter these genes from downstream interpretation; they are not "zero effect" but undetermined. If many genes show NaN, the screen depth is too low or many guides have failed -- audit with screen-qc.

### Every gene appears significant after RRA

**Trigger:** >40% of sgRNAs change direction (heavy selection screen).
**Mechanism:** Median normalization assumes most guides are non-changing; under heavy selection, median is biased and dispersion estimation breaks.
**Symptom:** Thousands of "hits" at FDR <0.05; volcano plot looks like a U with no separation.
**Fix:** Switch to `--norm-method control` with non-targeting controls; or use BAGEL2 / Chronos for essentiality analysis (`mageck` is not designed for screens with high-fraction true essentiality).

### MLE beta absorbs batch variance

**Trigger:** Multi-batch screen run without explicit batch covariate in design matrix.
**Mechanism:** Without a batch column, beta_condition includes both biological signal and batch difference between conditions.
**Symptom:** Beta scores differ between batches for the same biological condition; PCA shows samples cluster by batch not condition.
**Fix:** Add batch covariate columns to design matrix (e.g., one column per batch indicator); the biological betas are then estimated after batch adjustment. See [[batch-correction]].

### sgRNA-efficiency injection makes scores worse

**Trigger:** Supplied JACKS efficiency scores from a different cell line / chemistry than the current screen.
**Mechanism:** sgRNA efficacy is cell-line and modality dependent; HEK293T efficiency is not HCT116 efficiency; Cas9 efficiency is not CRISPRi efficiency.
**Symptom:** Some genes that were hits without efficiency become non-significant with efficiency.
**Fix:** Use efficiency derived from JACKS run on the matching cell line / library / chemistry, or use no efficiency (treat all sgRNAs as equal).

### Lack of cell-line covariate in multi-line screen

**Trigger:** Running `mageck mle` on multiple cell lines without a cell-line indicator column.
**Mechanism:** Each cell line has different essentiality profile; combining without indicator pools variance and dilutes per-line signal.
**Symptom:** Per-line known essentials don't show up; results look like a meta-analysis without effect sizes.
**Fix:** Add cell-line indicator columns; or move to Chronos which models cell-line and screen-quality jointly (see [[hit-calling]]).

## Reconciliation: When MAGeCK Disagrees With Other Tools

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MAGeCK FDR significant, BAGEL2 Bayes factor not | Different reference set; BAGEL2 trained on CEGv2/NEGv1 | Trust BAGEL2 for essentiality calling; MAGeCK for general LFC |
| MAGeCK significant, JACKS not | JACKS down-weighted noisy guides; MAGeCK trusts all 4 | Inspect per-sgRNA scores; if one sgRNA dominates, JACKS is right |
| MAGeCK and Chronos agree at top 100; disagree at 100-300 | Chronos accounts for CN and screen quality; MAGeCK does not | Trust Chronos for cancer-line screens; MAGeCK lacks CN correction |
| MAGeCK significant, drugZ not | drugZ uses bidirectional Z; MAGeCK uses RRA | For chemogenomic, trust drugZ (vehicle-anchored) |
| RRA and MLE disagree on same data | RRA more robust to outliers; MLE more sensitive | Higher-ranked hits in both = high confidence; only-one-method = orthogonal validate |

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Hit FDR (gene-level) | <0.05 | MAGeCK paper Li 2014; standard |
| Hit LFC magnitude | >1 (2-fold) | Biologically interpretable effect size; FDR-only hits at low LFC need validation |
| Normalization fraction-changing threshold | <40% guides change for median norm | Robson 2024 MAGeCK benchmark |
| sgRNAs needed for reliable MLE | ≥3 per gene with non-zero counts in each condition | MAGeCK 0.5+ behavior; below this, NaN |
| Time-course conditions for MLE | ≥3 (otherwise use test) | MLE statistical power |
| PR-AUC of ranked hits against CEGv2 | >0.7 for "passing" essentiality screen | Hart 2017; see [[screen-qc]] |
| Permutation iterations for RRA gene-level p-value | 1000 (default), 10000 for tight FDR | `--max-sgrnapergene-permutation`; trades runtime |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `mageck count` outputs mostly zero counts | Wrong `--trim-5` adapter | Check library plasmid; lentiGuide is CACCG, Brunello may differ |
| All genes significant | Median normalization break (heavy selection) | `--norm-method control` |
| NaN beta in MLE | Gene with insufficient non-zero counts | Exclude from interpretation |
| Two-condition MLE works but ranks differ from RRA | Different test statistic | Both correct; check direction and use the appropriate one |
| Hits include amplified genes | No CN correction | See [[copy-number-correction]] |
| Library not detected in screen.countsummary.txt | Library file format wrong | Check header is `sgRNA,Gene,Sequence` (no spaces) |

## References

- Li W et al. 2014. *Genome Biol* 15:554. MAGeCK; original alpha-RRA algorithm.
- Li W et al. 2015. *Genome Biol* 16:281. MAGeCK-VISPR; QC + visualization.
- Wang B et al. 2019. *Nat Protoc* 14:756. MAGeCKFlute pathway analysis.
- Robson J et al. 2024. *bioRxiv*. MAGeCK normalization benchmark.
- Joung J et al. 2017. *Nat Protoc* 12:828. Screen protocol.
- Hart T et al. 2017. *G3* 7:2719. CEGv2/NEGv1 reference sets for benchmarking.

## Related Skills

- crispr-screens/screen-qc - Pre-MAGeCK QC; gates whether to use median or control normalization
- crispr-screens/jacks-analysis - Joint efficacy + essentiality; alternative to MLE
- crispr-screens/bagel-essentiality - BAGEL2 Bayes-factor calling; alternative for essentiality
- crispr-screens/drugz-chemogenomic - drugZ for chemogenomic screens; alternative for drug screens
- crispr-screens/hit-calling - Cross-method decision tree
- crispr-screens/copy-number-correction - CRISPRcleanR / Chronos for cancer-line CN correction
- crispr-screens/batch-correction - Multi-batch design-matrix construction
- pathway-analysis/gsea - Downstream GSEA on ranked hits
- pathway-analysis/go-enrichment - GO enrichment of hit lists
