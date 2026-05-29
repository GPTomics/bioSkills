---
name: bio-crispr-screens-screen-qc
description: Quality control for pooled CRISPR screens covering library representation, Gini index, log-skew, replicate Pearson and Spearman concordance, essentialome precision-recall AUC against CEGv2 (Hart 2017), Cas9 cut-toxicity diagnostics, copy-number amplicon detection (Aguirre 2016 / Munoz 2016), bottleneck propagation through plasmid pool, infection, selection, and endpoint stages, MOI verification, and DepMap-style screen-quality scoring. Use when assessing screen quality before hit calling, deciding whether to repeat or rescue a screen, diagnosing low-confidence hits, choosing between MAGeCK / BAGEL2 / Chronos based on quality grade, picking a normalization strategy from QC signatures, or evaluating whether an in-vivo screen retained adequate library complexity.
tool_type: python
primary_tool: MAGeCK-VISPR
---

## Version Compatibility

Reference examples tested with: MAGeCK 0.5+ (count + VISPR), MAGeCKFlute 2.0+ (R), pandas 2.2+, numpy 1.26+, scikit-learn 1.4+, matplotlib 3.8+, seaborn 0.13+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show mageck` then `mageck count --help`; `pip show mageckflute`
- R: `packageVersion('MAGeCKFlute')` then `?BatchRemove` / `?FluteRRA`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## CRISPR Screen Quality Control

**"Audit my CRISPR screen quality before hit calling"** -> Assess library representation, replicate concordance, depth, drift, and biological signal recovery using DepMap-grade metrics, then decide whether the screen is usable, salvageable, or must be repeated.

- Python: `pandas` + `scikit-learn` for Gini, AUC, PCA; `MAGeCKFlute` (R) for one-shot QC dashboard
- CLI: `mageck count` `--gini`, `--mapping_summary` flags; MAGeCK-VISPR for interactive dashboard

## QC Stage Hierarchy

A pooled screen has six distinct bottlenecks where complexity can collapse. Audit each:

| Stage | Metric | Acceptable threshold | Failure consequence |
|-------|--------|----------------------|----------------------|
| Plasmid pool | Gini, skew, % zero-count guides | Gini <0.1, skew <2, zero <0.5% | Missing guides cannot be screened; dropout indistinguishable from non-coverage |
| Day-0 infection | Library coverage, MOI verification | ≥99% guide detection at 500x cells/sgRNA; MOI 0.3 | Founder effects; polyclonality with high MOI |
| Selection (puro/blast) | % cells surviving, time-course Gini | 30-40% survival at 5-7 days; Gini drift <0.05 | Selection artifact; fast-growers enriched |
| Endpoint | Replicate correlation, depth | Pearson >0.9 on log-counts (MAGeCK-VISPR), Spearman >0.7-0.8, ≥200 reads/sgRNA | Noise dominates; FDR inflates |
| Biological signal | CEGv2 PR-AUC, NEGv1 false-positive rate | PR-AUC >0.7 at FDR 5% (DepMap "passing"); CEGv2 enrichment in top 1k | Screen lacks essentiality signal; hits not credible |
| Copy-number artifact | Amplified-region enrichment, sgRNA-cut-count correlation | No correlation between sgRNA off-target count and depletion | False-positive essentiality at amplicons; ERBB2 in HER2+ etc. |

Each metric below quantifies one of these stages.

## Library Representation Metrics

**Goal:** Detect dropout, oversaturation, and library bottlenecks at each sequencing stage.

**Approach:** Compute per-sample zero-count fraction, low-count fraction (<30 reads, Joung 2017 convention), and percentile-based skew, then track how these change between plasmid -> Day-0 -> endpoint to localize the bottleneck.

```python
import pandas as pd
import numpy as np

def library_representation(counts_df):
    '''Per-sample library coverage diagnostics.
    counts_df: rows = sgRNAs, columns = samples (numeric counts).'''
    out = pd.DataFrame(index=counts_df.columns)
    out['n_sgrnas_detected'] = (counts_df > 0).sum()
    out['pct_zero'] = (counts_df == 0).sum() / len(counts_df) * 100
    out['pct_lowcount'] = (counts_df < 30).sum() / len(counts_df) * 100
    out['median_count'] = counts_df.median()
    out['p10_count'] = counts_df.quantile(0.10)
    out['p90_count'] = counts_df.quantile(0.90)
    out['skew_ratio'] = out['p90_count'] / out['p10_count'].replace(0, np.nan)
    return out

def stage_specific_thresholds():
    '''Joung 2017 + DepMap conventions per stage.'''
    return {
        'plasmid':  {'pct_zero_max': 0.5, 'skew_max': 2.0, 'gini_max': 0.10},
        'day_0':    {'pct_zero_max': 1.0, 'skew_max': 2.5, 'gini_max': 0.12},
        'endpoint': {'pct_zero_max': 5.0, 'skew_max': 10.0, 'gini_max': 0.30},
    }
```

**Interpretation:** Plasmid pool failing Gini <0.1 indicates synthesis or amplification bias; the screen is unfit for use. Endpoint Gini drifting above 0.30 indicates either heavy biological selection (acceptable for strong-phenotype drug screens) or a bottleneck (must be diagnosed). The Day-0 vs plasmid delta isolates whether the issue arose during infection (cloning is unlikely to lose specific guides between extraction and infection -- the change happens in cells).

## Gini Coefficient

**Goal:** Quantify how unevenly reads are distributed across sgRNAs in a single sample.

**Approach:** Sort non-zero counts ascending, compute Gini via the cumulative-fraction formula. Compare against stage-specific thresholds.

```python
def gini(x):
    '''Gini coefficient: 0 = perfect equality, 1 = maximal inequality.
    Uses non-zero counts only; zero-count sgRNAs handled separately by % zero.'''
    x = np.sort(x[x > 0].astype(float))
    if x.size == 0:
        return np.nan
    n = x.size
    cumx = np.cumsum(x)
    return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
```

**Stage-specific thresholds (Joung 2017 Nat Protoc; DepMap quality grades):**

| Stage | Excellent | Acceptable | Concerning | Failure |
|-------|-----------|------------|------------|---------|
| Plasmid pool | <0.10 | <0.15 | 0.15-0.20 | >0.20 |
| Day 0 (post-infection) | <0.12 | <0.18 | 0.18-0.25 | >0.25 |
| Endpoint (post-selection) | <0.30 | <0.40 | 0.40-0.55 | >0.55 |

A Gini that climbs from 0.10 (plasmid) to 0.45 (endpoint) is expected when the screen exerts strong selection (drug, lethal-condition). A Gini that climbs to 0.45 without any biological selection (e.g., a control-vs-control timepoint comparison) indicates technical drift.

## Replicate Concordance

**Goal:** Verify that biological/technical replicates agree before testing for between-condition differences.

**Approach:** Compute pairwise Pearson on log10(counts+1) (MAGeCK-VISPR convention) and Spearman ρ on raw rank, between every replicate pair within a condition. Flag any pair below 0.8 Pearson on log-scale.

```python
def replicate_concordance(counts_df, condition_map):
    '''condition_map: {condition_name: [sample_col1, sample_col2, ...]}.'''
    log_counts = np.log10(counts_df + 1)
    rows = []
    for cond, samples in condition_map.items():
        if len(samples) < 2:
            continue
        for i in range(len(samples)):
            for j in range(i+1, len(samples)):
                r_pearson = log_counts[[samples[i], samples[j]]].corr().iloc[0, 1]
                r_spearman = counts_df[[samples[i], samples[j]]].corr(method='spearman').iloc[0, 1]
                rows.append({'condition': cond, 'rep1': samples[i], 'rep2': samples[j],
                             'pearson_log': r_pearson, 'spearman': r_spearman})
    return pd.DataFrame(rows)
```

**Thresholds (MAGeCK-VISPR + DepMap):**

| Metric | Excellent | Acceptable | Failure |
|--------|-----------|------------|---------|
| Pearson on log10(counts+1) | >0.95 | >0.85 | <0.80 |
| Spearman on raw ranks | >0.85 | >0.70 | <0.60 |

**When Pearson is high but Spearman is low**, a few outlier sgRNAs are driving correlation (one extreme guide dominates). Inspect the scatterplot; typically caused by PCR jackpotting at a single guide. Hit calling should use a method that ranks (RRA, drugZ) rather than one that fits per-sgRNA fold change directly.

## Essentialome Recovery (CEGv2 PR-AUC)

**Goal:** Verify the screen has detectable biological essentiality signal by checking whether known essentials (Hart 2017 CEGv2) drop out faster than known non-essentials (NEGv1).

**Approach:** Compute precision-recall AUC where positives are CEGv2 genes and negatives are NEGv1; the screen "passes" if PR-AUC >0.7 (Hart 2017 calibration).

```python
from sklearn.metrics import precision_recall_curve, auc, roc_auc_score

def essentialome_recovery(gene_lfc_df, cegv2_set, negv1_set):
    '''gene_lfc_df: must have ["gene", "lfc"] columns (gene-level mean LFC, negative = depleted).
    cegv2_set, negv1_set: sets of gene symbols from Hart 2017.'''
    labeled = gene_lfc_df[gene_lfc_df['gene'].isin(cegv2_set | negv1_set)].copy()
    labeled['is_essential'] = labeled['gene'].isin(cegv2_set).astype(int)
    y_score = -labeled['lfc']  # negative LFC = depleted = more essential -> higher score
    precision, recall, _ = precision_recall_curve(labeled['is_essential'], y_score)
    return {
        'pr_auc': auc(recall, precision),
        'roc_auc': roc_auc_score(labeled['is_essential'], y_score),
        'n_essential_detected': labeled['is_essential'].sum(),
        'n_nonessential_detected': (1 - labeled['is_essential']).sum(),
    }
```

**Source / threshold:** Hart 2017 *G3* 7:2719 defines CEGv2 (~684 core essentials) and NEGv1 (~927 non-essentials); both lists at https://github.com/hart-lab/bagel/blob/master/CEGv2.txt and NEGv1.txt. DepMap convention: PR-AUC >0.7 at FDR 5% is the "passing" threshold; <0.5 means the screen has no essentiality signal and is not interpretable.

**When PR-AUC is low despite good Gini and Pearson**: cause is usually one of (a) Cas9 was not selected for before screen start (lots of Cas9-negative cells in the pool diluting signal), (b) puromycin selection truncated too aggressively (over-bottleneck), (c) the timepoint is too early (need 14-21 days for KO + decay + selection to manifest). Each has a different remediation.

## Copy-Number Amplicon Bias Diagnostic

**Goal:** Detect the Aguirre 2016 / Munoz 2016 copy-number artifact where sgRNAs targeting amplified loci appear "essential" purely from DNA-damage burden.

**Approach:** Bin genes by copy number (if known from matched WGS/SNP-array) and check whether mean LFC correlates with CN. Alternatively, count off-target cut sites per sgRNA and check correlation with depletion -- amplified loci share many identical cut sites.

```python
def cn_bias_diagnostic(gene_lfc_df, cn_df):
    '''cn_df: per-gene copy number (from WGS/SNP-array/matched ASCAT).
    Tests whether amplified genes show systematically lower LFC.'''
    merged = gene_lfc_df.merge(cn_df, on='gene')
    bins = pd.qcut(merged['copy_number'], q=5, duplicates='drop')
    bin_lfc = merged.groupby(bins, observed=True)['lfc'].agg(['mean', 'median', 'std', 'count'])
    from scipy.stats import spearmanr
    rho, p = spearmanr(merged['copy_number'], merged['lfc'])
    return {'cn_vs_lfc_rho': rho, 'cn_vs_lfc_p': p,
            'amplified_mean_lfc': merged[merged['copy_number'] > 4]['lfc'].mean(),
            'diploid_mean_lfc': merged[(merged['copy_number'] >= 1.5) & (merged['copy_number'] <= 2.5)]['lfc'].mean(),
            'per_bin': bin_lfc}
```

**Interpretation:** A Spearman ρ < -0.1 between copy number and LFC indicates copy-number artifact. The diagnostic threshold is conservative -- Aguirre 2016 showed that even 4-copy amplifications generate detectable artifact. Remediation: use CRISPRcleanR, CERES, or Chronos (see [[copy-number-correction]]) before hit calling.

## Sequencing Depth Audit

**Goal:** Verify that sequencing depth is sufficient to resolve fold changes at the smallest interesting effect size.

**Approach:** Compute reads/sgRNA per sample and the coefficient of variation (CV) of total reads across samples. Compare against the 200x minimum (Joung 2017) or 300x DepMap convention.

```python
def depth_audit(counts_df):
    '''Verify depth: 200x per sgRNA per sample is Joung 2017 minimum;
    300x is DepMap convention; 500x for low-confidence-hit recovery.'''
    total = counts_df.sum()
    n_sgrnas = len(counts_df)
    depth = total / n_sgrnas
    cv = total.std() / total.mean()
    return pd.DataFrame({'total_reads': total, 'reads_per_sgrna': depth,
                          'depth_grade': np.where(depth < 200, 'FAIL',
                                          np.where(depth < 300, 'CAUTION',
                                          np.where(depth < 500, 'OK', 'EXCELLENT')))}).assign(across_sample_cv=cv)
```

**CV interpretation:** CV >0.5 across samples in total reads indicates demultiplexing imbalance or library-pooling error; even if individual samples pass depth thresholds, the relative count is then biased.

## MOI Verification

**Goal:** Confirm that infection occurred at MOI 0.3-0.5 so that ≤1 sgRNA/cell predominates.

**Approach:** From titration plate (control wells with serial-diluted virus), compute infection efficiency, then verify by qPCR of integrated proviral copy number in the screen pool.

| MOI | P(≥1 sgRNA/cell) | P(≥2 sgRNAs/cell) | Cells with 2+ guides as fraction of infected |
|-----|------------------|--------------------|-----------------------------------------------|
| 0.3 | 26% | 4% | 14% |
| 0.5 | 39% | 9% | 23% |
| 1.0 | 63% | 26% | 41% |

**Decision rule:** Always titrate to 0.3. At 0.5, 14-23% of "perturbed" cells carry combinatorial perturbations that confound single-gene scoring. The Poisson math is non-negotiable -- there is no analytical correction for high-MOI confounding.

## PCA and Batch Effect Detection

**Goal:** Visualize whether samples cluster by biology or by batch.

**Approach:** PCA on log10(counts+1); samples should cluster by condition, not by batch/replicate-day/library-lot.

```python
from sklearn.decomposition import PCA

def screen_pca(counts_df, metadata_df, condition_col='condition'):
    '''metadata_df: rows = samples, columns include condition_col, batch (optional).'''
    log_counts = np.log10(counts_df + 1).T  # samples as rows for PCA
    pca = PCA(n_components=3)
    pcs = pca.fit_transform(log_counts)
    out = pd.DataFrame(pcs, columns=['PC1', 'PC2', 'PC3'], index=counts_df.columns)
    out = out.join(metadata_df)
    return out, pca.explained_variance_ratio_
```

**Interpretation:** If PC1 separates batches, see [[batch-correction]]. If PC1 separates conditions cleanly, the screen has interpretable biology. If neither separates anything, the screen has no signal (failed) or is dominated by technical noise.

## Composite DepMap-Style Quality Score

**Goal:** Generate a single quality grade combining all metrics for pipeline gating.

**Approach:** Z-normalize each metric against DepMap's distribution (Pacini 2024) and aggregate to a "screen quality score." Screens scoring <-1 SD are typically excluded from DepMap.

```python
def composite_qc_score(per_sample_qc):
    '''per_sample_qc: one row per sample with columns from library_representation()
    plus pearson, spearman, pr_auc, depth.'''
    metrics = {
        'gini_inv': 1 - per_sample_qc['gini'],
        'pearson': per_sample_qc['pearson_min_replicate'],
        'pr_auc': per_sample_qc['pr_auc'],
        'depth_log': np.log10(per_sample_qc['reads_per_sgrna']),
        'detected_frac': per_sample_qc['n_sgrnas_detected'] / per_sample_qc['n_sgrnas_total'],
    }
    return pd.DataFrame(metrics).mean(axis=1)
```

This is a pipeline gate, not a publication metric -- DepMap reports separate scores for `gene effect score quality` (Chronos-derived) and `screen quality` (Gini-based). See Pacini 2024 *Nat Commun* for the canonical implementation.

## Failure Modes

### High Gini in plasmid pool despite passing all design rules

**Trigger:** Library was cloned and amplified through too many PCR cycles (>20) or used a high-GC-bias polymerase.
**Mechanism:** Each PCR cycle compounds GC bias by ~5%; high-GC and low-GC guides become non-linear functions of starting abundance.
**Symptom:** Gini >0.15 in plasmid, GC-content stratification of dropout.
**Fix:** Cap PCR at 15 cycles for amplification; use Q5 / NEBNext Ultra II / KAPA HiFi (low-bias); re-sequence post-amp; if still bad, re-clone from glycerol stock.

### Falling PR-AUC across timepoints despite stable Gini

**Trigger:** Cas9 was not selected for before screen start; Cas9-negative cells in the pool dilute essentiality signal.
**Mechanism:** Each Cas9-negative cell carries a sgRNA but no editing; its sgRNA persists despite biological essentiality of the target.
**Symptom:** PR-AUC declines from 0.7 at week 1 to 0.4 at week 3; Gini and Pearson both pass.
**Fix:** Always select Cas9-positive cells (FACS or blast) before infection. For a salvage of an already-run screen, model Cas9-expression heterogeneity as a noise floor and accept reduced sensitivity.

### Apparent essentiality of amplified loci

**Trigger:** Cancer cell line with focal amplification (ERBB2 in SK-BR-3, MYC in colorectal, FGFR1 in head and neck).
**Mechanism:** Aguirre 2016 / Munoz 2016: many simultaneous Cas9 cuts trigger p53-dependent DNA-damage response and G2 arrest; sgRNAs at amplified loci appear depleted independently of target essentiality.
**Symptom:** Hits include genes within known amplicons; sgRNAs with more genome-wide cut sites are more depleted.
**Fix:** Apply CRISPRcleanR pre-hoc or use Chronos/CERES with matched CN profile (see [[copy-number-correction]]). Always required for cancer-cell-line screens, not optional.

### Outlier replicate dragging Pearson down

**Trigger:** One technical replicate had a library-prep failure (low input, PCR jackpot, sequencing-lane swap).
**Mechanism:** Outlier sample has different total reads or different per-sgRNA distribution but passes individual sample QC.
**Symptom:** Pearson between replicates 0.85-0.90 with one pair as outlier; condition-level means look fine.
**Fix:** Drop the outlier replicate; re-derive Pearson on the remaining pair. If only two replicates and one is outlier, the condition lacks replication and must be re-run.

### Low Day-0 coverage from high MOI

**Trigger:** Infection at MOI >0.5.
**Mechanism:** Poisson: at MOI 0.5, 23% of infected cells carry multiple sgRNAs; the "single-perturbation" assumption underlying every analysis method is violated.
**Symptom:** Apparent gene-gene interactions in single-gene screens; gene-level z-scores noisy; Pearson lower than expected for high-quality counts.
**Fix:** No analytical correction. Re-titrate, re-infect at MOI 0.3, re-run screen.

### CRISPRi/a screen with no signal on validated essentials

**Trigger:** Library targets wrong TSS (Ensembl canonical vs FANTOM5 highest-rank).
**Mechanism:** dCas9-KRAB knockdown is maximal within ±100 bp of the actual Pol II loading site; canonical annotation can be off 1-10 kb.
**Symptom:** RPS/RPL/EIF families dropping out as expected (these have clean canonical TSSs) but downstream genes failing; PR-AUC on broader CEGv2 panel drops.
**Fix:** Re-design library against FANTOM5 highest-CAGE-peak TSS (Sanson 2018); for tissue-specific lines, use matched CAGE / GRO-seq.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Plasmid Gini | <0.10 | Joung 2017 *Nat Protoc* 12:828 |
| Plasmid skew ratio (p90/p10) | <2 | Joung 2017 |
| % zero-count sgRNAs (plasmid) | <0.5% | DepMap convention |
| % zero-count sgRNAs (endpoint) | <5% | MAGeCK-VISPR; >5% impairs FDR |
| Replicate Pearson on log10(counts+1) | >0.85 acceptable, >0.95 ideal | Li et al 2014 MAGeCK paper; MAGeCK-VISPR |
| Replicate Spearman | >0.70 | MAGeCKFlute (Wang 2019 Nat Protoc) |
| CEGv2 PR-AUC at FDR 5% | >0.70 passing; >0.85 high quality | Hart 2017 *G3* 7:2719 |
| Reads per sgRNA per sample | ≥200 minimum, 300+ DepMap, 500+ low-effect screens | Joung 2017; Pacini 2024 DepMap |
| Library coverage at infection | 500x cells/sgRNA | Joung 2017; DepMap |
| In-vivo coverage | 200-1000x bottleneck-adjusted | Chen 2015 / Manguso 2017; see [[in-vivo-screens]] |
| MOI at infection | 0.3 strict | Poisson: P(≥2)=4% at 0.3 vs 9% at 0.5 |
| CN-bias Spearman ρ (LFC vs copy number) | abs(ρ) <0.10 | Aguirre 2016 Cancer Discov |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Plasmid Gini >0.2 | PCR over-amplification | Re-sequence; cap PCR at 15 cycles |
| Endpoint PR-AUC <0.5 | Cas9 not selected pre-screen | Select Cas9+; redo if mid-screen |
| Pearson high, Spearman low | A few outlier sgRNAs dominate | Use RRA / rank-based hit calling |
| Replicate Pearson <0.8 | Library-prep failure on one rep | Drop outlier; re-run if singleton |
| % zero increases dramatically Day-0 -> endpoint | Selection bottleneck | Reduce selection pressure; or increase coverage |
| Top hits include amplified-region genes | CN bias | CRISPRcleanR or Chronos |
| MOI verification shows 0.6+ | Over-infected | Re-run at lower MOI; no rescue |
| Detected fraction <90% in Day 0 | Coverage too low | Increase cells; expect drift |

## References

- Joung J et al. 2017. *Nat Protoc* 12:828. Genome-wide screen protocol; coverage and depth conventions.
- Li W et al. 2014. *Genome Biol* 15:554. MAGeCK; original Pearson cutoffs.
- Wang B et al. 2019. *Nat Protoc* 14:756. MAGeCKFlute; QC dashboard.
- Hart T et al. 2017. *G3* 7:2719. CEGv2 / NEGv1 reference essentiality gene sets.
- Aguirre AJ et al. 2016. *Cancer Discov* 6:914. Copy-number amplicon false-essentiality.
- Munoz DM et al. 2016. *Cancer Discov* 6:900. Copy-number gene-independent toxicity.
- Pacini C et al. 2024. *Nat Commun* 15:1230. DepMap screen-quality scoring.
- Meyers RM et al. 2017. *Nat Genet* 49:1779. CERES; mechanism of CN bias.
- Sanson KR et al. 2018. *Nat Commun* 9:5416. Dolcetto/Calabrese TSS rules.
- Dempster JM et al. 2021. *Genome Biol* 22:343. Chronos screen-quality model.

## Related Skills

- crispr-screens/library-design - Compose libraries that pass plasmid QC
- crispr-screens/mageck-analysis - Run MAGeCK count to generate QC inputs
- crispr-screens/copy-number-correction - Remediate Aguirre / Munoz CN artifact
- crispr-screens/batch-correction - Address inter-batch / cell-line confounding
- crispr-screens/hit-calling - Pick method by QC grade
- crispr-screens/in-vivo-screens - In-vivo-specific bottleneck QC
