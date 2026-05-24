---
name: bio-workflows-crispr-screen-pipeline
description: End-to-end pooled and single-cell CRISPR screen analysis from FASTQ to hit genes. Orchestrates library design QC, guide counting, six-stage screen QC (plasmid Gini, replicate Pearson, CEGv2 PR-AUC, copy-number artifact), method-appropriate hit calling across MAGeCK RRA/MLE, BAGEL2, drugZ, JACKS, and Chronos, cancer-cell-line copy-number correction (CRISPRcleanR / Chronos), batch correction for multi-batch screens, and the specialized branches for combinatorial paralog screens, single-cell Perturb-seq, base-editor variant-function screens, prime-editor screens, and in vivo bottleneck-aware screens. Use when analyzing any pooled CRISPR screen end-to-end, choosing the correct hit-calling method by experimental design, integrating copy-number correction into the pipeline, or branching the workflow for single-cell, combinatorial, base-editor, prime-editor, or in vivo variants.
tool_type: mixed
primary_tool: MAGeCK
workflow: true
depends_on:
  - crispr-screens/library-design
  - crispr-screens/screen-qc
  - crispr-screens/mageck-analysis
  - crispr-screens/bagel-essentiality
  - crispr-screens/drugz-chemogenomic
  - crispr-screens/jacks-analysis
  - crispr-screens/hit-calling
  - crispr-screens/copy-number-correction
  - crispr-screens/batch-correction
  - crispr-screens/crispresso-editing
  - crispr-screens/base-editing-analysis
  - crispr-screens/prime-editing-screens
  - crispr-screens/perturb-seq-analysis
  - crispr-screens/combinatorial-screens
  - crispr-screens/in-vivo-screens
qc_checkpoints:
  - after_counting: ">65% mapping rate; <0.5% zero-count in plasmid; Gini <0.1 on plasmid"
  - after_qc: "Replicate Pearson on log-counts >0.85; Spearman >0.7; CEGv2 PR-AUC >0.7"
  - after_cn_correction: "Spearman ρ between CN and gene LFC abs <0.05 post-correction"
  - after_hit_calling: "Tier-1 hits = 3-method consensus; Tier-2 = 2 of 3; Tier-3 = single-method exploratory"
---

## Version Compatibility

Reference examples tested with: MAGeCK 0.5.9+, BAGEL2 1.0.5+, drugZ Aug 2019+, JACKS 0.2.0+, Chronos 2.0+, CRISPRcleanR 3.0+ (R), Pertpy 0.6+, PRIDICT2, CRISPResso2 2.2.14+, MAGeCKFlute 2.0+, pandas 2.2+, numpy 1.26+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mageck --version`, `BAGEL.py fc --help`, `drugz -h`, `CRISPResso --version`
- Python: `pip show pertpy mageck-vispr jacks chronos-cn`
- R: `packageVersion('CRISPRcleanR')`, `packageVersion('MAGeCKFlute')`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## CRISPR Screen Pipeline

**"Analyze my pooled or single-cell CRISPR screen end-to-end"** -> Pick the screen design branch, run guide counting, audit six QC stages, apply copy-number and batch correction as needed, run the design-matched hit-calling method, and consolidate across methods for high-confidence hits.

## Pipeline Branches by Screen Design

```
                    Library Design ([[library-design]])
                              |
                              v
                FASTQ Files -> mageck count -> count matrix
                              |
                              v
                Six-Stage QC ([[screen-qc]])
                              |
        +---------------------+---------------------+
        |                                            |
        v                                            v
  Cancer cell line?                          Non-cancer?
  Apply CN correction                        No CN correction needed
  ([[copy-number-correction]])
        |                                            |
        +---------------------+---------------------+
                              v
                  Multi-batch? Apply batch covariate
                  ([[batch-correction]])
                              |
                              v
                 Pick hit-calling method by design ([[hit-calling]])
                              |
        +-----------+---------+---------+-----------+-----------+
        |           |         |         |           |           |
        v           v         v         v           v           v
    2-cond       Time      Drug      Essential   Multi-       Specialized
    MAGeCK RRA   MAGeCK    drugZ     BAGEL2      screen       (PE/BE/SC/
                 MLE                              JACKS or     in vivo/
                                                  Chronos      combinat)
        |           |         |         |           |           |
        +-----------+---------+---------+-----------+-----------+
                              v
                   Tier-based consensus
                              v
                Orthogonal validation
```

## Step 1: Library Design and Pre-Screen Validation

Reference [[library-design]] for full library composition. Verify before sequencing:
- Plasmid pool Gini <0.1 (Joung 2017 Nat Protoc 12:828)
- >=99% guides detected at >25 reads/guide
- Skew (p90/p10) <2
- NTCs comprise ~1% of library; CEGv2 reference essentials + NEGv1 non-essentials included

## Step 2: Guide Counting

```bash
mageck count \
    --list-seq library.csv \
    --sample-label Plasmid,Day0,Veh_r1,Veh_r2,Drug_r1,Drug_r2 \
    --fastq Plasmid.fq.gz Day0.fq.gz Veh_r1.fq.gz Veh_r2.fq.gz Drug_r1.fq.gz Drug_r2.fq.gz \
    --norm-method median \
    --output-prefix experiment \
    --trim-5 CACCG
```

For Cas12a libraries (Inzolia, in4mer): see [[combinatorial-screens]]. For 10X single-cell direct capture: use cellranger-arc or pertpy-aware counting; see [[perturb-seq-analysis]].

## Step 3: Six-Stage Quality Control

```python
import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_curve, auc

counts = pd.read_csv('experiment.count.txt', sep='\t', index_col=0)
genes = counts['Gene']
count_matrix = counts.drop('Gene', axis=1)

def gini(x):
    x = np.sort(x[x > 0].astype(float))
    if x.size == 0:
        return np.nan
    n = x.size
    cumx = np.cumsum(x)
    return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n

per_sample = pd.DataFrame({
    'pct_zero': (count_matrix == 0).sum() / len(count_matrix) * 100,
    'gini': count_matrix.apply(gini),
    'reads_per_sgrna': count_matrix.sum() / len(count_matrix),
})

log_counts = np.log10(count_matrix + 1)
pearson = log_counts.corr()
print(per_sample)
print('Replicate Pearson:', pearson.values[pearson.values < 1].mean())
```

Hard gates from [[screen-qc]]:
- Plasmid Gini <0.1; endpoint <0.3 (or <0.55 for heavy drug screens)
- Replicate Pearson on log-counts >0.85
- CEGv2 PR-AUC >0.7 against Hart 2017 reference essential gene set
- Reads per sgRNA per sample >=300 (DepMap convention)

## Step 4: Copy-Number Correction (Cancer Cell Lines Only)

If screening in a cancer cell line, apply CRISPRcleanR (unsupervised, no CN profile needed) or Chronos (joint with CN profile). Required to remove Aguirre 2016 / Munoz 2016 amplicon artifact.

```r
library(CRISPRcleanR)
data(KY_Library_v1.0)
norm <- ccr.NormfoldChanges(read.table('experiment.count.txt', header=TRUE, sep='\t'),
                              min_reads = 30, EXPname = 'screen',
                              libraryAnnotation = KY_Library_v1.0)
gw_lfc <- ccr.logFCs2chromPos(norm$norm_fold_changes, KY_Library_v1.0)
cleaned <- ccr.GWclean(gw_lfc, display = TRUE, label = 'screen')
corrected_counts <- ccr.correctCounts(my_screen = norm, correction = cleaned,
                                        outprefix = 'screen_cleanr',
                                        libraryAnnotation = KY_Library_v1.0)
# Feed corrected counts into MAGeCK / BAGEL2 / drugZ downstream
```

For DepMap-scale panels with longitudinal data + matched CN, use Chronos. See [[copy-number-correction]].

## Step 5: Batch Correction (Multi-Batch Screens)

For multi-batch screens, add batch as a covariate in MAGeCK MLE rather than pre-correcting with ComBat. See [[batch-correction]] for full decision tree.

## Step 6: Method-Matched Hit Calling

### 6a. Two-condition essentiality (MAGeCK RRA or BAGEL2)

```bash
mageck test \
    --count-table experiment.count.txt \
    --treatment-id Day14_r1,Day14_r2,Day14_r3 \
    --control-id Day0 \
    --norm-method median \
    --output-prefix essentiality_rra
```

```bash
BAGEL.py fc -i experiment.count.txt -o foldchange.txt -c Day0 --min-reads 30
BAGEL.py bf -i foldchange.txt -o bayes_factor.txt -e CEGv2.txt -n NEGv1.txt \
    -c Day14_r1,Day14_r2,Day14_r3 -k 1000
```

### 6b. Time-course / multi-condition (MAGeCK MLE)

```bash
mageck mle --count-table experiment.count.txt --design-matrix design.txt \
    --output-prefix timecourse_mle --norm-method median
```

### 6c. Drug-modifier (drugZ)

```bash
python drugz.py \
    -i experiment.count.txt \
    -o drugz_output.txt \
    -c Veh_r1,Veh_r2,Veh_r3 \
    -x Drug_r1,Drug_r2,Drug_r3 \
    -p 5
```

drugZ requires vehicle as control, not Day-0. See [[drugz-chemogenomic]].

### 6d. Multi-screen joint analysis (JACKS)

```bash
python run_JACKS.py experiment.count.txt replicatemap.txt guidemap.txt \
    --rep_hdr Replicate --sample_hdr Sample --ctrl_sample_hdr Control \
    --sgrna_hdr sgRNA --gene_hdr Gene --outprefix jacks_out --apply_w_hp
```

### 6e. Cancer cell-line panels (Chronos)

```python
import chronos
model = chronos.Chronos(sequence_map=sequence_map, guide_gene_map=guide_gene_map,
                          reads=counts_df, copy_number=cn_df)
model.train(n_steps=2000)
gene_effects = model.gene_effect()
```

DepMap quarterly standard; handles CN bias + screen quality + longitudinal jointly.

## Step 7: Tier-Based Consensus

```python
mageck = pd.read_csv('essentiality_rra.gene_summary.txt', sep='\t')[['id', 'neg|fdr']].rename(
    columns={'id': 'gene', 'neg|fdr': 'mageck_neg_fdr'})
bagel = pd.read_csv('bayes_factor.txt', sep='\t')[['GENE', 'BF']].rename(
    columns={'GENE': 'gene', 'BF': 'bagel_bf'})
drugz_df = pd.read_csv('drugz_output.txt', sep='\t')[['GENE', 'fdr_synth']].rename(
    columns={'GENE': 'gene', 'fdr_synth': 'drugz_synth_fdr'})

merged = mageck.merge(bagel, on='gene', how='outer').merge(drugz_df, on='gene', how='outer')
merged['mageck_hit'] = merged['mageck_neg_fdr'] < 0.05
merged['bagel_hit'] = merged['bagel_bf'] > 6
merged['drugz_hit'] = merged['drugz_synth_fdr'] < 0.05
merged['tier'] = merged[['mageck_hit', 'bagel_hit', 'drugz_hit']].astype(int).sum(axis=1)
tier1 = merged[merged['tier'] >= 3]
tier2 = merged[merged['tier'] == 2]
```

## Specialized Branches

| Screen design | Specialized workflow |
|---------------|----------------------|
| Single-cell Perturb-seq / CROP-seq / Multiome | [[perturb-seq-analysis]] -- Pertpy + Mixscape + SCEPTRE |
| Combinatorial paralog (Cas12a Inzolia / Big Papi) | [[combinatorial-screens]] -- GI scoring; synthetic-lethal identification |
| Base-editor variant-function (Hanna 2021 style) | [[base-editing-analysis]] + [[crispresso-editing]] |
| Prime-editor variant installation | [[prime-editing-screens]] -- PRIDICT2 pegRNA design |
| In vivo tumor / immune screens | [[in-vivo-screens]] -- focused library; per-animal meta-analysis |

## Visualization

```python
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10, 8))
gene_summary = pd.read_csv('essentiality_rra.gene_summary.txt', sep='\t')
sig = gene_summary['neg|fdr'] < 0.05
ax.scatter(gene_summary.loc[~sig, 'neg|lfc'],
            -np.log10(gene_summary.loc[~sig, 'neg|fdr'].clip(lower=1e-10)),
            c='lightgray', alpha=0.5, s=10)
ax.scatter(gene_summary.loc[sig, 'neg|lfc'],
            -np.log10(gene_summary.loc[sig, 'neg|fdr'].clip(lower=1e-10)),
            c='red', alpha=0.7, s=18)
ax.axhline(-np.log10(0.05), ls='--', c='black', lw=0.5)
ax.set_xlabel('Log2 Fold Change')
ax.set_ylabel('-Log10(FDR)')
plt.savefig('volcano.png', dpi=150)
```

MAGeCKFlute R package provides one-shot FluteRRA / FluteMLE dashboards with KEGG/Reactome enrichment.

## Output Files

| File | Source step | Description |
|------|-------------|-------------|
| experiment.count.txt | mageck count | Raw count matrix |
| experiment.countsummary.txt | mageck count | Per-sample Gini, mapping, % zero |
| screen_cleanr_corrected_counts.txt | CRISPRcleanR | CN-corrected counts (cancer lines) |
| essentiality_rra.gene_summary.txt | mageck test | Gene-level RRA scores |
| bayes_factor.txt | BAGEL2 | Per-gene Bayes factors |
| drugz_output.txt | drugZ | sumZ, normZ, per-direction FDR |
| jacks_out_gene_JACKS_results.txt | JACKS | Gene effect + sgRNA efficacy |
| tier_consensus.csv | Custom aggregation | Tier-1/2/3 hits across methods |

## Related Skills

- crispr-screens/library-design - Library composition and design rules
- crispr-screens/screen-qc - Six-stage QC + CEGv2 PR-AUC
- crispr-screens/mageck-analysis - MAGeCK RRA + MLE detail
- crispr-screens/bagel-essentiality - BAGEL2 Bayes factor essentiality
- crispr-screens/drugz-chemogenomic - drugZ for drug-modifier screens
- crispr-screens/jacks-analysis - Joint multi-screen analysis with shared efficacy
- crispr-screens/hit-calling - Cross-method decision tree + reconciliation
- crispr-screens/copy-number-correction - CRISPRcleanR / CERES / Chronos
- crispr-screens/batch-correction - Multi-batch design matrix
- crispr-screens/crispresso-editing - CRISPResso2 editing quantification
- crispr-screens/base-editing-analysis - Variant-function BE screens
- crispr-screens/prime-editing-screens - PRIDICT2 pegRNA design
- crispr-screens/perturb-seq-analysis - Single-cell screen analysis
- crispr-screens/combinatorial-screens - Cas12a multiplex + GI scoring
- crispr-screens/in-vivo-screens - Bottleneck-aware in vivo design
- pathway-analysis/go-enrichment - Functional enrichment of hits
- pathway-analysis/gsea - Pre-ranked GSEA on hit lists
