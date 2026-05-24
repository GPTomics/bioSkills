---
name: bio-crispr-screens-hit-calling
description: Cross-method decision tree for calling hits in pooled CRISPR screens. Catalogs statistical models (MAGeCK RRA, MAGeCK MLE, BAGEL2, drugZ, JACKS, Chronos, CERES), experimental designs each is built for, failure modes outside design domain, reconciliation when methods disagree, multiple-testing and effect-size thresholds, the order of operations (count -> QC -> CN-correct -> hit-call -> validate), the second-best-sgRNA conservative rule, and consensus-hit strategy. Use when choosing among MAGeCK / BAGEL2 / drugZ / JACKS / Chronos for a given design, reconciling disagreement across two or three methods on the same screen, deciding whether to require consensus, gating downstream validation by hit-confidence tier, or interpreting unstable hit lists across reruns.
tool_type: mixed
primary_tool: MAGeCK
---

## Version Compatibility

Reference examples tested with: MAGeCK 0.5.9+, BAGEL2 1.0+, drugZ Aug 2019+, JACKS 0.2.0+, Chronos 2.0+ (DepMap), CERES 1.0+, pandas 2.2+, numpy 1.26+, scipy 1.12+, statsmodels 0.14+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mageck --version`, `BAGEL.py --version`, `drugz --help`
- Python: `pip show jacks chronos-cn`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Hit Calling Decision Tree

**"Identify significant hits in my CRISPR screen"** -> Choose the analysis method that matches the experimental design, statistical assumptions, and quality grade of the screen. Reconcile across methods when high-stakes hits must be validated.

The five primary hit-calling methods cover non-overlapping niches; the decision is not "which is best" but "which matches the design."

| Design / question | Primary method | Why | Secondary check |
|--------------------|----------------|-----|------------------|
| Two-condition essentiality, one cell line, no CN concerns | MAGeCK RRA | Robust, fast, gold-standard for ranked analysis | BAGEL2 (Bayes factor on same data) |
| Time course (3+ timepoints) | MAGeCK MLE | RRA cannot model multi-condition | JACKS (efficacy-aware) |
| Multi-cell-line panel (cancer dependency) | Chronos | Models CN bias + screen quality jointly | MAGeCK MLE per line + meta-analysis |
| Drug screen (vehicle vs drug) | drugZ | Bidirectional Z; vehicle-anchored | MAGeCK MLE with dose covariate |
| Multi-screen joint, same library | JACKS | Joint efficacy + 2.5x sample-size reduction | MAGeCK MLE; results should converge |
| Essentiality classification with reference sets | BAGEL2 | Bayes factor with CEGv2/NEGv1 calibration | MAGeCK RRA |
| Combinatorial / paired guide | MAGeCK MLE with GI scoring | Models interaction term; see [[combinatorial-screens]] | Custom GI scoring |
| Single-cell perturbation (Perturb-seq) | SCEPTRE | NB GLM + permutation; see [[perturb-seq-analysis]] | Mixscape pre-filter |
| Cancer-line copy-number screen | Chronos (preferred) or CERES | Joint CN-bias + gene-effect modeling; see [[copy-number-correction]] | CRISPRcleanR pre-hoc + MAGeCK |

## Statistical Models Compared

| Method | Year | Statistical model | Tests | Best for | Fails when |
|--------|------|-------------------|-------|----------|------------|
| MAGeCK RRA | 2014 | NB per-sgRNA -> alpha-RRA per gene | Two-sided | General two-condition | >40% guides change (median norm breaks); time course; cancer-line CN |
| MAGeCK MLE | 2015 | NB GLM with design matrix; per-gene beta | Wald per condition | Multi-condition / time course | Cell-line specific essentiality; CN bias |
| BAGEL2 | 2021 | Bayes factor from log-likelihood ratio | Essential vs non-essential | Essentiality classification | Non-essentiality screens; drug screens |
| drugZ | 2019 | Bidirectional Z-score on guide-level LFC | Sensitizer vs suppressor | Drug-modifier / chemogenomic | Essentiality (no biological prior); time-course |
| JACKS | 2019 | Variational Bayes: LFC = gene * efficacy | Per-gene posterior | Multi-screen joint, library calibration | Single screen; cross-chemistry |
| Chronos | 2021 | Cell-population dynamics ODE + NB | Gene effect adjusted for screen quality | Cancer-line panels, longitudinal | Single screen; non-cancer applications |
| CERES | 2017 | Nonlinear model decoupling CN-bias from gene effect | Per-gene effect | Cancer-line panel with CN profile | Superseded by Chronos at DepMap |

## RRA vs MLE Within MAGeCK

| Property | RRA (`mageck test`) | MLE (`mageck mle`) |
|----------|----------------------|---------------------|
| Conditions supported | 2 | Multiple (design matrix) |
| Statistical test | Robust rank aggregation | Wald on beta from NB GLM |
| Output | neg/pos score, FDR per direction | beta per condition |
| sgRNA efficiency | Not modeled (optional fixed input) | Modeled via `--sgrna-efficiency` |
| Outlier robustness | High (rank-based) | Lower (likelihood-based) |
| Best for | Standard 2-condition screen | Time course, drug screen, multi-cell-line, paired |
| Speed | Fast | Slow (per-gene optimization) |

## Algorithmic Taxonomy: Why Each Was Built

| Method | Designed to solve |
|--------|--------------------|
| MAGeCK RRA | First robust statistical framework for CRISPR-screen ranking; alpha-RRA borrowed from RRA in microarray meta-analysis |
| MAGeCK MLE | Extend MAGeCK to multi-condition; explicit beta scores allow direct LFC interpretation |
| BAGEL2 | Reference-set-anchored Bayesian classification; precision-recall calibrated; tumor-suppressor sensitivity (BAGEL1 was uni-directional) |
| drugZ | Drug-modifier screens have low effect sizes and need bidirectional sensitivity; STARS/MAGeCK miss synthetic-lethal hits |
| JACKS | Sample-size reduction via library-shared efficacy; library calibration as side product |
| Chronos | DepMap-scale (1000+ cell lines, billions of cell-divisions) needs population-dynamics model; CN bias + screen quality first-class |
| CERES | First to formally decouple CN from gene effect at DepMap scale; superseded but historically important |

## Run All Five on the Same Data (Consensus Strategy)

**Goal:** For high-stakes hits (drug-target nomination, paper-level claims), require agreement across 2-3 orthogonal methods.

**Approach:** Run MAGeCK + BAGEL2 + (drugZ or JACKS) on the same count matrix; rank by each; classify hits as called by 1, 2, or 3 methods.

```python
import pandas as pd

def consensus_hits(mageck_path, bagel_path, drugz_path,
                   mageck_fdr_thresh=0.05, bagel_bf_thresh=5, drugz_fdr_thresh=0.05):
    '''Build consensus across MAGeCK / BAGEL2 / drugZ on the same screen.
    Each hit gets a count of supporting methods.'''
    mageck = pd.read_csv(mageck_path, sep='\t')[['id', 'neg|fdr']].rename(columns={'id': 'gene', 'neg|fdr': 'mageck_neg_fdr'})
    bagel = pd.read_csv(bagel_path, sep='\t')[['GENE', 'BF']].rename(columns={'GENE': 'gene', 'BF': 'bagel_bf'})
    drugz = pd.read_csv(drugz_path, sep='\t')[['GENE', 'fdr_synth']].rename(columns={'GENE': 'gene', 'fdr_synth': 'drugz_synth_fdr'})
    merged = mageck.merge(bagel, on='gene', how='outer').merge(drugz, on='gene', how='outer')
    merged['mageck_hit'] = merged['mageck_neg_fdr'] < mageck_fdr_thresh
    merged['bagel_hit'] = merged['bagel_bf'] > bagel_bf_thresh
    merged['drugz_hit'] = merged['drugz_synth_fdr'] < drugz_fdr_thresh
    merged['consensus_count'] = (merged[['mageck_hit', 'bagel_hit', 'drugz_hit']].astype(int)).sum(axis=1)
    return merged.sort_values('consensus_count', ascending=False)
```

**Confidence tiers:**

| Tier | Definition | Validation requirement |
|------|------------|-------------------------|
| Tier 1 (high) | Called by 3/3 methods | Arrayed validation; orthogonal modality (CRISPRi if originally Cas9) |
| Tier 2 (medium) | Called by 2/3 methods | Arrayed validation in matched line |
| Tier 3 (exploratory) | Called by 1/3 methods | Treat as hypothesis; further screens before publication |

## Reconciliation: When Two Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MAGeCK significant, BAGEL2 not | BAGEL2 trained on CEGv2/NEGv1; gene is essential but not in reference | Trust MAGeCK; flag for follow-up |
| BAGEL2 significant, MAGeCK not | BAGEL2 has tumor-suppressor sensitivity MAGeCK lacks | Investigate sgrna_summary for one weak guide |
| MAGeCK significant, JACKS not | JACKS down-weighted one outlier guide | Trust JACKS if guides agree; outlier may be off-target |
| Chronos and MAGeCK disagree on cancer line | Chronos accounts for CN; MAGeCK does not | Trust Chronos; apply [[copy-number-correction]] |
| drugZ significant, MAGeCK not on drug screen | drugZ bidirectional Z is more sensitive | Trust drugZ for chemogenomic; MAGeCK may miss small effects |
| MAGeCK MLE significant, MAGeCK RRA not in 2-condition | Beta-score effect size is significant but rank-based not | Trust MLE if guides consistent; RRA may be over-conservative |
| All methods disagree | Either no real biology or all methods are mis-applied | Stop. Re-audit QC; check chemistry / library / design matrix |

## Second-Best sgRNA Conservative Rule

**Goal:** Reduce false positives from single outlier sgRNAs by requiring the second-most-extreme guide per gene to also be a hit.

**Approach:** For each gene, sort sgRNAs by LFC; require the second-best LFC to exceed a threshold. Rejects genes that depend on one extreme guide.

```python
def second_best_lfc(sgrna_lfc_df, genes_series, direction='neg'):
    '''Return per-gene LFC of the second-best sgRNA in the direction of interest.
    For dropout (direction="neg"), second-most-negative LFC.'''
    results = []
    for gene in genes_series.unique():
        gene_lfc = sgrna_lfc_df[genes_series == gene].sort_values()
        if direction == 'neg':
            second = gene_lfc.iloc[1] if len(gene_lfc) >= 2 else gene_lfc.iloc[0]
        else:
            second = gene_lfc.iloc[-2] if len(gene_lfc) >= 2 else gene_lfc.iloc[-1]
        results.append({'gene': gene, 'second_best_lfc': second})
    return pd.DataFrame(results)
```

**Rule:** A high-confidence hit has second-best LFC also passing the threshold. A guide-of-one hit has only one extreme guide and should be flagged for orthogonal validation. This rule predates JACKS and is implicit in MAGeCK RRA but explicit elsewhere.

## Multiple-Testing Correction Conventions

| Method | Native correction | Cross-method comparison |
|--------|---------------------|---------------------------|
| MAGeCK RRA | BH per direction | `neg|fdr`, `pos|fdr` |
| MAGeCK MLE | BH per condition | `<cond>|fdr` |
| BAGEL2 | Bootstrap BF; reports BF threshold | BF > 6 ≈ FDR 0.05 (Hart 2017 calibration) |
| drugZ | BH per direction | `fdr_synth`, `fdr_supp` |
| JACKS | Posterior probability + BH | `fdr_log10` (log10 FDR) |
| Chronos | DepMap gene-effect probability | `effect_probability` |

**Reconciliation:** BF >6 in BAGEL2 corresponds roughly to FDR 0.05 in MAGeCK; this calibration was confirmed in Hart 2017 G3 by overlap with CEGv2. drugZ FDR is per-direction; the `fdr_synth` and `fdr_supp` columns are independent BH corrections.

## Order of Operations

```
1. Library design (see library-design)         <- design quality dictates hit calling
2. Plasmid pool sequencing                     <- baseline; non-negotiable
3. Run screen at MOI 0.3, 500x coverage
4. Sequence endpoint
5. Run mageck count                            <- generates raw + normalized counts
6. Screen QC (see screen-qc)                   <- gates downstream method choice
7. Copy-number correction if cancer line       <- CRISPRcleanR or Chronos; see copy-number-correction
8. Batch correction if multi-batch             <- see batch-correction
9. Hit calling (this skill)                    <- choose method by design
10. Consensus across 2-3 methods               <- for high-stakes hits
11. Orthogonal validation                      <- arrayed; different chemistry
12. Pathway analysis                           <- see pathway-analysis/gsea
```

## Custom z-score Hit Calling (when standard tools don't fit)

**Goal:** Compute gene-level z-scores when neither MAGeCK nor BAGEL2 fits the experimental design.

**Approach:** RPM-normalize, compute per-sgRNA log2 fold-changes, aggregate to gene level, derive z-score from the null distribution of non-targeting controls (cleanest) or all genes (assumes <40% changing), apply BH correction.

```python
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def custom_zscore_hit_calling(counts_df, ctrl_cols, treat_cols, genes_series, ntc_genes=None):
    '''Z-score gene-level hit calling. If ntc_genes provided, null derived from NTCs only;
    otherwise from all genes (assumes <40% changing).'''
    def rpm(df):
        return df.div(df.sum(axis=0), axis=1) * 1e6
    ctrl_rpm = rpm(counts_df[ctrl_cols])
    treat_rpm = rpm(counts_df[treat_cols])
    lfc_per_sgrna = np.log2((treat_rpm.mean(axis=1) + 1) / (ctrl_rpm.mean(axis=1) + 1))
    gene_lfc = pd.DataFrame({'gene': genes_series, 'lfc': lfc_per_sgrna}).groupby('gene')['lfc'].agg(['mean', 'std', 'count'])
    gene_lfc.columns = ['mean_lfc', 'std_lfc', 'n_sgrnas']
    if ntc_genes is not None:
        null = gene_lfc.loc[gene_lfc.index.isin(ntc_genes), 'mean_lfc']
        null_mean, null_std = null.median(), null.std()
    else:
        null_mean = gene_lfc['mean_lfc'].median()
        null_std = gene_lfc['mean_lfc'].std()
    gene_lfc['z'] = (gene_lfc['mean_lfc'] - null_mean) / null_std
    gene_lfc['p'] = 2 * stats.norm.sf(np.abs(gene_lfc['z']))
    gene_lfc['fdr'] = multipletests(gene_lfc['p'], method='fdr_bh')[1]
    return gene_lfc.sort_values('z')
```

## Failure Modes

### MAGeCK and BAGEL2 disagree by 200+ hits at FDR 0.05

**Trigger:** Heavy-selection screen (>40% guides change) or cancer-line CN bias.
**Mechanism:** MAGeCK median normalization breaks; BAGEL2 is robust due to reference-set anchoring.
**Symptom:** MAGeCK hit list inflated; BAGEL2 list closer to expected size.
**Fix:** Run MAGeCK with `--norm-method control`; apply CN correction; trust BAGEL2 for essentiality.

### Chronos and MAGeCK disagree at the top 10 in a cancer line

**Trigger:** Top hits are at amplified loci.
**Mechanism:** Chronos models CN bias; MAGeCK does not.
**Symptom:** ERBB2 in HER2+, MYC in MYC-amplified, etc. are top hits in MAGeCK but not Chronos.
**Fix:** Apply [[copy-number-correction]] before MAGeCK or switch to Chronos.

### drugZ and MAGeCK disagree on small-effect drug-modifier screen

**Trigger:** Effect size is small; MAGeCK rank-based test is less sensitive than drugZ bidirectional Z.
**Mechanism:** drugZ specifically optimized for small effects in drug screens (Li & Hart 2019); MAGeCK RRA loses sensitivity at small effects.
**Symptom:** drugZ has ~2x more hits at the same FDR threshold; MAGeCK misses real chemogenomic interactions.
**Fix:** Use drugZ as primary for chemogenomic; MAGeCK as confirmatory. See [[drugz-chemogenomic]].

### JACKS down-weights efficiency, MAGeCK doesn't, disagreement

**Trigger:** A gene has one or two strong sgRNAs and 2-3 weak ones; MAGeCK averages them, JACKS down-weights the weak.
**Mechanism:** JACKS variational Bayes correctly identifies low-efficacy guides; MAGeCK aggregates without this prior.
**Symptom:** Gene is JACKS hit but not MAGeCK.
**Fix:** Inspect per-sgRNA LFC; if strong guides are consistent, JACKS is correct. Validate gene orthogonally.

### Consensus across 3 methods is empty (no hits)

**Trigger:** Either no real biology, or each method has different failure mode being triggered.
**Mechanism:** Screen quality is low; signal-to-noise across all methods is poor.
**Symptom:** Tier 1 consensus list is empty.
**Fix:** Re-audit QC. Check Cas9 selection, MOI, timepoint, library positioning. Re-run screen if QC fails.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| MAGeCK RRA FDR (gene-level) | <0.05 | Li 2014; standard publication |
| MAGeCK RRA LFC | abs(LFC) >1 | 2-fold; biological |
| BAGEL2 Bayes Factor | >6 (≈ FDR 0.05); >12 (≈ FDR 0.005) | Kim & Hart 2021; calibrated against CEGv2 |
| drugZ FDR | <0.05 per direction | Li & Hart 2019 |
| JACKS fdr_log10 | <-1 (FDR <0.1); <-2 (FDR <0.01) | Allen 2019 |
| Chronos gene-effect probability | >0.7 | DepMap convention |
| Tier 1 consensus (3 methods) | 100% agreement | High confidence; minimal validation needed |
| Tier 2 consensus (2 of 3) | 67% agreement | Arrayed validation required |
| Tier 3 (1 method only) | Hypothesis; flag for follow-up | Multiple screens or arrayed required |
| Second-best sgRNA rule | Second-best LFC also passes threshold | Reduces single-guide outliers |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| All genes significant in MAGeCK RRA | Heavy selection breaks median norm | `--norm-method control`; or use BAGEL2 |
| BAGEL2 returns no hits despite known essentials | Wrong reference gene set | Verify CEGv2/NEGv1 files match library |
| drugZ output empty | Used Day 0 as control instead of vehicle | Re-run with vehicle as control |
| Chronos errors out | Missing CN profile for cell line | Use CRISPRcleanR (unsupervised) instead |
| Methods disagree by orders of magnitude | Quality issue or design mismatch | Re-audit QC; reconcile via tier consensus |
| Empty tier 1 consensus | No real biology OR QC failure | Re-audit QC |
| Single-guide-driven hits | Outlier sgRNA | Apply second-best rule; orthogonal validate |

## References

- Li W et al. 2014. *Genome Biol* 15:554. MAGeCK alpha-RRA.
- Li W et al. 2015. *Genome Biol* 16:281. MAGeCK MLE.
- Kim E & Hart T. 2021. *Genome Med* 13:2. BAGEL2.
- Li G et al. 2019. *Genome Med* 11:52. drugZ.
- Allen F et al. 2019. *Genome Res* 29:464. JACKS.
- Dempster J et al. 2021. *Genome Biol* 22:343. Chronos.
- Meyers R et al. 2017. *Nat Genet* 49:1779. CERES.
- Hart T & Moffat J. 2016. *BMC Bioinformatics* 17:164. BAGEL Bayes factor framework.
- Hart T et al. 2017. *G3* 7:2719. CEGv2/NEGv1 calibration.

## Related Skills

- crispr-screens/mageck-analysis - Full MAGeCK RRA + MLE detail
- crispr-screens/bagel-essentiality - Full BAGEL2 detail
- crispr-screens/drugz-chemogenomic - Full drugZ detail for drug screens
- crispr-screens/jacks-analysis - Full JACKS detail and library calibration
- crispr-screens/copy-number-correction - Chronos, CERES, CRISPRcleanR
- crispr-screens/screen-qc - Quality gates that drive method choice
- crispr-screens/library-design - Library type dictates analysis method
- crispr-screens/combinatorial-screens - GI scoring (synthetic lethality)
- crispr-screens/perturb-seq-analysis - SCEPTRE for single-cell screens
- crispr-screens/batch-correction - Multi-batch normalization upstream of hit calling
- pathway-analysis/gsea - Downstream pathway enrichment
- pathway-analysis/go-enrichment - GO enrichment of hit lists
